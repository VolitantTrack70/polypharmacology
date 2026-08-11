//! The cascade query: structure -> similar compounds -> targets -> pathways.
//!
//! Runs in two phases because similarity is not expressible in SQL against our
//! storage (see docs/decisions/0001):
//!
//!   1. chemworker scans the fingerprint index and returns (chembl_id, tanimoto).
//!   2. Postgres expands that set through `binds_to` and `target_pathway`.
//!
//! Phase 2 is a single round trip -- the hit set is passed as an array and
//! UNNESTed, rather than issuing one query per hit.

use std::collections::{HashMap, HashSet};

use axum::extract::State;
use axum::Json;
use serde::{Deserialize, Serialize};

use crate::error::{ApiError, ApiResult};
use crate::AppState;

#[derive(Deserialize)]
pub struct OffTargetRequest {
    /// One of `smiles` or `chembl_id` is required.
    pub smiles: Option<String>,
    pub chembl_id: Option<String>,
    /// Morgan/ECFP4 Tanimoto floor. Defaults to 0.40 -- see docs/decisions/0002.
    pub tanimoto: Option<f32>,
    /// pChEMBL floor. 6.0 == 1 uM.
    pub pchembl: Option<f32>,
    /// Cap on similar compounds pulled through to phase 2.
    pub limit: Option<u32>,
    /// Restrict to a single organism, e.g. "Homo sapiens".
    pub organism: Option<String>,
}

#[derive(Serialize)]
pub struct OffTargetResponse {
    pub query: QueryEcho,
    pub stats: Stats,
    pub nodes: Vec<Node>,
    pub edges: Vec<Edge>,
}

#[derive(Serialize)]
pub struct QueryEcho {
    pub canonical_smiles: String,
    pub tanimoto_cutoff: f32,
    pub pchembl_cutoff: f32,
}

#[derive(Serialize)]
pub struct Stats {
    pub similar_compounds: usize,
    pub targets: usize,
    pub pathways: usize,
    pub compounds_scanned: u32,
    pub search_ms: f32,
}

/// Flat node list, shaped for cytoscape.js. `kind` drives styling client-side.
#[derive(Serialize, Debug, Clone, PartialEq)]
pub struct Node {
    pub id: String,
    pub label: String,
    pub kind: &'static str, // "query" | "compound" | "target" | "pathway"
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tanimoto: Option<f32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub organism: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub biological_domain: Option<String>,
}

#[derive(Serialize, Debug, Clone, PartialEq)]
pub struct Edge {
    pub source: String,
    pub target: String,
    pub kind: &'static str, // "similar_to" | "binds_to" | "participates_in"
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tanimoto: Option<f32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pchembl: Option<f32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub n_measurements: Option<i64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub activity_types: Option<Vec<String>>,
}

#[derive(sqlx::FromRow, Debug, Clone)]
pub struct CascadeRow {
    pub chembl_id: String,
    pub target_chembl_id: String,
    pub pref_name: Option<String>,
    pub organism: Option<String>,
    pub max_pchembl: Option<f32>,
    pub n_measurements: Option<i64>,
    pub activity_types: Option<Vec<String>>,
    pub reactome_id: Option<String>,
    pub pathway_name: Option<String>,
    pub biological_domain: Option<String>,
}

pub const QUERY_NODE_ID: &str = "__query__";

const CASCADE_SQL: &str = r#"
WITH hits AS (
    SELECT * FROM UNNEST($1::text[]) AS t(chembl_id)
)
SELECT
    b.chembl_id,
    b.target_chembl_id,
    t.pref_name,
    t.organism,
    b.max_pchembl,
    b.n_measurements,
    b.activity_types,
    tp.reactome_id,
    tp.pathway_name,
    tp.biological_domain
FROM hits h
JOIN chem.binds_to b       ON b.chembl_id = h.chembl_id
JOIN chem.target   t       ON t.target_chembl_id = b.target_chembl_id
LEFT JOIN chem.target_pathway tp ON tp.target_chembl_id = b.target_chembl_id
WHERE b.max_pchembl >= $2
  AND ($3::text IS NULL OR t.organism = $3)
"#;

/// Result of folding the fanned-out SQL rows into a deduplicated graph.
pub struct Graph {
    pub nodes: Vec<Node>,
    pub edges: Vec<Edge>,
    pub n_compounds: usize,
    pub n_targets: usize,
    pub n_pathways: usize,
}

/// Fold `compound x target x pathway` rows into a node/edge graph.
///
/// The SQL join fans out: one row per (compound, target, pathway) combination,
/// so a compound binding 3 targets that sit in 20 pathways each arrives as 60
/// rows describing 1 compound, 3 targets and up to 60 pathways. Everything
/// here is deduplication.
///
/// Kept separate from the handler so it can be tested without a database.
pub fn build_graph(rows: &[CascadeRow], similarity: &HashMap<String, f32>) -> Graph {
    let mut nodes = vec![Node {
        id: QUERY_NODE_ID.to_string(),
        label: "query structure".into(),
        kind: "query",
        tanimoto: None,
        organism: None,
        biological_domain: None,
    }];
    let mut edges: Vec<Edge> = Vec::new();

    let mut seen_compounds: HashSet<&str> = HashSet::new();
    let mut seen_targets: HashSet<&str> = HashSet::new();
    let mut seen_pathways: HashSet<&str> = HashSet::new();
    let mut seen_binds: HashSet<(&str, &str)> = HashSet::new();
    let mut seen_participates: HashSet<(&str, &str)> = HashSet::new();

    for row in rows {
        let ts = similarity.get(&row.chembl_id).copied();

        if seen_compounds.insert(&row.chembl_id) {
            nodes.push(Node {
                id: row.chembl_id.clone(),
                label: row.chembl_id.clone(),
                kind: "compound",
                tanimoto: ts,
                organism: None,
                biological_domain: None,
            });
            edges.push(Edge {
                source: QUERY_NODE_ID.to_string(),
                target: row.chembl_id.clone(),
                kind: "similar_to",
                tanimoto: ts,
                pchembl: None,
                n_measurements: None,
                activity_types: None,
            });
        }

        if seen_targets.insert(&row.target_chembl_id) {
            nodes.push(Node {
                id: row.target_chembl_id.clone(),
                label: row
                    .pref_name
                    .clone()
                    .unwrap_or_else(|| row.target_chembl_id.clone()),
                kind: "target",
                tanimoto: None,
                organism: row.organism.clone(),
                biological_domain: None,
            });
        }

        if seen_binds.insert((&row.chembl_id, &row.target_chembl_id)) {
            edges.push(Edge {
                source: row.chembl_id.clone(),
                target: row.target_chembl_id.clone(),
                kind: "binds_to",
                tanimoto: None,
                pchembl: row.max_pchembl,
                n_measurements: row.n_measurements,
                activity_types: row.activity_types.clone(),
            });
        }

        // LEFT JOIN: a target with no pathway annotation still yields a row,
        // with the pathway columns null. It must appear as a target node with
        // no outgoing edge, not be dropped.
        if let Some(reactome_id) = &row.reactome_id {
            if seen_pathways.insert(reactome_id) {
                nodes.push(Node {
                    id: reactome_id.clone(),
                    label: row
                        .pathway_name
                        .clone()
                        .unwrap_or_else(|| reactome_id.clone()),
                    kind: "pathway",
                    tanimoto: None,
                    organism: None,
                    biological_domain: row.biological_domain.clone(),
                });
            }
            if seen_participates.insert((&row.target_chembl_id, reactome_id)) {
                edges.push(Edge {
                    source: row.target_chembl_id.clone(),
                    target: reactome_id.clone(),
                    kind: "participates_in",
                    tanimoto: None,
                    pchembl: None,
                    n_measurements: None,
                    activity_types: None,
                });
            }
        }
    }

    Graph {
        n_compounds: seen_compounds.len(),
        n_targets: seen_targets.len(),
        n_pathways: seen_pathways.len(),
        nodes,
        edges,
    }
}

pub async fn off_targets(
    State(state): State<AppState>,
    Json(req): Json<OffTargetRequest>,
) -> ApiResult<Json<OffTargetResponse>> {
    let tanimoto = req.tanimoto.unwrap_or(state.defaults.tanimoto);
    let pchembl = req.pchembl.unwrap_or(state.defaults.pchembl);
    let limit = req.limit.unwrap_or(state.defaults.similar_limit);

    if !(0.0..=1.0).contains(&tanimoto) || tanimoto == 0.0 {
        return Err(ApiError::BadRequest("tanimoto must be in (0, 1]".into()));
    }

    // ---- Phase 1: similarity ------------------------------------------------
    let search = match (&req.smiles, &req.chembl_id) {
        (Some(s), _) if !s.trim().is_empty() => {
            state.chem.search_by_smiles(s, tanimoto, limit).await?
        }
        (_, Some(id)) if !id.trim().is_empty() => {
            state.chem.search_by_chembl_id(id, tanimoto, limit).await?
        }
        _ => {
            return Err(ApiError::BadRequest(
                "provide either `smiles` or `chembl_id`".into(),
            ))
        }
    };

    let similarity: HashMap<String, f32> = search
        .hits
        .iter()
        .map(|h| (h.chembl_id.clone(), h.tanimoto))
        .collect();
    let hit_ids: Vec<String> = search.hits.iter().map(|h| h.chembl_id.clone()).collect();

    let query_echo = QueryEcho {
        canonical_smiles: search.canonical_smiles.clone(),
        tanimoto_cutoff: tanimoto,
        pchembl_cutoff: pchembl,
    };

    // Not an error. Returning the full shape with zero results lets the UI
    // suggest lowering the threshold rather than showing a failure.
    if hit_ids.is_empty() {
        let graph = build_graph(&[], &similarity);
        return Ok(Json(OffTargetResponse {
            query: query_echo,
            stats: Stats {
                similar_compounds: 0,
                targets: 0,
                pathways: 0,
                compounds_scanned: search.searched,
                search_ms: search.elapsed_ms,
            },
            nodes: graph.nodes,
            edges: graph.edges,
        }));
    }

    // ---- Phase 2: graph expansion ------------------------------------------
    let rows: Vec<CascadeRow> = sqlx::query_as(CASCADE_SQL)
        .bind(&hit_ids)
        .bind(pchembl)
        .bind(req.organism.as_deref())
        .fetch_all(&state.db)
        .await?;

    let graph = build_graph(&rows, &similarity);

    Ok(Json(OffTargetResponse {
        query: query_echo,
        stats: Stats {
            similar_compounds: graph.n_compounds,
            targets: graph.n_targets,
            pathways: graph.n_pathways,
            compounds_scanned: search.searched,
            search_ms: search.elapsed_ms,
        },
        nodes: graph.nodes,
        edges: graph.edges,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(
        compound: &str,
        target: &str,
        pathway: Option<&str>,
        pchembl: Option<f32>,
    ) -> CascadeRow {
        CascadeRow {
            chembl_id: compound.into(),
            target_chembl_id: target.into(),
            pref_name: Some(format!("{target} name")),
            organism: Some("Homo sapiens".into()),
            max_pchembl: pchembl,
            n_measurements: Some(3),
            activity_types: Some(vec!["IC50".into()]),
            reactome_id: pathway.map(Into::into),
            pathway_name: pathway.map(|p| format!("{p} pathway")),
            biological_domain: pathway.map(|_| "Signal Transduction".into()),
        }
    }

    fn sim(pairs: &[(&str, f32)]) -> HashMap<String, f32> {
        pairs.iter().map(|(k, v)| (k.to_string(), *v)).collect()
    }

    fn ids_of(g: &Graph, kind: &str) -> Vec<String> {
        let mut v: Vec<String> = g
            .nodes
            .iter()
            .filter(|n| n.kind == kind)
            .map(|n| n.id.clone())
            .collect();
        v.sort();
        v
    }

    fn edges_of<'a>(g: &'a Graph, kind: &str) -> Vec<&'a Edge> {
        g.edges.iter().filter(|e| e.kind == kind).collect()
    }

    #[test]
    fn empty_input_still_yields_the_query_node() {
        let g = build_graph(&[], &sim(&[]));
        assert_eq!(g.nodes.len(), 1);
        assert_eq!(g.nodes[0].kind, "query");
        assert!(g.edges.is_empty());
        assert_eq!(g.n_compounds, 0);
    }

    #[test]
    fn fanned_out_rows_collapse_to_unique_nodes() {
        // One compound, one target, three pathways -> 3 rows, but 1 target node.
        let rows = vec![
            row("C1", "T1", Some("P1"), Some(7.0)),
            row("C1", "T1", Some("P2"), Some(7.0)),
            row("C1", "T1", Some("P3"), Some(7.0)),
        ];
        let g = build_graph(&rows, &sim(&[("C1", 0.9)]));

        assert_eq!(ids_of(&g, "compound"), vec!["C1"]);
        assert_eq!(ids_of(&g, "target"), vec!["T1"]);
        assert_eq!(g.n_pathways, 3);
        // The compound must not get one similar_to edge per row.
        assert_eq!(edges_of(&g, "similar_to").len(), 1);
        assert_eq!(edges_of(&g, "binds_to").len(), 1);
        assert_eq!(edges_of(&g, "participates_in").len(), 3);
    }

    #[test]
    fn shared_target_across_compounds_is_one_node_two_edges() {
        let rows = vec![
            row("C1", "T_SHARED", Some("P1"), Some(8.0)),
            row("C2", "T_SHARED", Some("P1"), Some(6.5)),
        ];
        let g = build_graph(&rows, &sim(&[("C1", 1.0), ("C2", 0.5)]));

        assert_eq!(ids_of(&g, "target"), vec!["T_SHARED"]);
        assert_eq!(edges_of(&g, "binds_to").len(), 2);
        // The shared pathway must not get a duplicate participates_in edge.
        assert_eq!(edges_of(&g, "participates_in").len(), 1);
    }

    #[test]
    fn shared_pathway_across_targets_gets_one_node_two_edges() {
        let rows = vec![
            row("C1", "T1", Some("P_SHARED"), Some(7.0)),
            row("C1", "T2", Some("P_SHARED"), Some(7.0)),
        ];
        let g = build_graph(&rows, &sim(&[("C1", 0.8)]));

        assert_eq!(g.n_pathways, 1);
        assert_eq!(edges_of(&g, "participates_in").len(), 2);
    }

    #[test]
    fn target_without_pathway_survives_the_left_join() {
        let rows = vec![row("C1", "T_ORPHAN", None, Some(7.0))];
        let g = build_graph(&rows, &sim(&[("C1", 0.7)]));

        assert_eq!(ids_of(&g, "target"), vec!["T_ORPHAN"]);
        assert_eq!(g.n_pathways, 0);
        assert!(edges_of(&g, "participates_in").is_empty());
        assert_eq!(edges_of(&g, "binds_to").len(), 1);
    }

    #[test]
    fn tanimoto_is_attached_to_compound_nodes_and_similarity_edges() {
        let rows = vec![row("C1", "T1", Some("P1"), Some(7.0))];
        let g = build_graph(&rows, &sim(&[("C1", 0.517)]));

        let c = g.nodes.iter().find(|n| n.id == "C1").unwrap();
        assert_eq!(c.tanimoto, Some(0.517));
        assert_eq!(edges_of(&g, "similar_to")[0].tanimoto, Some(0.517));
        // Affinity belongs on the binding edge, not the similarity edge.
        assert_eq!(edges_of(&g, "similar_to")[0].pchembl, None);
        assert_eq!(edges_of(&g, "binds_to")[0].pchembl, Some(7.0));
    }

    #[test]
    fn missing_similarity_score_does_not_drop_the_compound() {
        // Defensive: chemworker returned a hit the map somehow lacks.
        let rows = vec![row("C_UNKNOWN", "T1", Some("P1"), Some(7.0))];
        let g = build_graph(&rows, &sim(&[]));

        assert_eq!(ids_of(&g, "compound"), vec!["C_UNKNOWN"]);
        assert_eq!(g.nodes.iter().find(|n| n.id == "C_UNKNOWN").unwrap().tanimoto, None);
    }

    #[test]
    fn target_label_falls_back_to_id_when_name_is_null() {
        let mut r = row("C1", "T1", None, Some(7.0));
        r.pref_name = None;
        let g = build_graph(&[r], &sim(&[("C1", 0.5)]));

        let t = g.nodes.iter().find(|n| n.id == "T1").unwrap();
        assert_eq!(t.label, "T1");
    }

    #[test]
    fn every_edge_endpoint_refers_to_a_real_node() {
        // Cytoscape silently drops edges with dangling endpoints, which would
        // show up as a mysteriously sparse graph rather than an error.
        let rows = vec![
            row("C1", "T1", Some("P1"), Some(7.0)),
            row("C1", "T2", Some("P2"), Some(6.0)),
            row("C2", "T1", None, Some(8.0)),
        ];
        let g = build_graph(&rows, &sim(&[("C1", 0.9), ("C2", 0.6)]));

        let ids: HashSet<&str> = g.nodes.iter().map(|n| n.id.as_str()).collect();
        for e in &g.edges {
            assert!(ids.contains(e.source.as_str()), "dangling source {}", e.source);
            assert!(ids.contains(e.target.as_str()), "dangling target {}", e.target);
        }
    }

    #[test]
    fn node_ids_are_unique() {
        let rows = vec![
            row("C1", "T1", Some("P1"), Some(7.0)),
            row("C1", "T1", Some("P1"), Some(7.0)),
            row("C2", "T1", Some("P1"), Some(7.0)),
        ];
        let g = build_graph(&rows, &sim(&[("C1", 0.9), ("C2", 0.4)]));

        let mut ids: Vec<&str> = g.nodes.iter().map(|n| n.id.as_str()).collect();
        let before = ids.len();
        ids.sort_unstable();
        ids.dedup();
        assert_eq!(ids.len(), before, "duplicate node ids emitted");
    }

    #[test]
    fn stats_match_the_emitted_nodes() {
        let rows = vec![
            row("C1", "T1", Some("P1"), Some(7.0)),
            row("C1", "T2", Some("P2"), Some(7.0)),
            row("C2", "T1", Some("P1"), Some(6.0)),
        ];
        let g = build_graph(&rows, &sim(&[("C1", 0.9), ("C2", 0.5)]));

        assert_eq!(g.n_compounds, ids_of(&g, "compound").len());
        assert_eq!(g.n_targets, ids_of(&g, "target").len());
        assert_eq!(g.n_pathways, ids_of(&g, "pathway").len());
    }
}
