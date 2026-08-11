//! The cascade query: structure -> similar compounds -> targets -> pathways.
//!
//! Runs in two phases because similarity is not expressible in SQL against our
//! storage (see docs/decisions/0001):
//!
//!   1. chemworker scans the fingerprint index and returns (chembl_id, tanimoto).
//!   2. Postgres expands that set through `binds_to` and `target_pathway`.
//!
//! Phase 2 is a single round trip -- the hit set is passed as arrays and
//! UNNESTed, rather than issuing one query per hit.

use std::collections::HashMap;

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
#[derive(Serialize)]
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

#[derive(Serialize)]
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

#[derive(sqlx::FromRow)]
struct CascadeRow {
    chembl_id: String,
    target_chembl_id: String,
    pref_name: Option<String>,
    organism: Option<String>,
    max_pchembl: Option<f32>,
    n_measurements: Option<i64>,
    activity_types: Option<Vec<String>>,
    reactome_id: Option<String>,
    pathway_name: Option<String>,
    biological_domain: Option<String>,
}

const QUERY_NODE_ID: &str = "__query__";

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

pub async fn off_targets(
    State(state): State<AppState>,
    Json(req): Json<OffTargetRequest>,
) -> ApiResult<Json<OffTargetResponse>> {
    let tanimoto = req.tanimoto.unwrap_or(state.defaults.tanimoto);
    let pchembl = req.pchembl.unwrap_or(state.defaults.pchembl);
    let limit = req.limit.unwrap_or(state.defaults.similar_limit);

    if !(0.0..=1.0).contains(&tanimoto) || tanimoto == 0.0 {
        return Err(ApiError::BadRequest(
            "tanimoto must be in (0, 1]".into(),
        ));
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

    let mut nodes = vec![Node {
        id: QUERY_NODE_ID.to_string(),
        label: "query structure".into(),
        kind: "query",
        tanimoto: None,
        organism: None,
        biological_domain: None,
    }];
    let mut edges: Vec<Edge> = Vec::new();

    if hit_ids.is_empty() {
        // Not an error. Returning the shape with zero results lets the UI
        // suggest lowering the threshold rather than showing a failure.
        return Ok(Json(OffTargetResponse {
            query: QueryEcho {
                canonical_smiles: search.canonical_smiles,
                tanimoto_cutoff: tanimoto,
                pchembl_cutoff: pchembl,
            },
            stats: Stats {
                similar_compounds: 0,
                targets: 0,
                pathways: 0,
                compounds_scanned: search.searched,
                search_ms: search.elapsed_ms,
            },
            nodes,
            edges,
        }));
    }

    // ---- Phase 2: graph expansion ------------------------------------------
    let rows: Vec<CascadeRow> = sqlx::query_as(CASCADE_SQL)
        .bind(&hit_ids)
        .bind(pchembl)
        .bind(req.organism.as_deref())
        .fetch_all(&state.db)
        .await?;

    // The join fans out (one row per compound x target x pathway), so
    // deduplicate into node and edge sets as we go.
    let mut seen_compounds = HashMap::new();
    let mut seen_targets = HashMap::new();
    let mut seen_pathways = HashMap::new();
    let mut seen_binds = std::collections::HashSet::new();
    let mut seen_participates = std::collections::HashSet::new();

    for row in &rows {
        let ts = similarity.get(&row.chembl_id).copied();

        if seen_compounds.insert(row.chembl_id.clone(), ()).is_none() {
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

        if seen_targets.insert(row.target_chembl_id.clone(), ()).is_none() {
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

        if seen_binds.insert((row.chembl_id.clone(), row.target_chembl_id.clone())) {
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

        if let Some(reactome_id) = &row.reactome_id {
            if seen_pathways.insert(reactome_id.clone(), ()).is_none() {
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
            if seen_participates.insert((row.target_chembl_id.clone(), reactome_id.clone())) {
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

    Ok(Json(OffTargetResponse {
        query: QueryEcho {
            canonical_smiles: search.canonical_smiles,
            tanimoto_cutoff: tanimoto,
            pchembl_cutoff: pchembl,
        },
        stats: Stats {
            similar_compounds: seen_compounds.len(),
            targets: seen_targets.len(),
            pathways: seen_pathways.len(),
            compounds_scanned: search.searched,
            search_ms: search.elapsed_ms,
        },
        nodes,
        edges,
    }))
}
