-- ============================================================================
-- 002_indexes_and_views.sql
-- Run AFTER bulk load, not before -- building these indexes up front makes a
-- 21M-row COPY several times slower. See services/ingest README.
-- ============================================================================

BEGIN;
SET search_path TO chem, public;

-- --- Lookup paths the API actually uses -------------------------------------
CREATE INDEX IF NOT EXISTS idx_compound_inchikey    ON compound (standard_inchi_key);
CREATE INDEX IF NOT EXISTS idx_compound_max_phase   ON compound (max_phase) WHERE max_phase IS NOT NULL;

CREATE INDEX IF NOT EXISTS idx_protein_gene         ON protein (gene_symbol);
CREATE INDEX IF NOT EXISTS idx_target_tax           ON target (tax_id);
CREATE INDEX IF NOT EXISTS idx_target_component_uni ON target_component (uniprot_accession);

-- The hot path: "given these N compound ids, give me their high-affinity
-- targets". Composite + INCLUDE makes this index-only for the common query.
CREATE INDEX IF NOT EXISTS idx_activity_compound_pchembl
    ON activity (chembl_id, pchembl_value DESC)
    INCLUDE (target_chembl_id, standard_type)
    WHERE pchembl_value IS NOT NULL;

CREATE INDEX IF NOT EXISTS idx_activity_target
    ON activity (target_chembl_id)
    WHERE pchembl_value IS NOT NULL;

CREATE INDEX IF NOT EXISTS idx_protein_pathway_reactome ON protein_pathway (reactome_id);
CREATE INDEX IF NOT EXISTS idx_pathway_hierarchy_child  ON pathway_hierarchy (child_reactome_id);

-- ---------------------------------------------------------------------------
-- binds_to: the aggregated compound->target edge.
--
-- Raw `activity` has many measurements per (compound, target) pair from
-- different assays and papers. The graph edge needs ONE row per pair, so we
-- aggregate -- keeping max pchembl (best observed affinity), the measurement
-- count, and which assay types contributed.
--
-- We deliberately exclude:
--   * confidence_score < 7  -- assay not confidently tied to a single target
--   * standard_relation '>' -- "activity was worse than X", i.e. a non-binder
--   * rows flagged with a data_validity_comment
--
-- 7 is a hard floor baked into the view. best_confidence is retained so
-- callers can raise it per query (the API's min_confidence parameter); it
-- cannot be lowered without rebuilding this view.
-- ---------------------------------------------------------------------------
DROP MATERIALIZED VIEW IF EXISTS binds_to CASCADE;
CREATE MATERIALIZED VIEW binds_to AS
SELECT
    a.chembl_id,
    a.target_chembl_id,
    MAX(a.pchembl_value)::REAL              AS max_pchembl,
    AVG(a.pchembl_value)::REAL              AS mean_pchembl,
    COUNT(*)                                AS n_measurements,
    COUNT(DISTINCT a.doc_chembl_id)         AS n_documents,
    ARRAY_AGG(DISTINCT a.standard_type)     AS activity_types,
    MAX(a.confidence_score)                 AS best_confidence
FROM activity a
WHERE a.pchembl_value IS NOT NULL
  AND a.confidence_score >= 7
  AND (a.standard_relation IS NULL OR a.standard_relation IN ('=', '<', '<='))
  AND a.data_validity_comment IS NULL
GROUP BY a.chembl_id, a.target_chembl_id;

CREATE UNIQUE INDEX IF NOT EXISTS idx_binds_to_pk
    ON binds_to (chembl_id, target_chembl_id);
CREATE INDEX IF NOT EXISTS idx_binds_to_compound
    ON binds_to (chembl_id, max_pchembl DESC);
CREATE INDEX IF NOT EXISTS idx_binds_to_target
    ON binds_to (target_chembl_id, max_pchembl DESC);

-- ---------------------------------------------------------------------------
-- target_pathway: collapses target -> component protein -> pathway so the
-- traversal doesn't have to hop through target_component every query.
-- ---------------------------------------------------------------------------
DROP MATERIALIZED VIEW IF EXISTS target_pathway CASCADE;
CREATE MATERIALIZED VIEW target_pathway AS
SELECT DISTINCT
    tc.target_chembl_id,
    pp.reactome_id,
    p.pathway_name,
    p.biological_domain
FROM target_component tc
JOIN protein_pathway pp ON pp.uniprot_accession = tc.uniprot_accession
JOIN pathway p          ON p.reactome_id        = pp.reactome_id;

CREATE UNIQUE INDEX IF NOT EXISTS idx_target_pathway_pk
    ON target_pathway (target_chembl_id, reactome_id);
CREATE INDEX IF NOT EXISTS idx_target_pathway_reactome
    ON target_pathway (reactome_id);

COMMIT;
