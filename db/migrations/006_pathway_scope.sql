-- ============================================================================
-- 006_pathway_scope.sql
-- Give `target_pathway` the hierarchy metadata needed to control how much of
-- Reactome a query pulls back.
--
-- WHY: we ingest UniProt2Reactome_All_Levels.txt, which annotates a protein to
-- every ancestor pathway, not just the leaf. That is the right thing to store
-- -- it makes roll-up possible -- but it means a single well-studied target is
-- genuinely a member of dozens of pathways. ABL1 alone returns ~52.
--
-- Returning all of them is correct and useless: the cascade becomes a hairball
-- in which "Disease" and "Signal Transduction" sit beside the one specific
-- pathway a researcher cares about. Callers need to choose a level.
-- ============================================================================

BEGIN;
SET search_path TO chem, public;

DROP MATERIALIZED VIEW IF EXISTS target_pathway CASCADE;
CREATE MATERIALIZED VIEW target_pathway AS
SELECT DISTINCT
    tc.target_chembl_id,
    pp.reactome_id,
    p.pathway_name,
    p.biological_domain,
    p.is_top_level,
    -- Distance from the top-level ancestor. 0 = a root pathway; larger is more
    -- specific. NULL when the pathway is not reachable from any root.
    COALESCE(d.depth_from_root, 0) AS depth_from_root
FROM target_component tc
JOIN protein_pathway pp ON pp.uniprot_accession = tc.uniprot_accession
JOIN pathway p          ON p.reactome_id        = pp.reactome_id
LEFT JOIN pathway_domain d ON d.reactome_id     = p.reactome_id;

CREATE UNIQUE INDEX IF NOT EXISTS idx_target_pathway_pk
    ON target_pathway (target_chembl_id, reactome_id);
CREATE INDEX IF NOT EXISTS idx_target_pathway_reactome
    ON target_pathway (reactome_id);
-- Supports the `scope` filter without re-scanning the whole view.
CREATE INDEX IF NOT EXISTS idx_target_pathway_scope
    ON target_pathway (target_chembl_id, is_top_level, depth_from_root);

COMMIT;
