-- ============================================================================
-- 003_pathway_rollup.sql
-- Resolve each pathway's top-level ancestor ("biological domain") by walking
-- the Reactome hierarchy upward.
--
-- This is the one genuinely variable-depth traversal in the relational layer.
-- The hierarchy is a DAG, not a tree -- a pathway can have multiple parents --
-- so we pick a deterministic representative domain (alphabetically first) and
-- keep the full multi-domain set alongside it.
-- ============================================================================

BEGIN;
SET search_path TO chem, public;

DROP MATERIALIZED VIEW IF EXISTS pathway_domain CASCADE;
CREATE MATERIALIZED VIEW pathway_domain AS
WITH RECURSIVE up AS (
    -- Seed: every pathway is its own ancestor at depth 0.
    SELECT
        p.reactome_id AS reactome_id,
        p.reactome_id AS ancestor_id,
        0             AS depth
    FROM pathway p

    UNION

    -- Step: climb to parents. UNION (not UNION ALL) terminates cleanly on the
    -- DAG by discarding already-seen (node, ancestor) pairs.
    SELECT
        u.reactome_id,
        h.parent_reactome_id,
        u.depth + 1
    FROM up u
    JOIN pathway_hierarchy h ON h.child_reactome_id = u.ancestor_id
    WHERE u.depth < 25          -- guard against pathological cycles
)
SELECT
    u.reactome_id,
    MIN(anc.pathway_name)                              AS biological_domain,
    ARRAY_AGG(DISTINCT anc.pathway_name ORDER BY anc.pathway_name) AS all_domains,
    MAX(u.depth)                                       AS depth_from_root
FROM up u
JOIN pathway anc ON anc.reactome_id = u.ancestor_id
WHERE anc.is_top_level
GROUP BY u.reactome_id;

CREATE UNIQUE INDEX IF NOT EXISTS idx_pathway_domain_pk
    ON pathway_domain (reactome_id);

-- Denormalise back onto `pathway` so the graph projection is a flat read.
UPDATE pathway p
SET biological_domain = d.biological_domain
FROM pathway_domain d
WHERE d.reactome_id = p.reactome_id
  AND p.biological_domain IS DISTINCT FROM d.biological_domain;

COMMIT;
