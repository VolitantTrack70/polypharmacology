-- requires-extension: age
-- ============================================================================
-- 004_age_graph.sql
-- Apache AGE property-graph overlay.
--
-- SKIPPED AUTOMATICALLY when the `age` extension isn't installed (e.g. a stock
-- Postgres rather than the apache/age image). Nothing downstream depends on
-- it: the relational tables are the system of record and answer every query
-- the API currently makes.
--
-- IMPORTANT -- what this graph is and isn't:
--   * It is a PROJECTION of the relational tables. `chem.*` is the source of
--     truth; drop and rebuild this graph freely.
--   * It contains NO `IS_STRUCTURALLY_SIMILAR_TO` edges. Chemical similarity
--     is computed at query time against the FPSim2 index -- materialising
--     all-pairs Tanimoto over 2.4M compounds is ~3e12 comparisons and is not
--     something you want on disk. See docs/decisions/0001.
--   * Only compounds that actually have a qualifying bioactivity edge are
--     projected. Isolated compound nodes add millions of vertices and answer
--     no question the relational layer can't answer faster.
--
-- Population is done by the ingest loader (`chemmed-ingest project-graph`),
-- which writes CSVs and uses AGE's bulk file loaders. Creating 3M edges via
-- individual cypher() CREATE statements is prohibitively slow.
-- ============================================================================

CREATE EXTENSION IF NOT EXISTS age;
LOAD 'age';
SET search_path = ag_catalog, "$user", chem, public;

-- Idempotent graph creation.
DO $$
BEGIN
    IF NOT EXISTS (SELECT 1 FROM ag_catalog.ag_graph WHERE name = 'polypharm') THEN
        PERFORM ag_catalog.create_graph('polypharm');
    END IF;
END
$$;

-- Vertex labels
DO $$
DECLARE
    lbl TEXT;
BEGIN
    FOREACH lbl IN ARRAY ARRAY['Compound', 'Target', 'Pathway'] LOOP
        IF NOT EXISTS (
            SELECT 1 FROM ag_catalog.ag_label l
            JOIN ag_catalog.ag_graph g ON g.oid = l.graph
            WHERE g.name = 'polypharm' AND l.name = lbl
        ) THEN
            PERFORM ag_catalog.create_vlabel('polypharm', lbl);
        END IF;
    END LOOP;
END
$$;

-- Edge labels
DO $$
DECLARE
    lbl TEXT;
BEGIN
    FOREACH lbl IN ARRAY ARRAY['BINDS_TO', 'PARTICIPATES_IN', 'PARENT_OF'] LOOP
        IF NOT EXISTS (
            SELECT 1 FROM ag_catalog.ag_label l
            JOIN ag_catalog.ag_graph g ON g.oid = l.graph
            WHERE g.name = 'polypharm' AND l.name = lbl
        ) THEN
            PERFORM ag_catalog.create_elabel('polypharm', lbl);
        END IF;
    END LOOP;
END
$$;

-- Property indexes on the AGE-backing tables. AGE stores properties as
-- agtype JSONB, so these are expression indexes on the extracted key.
CREATE INDEX IF NOT EXISTS idx_age_compound_chembl
    ON polypharm."Compound" (ag_catalog.agtype_access_operator(properties, '"chembl_id"'::agtype));

CREATE INDEX IF NOT EXISTS idx_age_target_chembl
    ON polypharm."Target" (ag_catalog.agtype_access_operator(properties, '"target_chembl_id"'::agtype));

CREATE INDEX IF NOT EXISTS idx_age_pathway_reactome
    ON polypharm."Pathway" (ag_catalog.agtype_access_operator(properties, '"reactome_id"'::agtype));
