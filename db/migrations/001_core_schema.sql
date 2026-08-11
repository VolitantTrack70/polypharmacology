-- ============================================================================
-- 001_core_schema.sql
-- Relational system of record for compounds, proteins, targets, bioactivity,
-- and pathways. The AGE property graph (003) is a *projection* of these
-- tables, not a separate source of truth.
--
-- Design note: ChEMBL is natively relational and we load it losslessly here.
-- Numeric-range filtering over ~21M activity rows is a relational workload;
-- the graph overlay exists for traversal ergonomics, not for storage.
-- ============================================================================

BEGIN;

CREATE SCHEMA IF NOT EXISTS chem;
SET search_path TO chem, public;

-- ---------------------------------------------------------------------------
-- Provenance. Every bulk load stamps a row here so you can tell which ChEMBL /
-- Reactome release a given fact came from, and roll back a bad ingest.
-- ---------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS data_release (
    release_id      SERIAL PRIMARY KEY,
    source          TEXT        NOT NULL,   -- 'chembl' | 'reactome' | 'uniprot'
    version         TEXT        NOT NULL,   -- e.g. '35', '2025-06'
    source_url      TEXT,
    downloaded_at   TIMESTAMPTZ,
    loaded_at       TIMESTAMPTZ NOT NULL DEFAULT now(),
    row_counts      JSONB       NOT NULL DEFAULT '{}'::jsonb,
    notes           TEXT,
    UNIQUE (source, version)
);

-- ---------------------------------------------------------------------------
-- Compounds.
-- MW / logP / HBD / HBA come straight from ChEMBL's compound_properties.
-- We do NOT recompute them with RDKit -- that would burn hours and introduce
-- silent disagreement with the published values.
-- ---------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS compound (
    chembl_id            TEXT        PRIMARY KEY,
    canonical_smiles     TEXT,
    standard_inchi_key   TEXT,
    molformula           TEXT,
    mw_freebase          REAL,
    alogp                REAL,
    hba                  SMALLINT,
    hbd                  SMALLINT,
    psa                  REAL,
    rtb                  SMALLINT,
    aromatic_rings       SMALLINT,
    heavy_atoms          SMALLINT,
    num_ro5_violations   SMALLINT,
    -- Highest clinical phase reached (0=preclinical .. 4=approved).
    -- Very useful for ranking off-target hits by real-world relevance.
    max_phase            REAL,
    first_approval       SMALLINT,
    withdrawn_flag       BOOLEAN     NOT NULL DEFAULT FALSE,
    release_id           INTEGER     REFERENCES data_release(release_id)
);

-- Morgan fingerprints, stored as raw packed bits.
-- The *search* index lives in FPSim2 (see docs/decisions/0001), not here --
-- this table is for exact lookup, reproducibility, and rebuilding that index.
CREATE TABLE IF NOT EXISTS compound_fingerprint (
    chembl_id   TEXT     PRIMARY KEY REFERENCES compound(chembl_id) ON DELETE CASCADE,
    radius      SMALLINT NOT NULL,
    n_bits      SMALLINT NOT NULL,
    -- 2048 bits = 256 bytes.
    fp          BYTEA    NOT NULL,
    -- Precomputed popcount: lets you cheaply bound Tanimoto without unpacking,
    -- since Ts(A,B) <= min(|A|,|B|) / max(|A|,|B|).
    popcount    SMALLINT NOT NULL,
    CONSTRAINT fp_config_consistent CHECK (n_bits % 8 = 0)
);

-- ---------------------------------------------------------------------------
-- Proteins (UniProt) and ChEMBL targets.
-- A ChEMBL target is NOT 1:1 with a protein -- it may be a complex, a protein
-- family, or a whole organism. Modelling the join table explicitly avoids
-- silently collapsing complexes into single proteins.
-- ---------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS protein (
    uniprot_accession  TEXT PRIMARY KEY,
    gene_symbol        TEXT,
    protein_name       TEXT,
    organism           TEXT,
    tax_id             INTEGER,
    sequence_length    INTEGER,
    is_reviewed        BOOLEAN,          -- SwissProt (true) vs TrEMBL (false)
    release_id         INTEGER REFERENCES data_release(release_id)
);

CREATE TABLE IF NOT EXISTS target (
    target_chembl_id   TEXT PRIMARY KEY,
    pref_name          TEXT,
    target_type        TEXT,             -- 'SINGLE PROTEIN', 'PROTEIN COMPLEX', ...
    organism           TEXT,
    tax_id             INTEGER,
    release_id         INTEGER REFERENCES data_release(release_id)
);

CREATE TABLE IF NOT EXISTS target_component (
    target_chembl_id   TEXT NOT NULL REFERENCES target(target_chembl_id) ON DELETE CASCADE,
    uniprot_accession  TEXT NOT NULL,
    PRIMARY KEY (target_chembl_id, uniprot_accession)
);

-- ---------------------------------------------------------------------------
-- Bioactivity. One row per measurement, kept at full granularity.
--
-- pchembl_value is ChEMBL's -log10 normalisation across IC50/Ki/Kd/EC50.
-- It is a *rough* affinity proxy: IC50 is assay-condition dependent, so
-- standard_type is retained on every row and surfaced in the API to let
-- callers filter to a single assay type when rigour matters.
-- ---------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS activity (
    activity_id        BIGINT PRIMARY KEY,
    chembl_id          TEXT   NOT NULL REFERENCES compound(chembl_id) ON DELETE CASCADE,
    target_chembl_id   TEXT   NOT NULL REFERENCES target(target_chembl_id) ON DELETE CASCADE,
    assay_chembl_id    TEXT,
    doc_chembl_id      TEXT,
    standard_type      TEXT,             -- 'IC50' | 'Ki' | 'Kd' | 'EC50' | ...
    standard_relation  TEXT,             -- '=' | '>' | '<' ...  ('>' means "no better than")
    standard_value     DOUBLE PRECISION,
    standard_units     TEXT,
    pchembl_value      REAL,
    -- ChEMBL assay-to-target confidence, 0-9. 9 = direct single-protein assay.
    -- Anything below ~7 is a weak basis for an off-target claim.
    confidence_score   SMALLINT,
    data_validity_comment TEXT,          -- non-null often means "suspect"
    potential_duplicate   BOOLEAN,
    release_id         INTEGER REFERENCES data_release(release_id)
);

-- ---------------------------------------------------------------------------
-- Reactome pathways.
-- Hierarchy is genuinely recursive (arbitrary depth) but small (~2.7k human
-- pathways), so a recursive CTE over this table is cheap.
-- ---------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS pathway (
    reactome_id        TEXT PRIMARY KEY,
    pathway_name       TEXT NOT NULL,
    species            TEXT,
    is_top_level       BOOLEAN NOT NULL DEFAULT FALSE,
    -- Denormalised top-level ancestor name; filled by 004 rollup.
    biological_domain  TEXT,
    release_id         INTEGER REFERENCES data_release(release_id)
);

CREATE TABLE IF NOT EXISTS pathway_hierarchy (
    parent_reactome_id TEXT NOT NULL REFERENCES pathway(reactome_id) ON DELETE CASCADE,
    child_reactome_id  TEXT NOT NULL REFERENCES pathway(reactome_id) ON DELETE CASCADE,
    PRIMARY KEY (parent_reactome_id, child_reactome_id)
);

CREATE TABLE IF NOT EXISTS protein_pathway (
    uniprot_accession  TEXT NOT NULL,
    reactome_id        TEXT NOT NULL REFERENCES pathway(reactome_id) ON DELETE CASCADE,
    evidence_code      TEXT,             -- 'TAS' (traceable) | 'IEA' (inferred)
    PRIMARY KEY (uniprot_accession, reactome_id)
);

COMMIT;
