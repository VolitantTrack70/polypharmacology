-- ============================================================================
-- 005_compound_pref_name.sql
-- Carry ChEMBL's preferred compound name through to `chem.compound`.
--
-- Without this, /api/resolve can tell you a structure is CHEMBL941 but not
-- that CHEMBL941 is imatinib -- which is the part a human actually reads.
--
-- Added as its own migration rather than edited into 001: migrations that have
-- already been applied are history, and rewriting them means existing
-- databases silently diverge from fresh ones.
-- ============================================================================

BEGIN;
SET search_path TO chem, public;

ALTER TABLE compound ADD COLUMN IF NOT EXISTS pref_name TEXT;

-- Names are searched case-insensitively from the UI ("imatinib", "Imatinib").
CREATE INDEX IF NOT EXISTS idx_compound_pref_name_lower
    ON compound (LOWER(pref_name))
    WHERE pref_name IS NOT NULL;

COMMIT;
