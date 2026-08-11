"""Builds a miniature SQLite database with ChEMBL's real table structure.

The full ChEMBL SQLite dump is ~5 GB compressed and tens of GB expanded, which
makes it useless as a test dependency. This reproduces the subset of the schema
that `sources/chembl.py` actually queries, populated with real drugs and real
targets, so the parser's joins and filters can be verified in milliseconds.

Column names and types mirror ChEMBL. If a query in `sources/chembl.py` starts
failing against the real dump but passes here, this fixture has drifted --
check it against the official schema documentation first.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

DDL = """
CREATE TABLE molecule_dictionary (
    molregno       INTEGER PRIMARY KEY,
    chembl_id      TEXT NOT NULL UNIQUE,
    pref_name      TEXT,
    max_phase      NUMERIC,
    first_approval INTEGER,
    withdrawn_flag INTEGER DEFAULT 0
);

CREATE TABLE compound_structures (
    molregno            INTEGER PRIMARY KEY,
    canonical_smiles    TEXT,
    standard_inchi      TEXT,
    standard_inchi_key  TEXT
);

CREATE TABLE compound_properties (
    molregno            INTEGER PRIMARY KEY,
    full_molformula     TEXT,
    mw_freebase         NUMERIC,
    alogp               NUMERIC,
    hba                 INTEGER,
    hbd                 INTEGER,
    psa                 NUMERIC,
    rtb                 INTEGER,
    aromatic_rings      INTEGER,
    heavy_atoms         INTEGER,
    num_ro5_violations  INTEGER
);

CREATE TABLE target_dictionary (
    tid          INTEGER PRIMARY KEY,
    chembl_id    TEXT NOT NULL UNIQUE,
    pref_name    TEXT,
    target_type  TEXT,
    organism     TEXT,
    tax_id       INTEGER
);

CREATE TABLE component_sequences (
    component_id INTEGER PRIMARY KEY,
    accession    TEXT,
    sequence     TEXT,
    organism     TEXT,
    tax_id       INTEGER,
    description  TEXT,
    db_source    TEXT
);

CREATE TABLE target_components (
    tid          INTEGER NOT NULL,
    component_id INTEGER NOT NULL,
    PRIMARY KEY (tid, component_id)
);

CREATE TABLE docs (
    doc_id     INTEGER PRIMARY KEY,
    chembl_id  TEXT
);

CREATE TABLE assays (
    assay_id          INTEGER PRIMARY KEY,
    chembl_id         TEXT,
    tid               INTEGER,
    confidence_score  INTEGER,
    assay_type        TEXT
);

CREATE TABLE activities (
    activity_id            INTEGER PRIMARY KEY,
    assay_id               INTEGER NOT NULL,
    molregno               INTEGER NOT NULL,
    doc_id                 INTEGER,
    standard_type          TEXT,
    standard_relation      TEXT,
    standard_value         NUMERIC,
    standard_units         TEXT,
    pchembl_value          NUMERIC,
    data_validity_comment  TEXT,
    potential_duplicate    INTEGER DEFAULT 0
);
"""

# (molregno, chembl_id, name, smiles, max_phase, first_approval)
COMPOUNDS = [
    (1, "CHEMBL941", "IMATINIB",
     "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1", 4.0, 2001),
    (2, "CHEMBL255863", "NILOTINIB",
     "Cc1cn(-c2cc(NC(=O)c3ccc(C)c(Nc4nccc(-c5cccnc5)n4)c3)cc(C(F)(F)F)c2)cn1", 4.0, 2007),
    (3, "CHEMBL25", "ASPIRIN", "CC(=O)Oc1ccccc1C(=O)O", 4.0, 1950),
    (4, "CHEMBL521", "IBUPROFEN", "CC(C)Cc1ccc(C(C)C(=O)O)cc1", 4.0, 1974),
    (5, "CHEMBL112", "PARACETAMOL", "CC(=O)Nc1ccc(O)cc1", 4.0, 1950),
    # No structure row -- must be excluded by the parser's INNER JOIN.
    (6, "CHEMBL9999", "STRUCTURELESS", None, None, None),
]

# (molregno, formula, mw, alogp, hba, hbd, psa, rtb, arom, heavy, ro5)
PROPERTIES = [
    (1, "C29H31N7O", 493.6, 3.0, 7, 2, 86.3, 7, 4, 37, 0),
    (2, "C28H22F3N7O", 529.5, 4.6, 7, 2, 97.6, 6, 5, 38, 1),
    (3, "C9H8O4", 180.2, 1.3, 4, 1, 63.6, 3, 1, 13, 0),
    (4, "C13H18O2", 206.3, 3.5, 2, 1, 37.3, 4, 1, 15, 0),
    # molregno 5 deliberately has NO properties row -- LEFT JOIN must survive it.
]

# (tid, chembl_id, pref_name, target_type, organism, tax_id)
TARGETS = [
    (10, "CHEMBL1862", "Tyrosine-protein kinase ABL1", "SINGLE PROTEIN", "Homo sapiens", 9606),
    (11, "CHEMBL1936", "Stem cell growth factor receptor", "SINGLE PROTEIN", "Homo sapiens", 9606),
    (12, "CHEMBL2007", "PDGF receptor alpha", "SINGLE PROTEIN", "Homo sapiens", 9606),
    (13, "CHEMBL221", "Cyclooxygenase-1", "SINGLE PROTEIN", "Homo sapiens", 9606),
    (14, "CHEMBL230", "Cyclooxygenase-2", "SINGLE PROTEIN", "Homo sapiens", 9606),
    # A multi-protein complex: exercises the many-to-many target_components join.
    (15, "CHEMBL2111445", "BCR/ABL fusion protein", "PROTEIN COMPLEX", "Homo sapiens", 9606),
    (16, "CHEMBL3468", "Rat kinase (non-human)", "SINGLE PROTEIN", "Rattus norvegicus", 10116),
]

# (component_id, accession, description, organism, tax_id, db_source)
COMPONENTS = [
    (100, "P00519", "Tyrosine-protein kinase ABL1", "Homo sapiens", 9606, "SWISS-PROT"),
    (101, "P10721", "Mast/stem cell growth factor receptor Kit",
     "Homo sapiens", 9606, "SWISS-PROT"),
    (102, "P16234", "Platelet-derived growth factor receptor alpha",
     "Homo sapiens", 9606, "SWISS-PROT"),
    (103, "P23219", "Prostaglandin G/H synthase 1", "Homo sapiens", 9606, "SWISS-PROT"),
    (104, "P35354", "Prostaglandin G/H synthase 2", "Homo sapiens", 9606, "SWISS-PROT"),
    (105, "P11274", "Breakpoint cluster region protein", "Homo sapiens", 9606, "SWISS-PROT"),
    (106, "Q00000", "Unreviewed rat kinase", "Rattus norvegicus", 10116, "TREMBL"),
]

TARGET_COMPONENTS = [
    (10, 100), (11, 101), (12, 102), (13, 103), (14, 104),
    (15, 100), (15, 105),   # complex -> two proteins
    (16, 106),
]

DOCS = [(500, "CHEMBL_DOC_1"), (501, "CHEMBL_DOC_2")]

# (assay_id, chembl_id, tid, confidence_score)
ASSAYS = [
    (200, "CHEMBL_A1", 10, 9),
    (201, "CHEMBL_A2", 11, 9),
    (202, "CHEMBL_A3", 12, 8),
    (203, "CHEMBL_A4", 13, 9),
    (204, "CHEMBL_A5", 14, 9),
    (205, "CHEMBL_A6", 10, 4),   # low confidence -- excluded from binds_to
    (206, "CHEMBL_A7", 16, 9),   # rat target
]

# (activity_id, assay_id, molregno, doc_id, type, relation, value, units,
#  pchembl, validity_comment, potential_duplicate)
ACTIVITIES = [
    (1000, 200, 1, 500, "IC50", "=", 25.0, "nM", 7.6, None, 0),
    (1001, 200, 1, 501, "Kd", "=", 10.0, "nM", 8.0, None, 0),      # 2nd measurement, same pair
    (1002, 201, 1, 500, "IC50", "=", 100.0, "nM", 7.0, None, 0),
    (1003, 202, 1, 500, "IC50", "=", 200.0, "nM", 6.7, None, 0),
    (1004, 200, 2, 500, "IC50", "=", 20.0, "nM", 7.7, None, 0),
    (1005, 201, 2, 501, "IC50", "=", 60.0, "nM", 7.2, None, 0),
    (1006, 203, 3, 500, "IC50", "=", 5000.0, "nM", 5.3, None, 0),
    (1007, 204, 3, 500, "IC50", "=", 8000.0, "nM", 5.1, None, 0),
    (1008, 203, 4, 500, "IC50", "=", 3000.0, "nM", 5.5, None, 0),
    # --- rows the pipeline must handle specially ---
    (1009, 200, 1, 500, "IC50", ">", 10000.0, "nM", 5.0, None, 0),      # '>' = non-binder
    (1010, 201, 1, 500, "IC50", "=", 1.0, "nM", 9.0, "Outside typical range", 0),  # suspect
    (1011, 205, 1, 500, "IC50", "=", 50.0, "nM", 7.3, None, 0),         # low-confidence assay
    (1012, 206, 1, 500, "IC50", "=", 30.0, "nM", 7.5, None, 0),         # rat target
    (1013, 204, 4, 501, "IC50", None, None, None, None, None, 0),       # no pchembl -> excluded
]


def build_fixture(path: Path) -> Path:
    """Create the fixture database at `path`, replacing any existing file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()

    conn = sqlite3.connect(path)
    try:
        conn.executescript(DDL)

        conn.executemany(
            "INSERT INTO molecule_dictionary "
            "(molregno, chembl_id, pref_name, max_phase, first_approval) VALUES (?,?,?,?,?)",
            [(m[0], m[1], m[2], m[4], m[5]) for m in COMPOUNDS],
        )
        conn.executemany(
            "INSERT INTO compound_structures "
            "(molregno, canonical_smiles, standard_inchi_key) VALUES (?,?,?)",
            [
                (m[0], m[3], f"FAKEKEY{m[0]:022d}"[:27])
                for m in COMPOUNDS
                if m[3] is not None
            ],
        )
        conn.executemany(
            "INSERT INTO compound_properties (molregno, full_molformula, mw_freebase, alogp, "
            "hba, hbd, psa, rtb, aromatic_rings, heavy_atoms, num_ro5_violations) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?)",
            PROPERTIES,
        )
        conn.executemany(
            "INSERT INTO target_dictionary "
            "(tid, chembl_id, pref_name, target_type, organism, tax_id) VALUES (?,?,?,?,?,?)",
            TARGETS,
        )
        conn.executemany(
            "INSERT INTO component_sequences "
            "(component_id, accession, description, organism, tax_id, db_source, sequence) "
            "VALUES (?,?,?,?,?,?,?)",
            [(c[0], c[1], c[2], c[3], c[4], c[5], "MOCKSEQUENCE") for c in COMPONENTS],
        )
        conn.executemany(
            "INSERT INTO target_components (tid, component_id) VALUES (?,?)",
            TARGET_COMPONENTS,
        )
        conn.executemany("INSERT INTO docs (doc_id, chembl_id) VALUES (?,?)", DOCS)
        conn.executemany(
            "INSERT INTO assays (assay_id, chembl_id, tid, confidence_score) VALUES (?,?,?,?)",
            ASSAYS,
        )
        conn.executemany(
            "INSERT INTO activities (activity_id, assay_id, molregno, doc_id, standard_type, "
            "standard_relation, standard_value, standard_units, pchembl_value, "
            "data_validity_comment, potential_duplicate) VALUES (?,?,?,?,?,?,?,?,?,?,?)",
            ACTIVITIES,
        )
        conn.commit()
    finally:
        conn.close()

    return path


if __name__ == "__main__":
    import sys

    out = Path(sys.argv[1] if len(sys.argv) > 1 else "chembl_fixture.db")
    print(f"wrote {build_fixture(out)}")
