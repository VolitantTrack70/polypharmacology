"""ChEMBL SQLite -> normalised Parquet.

ChEMBL ships as a full relational database, so this module is mostly careful
SQL rather than parsing. Everything streams in batches: `activities` alone is
~21M rows and will not fit in memory as a single frame.

Two things worth knowing about the ChEMBL model:

1. `molregno` is the internal PK; `chembl_id` is the public one. Joins happen
   on molregno, but only chembl_id is exported.
2. A ChEMBL *target* is not a protein. It may be a complex, a family, or a cell
   line. `target_components` -> `component_sequences.accession` is what maps a
   target onto UniProt proteins, and it is genuinely many-to-many. Collapsing
   it would silently merge distinct protein complexes.
"""

from __future__ import annotations

import logging
import sqlite3
from collections.abc import Iterator
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

log = logging.getLogger(__name__)

BATCH_ROWS = 250_000


# ---------------------------------------------------------------------------
# Queries. Kept as module constants so they can be inspected and tested.
# ---------------------------------------------------------------------------

# Compound columns as (output_name, sql_expression, required).
#
# ChEMBL's schema drifts between releases -- columns are added, moved between
# tables, and occasionally dropped. Anything marked optional is substituted
# with NULL when the running database does not have it, so a schema change
# degrades a single field instead of killing a multi-hour ingest.
_COMPOUND_COLUMNS: list[tuple[str, str, bool]] = [
    ("chembl_id",          "md.chembl_id",           True),
    ("pref_name",          "md.pref_name",           False),
    ("canonical_smiles",   "cs.canonical_smiles",    True),
    ("standard_inchi_key", "cs.standard_inchi_key",  True),
    ("molformula",         "cp.full_molformula",     False),
    ("mw_freebase",        "cp.mw_freebase",         False),
    ("alogp",              "cp.alogp",               False),
    ("hba",                "cp.hba",                 False),
    ("hbd",                "cp.hbd",                 False),
    ("psa",                "cp.psa",                 False),
    ("rtb",                "cp.rtb",                 False),
    ("aromatic_rings",     "cp.aromatic_rings",      False),
    ("heavy_atoms",        "cp.heavy_atoms",         False),
    ("num_ro5_violations", "cp.num_ro5_violations",  False),
    ("max_phase",          "md.max_phase",           False),
    ("first_approval",     "md.first_approval",      False),
    ("withdrawn_flag",     "COALESCE(md.withdrawn_flag, 0)", False),
]

_ALIAS_TO_TABLE = {
    "md": "molecule_dictionary",
    "cs": "compound_structures",
    "cp": "compound_properties",
}

Q_COMPOUNDS_TEMPLATE = """
SELECT
{columns}
FROM molecule_dictionary md
JOIN compound_structures  cs ON cs.molregno = md.molregno
LEFT JOIN compound_properties cp ON cp.molregno = md.molregno
WHERE cs.canonical_smiles IS NOT NULL
"""

Q_TARGETS = """
SELECT
    td.chembl_id   AS target_chembl_id,
    td.pref_name   AS pref_name,
    td.target_type AS target_type,
    td.organism    AS organism,
    td.tax_id      AS tax_id
FROM target_dictionary td
"""

Q_TARGET_COMPONENTS = """
SELECT DISTINCT
    td.chembl_id     AS target_chembl_id,
    cseq.accession   AS uniprot_accession
FROM target_dictionary td
JOIN target_components   tc   ON tc.tid = td.tid
JOIN component_sequences cseq ON cseq.component_id = tc.component_id
WHERE cseq.accession IS NOT NULL
"""

# ChEMBL already carries enough protein annotation that a separate UniProt
# download is optional. `db_source` distinguishes SWISS-PROT from TREMBL.
Q_PROTEINS = """
SELECT DISTINCT
    cseq.accession                        AS uniprot_accession,
    cseq.description                      AS protein_name,
    cseq.organism                         AS organism,
    cseq.tax_id                           AS tax_id,
    LENGTH(cseq.sequence)                 AS sequence_length,
    CASE WHEN cseq.db_source = 'SWISS-PROT' THEN 1 ELSE 0 END AS is_reviewed
FROM component_sequences cseq
WHERE cseq.accession IS NOT NULL
"""

# Only rows with a pchembl_value are exported. Rows without one carry no
# comparable affinity and would just bloat a 21M-row table.
Q_ACTIVITIES = """
SELECT
    act.activity_id            AS activity_id,
    md.chembl_id               AS chembl_id,
    td.chembl_id               AS target_chembl_id,
    ass.chembl_id              AS assay_chembl_id,
    d.chembl_id                AS doc_chembl_id,
    act.standard_type          AS standard_type,
    act.standard_relation      AS standard_relation,
    act.standard_value         AS standard_value,
    act.standard_units         AS standard_units,
    act.pchembl_value          AS pchembl_value,
    ass.confidence_score       AS confidence_score,
    act.data_validity_comment  AS data_validity_comment,
    COALESCE(act.potential_duplicate, 0) AS potential_duplicate
FROM activities act
JOIN assays              ass ON ass.assay_id = act.assay_id
JOIN target_dictionary   td  ON td.tid       = ass.tid
JOIN molecule_dictionary md  ON md.molregno  = act.molregno
LEFT JOIN docs           d   ON d.doc_id     = act.doc_id
WHERE act.pchembl_value IS NOT NULL
"""


# Explicit Arrow schemas. Without these, a batch of all-NULL values in a
# nullable column infers as null-typed and the Parquet write fails partway
# through a multi-hour job.
SCHEMAS: dict[str, pa.Schema] = {
    "compound": pa.schema([
        ("chembl_id", pa.string()),
        ("pref_name", pa.string()),
        ("canonical_smiles", pa.string()),
        ("standard_inchi_key", pa.string()),
        ("molformula", pa.string()),
        ("mw_freebase", pa.float32()),
        ("alogp", pa.float32()),
        ("hba", pa.int16()),
        ("hbd", pa.int16()),
        ("psa", pa.float32()),
        ("rtb", pa.int16()),
        ("aromatic_rings", pa.int16()),
        ("heavy_atoms", pa.int16()),
        ("num_ro5_violations", pa.int16()),
        ("max_phase", pa.float32()),
        ("first_approval", pa.int16()),
        ("withdrawn_flag", pa.bool_()),
    ]),
    "target": pa.schema([
        ("target_chembl_id", pa.string()),
        ("pref_name", pa.string()),
        ("target_type", pa.string()),
        ("organism", pa.string()),
        ("tax_id", pa.int32()),
    ]),
    "target_component": pa.schema([
        ("target_chembl_id", pa.string()),
        ("uniprot_accession", pa.string()),
    ]),
    "protein": pa.schema([
        ("uniprot_accession", pa.string()),
        ("protein_name", pa.string()),
        ("organism", pa.string()),
        ("tax_id", pa.int32()),
        ("sequence_length", pa.int32()),
        ("is_reviewed", pa.bool_()),
    ]),
    "activity": pa.schema([
        ("activity_id", pa.int64()),
        ("chembl_id", pa.string()),
        ("target_chembl_id", pa.string()),
        ("assay_chembl_id", pa.string()),
        ("doc_chembl_id", pa.string()),
        ("standard_type", pa.string()),
        ("standard_relation", pa.string()),
        ("standard_value", pa.float64()),
        ("standard_units", pa.string()),
        ("pchembl_value", pa.float32()),
        ("confidence_score", pa.int16()),
        ("data_validity_comment", pa.string()),
        ("potential_duplicate", pa.bool_()),
    ]),
}

_QUERIES: dict[str, str] = {
    "target": Q_TARGETS,
    "target_component": Q_TARGET_COMPONENTS,
    "protein": Q_PROTEINS,
    "activity": Q_ACTIVITIES,
}

# Load order matters downstream: parents before the tables that reference them.
ENTITIES: list[str] = ["compound", "target", "protein", "target_component", "activity"]

# Tables the pipeline cannot work without, and the columns it joins on. Checked
# up front so a schema mismatch surfaces immediately rather than several hours
# into a run.
REQUIRED_SCHEMA: dict[str, set[str]] = {
    "molecule_dictionary": {"molregno", "chembl_id"},
    "compound_structures": {"molregno", "canonical_smiles", "standard_inchi_key"},
    "target_dictionary": {"tid", "chembl_id"},
    "target_components": {"tid", "component_id"},
    "component_sequences": {"component_id", "accession"},
    "assays": {"assay_id", "tid", "confidence_score"},
    "activities": {"activity_id", "assay_id", "molregno", "pchembl_value"},
}


class SchemaMismatch(RuntimeError):
    """The ChEMBL database does not have the structure the queries expect."""


def table_columns(conn: sqlite3.Connection, table: str) -> set[str]:
    """Column names of `table`, or an empty set if the table is absent."""
    rows = conn.execute(f"PRAGMA table_info({table})").fetchall()
    return {r[1] for r in rows}


def check_schema(conn: sqlite3.Connection) -> dict[str, set[str]]:
    """Verify the required tables and columns exist.

    Raises SchemaMismatch listing everything missing, rather than failing on
    the first problem -- if the schema has moved on, you want the whole picture
    in one go, not a fix-and-retry loop measured in hours.

    Returns the full {table: columns} map, which the query builders use to
    decide which optional columns are available.
    """
    available = {
        table: table_columns(conn, table)
        for table in set(REQUIRED_SCHEMA) | set(_ALIAS_TO_TABLE.values()) | {"docs"}
    }

    problems: list[str] = []
    for table, needed in REQUIRED_SCHEMA.items():
        present = available.get(table, set())
        if not present:
            problems.append(f"  table {table!r} is missing entirely")
            continue
        for column in sorted(needed - present):
            problems.append(f"  {table}.{column} is missing")

    if problems:
        raise SchemaMismatch(
            "This ChEMBL database does not match the expected schema:\n"
            + "\n".join(problems)
            + "\n\nChEMBL's schema changes between releases. Check the release "
            "notes and update sources/chembl.py accordingly."
        )

    return available


def build_compound_query(available: dict[str, set[str]]) -> tuple[str, list[str]]:
    """Assemble the compound SELECT from whatever columns actually exist.

    Returns (sql, dropped_columns).
    """
    selected: list[str] = []
    dropped: list[str] = []

    for name, expr, required in _COMPOUND_COLUMNS:
        alias = expr.split(".")[0].split("(")[-1].strip()
        table = _ALIAS_TO_TABLE.get(alias)
        column = name if table is None else expr.split(".")[-1].rstrip(", 0)")

        present = table is not None and column in available.get(table, set())
        if present or table is None:
            selected.append(f"    {expr} AS {name}")
        elif required:
            raise SchemaMismatch(f"required compound column {table}.{column} is missing")
        else:
            selected.append(f"    NULL AS {name}")
            dropped.append(f"{table}.{column}")

    return Q_COMPOUNDS_TEMPLATE.format(columns=",\n".join(selected)), dropped


def _coerce(value: object, field: pa.Field) -> object:
    """SQLite is dynamically typed and hands back ints where the schema wants
    bools. Normalise the few cases that actually occur."""
    if value is None:
        return None
    if pa.types.is_boolean(field.type):
        return bool(value)
    return value


def _iter_batches(
    conn: sqlite3.Connection,
    sql: str,
    schema: pa.Schema,
    limit: int | None,
) -> Iterator[pa.RecordBatch]:
    stmt = sql if limit is None else f"{sql}\nLIMIT {int(limit)}"
    cursor = conn.execute(stmt)
    names = [f.name for f in schema]

    while True:
        rows = cursor.fetchmany(BATCH_ROWS)
        if not rows:
            break
        columns = [
            pa.array([_coerce(row[i], schema.field(i)) for row in rows], type=schema.field(i).type)
            for i in range(len(names))
        ]
        yield pa.RecordBatch.from_arrays(columns, names=names)


def export_table(
    sqlite_path: Path,
    entity: str,
    out_dir: Path,
    limit: int | None = None,
) -> tuple[Path, int]:
    """Stream one entity out of the ChEMBL SQLite into a Parquet file.

    Returns (path, row_count).
    """
    if entity not in ENTITIES:
        raise KeyError(f"unknown entity {entity!r}; known: {ENTITIES}")

    schema = SCHEMAS[entity]
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{entity}.parquet"

    conn = sqlite3.connect(f"file:{sqlite_path}?mode=ro", uri=True)
    # Read-only bulk scan: durability settings are irrelevant, cache is not.
    conn.execute("PRAGMA query_only = ON")
    conn.execute("PRAGMA cache_size = -262144")  # ~256MB
    conn.execute("PRAGMA temp_store = MEMORY")

    total = 0
    try:
        available = check_schema(conn)
        if entity == "compound":
            sql, dropped = build_compound_query(available)
            for col in dropped:
                log.warning("compound.%s not present in this release; exporting NULL", col)
        else:
            sql = _QUERIES[entity]

        writer = pq.ParquetWriter(out_path, schema, compression="zstd")
        try:
            for batch in _iter_batches(conn, sql, schema, limit):
                writer.write_batch(batch)
                total += batch.num_rows
                log.info("%s: %d rows", entity, total)
        finally:
            writer.close()
    finally:
        conn.close()

    log.info("wrote %s (%d rows)", out_path, total)
    return out_path, total


def export_all(
    sqlite_path: Path,
    out_dir: Path,
    limit: int | None = None,
) -> dict[str, int]:
    """Export every entity. `limit` applies per-table and is for smoke tests.

    Note that a `--limit`ed run produces a deliberately inconsistent slice:
    activities will reference compounds outside the truncated compound set.
    The Postgres loader filters dangling references rather than failing.
    """
    if not sqlite_path.exists():
        raise FileNotFoundError(
            f"ChEMBL SQLite not found at {sqlite_path}. "
            f"Run `chemmed-ingest download --source chembl` first."
        )

    # Preflight once, before any work, so a schema mismatch is reported in
    # seconds rather than after the first table has already been written.
    conn = sqlite3.connect(f"file:{sqlite_path}?mode=ro", uri=True)
    try:
        check_schema(conn)
    finally:
        conn.close()

    counts: dict[str, int] = {}
    for entity in ENTITIES:
        _, n = export_table(sqlite_path, entity, out_dir, limit=limit)
        counts[entity] = n
    return counts
