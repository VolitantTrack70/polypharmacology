"""Bulk load Parquet into Postgres via COPY.

Everything routes through an UNLOGGED staging table before landing in the real
one. That costs an extra pass but buys three things worth more than the pass:

  * Dangling foreign keys are filtered instead of aborting the load. A
    `--limit`ed smoke run produces activities referencing compounds outside the
    truncated set; failing on that would make small runs impossible.
  * Re-running a stage is idempotent (`ON CONFLICT DO NOTHING`).
  * A malformed row kills the staging load, not a half-populated real table.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import psycopg
import pyarrow.parquet as pq

log = logging.getLogger(__name__)

COPY_BATCH = 50_000

# table -> (conflict target, referential guards applied when moving staging->real)
#   guard = (column_in_staging, referenced_table, referenced_column)
_LOAD_SPEC: dict[str, dict[str, Any]] = {
    "compound": {"pk": ["chembl_id"], "guards": []},
    "target": {"pk": ["target_chembl_id"], "guards": []},
    "protein": {"pk": ["uniprot_accession"], "guards": []},
    "target_component": {
        "pk": ["target_chembl_id", "uniprot_accession"],
        "guards": [("target_chembl_id", "target", "target_chembl_id")],
    },
    "activity": {
        "pk": ["activity_id"],
        "guards": [
            ("chembl_id", "compound", "chembl_id"),
            ("target_chembl_id", "target", "target_chembl_id"),
        ],
    },
    "pathway": {"pk": ["reactome_id"], "guards": []},
    "pathway_hierarchy": {
        "pk": ["parent_reactome_id", "child_reactome_id"],
        "guards": [
            ("parent_reactome_id", "pathway", "reactome_id"),
            ("child_reactome_id", "pathway", "reactome_id"),
        ],
    },
    "protein_pathway": {
        "pk": ["uniprot_accession", "reactome_id"],
        "guards": [("reactome_id", "pathway", "reactome_id")],
    },
    "compound_fingerprint": {
        "pk": ["chembl_id"],
        "guards": [("chembl_id", "compound", "chembl_id")],
    },
}


def _iter_rows(path: Path) -> Iterator[tuple[list[str], list[tuple]]]:
    """Stream a Parquet file as (column_names, row_batch)."""
    parquet = pq.ParquetFile(path)
    columns = [f.name for f in parquet.schema_arrow]
    for batch in parquet.iter_batches(batch_size=COPY_BATCH):
        table = batch.to_pydict()
        rows = list(zip(*(table[c] for c in columns), strict=True))
        yield columns, rows


def load_table(
    conn: psycopg.Connection,
    table: str,
    parquet_path: Path,
    release_id: int | None = None,
) -> int:
    """COPY one Parquet file into `chem.<table>`. Returns rows actually inserted."""
    if table not in _LOAD_SPEC:
        raise KeyError(f"no load spec for {table!r}")
    if not parquet_path.exists():
        log.warning("skipping %s: %s not found", table, parquet_path)
        return 0

    spec = _LOAD_SPEC[table]
    staging = f"_stage_{table}"

    with conn.cursor() as cur:
        cur.execute(f"DROP TABLE IF EXISTS {staging}")
        # AS SELECT, not LIKE: LIKE copies NOT NULL but not defaults, so any
        # Parquet omitting a NOT NULL DEFAULT column would fail the COPY.
        # Staging takes anything; the real table's defaults apply on insert.
        cur.execute(
            f"CREATE UNLOGGED TABLE {staging} AS "
            f"SELECT * FROM chem.{table} WITH NO DATA"
        )

        copied = 0
        first = True
        columns: list[str] = []
        for cols, rows in _iter_rows(parquet_path):
            if first:
                columns = cols
                first = False
            collist = ", ".join(f'"{c}"' for c in columns)
            with cur.copy(f"COPY {staging} ({collist}) FROM STDIN") as copy:
                for row in rows:
                    copy.write_row(row)
            copied += len(rows)
            log.info("%s: staged %d rows", table, copied)

        if copied == 0:
            cur.execute(f"DROP TABLE IF EXISTS {staging}")
            return 0

        # Move staging -> real, dropping rows that would violate an FK.
        where = " AND ".join(
            f"EXISTS (SELECT 1 FROM chem.{ref_t} r WHERE r.{ref_c} = s.{col})"
            for col, ref_t, ref_c in spec["guards"]
        )
        where_clause = f"WHERE {where}" if where else ""

        collist = ", ".join(f'"{c}"' for c in columns)
        select_cols = ", ".join(f"s.{c}" for c in columns)
        conflict = ", ".join(spec["pk"])

        # Only stamp provenance on tables that actually carry it -- the pure
        # join tables (target_component, pathway_hierarchy, protein_pathway,
        # compound_fingerprint) have no release_id column.
        stamped_release = False
        if release_id is not None and "release_id" not in columns:
            cur.execute(
                """
                SELECT 1 FROM information_schema.columns
                WHERE table_schema = 'chem' AND table_name = %s
                  AND column_name = 'release_id'
                """,
                (table,),
            )
            if cur.fetchone() is not None:
                collist += ", release_id"
                select_cols += f", {int(release_id)}"
                stamped_release = True

        # Upsert, so re-ingesting a newer release refreshes existing rows.
        # Join tables have no non-key columns and fall back to DO NOTHING
        # (an empty SET clause is a syntax error).
        updatable = [c for c in columns if c not in spec["pk"]]
        # release_id isn't a Parquet column, so add it explicitly or a row
        # keeps whichever release first loaded it.
        if stamped_release:
            updatable.append("release_id")
        if updatable:
            assignments = ", ".join(f'"{c}" = EXCLUDED."{c}"' for c in updatable)
            action = f"DO UPDATE SET {assignments}"
        else:
            action = "DO NOTHING"

        cur.execute(
            f"INSERT INTO chem.{table} ({collist}) "
            f"SELECT {select_cols} FROM {staging} s {where_clause} "
            f"ON CONFLICT ({conflict}) {action}"
        )
        inserted = cur.rowcount
        cur.execute(f"DROP TABLE IF EXISTS {staging}")

    conn.commit()
    dropped = copied - inserted
    if dropped:
        log.info("%s: inserted %d, skipped %d (dangling FK or duplicate)",
                 table, inserted, dropped)
    else:
        log.info("%s: inserted %d", table, inserted)
    return inserted


# Every table holding ingested source data, children first.
LOADABLE_TABLES = [
    "activity",
    "compound_fingerprint",
    "target_component",
    "protein_pathway",
    "pathway_hierarchy",
    "compound",
    "target",
    "protein",
    "pathway",
]


@contextmanager
def deferred_indexes(conn: psycopg.Connection, tables: list[str]) -> Iterator[int]:
    """Drop secondary indexes for the duration of a bulk load, then rebuild.

    Maintaining an index row-by-row during a multi-million-row INSERT is much
    slower than building it once at the end.

    Unique and constraint-backed indexes are left alone: ON CONFLICT needs the
    unique index on its conflict target, and dropping a PK would cascade.
    Definitions come from pg_get_indexdef so there is no DDL to keep in sync.
    """
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT c.relname, pg_get_indexdef(c.oid)
            FROM pg_index x
            JOIN pg_class c      ON c.oid = x.indexrelid
            JOIN pg_class t      ON t.oid = x.indrelid
            JOIN pg_namespace n  ON n.oid = t.relnamespace
            LEFT JOIN pg_constraint con ON con.conindid = c.oid
            WHERE n.nspname = 'chem'
              AND t.relname = ANY(%s)
              AND con.oid IS NULL
              AND NOT x.indisunique
            """,
            (tables,),
        )
        saved = cur.fetchall()

        for name, _ in saved:
            log.info("dropping index %s", name)
            cur.execute(f"DROP INDEX IF EXISTS chem.{name}")
    conn.commit()

    try:
        yield len(saved)
    finally:
        with conn.cursor() as cur:
            for name, definition in saved:
                log.info("rebuilding index %s", name)
                cur.execute(definition)
        conn.commit()


def truncate_all(conn: psycopg.Connection) -> None:
    """Empty every ingested table. Used before a load that must not merge with
    what is already there -- e.g. replacing fixture data with a real release,
    where reused identifiers would otherwise silently mean different things."""
    with conn.cursor() as cur:
        tables = ", ".join(f"chem.{t}" for t in LOADABLE_TABLES)
        log.info("truncating %s", tables)
        cur.execute(f"TRUNCATE {tables} RESTART IDENTITY CASCADE")
    conn.commit()


def record_release(
    conn: psycopg.Connection,
    source: str,
    version: str,
    source_url: str | None = None,
    row_counts: dict[str, int] | None = None,
) -> int:
    """Upsert a provenance row and return its release_id."""
    import json

    with conn.cursor() as cur:
        cur.execute(
            """
            INSERT INTO chem.data_release (source, version, source_url, row_counts)
            VALUES (%s, %s, %s, %s::jsonb)
            ON CONFLICT (source, version) DO UPDATE
              SET loaded_at = now(),
                  row_counts = EXCLUDED.row_counts,
                  source_url = COALESCE(EXCLUDED.source_url, chem.data_release.source_url)
            RETURNING release_id
            """,
            (source, version, source_url, json.dumps(row_counts or {})),
        )
        row = cur.fetchone()
    conn.commit()
    assert row is not None
    return int(row[0])


def refresh_derived(conn: psycopg.Connection) -> None:
    """Rebuild the materialised views after a load.

    ORDER MATTERS and is not alphabetical:

      1. `pathway_domain` walks the Reactome hierarchy upward to find each
         pathway's top-level ancestor.
      2. That result is denormalised onto `pathway.biological_domain`. This
         UPDATE also lives in migration 003, but a migration runs once against
         a table that is usually still empty -- so it has to happen here too,
         after every load, or the column stays NULL forever.
      3. `target_pathway` reads `pathway.biological_domain`, so it must refresh
         only after step 2.

    CONCURRENTLY is deliberately not used: it requires the view to be already
    populated and is far slower on a first build.
    """
    with conn.cursor() as cur:
        log.info("refreshing chem.binds_to")
        cur.execute("REFRESH MATERIALIZED VIEW chem.binds_to")

        log.info("refreshing chem.pathway_domain")
        cur.execute("REFRESH MATERIALIZED VIEW chem.pathway_domain")

        log.info("denormalising biological_domain onto chem.pathway")
        cur.execute(
            """
            UPDATE chem.pathway p
            SET biological_domain = d.biological_domain
            FROM chem.pathway_domain d
            WHERE d.reactome_id = p.reactome_id
              AND p.biological_domain IS DISTINCT FROM d.biological_domain
            """
        )
        log.info("  %d pathways updated", cur.rowcount)

        log.info("refreshing chem.target_pathway")
        cur.execute("REFRESH MATERIALIZED VIEW chem.target_pathway")
    conn.commit()


REQUIRES_MARKER = "-- requires-extension:"


def _required_extension(sql: str) -> str | None:
    """Read an optional `-- requires-extension: <name>` marker from a migration."""
    for line in sql.splitlines()[:20]:
        if line.strip().startswith(REQUIRES_MARKER):
            return line.split(":", 1)[1].strip()
    return None


def _extension_available(conn: psycopg.Connection, name: str) -> bool:
    with conn.cursor() as cur:
        cur.execute("SELECT 1 FROM pg_available_extensions WHERE name = %s", (name,))
        return cur.fetchone() is not None


def apply_migrations(
    conn: psycopg.Connection, migrations_dir: Path
) -> tuple[list[str], list[tuple[str, str]]]:
    """Run every .sql in `db/migrations` in filename order.

    Migrations are idempotent, so re-running is safe.

    A migration may declare `-- requires-extension: <name>`. If that extension
    isn't available on the server it is SKIPPED rather than failing the run.
    This is what lets the pipeline work against a stock Postgres without Apache
    AGE -- the relational layer is the system of record and answers every
    current query, so the graph overlay is genuinely optional.

    Returns (applied, skipped) where skipped is [(filename, reason)].
    """
    applied: list[str] = []
    skipped: list[tuple[str, str]] = []

    for path in sorted(migrations_dir.glob("*.sql")):
        sql = path.read_text(encoding="utf-8")

        required = _required_extension(sql)
        if required and not _extension_available(conn, required):
            reason = f"extension {required!r} not available on this server"
            log.warning("skipping %s: %s", path.name, reason)
            skipped.append((path.name, reason))
            continue

        log.info("applying %s", path.name)
        try:
            with conn.cursor() as cur:
                cur.execute(sql)
            conn.commit()
        except psycopg.Error:
            # Leave the connection usable so later migrations can still run,
            # then re-raise -- a failed migration is not something to swallow.
            conn.rollback()
            log.error("migration %s failed", path.name)
            raise
        applied.append(path.name)

    return applied, skipped
