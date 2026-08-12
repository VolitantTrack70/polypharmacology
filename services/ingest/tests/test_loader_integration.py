"""Integration tests for the Postgres loader.

These run against a real database because the behaviour under test IS the
database behaviour -- upsert semantics, foreign-key guards, and which tables
carry provenance. Mocking psycopg would test the mock.

A scratch database is created and dropped per session, so the dev database is
never touched. Skipped entirely if Postgres isn't reachable.
"""

from __future__ import annotations

import os
from pathlib import Path
from urllib.parse import urlparse, urlunparse

import polars as pl
import pytest

psycopg = pytest.importorskip("psycopg")

from chemmed_ingest.config import DATABASE_URL, REPO_ROOT  # noqa: E402
from chemmed_ingest.load.postgres import (  # noqa: E402
    apply_migrations,
    deferred_indexes,
    load_table,
    record_release,
)

TEST_DB = "chemmed_pytest"


def _url_for(database: str) -> str:
    parts = urlparse(DATABASE_URL)
    return urlunparse(parts._replace(path=f"/{database}"))


def _server_reachable() -> bool:
    try:
        with psycopg.connect(_url_for("postgres"), connect_timeout=3):
            return True
    except Exception:
        return False


pytestmark = pytest.mark.skipif(
    not _server_reachable(),
    reason="PostgreSQL not reachable; set DATABASE_URL to run loader integration tests",
)


@pytest.fixture(scope="module")
def db():
    """Create a scratch database, migrate it, yield a connection, then drop it."""
    admin = psycopg.connect(_url_for("postgres"), autocommit=True)
    try:
        admin.execute(f"DROP DATABASE IF EXISTS {TEST_DB} WITH (FORCE)")
        admin.execute(f"CREATE DATABASE {TEST_DB}")
    finally:
        admin.close()

    conn = psycopg.connect(_url_for(TEST_DB))
    try:
        applied, _ = apply_migrations(conn, REPO_ROOT / "db" / "migrations")
        assert "001_core_schema.sql" in applied
        yield conn
    finally:
        conn.close()
        admin = psycopg.connect(_url_for("postgres"), autocommit=True)
        try:
            admin.execute(f"DROP DATABASE IF EXISTS {TEST_DB} WITH (FORCE)")
        finally:
            admin.close()


def write_parquet(tmp_path: Path, name: str, rows: list[dict]) -> Path:
    path = tmp_path / f"{name}.parquet"
    pl.DataFrame(rows).write_parquet(path)
    return path


def count(conn, table: str) -> int:
    with conn.cursor() as cur:
        cur.execute(f"SELECT count(*) FROM chem.{table}")
        return cur.fetchone()[0]


class TestForeignKeyGuards:
    def test_dangling_activity_rows_are_dropped_not_fatal(self, db, tmp_path):
        """A --limit'ed run yields activities referencing compounds outside the
        truncated set. That must filter, not abort -- otherwise small runs are
        impossible."""
        load_table(db, "compound", write_parquet(tmp_path, "compound", [
            {"chembl_id": "C_REAL", "canonical_smiles": "CCO"},
        ]))
        load_table(db, "target", write_parquet(tmp_path, "target", [
            {"target_chembl_id": "T_REAL", "pref_name": "Real target"},
        ]))

        path = write_parquet(tmp_path, "activity", [
            {"activity_id": 1, "chembl_id": "C_REAL", "target_chembl_id": "T_REAL",
             "pchembl_value": 7.0},
            # references a compound that was never loaded
            {"activity_id": 2, "chembl_id": "C_MISSING", "target_chembl_id": "T_REAL",
             "pchembl_value": 8.0},
            # references a target that was never loaded
            {"activity_id": 3, "chembl_id": "C_REAL", "target_chembl_id": "T_MISSING",
             "pchembl_value": 9.0},
        ])
        inserted = load_table(db, "activity", path)

        assert inserted == 1, "only the fully-resolvable row should land"
        with db.cursor() as cur:
            cur.execute("SELECT activity_id FROM chem.activity ORDER BY activity_id")
            assert [r[0] for r in cur.fetchall()] == [1]


class TestUpsert:
    def test_reload_updates_existing_rows(self, db, tmp_path):
        """Re-ingesting a newer release must refresh existing rows. With
        ON CONFLICT DO NOTHING a corrected structure would be ignored forever."""
        first = write_parquet(tmp_path, "compound", [
            {"chembl_id": "C_UPSERT", "canonical_smiles": "CCO", "max_phase": 1.0},
        ])
        load_table(db, "compound", first)

        # Destination table is a separate argument from the file path, so the
        # "newer release" can just be a differently-named file.
        second = write_parquet(tmp_path, "compound_v2", [
            {"chembl_id": "C_UPSERT", "canonical_smiles": "CCOCC", "max_phase": 4.0},
        ])
        load_table(db, "compound", second)

        with db.cursor() as cur:
            cur.execute(
                "SELECT canonical_smiles, max_phase FROM chem.compound "
                "WHERE chembl_id = 'C_UPSERT'"
            )
            smiles, phase = cur.fetchone()
        assert smiles == "CCOCC"
        assert phase == pytest.approx(4.0)

    def test_reload_does_not_duplicate(self, db, tmp_path):
        before = count(db, "compound")
        path = write_parquet(tmp_path, "compound", [
            {"chembl_id": "C_DUP", "canonical_smiles": "CCC"},
        ])
        load_table(db, "compound", path)
        after_first = count(db, "compound")
        load_table(db, "compound", path)
        after_second = count(db, "compound")

        assert after_first == before + 1
        assert after_second == after_first

    def test_join_table_reload_is_safe(self, db, tmp_path):
        """Pure join tables have no non-key columns, so the upsert must fall
        back to DO NOTHING -- an empty SET clause is a syntax error."""
        load_table(db, "target", write_parquet(tmp_path, "target", [
            {"target_chembl_id": "T_JOIN", "pref_name": "T"},
        ]))
        path = write_parquet(tmp_path, "target_component", [
            {"target_chembl_id": "T_JOIN", "uniprot_accession": "P00001"},
        ])
        load_table(db, "target_component", path)
        load_table(db, "target_component", path)  # must not raise

        with db.cursor() as cur:
            cur.execute(
                "SELECT count(*) FROM chem.target_component WHERE target_chembl_id='T_JOIN'"
            )
            assert cur.fetchone()[0] == 1


class TestProvenance:
    def test_release_id_stamped_only_where_column_exists(self, db, tmp_path):
        """target_component has no release_id column. Blindly appending it was
        a real bug that aborted the load."""
        rid = record_release(db, "chembl", "test-provenance", row_counts={"n": 1})

        load_table(db, "compound", write_parquet(tmp_path, "compound", [
            {"chembl_id": "C_PROV", "canonical_smiles": "CO"},
        ]), release_id=rid)
        load_table(db, "target", write_parquet(tmp_path, "target", [
            {"target_chembl_id": "T_PROV", "pref_name": "T"},
        ]), release_id=rid)
        # must not raise despite having no release_id column
        load_table(db, "target_component", write_parquet(tmp_path, "target_component", [
            {"target_chembl_id": "T_PROV", "uniprot_accession": "P00002"},
        ]), release_id=rid)

        with db.cursor() as cur:
            cur.execute("SELECT release_id FROM chem.compound WHERE chembl_id='C_PROV'")
            assert cur.fetchone()[0] == rid

    def test_record_release_is_idempotent(self, db):
        a = record_release(db, "reactome", "v-same")
        b = record_release(db, "reactome", "v-same")
        assert a == b

    def test_release_id_is_refreshed_on_reload(self, db, tmp_path):
        """release_id is not a Parquet column, so it has to be added to the
        upsert's SET clause explicitly. Without that a row keeps whichever
        release first loaded it forever -- and stale provenance is worse than
        none, because it still looks authoritative."""
        first = record_release(db, "chembl", "rel-1")
        second = record_release(db, "chembl", "rel-2")
        assert first != second

        path = write_parquet(tmp_path, "compound", [
            {"chembl_id": "C_RELOAD", "canonical_smiles": "CCN"},
        ])
        load_table(db, "compound", path, release_id=first)
        load_table(db, "compound", path, release_id=second)

        with db.cursor() as cur:
            cur.execute("SELECT release_id FROM chem.compound WHERE chembl_id='C_RELOAD'")
            assert cur.fetchone()[0] == second

    def test_row_counts_are_recorded(self, db):
        rid = record_release(db, "chembl", "with-counts", row_counts={"compound": 42})
        with db.cursor() as cur:
            cur.execute("SELECT row_counts FROM chem.data_release WHERE release_id = %s", (rid,))
            assert cur.fetchone()[0] == {"compound": 42}


class TestMigrations:
    def test_age_migration_skipped_without_extension(self, db):
        """004 declares `-- requires-extension: age`. On a stock Postgres it
        must be skipped with a reason, not blow up the run."""
        _, skipped = apply_migrations(db, REPO_ROOT / "db" / "migrations")
        names = [n for n, _ in skipped]
        if not _extension_present(db, "age"):
            assert "004_age_graph.sql" in names

    def test_migrations_are_idempotent(self, db):
        applied_a, _ = apply_migrations(db, REPO_ROOT / "db" / "migrations")
        applied_b, _ = apply_migrations(db, REPO_ROOT / "db" / "migrations")
        assert applied_a == applied_b


def _extension_present(conn, name: str) -> bool:
    with conn.cursor() as cur:
        cur.execute("SELECT 1 FROM pg_available_extensions WHERE name = %s", (name,))
        return cur.fetchone() is not None


class TestDeferredIndexes:
    def _indexes(self, db, table: str) -> dict[str, bool]:
        """{index_name: is_unique} for one table."""
        with db.cursor() as cur:
            cur.execute(
                """
                SELECT c.relname, x.indisunique
                FROM pg_index x
                JOIN pg_class c     ON c.oid = x.indexrelid
                JOIN pg_class t     ON t.oid = x.indrelid
                JOIN pg_namespace n ON n.oid = t.relnamespace
                WHERE n.nspname = 'chem' AND t.relname = %s
                """,
                (table,),
            )
            return dict(cur.fetchall())

    def test_secondary_indexes_are_dropped_and_rebuilt(self, db):
        before = self._indexes(db, "activity")
        secondary = {n for n, uniq in before.items() if not uniq}
        assert secondary, "fixture precondition: activity should have secondary indexes"

        with deferred_indexes(db, ["activity"]) as n:
            assert n == len(secondary)
            during = self._indexes(db, "activity")
            assert not (secondary & set(during)), "secondary indexes should be gone"

        assert self._indexes(db, "activity") == before

    def test_unique_indexes_are_never_dropped(self, db):
        """ON CONFLICT needs the unique index on its conflict target."""
        before = self._indexes(db, "compound")
        unique = {n for n, uniq in before.items() if uniq}

        with deferred_indexes(db, ["compound"]):
            during = set(self._indexes(db, "compound"))
            assert unique <= during

    def test_indexes_are_rebuilt_even_if_the_load_fails(self, db):
        before = self._indexes(db, "activity")
        with pytest.raises(RuntimeError), deferred_indexes(db, ["activity"]):
            raise RuntimeError("simulated load failure")
        assert self._indexes(db, "activity") == before

    def test_load_still_works_with_indexes_deferred(self, db, tmp_path):
        load_table(db, "compound", write_parquet(tmp_path, "compound", [
            {"chembl_id": "C_DEFER", "canonical_smiles": "CCCC"},
        ]))
        load_table(db, "target", write_parquet(tmp_path, "target", [
            {"target_chembl_id": "T_DEFER", "pref_name": "T"},
        ]))
        path = write_parquet(tmp_path, "activity", [
            {"activity_id": 9001, "chembl_id": "C_DEFER",
             "target_chembl_id": "T_DEFER", "pchembl_value": 7.5},
        ])
        with deferred_indexes(db, ["activity"]):
            assert load_table(db, "activity", path) == 1

        with db.cursor() as cur:
            cur.execute("SELECT count(*) FROM chem.activity WHERE activity_id = 9001")
            assert cur.fetchone()[0] == 1


class TestMissingFiles:
    def test_absent_parquet_is_skipped_quietly(self, db, tmp_path):
        assert load_table(db, "pathway", tmp_path / "nope.parquet") == 0


def test_env_is_not_the_dev_database():
    """Guard: these tests drop their database. Make sure that is never the one
    the developer is actually using."""
    assert TEST_DB not in DATABASE_URL, (
        f"DATABASE_URL points at {TEST_DB!r}, which this module drops."
    )
    assert os.environ.get("DATABASE_URL", DATABASE_URL) != _url_for(TEST_DB)
