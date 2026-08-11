"""Verifies the ChEMBL SQL against a schema-accurate fixture.

The queries in `sources/chembl.py` encode real decisions -- which joins are
INNER vs LEFT, which rows get dropped -- and those decisions are invisible
until something silently loses data. These tests pin them down.
"""

from __future__ import annotations

import polars as pl
import pytest

from chemmed_ingest.sources import chembl

from .chembl_fixture import build_fixture


@pytest.fixture(scope="module")
def fixture_db(tmp_path_factory):
    return build_fixture(tmp_path_factory.mktemp("chembl") / "chembl_test.db")


@pytest.fixture(scope="module")
def exported(fixture_db, tmp_path_factory):
    out = tmp_path_factory.mktemp("parquet")
    counts = chembl.export_all(fixture_db, out)
    return out, counts


def read(out_dir, name: str) -> pl.DataFrame:
    return pl.read_parquet(out_dir / f"{name}.parquet")


class TestCompounds:
    def test_structureless_compound_is_dropped(self, exported):
        out, _ = exported
        df = read(out, "compound")
        # CHEMBL9999 has no compound_structures row. The INNER JOIN must drop
        # it -- a compound with no structure cannot be fingerprinted.
        assert "CHEMBL9999" not in df["chembl_id"].to_list()
        assert df.height == 5

    def test_compound_without_properties_row_survives(self, exported):
        """molregno 5 has no compound_properties row. The LEFT JOIN must keep
        the compound and null the descriptors, not drop it."""
        out, _ = exported
        df = read(out, "compound").filter(pl.col("chembl_id") == "CHEMBL112")
        assert df.height == 1
        assert df["mw_freebase"][0] is None

    def test_properties_come_through_unmodified(self, exported):
        """MW/logP are taken from ChEMBL, never recomputed."""
        out, _ = exported
        row = read(out, "compound").filter(pl.col("chembl_id") == "CHEMBL941")
        assert row["mw_freebase"][0] == pytest.approx(493.6, abs=0.1)
        assert row["alogp"][0] == pytest.approx(3.0, abs=0.1)
        assert row["num_ro5_violations"][0] == 0

    def test_withdrawn_flag_is_boolean(self, exported):
        out, _ = exported
        assert read(out, "compound")["withdrawn_flag"].dtype == pl.Boolean


class TestTargets:
    def test_all_targets_exported_including_non_human(self, exported):
        out, _ = exported
        assert read(out, "target").height == len(chembl_targets())

    def test_protein_complex_maps_to_multiple_proteins(self, exported):
        """A ChEMBL target is not 1:1 with a protein. CHEMBL2111445 (BCR/ABL)
        has two components; collapsing that would lose real information."""
        out, _ = exported
        df = read(out, "target_component").filter(
            pl.col("target_chembl_id") == "CHEMBL2111445"
        )
        assert sorted(df["uniprot_accession"].to_list()) == ["P00519", "P11274"]

    def test_one_protein_can_serve_multiple_targets(self, exported):
        """P00519 (ABL1) belongs to both the single-protein target and the
        fusion complex."""
        out, _ = exported
        df = read(out, "target_component").filter(
            pl.col("uniprot_accession") == "P00519"
        )
        assert sorted(df["target_chembl_id"].to_list()) == ["CHEMBL1862", "CHEMBL2111445"]

    def test_swissprot_flag_distinguishes_reviewed(self, exported):
        out, _ = exported
        df = read(out, "protein")
        reviewed = dict(zip(df["uniprot_accession"], df["is_reviewed"], strict=True))
        assert reviewed["P00519"] is True
        assert reviewed["Q00000"] is False


class TestActivities:
    def test_rows_without_pchembl_are_excluded(self, exported):
        """activity_id 1013 has a null pchembl_value. Without a comparable
        affinity it carries no signal and would just bloat the table."""
        out, _ = exported
        assert 1013 not in read(out, "activity")["activity_id"].to_list()

    def test_parser_does_not_apply_quality_filters(self, exported):
        """Confidence, relation, and validity filtering belong to the
        `binds_to` view, not the parser -- the raw table stays complete so
        those thresholds can be revisited without re-ingesting."""
        out, _ = exported
        ids = read(out, "activity")["activity_id"].to_list()
        assert 1009 in ids, "'>' relation row should reach the raw table"
        assert 1010 in ids, "suspect-validity row should reach the raw table"
        assert 1011 in ids, "low-confidence row should reach the raw table"

    def test_confidence_score_joined_from_assay(self, exported):
        out, _ = exported
        df = read(out, "activity").filter(pl.col("activity_id") == 1011)
        assert df["confidence_score"][0] == 4

    def test_target_and_compound_ids_are_public_not_internal(self, exported):
        """Joins happen on molregno/tid but only chembl_id is exported."""
        out, _ = exported
        df = read(out, "activity")
        assert df["chembl_id"].str.starts_with("CHEMBL").all()
        assert df["target_chembl_id"].str.starts_with("CHEMBL").all()

    def test_multiple_measurements_per_pair_are_preserved(self, exported):
        """1000 and 1001 are the same (compound, target) from different docs.
        Aggregation is the view's job; the parser must keep both."""
        out, _ = exported
        df = read(out, "activity").filter(
            (pl.col("chembl_id") == "CHEMBL941")
            & (pl.col("target_chembl_id") == "CHEMBL1862")
        )
        assert df.height >= 2


class TestLimit:
    def test_limit_caps_rows(self, fixture_db, tmp_path):
        _, n = chembl.export_table(fixture_db, "activity", tmp_path, limit=3)
        assert n == 3

    def test_schema_is_stable_under_limit(self, fixture_db, tmp_path):
        chembl.export_table(fixture_db, "compound", tmp_path, limit=1)
        df = pl.read_parquet(tmp_path / "compound.parquet")
        assert df.height == 1
        assert set(df.columns) == {f.name for f in chembl.SCHEMAS["compound"]}


class TestErrors:
    def test_missing_database_is_a_clear_error(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="download"):
            chembl.export_all(tmp_path / "nope.db", tmp_path)

    def test_unknown_entity_rejected(self, fixture_db, tmp_path):
        with pytest.raises(KeyError, match="unknown entity"):
            chembl.export_table(fixture_db, "not_a_table", tmp_path)


def chembl_targets():
    from .chembl_fixture import TARGETS

    return TARGETS
