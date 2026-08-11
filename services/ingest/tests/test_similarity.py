"""Tests for the Tanimoto search engine.

These deliberately avoid RDKit so the core maths can be verified without the
heavy dependency installed.
"""

from __future__ import annotations

import numpy as np
import pytest

from chemmed_ingest.chem.similarity import FingerprintIndex, tanimoto
from chemmed_ingest.config import FP


def make_fp(set_bits: set[int]) -> bytes:
    """Build a packed fingerprint with the given bit positions set."""
    bits = np.zeros(FP.n_bits, dtype=np.uint8)
    for b in set_bits:
        bits[b] = 1
    return np.packbits(bits).tobytes()


class TestTanimoto:
    def test_identical_vectors_score_one(self):
        fp = make_fp({1, 5, 9, 100, 2000})
        assert tanimoto(fp, fp) == pytest.approx(1.0)

    def test_disjoint_vectors_score_zero(self):
        assert tanimoto(make_fp({1, 2, 3}), make_fp({4, 5, 6})) == pytest.approx(0.0)

    def test_known_value(self):
        # |A|=4, |B|=4, |A and B|=2  ->  2 / (4 + 4 - 2) = 1/3
        a = make_fp({1, 2, 3, 4})
        b = make_fp({3, 4, 5, 6})
        assert tanimoto(a, b) == pytest.approx(1 / 3)

    def test_both_empty_is_defined_as_one(self):
        empty = make_fp(set())
        assert tanimoto(empty, empty) == pytest.approx(1.0)

    def test_symmetry(self):
        a, b = make_fp({1, 7, 19, 512}), make_fp({7, 19, 900})
        assert tanimoto(a, b) == pytest.approx(tanimoto(b, a))


class TestFingerprintIndex:
    @pytest.fixture
    def index(self) -> FingerprintIndex:
        records = [
            ("CHEMBL_EXACT", make_fp({1, 2, 3, 4})),
            ("CHEMBL_HALF", make_fp({3, 4, 5, 6})),
            ("CHEMBL_NONE", make_fp({500, 501, 502, 503})),
            ("CHEMBL_SUPER", make_fp({1, 2, 3, 4, 5, 6})),
        ]
        return FingerprintIndex.from_records(records)

    def test_finds_exact_match_first(self, index):
        hits = index.search(make_fp({1, 2, 3, 4}), threshold=0.1)
        assert hits[0].chembl_id == "CHEMBL_EXACT"
        assert hits[0].tanimoto == pytest.approx(1.0)

    def test_results_are_descending(self, index):
        hits = index.search(make_fp({1, 2, 3, 4}), threshold=0.0001)
        scores = [h.tanimoto for h in hits]
        assert scores == sorted(scores, reverse=True)

    def test_threshold_excludes_weak_matches(self, index):
        hits = index.search(make_fp({1, 2, 3, 4}), threshold=0.9)
        assert [h.chembl_id for h in hits] == ["CHEMBL_EXACT"]

    def test_disjoint_compound_never_returned(self, index):
        hits = index.search(make_fp({1, 2, 3, 4}), threshold=0.0001)
        assert "CHEMBL_NONE" not in {h.chembl_id for h in hits}

    def test_limit_truncates_to_best(self, index):
        hits = index.search(make_fp({1, 2, 3, 4}), threshold=0.0001, limit=2)
        assert len(hits) == 2
        assert hits[0].chembl_id == "CHEMBL_EXACT"

    def test_popcount_bound_does_not_drop_valid_hits(self, index):
        """The pruning bound must never discard a compound that would have
        passed. Brute-force every pair and compare against the index."""
        query = make_fp({1, 2, 3, 4})
        threshold = 0.3
        raw = {
            "CHEMBL_EXACT": make_fp({1, 2, 3, 4}),
            "CHEMBL_HALF": make_fp({3, 4, 5, 6}),
            "CHEMBL_NONE": make_fp({500, 501, 502, 503}),
            "CHEMBL_SUPER": make_fp({1, 2, 3, 4, 5, 6}),
        }
        expected = {k for k, v in raw.items() if tanimoto(query, v) >= threshold}
        actual = {h.chembl_id for h in index.search(query, threshold=threshold)}
        assert actual == expected

    def test_empty_query_returns_nothing(self, index):
        assert index.search(make_fp(set()), threshold=0.1) == []

    def test_rejects_wrong_width_query(self, index):
        with pytest.raises(ValueError, match="bytes"):
            index.search(b"\x00" * 8, threshold=0.5)

    def test_rejects_invalid_threshold(self, index):
        with pytest.raises(ValueError, match="threshold"):
            index.search(make_fp({1}), threshold=0.0)

    def test_roundtrip_through_disk(self, index, tmp_path):
        path = tmp_path / "idx.npz"
        index.save(path)
        reloaded = FingerprintIndex.load(path)
        assert len(reloaded) == len(index)
        a = index.search(make_fp({1, 2, 3, 4}), threshold=0.1)
        b = reloaded.search(make_fp({1, 2, 3, 4}), threshold=0.1)
        assert [(h.chembl_id, h.tanimoto) for h in a] == [(h.chembl_id, h.tanimoto) for h in b]
