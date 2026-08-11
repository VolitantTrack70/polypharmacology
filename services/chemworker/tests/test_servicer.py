"""Tests for the gRPC servicer.

These exercise the servicer directly rather than over a socket -- the logic
under test is the request handling and the error contract, not gRPC's transport.

The error contract is the important part and is easy to regress:
  * a malformed SMILES to Standardize is a CLIENT condition, reported as
    ok=false with a populated error field, because the UI shows it inline;
  * a malformed SMILES to SimilaritySearch aborts with INVALID_ARGUMENT,
    because there is no partial result to hand back.
"""

from __future__ import annotations

import numpy as np
import pytest

grpc = pytest.importorskip("grpc")

from chemmed_ingest.chem.fingerprint import fingerprint_bytes, standardize  # noqa: E402
from chemmed_ingest.chem.similarity import FingerprintIndex  # noqa: E402
from chemmed_ingest.config import FP  # noqa: E402

from chemworker import chemworker_pb2 as pb  # noqa: E402
from chemworker.server import ChemWorkerServicer  # noqa: E402

IMATINIB = "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"
NILOTINIB = "Cc1cn(-c2cc(NC(=O)c3ccc(C)c(Nc4nccc(-c5cccnc5)n4)c3)cc(C(F)(F)F)c2)cn1"
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"


class Aborted(Exception):
    """Raised by the fake context so tests can assert on abort()."""

    def __init__(self, code, details):
        self.code = code
        self.details = details
        super().__init__(f"{code}: {details}")


class FakeContext:
    """Minimal stand-in for grpc.ServicerContext.

    Real gRPC raises out of abort(); this mirrors that so a handler which
    aborts and then keeps going is caught rather than silently continuing.
    """

    def __init__(self):
        self.code = None
        self.details = None

    def abort(self, code, details):
        self.code = code
        self.details = details
        raise Aborted(code, details)


@pytest.fixture(scope="module")
def servicer() -> ChemWorkerServicer:
    records = [
        ("CHEMBL941", fingerprint_bytes(standardize(IMATINIB).mol)),
        ("CHEMBL255863", fingerprint_bytes(standardize(NILOTINIB).mol)),
        ("CHEMBL25", fingerprint_bytes(standardize(ASPIRIN).mol)),
    ]
    return ChemWorkerServicer(FingerprintIndex.from_records(records))


class TestStandardize:
    def test_valid_smiles(self, servicer):
        resp = servicer.Standardize(pb.StandardizeRequest(smiles=ASPIRIN), FakeContext())
        assert resp.ok
        assert resp.canonical_smiles
        assert resp.inchikey.startswith("BSYNRYMUTXBXSQ")
        assert resp.parent_inchikey == resp.inchikey.split("-")[0]
        assert resp.error == ""

    def test_invalid_smiles_is_reported_not_aborted(self, servicer):
        """A bad structure is a client condition the UI shows inline. It must
        NOT come back as a gRPC error status."""
        ctx = FakeContext()
        resp = servicer.Standardize(pb.StandardizeRequest(smiles="not a molecule"), ctx)
        assert resp.ok is False
        assert resp.error
        assert ctx.code is None, "Standardize must not abort on bad input"

    def test_empty_smiles_is_reported_not_aborted(self, servicer):
        ctx = FakeContext()
        resp = servicer.Standardize(pb.StandardizeRequest(smiles=""), ctx)
        assert resp.ok is False
        assert ctx.code is None

    def test_salt_is_stripped_to_parent(self, servicer):
        """Imatinib mesylate must standardise to the same structure as the
        free base, or the salt and parent index as different compounds."""
        salt = servicer.Standardize(
            pb.StandardizeRequest(smiles="CS(=O)(=O)O." + IMATINIB), FakeContext()
        )
        base = servicer.Standardize(pb.StandardizeRequest(smiles=IMATINIB), FakeContext())
        assert salt.ok and base.ok
        assert salt.canonical_smiles == base.canonical_smiles
        assert salt.inchikey == base.inchikey


class TestSimilaritySearch:
    def test_finds_self_and_analogue(self, servicer):
        resp = servicer.SimilaritySearch(
            pb.SimilaritySearchRequest(smiles=IMATINIB, threshold=0.35, limit=10),
            FakeContext(),
        )
        hits = {h.chembl_id: h.tanimoto for h in resp.hits}
        assert hits["CHEMBL941"] == pytest.approx(1.0)
        assert "CHEMBL255863" in hits
        assert 0.45 < hits["CHEMBL255863"] < 0.60

    def test_the_blueprint_threshold_finds_only_the_query(self, servicer):
        """Regression guard for docs/decisions/0002. At 0.85 the canonical
        polypharmacology pair disappears -- if this ever starts returning
        nilotinib, the fingerprint configuration has changed."""
        resp = servicer.SimilaritySearch(
            pb.SimilaritySearchRequest(smiles=IMATINIB, threshold=0.85, limit=10),
            FakeContext(),
        )
        assert [h.chembl_id for h in resp.hits] == ["CHEMBL941"]

    def test_results_are_descending(self, servicer):
        resp = servicer.SimilaritySearch(
            pb.SimilaritySearchRequest(smiles=IMATINIB, threshold=0.01, limit=10),
            FakeContext(),
        )
        scores = [h.tanimoto for h in resp.hits]
        assert scores == sorted(scores, reverse=True)

    def test_limit_is_respected(self, servicer):
        resp = servicer.SimilaritySearch(
            pb.SimilaritySearchRequest(smiles=IMATINIB, threshold=0.01, limit=1),
            FakeContext(),
        )
        assert len(resp.hits) == 1

    def test_canonical_smiles_is_echoed(self, servicer):
        resp = servicer.SimilaritySearch(
            pb.SimilaritySearchRequest(smiles=ASPIRIN, threshold=0.5, limit=5),
            FakeContext(),
        )
        assert resp.canonical_smiles == standardize(ASPIRIN).canonical_smiles

    def test_searched_count_reflects_index_size(self, servicer):
        resp = servicer.SimilaritySearch(
            pb.SimilaritySearchRequest(smiles=ASPIRIN, threshold=0.5, limit=5),
            FakeContext(),
        )
        assert resp.searched == 3

    def test_bad_smiles_aborts(self, servicer):
        """Unlike Standardize -- there is no partial result to return."""
        with pytest.raises(Aborted) as exc:
            servicer.SimilaritySearch(
                pb.SimilaritySearchRequest(smiles="???", threshold=0.4), FakeContext()
            )
        assert exc.value.code == grpc.StatusCode.INVALID_ARGUMENT

    def test_missing_query_aborts(self, servicer):
        with pytest.raises(Aborted) as exc:
            servicer.SimilaritySearch(
                pb.SimilaritySearchRequest(threshold=0.4), FakeContext()
            )
        assert exc.value.code == grpc.StatusCode.INVALID_ARGUMENT

    def test_threshold_above_one_aborts(self, servicer):
        with pytest.raises(Aborted) as exc:
            servicer.SimilaritySearch(
                pb.SimilaritySearchRequest(smiles=ASPIRIN, threshold=1.5), FakeContext()
            )
        assert exc.value.code == grpc.StatusCode.INVALID_ARGUMENT

    def test_zero_threshold_defaults_rather_than_aborting(self, servicer):
        """proto3 scalars default to 0 when omitted, so an unset threshold is
        indistinguishable from 0. It must fall back to the sane default, not
        be treated as an invalid value."""
        resp = servicer.SimilaritySearch(
            pb.SimilaritySearchRequest(smiles=IMATINIB, limit=10), FakeContext()
        )
        assert any(h.chembl_id == "CHEMBL941" for h in resp.hits)


class TestStatus:
    def test_reports_index_and_signature(self, servicer):
        resp = servicer.Status(pb.StatusRequest(), FakeContext())
        assert resp.index_loaded is True
        assert resp.compound_count == 3
        assert resp.fp_signature == FP.signature
        assert resp.version


class TestIndexGuards:
    def test_index_refuses_mismatched_fingerprint_width(self):
        """A width mismatch means the index and the query disagree about the
        fingerprint configuration, which yields silently wrong scores."""
        with pytest.raises(ValueError, match="width mismatch"):
            FingerprintIndex.from_records([("X", b"\x00" * 8)])

    def test_index_rejects_wrong_dtype(self):
        with pytest.raises(TypeError, match="uint64"):
            FingerprintIndex(
                np.array(["A"], dtype=object),
                np.zeros((1, 32), dtype=np.uint32),
                np.zeros(1, dtype=np.int32),
            )
