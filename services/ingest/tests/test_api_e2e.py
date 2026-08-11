"""End-to-end tests against the running stack.

Everything else in this suite tests one component in isolation. These verify
that the API, the gRPC worker and Postgres actually work *together* -- the
failures that only appear at the seams: a type mismatch between sqlx and a
Postgres column, a gRPC status mapped to the wrong HTTP code, a graph whose
edges point at nodes that were never emitted.

Skipped unless the API is reachable, so a normal `pytest` run does not fail
when the stack happens to be down.

Start the stack with:
    scripts\\dev-chemworker.cmd
    scripts\\dev-api.cmd
"""

from __future__ import annotations

import os

import pytest
import requests

API = os.environ.get("CHEMMED_API_BASE", "http://localhost:8080/api")
TIMEOUT = 15

IMATINIB = "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"


def _api_up() -> bool:
    try:
        return requests.get(f"{API}/health", timeout=3).ok
    except requests.RequestException:
        return False


pytestmark = pytest.mark.skipif(
    not _api_up(), reason=f"API not reachable at {API}; start it with scripts/dev-api.cmd"
)


def post(path: str, payload: dict) -> requests.Response:
    return requests.post(f"{API}{path}", json=payload, timeout=TIMEOUT)


@pytest.fixture(scope="module")
def cascade() -> dict:
    r = post("/offtargets", {"smiles": IMATINIB, "tanimoto": 0.35, "organism": "Homo sapiens"})
    r.raise_for_status()
    return r.json()


class TestHealth:
    def test_health(self):
        assert requests.get(f"{API}/health", timeout=TIMEOUT).json()["status"] == "ok"

    def test_both_dependencies_report_up(self):
        """Regression guard: `SELECT 1` is INT4, and decoding it as i64 made a
        perfectly healthy database report itself as down."""
        body = requests.get(f"{API}/status", timeout=TIMEOUT).json()
        assert body["database"]["state"] == "up", body["database"]
        assert body["chemworker"]["state"] == "up", body["chemworker"]

    def test_index_signature_is_reported(self):
        body = requests.get(f"{API}/status", timeout=TIMEOUT).json()
        assert body["fp_signature"] == "morgan-r2-2048"
        assert body["compounds_indexed"] >= 1


class TestResolve:
    def test_identifies_a_known_drug(self):
        body = post("/resolve", {"smiles": IMATINIB}).json()
        assert body["inchikey"].startswith("KTUFNOKKBVMGRW")
        assert body["chembl_id"] == "CHEMBL941"
        assert body["known_name"]

    def test_salt_resolves_to_the_parent(self):
        salt = post("/resolve", {"smiles": "CS(=O)(=O)O." + IMATINIB}).json()
        assert salt["chembl_id"] == "CHEMBL941"

    def test_bad_structure_is_422_not_500(self):
        r = post("/resolve", {"smiles": "not-a-molecule"})
        assert r.status_code == 422
        assert r.json()["error"]["code"] == "invalid_structure"


class TestCascade:
    def test_returns_the_expected_shape(self, cascade):
        assert cascade["stats"]["similar_compounds"] >= 1
        assert cascade["stats"]["targets"] >= 1
        assert {"query", "stats", "nodes", "edges"} <= cascade.keys()

    def test_finds_the_query_compound_itself(self, cascade):
        compounds = [n for n in cascade["nodes"] if n["kind"] == "compound"]
        assert any(c["tanimoto"] == pytest.approx(1.0) for c in compounds)

    def test_every_edge_endpoint_is_a_real_node(self, cascade):
        """cytoscape silently discards edges with dangling endpoints, which
        reads as a mysteriously sparse graph rather than an error."""
        ids = {n["id"] for n in cascade["nodes"]}
        for e in cascade["edges"]:
            assert e["source"] in ids, f"dangling source {e['source']}"
            assert e["target"] in ids, f"dangling target {e['target']}"

    def test_node_ids_are_unique(self, cascade):
        ids = [n["id"] for n in cascade["nodes"]]
        assert len(ids) == len(set(ids))

    def test_stats_agree_with_the_payload(self, cascade):
        kinds = [n["kind"] for n in cascade["nodes"]]
        assert cascade["stats"]["similar_compounds"] == kinds.count("compound")
        assert cascade["stats"]["targets"] == kinds.count("target")
        assert cascade["stats"]["pathways"] == kinds.count("pathway")

    def test_organism_filter_excludes_other_species(self, cascade):
        for n in cascade["nodes"]:
            if n["kind"] == "target" and n.get("organism"):
                assert n["organism"] == "Homo sapiens"

    def test_affinity_is_on_binding_edges_only(self, cascade):
        for e in cascade["edges"]:
            if e["kind"] == "binds_to":
                assert e.get("pchembl") is not None
            if e["kind"] == "similar_to":
                assert e.get("pchembl") is None
                assert e.get("tanimoto") is not None


class TestThresholdBehaviour:
    def test_higher_cutoff_never_returns_more(self):
        low = post("/offtargets", {"smiles": IMATINIB, "tanimoto": 0.35}).json()
        high = post("/offtargets", {"smiles": IMATINIB, "tanimoto": 0.85}).json()
        assert high["stats"]["similar_compounds"] <= low["stats"]["similar_compounds"]

    def test_blueprint_threshold_collapses_to_the_query_alone(self):
        """docs/decisions/0002 in executable form: at 0.85 the canonical
        polypharmacology pair disappears."""
        body = post("/offtargets", {"smiles": IMATINIB, "tanimoto": 0.85}).json()
        assert body["stats"]["similar_compounds"] == 1

    def test_no_matches_is_a_result_not_an_error(self):
        """An empty result must keep the response shape so the UI can suggest
        lowering the cutoff rather than showing a failure."""
        r = post("/offtargets", {"smiles": "C", "tanimoto": 0.99})
        assert r.status_code == 200
        body = r.json()
        assert body["stats"]["similar_compounds"] == 0
        assert body["nodes"]  # the query node still comes back


class TestPathwayScope:
    def test_domain_scope_is_smaller_than_all(self):
        args = {"smiles": IMATINIB, "tanimoto": 0.35, "max_pathways_per_target": 500}
        every = post("/offtargets", {**args, "pathway_scope": "all"}).json()
        domain = post("/offtargets", {**args, "pathway_scope": "domain"}).json()
        assert domain["stats"]["pathways"] < every["stats"]["pathways"]

    def test_cap_limits_pathways_per_target(self):
        body = post(
            "/offtargets",
            {"smiles": IMATINIB, "tanimoto": 0.35, "max_pathways_per_target": 2},
        ).json()
        per_target: dict[str, int] = {}
        for e in body["edges"]:
            if e["kind"] == "participates_in":
                per_target[e["source"]] = per_target.get(e["source"], 0) + 1
        assert all(n <= 2 for n in per_target.values()), per_target

    def test_targets_survive_an_aggressive_cap(self):
        """The cap must not delete targets whose pathways were all filtered."""
        loose = post("/offtargets", {"smiles": IMATINIB, "tanimoto": 0.35}).json()
        tight = post(
            "/offtargets",
            {"smiles": IMATINIB, "tanimoto": 0.35, "pathway_scope": "domain",
             "max_pathways_per_target": 1},
        ).json()
        assert tight["stats"]["targets"] == loose["stats"]["targets"]

    def test_invalid_scope_is_rejected(self):
        r = post("/offtargets", {"smiles": IMATINIB, "pathway_scope": "nonsense"})
        assert r.status_code == 400
        assert r.json()["error"]["code"] == "bad_request"


class TestErrors:
    def test_missing_query_is_400(self):
        r = post("/offtargets", {"tanimoto": 0.4})
        assert r.status_code == 400

    def test_out_of_range_tanimoto_is_400(self):
        assert post("/offtargets", {"smiles": IMATINIB, "tanimoto": 1.5}).status_code == 400

    def test_unknown_chembl_id_is_404(self):
        r = post("/offtargets", {"chembl_id": "CHEMBL_NOPE"})
        assert r.status_code == 404

    def test_errors_have_a_consistent_envelope(self):
        r = post("/offtargets", {"smiles": "!!!"})
        body = r.json()
        assert set(body["error"]) >= {"code", "message"}
