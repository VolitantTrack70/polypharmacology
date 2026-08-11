"""Query-time Tanimoto similarity search over the full compound set.

WHY THIS EXISTS INSTEAD OF PRECOMPUTED EDGES
--------------------------------------------
The original design called for materialised `IS_STRUCTURALLY_SIMILAR_TO`
edges. ChEMBL holds ~2.4M distinct structures, so all-pairs is ~2.9e12
comparisons -- infeasible to compute and pointless to store.

It is also unnecessary. 2.4M x 2048 bits packs to ~600MB of uint64, which sits
comfortably in RAM. A vectorised popcount scan over that runs in well under a
second, which means similarity becomes a *live* query parameter: the UI's
threshold slider re-searches instead of being locked to whatever cutoff the
ingest happened to bake in.

ALGORITHM
---------
Tanimoto over bit vectors:   Ts(A,B) = c / (a + b - c)
where a = popcount(A), b = popcount(B), c = popcount(A AND B).

The scan is preceded by a cheap exact bound. Since c <= min(a, b) and
(a + b - c) >= max(a, b):

    Ts(A,B) <= min(a, b) / max(a, b)

so any compound whose popcount falls outside [t*a, a/t] cannot reach threshold
t and is skipped without touching its bits.

How much that prune actually saves depends entirely on the threshold, and it is
easy to overestimate. Measured over 2.4M drug-like fingerprints it discards
68.6% at t=0.85 but *nothing* at t=0.40, because drug-like molecules have
similar bit counts and min/max rarely drops below ~0.44. The bound is free and
exact so it stays, but it is not what makes this fast -- a vectorised popcount
over 600 MB is just quick (~400 ms at the default threshold).
See docs/decisions/0001 and benchmarks/bench_similarity.py.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from chemmed_ingest.config import FP

# NumPy 2.0 exposes a vectorised popcount intrinsic. Fall back to a nibble
# lookup table on older builds rather than failing at import.
_HAS_BITWISE_COUNT = hasattr(np, "bitwise_count")
_POPCOUNT_TABLE = np.array([bin(i).count("1") for i in range(256)], dtype=np.uint8)


def _popcount_u64(arr: np.ndarray) -> np.ndarray:
    """Population count over a uint64 array, summing along the last axis."""
    if _HAS_BITWISE_COUNT:
        return np.bitwise_count(arr).sum(axis=-1, dtype=np.int32)
    view = arr.view(np.uint8).reshape(*arr.shape[:-1], -1)
    return _POPCOUNT_TABLE[view].sum(axis=-1, dtype=np.int32)


@dataclass(frozen=True)
class SimilarityHit:
    chembl_id: str
    tanimoto: float


class FingerprintIndex:
    """In-memory packed fingerprint matrix with exact Tanimoto search.

    Exact, not approximate -- no LSH, no ANN. At this scale a linear scan is
    fast enough that giving up exactness would buy nothing.
    """

    __slots__ = ("ids", "matrix", "popcounts", "n_words")

    def __init__(self, ids: np.ndarray, matrix: np.ndarray, popcounts: np.ndarray) -> None:
        if matrix.dtype != np.uint64:
            raise TypeError(f"matrix must be uint64, got {matrix.dtype}")
        if not (len(ids) == matrix.shape[0] == len(popcounts)):
            raise ValueError(
                f"ids, matrix, and popcounts must have equal length; got "
                f"{len(ids)}, {matrix.shape[0]}, {len(popcounts)}"
            )
        self.ids = ids
        self.matrix = matrix
        self.popcounts = popcounts
        self.n_words = matrix.shape[1]

    # -- construction --------------------------------------------------------

    @classmethod
    def from_records(cls, records: list[tuple[str, bytes]]) -> FingerprintIndex:
        """Build from (chembl_id, packed_fp_bytes) pairs."""
        if not records:
            raise ValueError("cannot build an index from zero records")

        n_bytes = FP.n_bytes
        ids = np.array([r[0] for r in records], dtype=object)
        blob = b"".join(r[1] for r in records)
        expected = len(records) * n_bytes
        if len(blob) != expected:
            raise ValueError(
                f"fingerprint width mismatch: got {len(blob)} bytes, "
                f"expected {expected} ({len(records)} x {n_bytes}). "
                f"Index config is {FP.signature}."
            )

        matrix = np.frombuffer(blob, dtype=np.uint8).reshape(len(records), n_bytes)
        matrix = np.ascontiguousarray(matrix).view(np.uint64)
        popcounts = _popcount_u64(matrix)
        return cls(ids, matrix, popcounts)

    def save(self, path: Path) -> None:
        """Persist to a single .npz so restarts don't re-read Postgres."""
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            path,
            ids=self.ids,
            matrix=self.matrix,
            popcounts=self.popcounts,
            signature=np.array(FP.signature),
        )

    @classmethod
    def load(cls, path: Path) -> FingerprintIndex:
        data = np.load(path, allow_pickle=True)
        stored = str(data["signature"])
        if stored != FP.signature:
            raise ValueError(
                f"fingerprint index was built with {stored!r} but the current "
                f"config is {FP.signature!r}. Rebuild the index -- comparing "
                f"across configurations yields silently wrong scores."
            )
        return cls(data["ids"], data["matrix"], data["popcounts"])

    # -- search --------------------------------------------------------------

    def search(
        self,
        query_fp: bytes,
        threshold: float = 0.40,
        limit: int | None = 500,
    ) -> list[SimilarityHit]:
        """Return every compound scoring >= `threshold`, best first.

        `limit` truncates the *result*, not the scan, so results are always the
        globally best matches rather than the first ones encountered.
        """
        if not 0.0 < threshold <= 1.0:
            raise ValueError(f"threshold must be in (0, 1], got {threshold}")

        q = np.frombuffer(query_fp, dtype=np.uint8)
        if q.size != FP.n_bytes:
            raise ValueError(f"query is {q.size} bytes, index expects {FP.n_bytes}")
        q = np.ascontiguousarray(q).view(np.uint64)
        q_pop = int(_popcount_u64(q))

        if q_pop == 0:
            return []  # e.g. a single bare atom; no bits set, no meaningful match

        # Exact popcount bound -- prune before touching any bits.
        lo = q_pop * threshold
        hi = q_pop / threshold
        candidates = np.flatnonzero((self.popcounts >= lo) & (self.popcounts <= hi))
        if candidates.size == 0:
            return []

        inter = _popcount_u64(self.matrix[candidates] & q)
        union = self.popcounts[candidates].astype(np.int32) + q_pop - inter
        # union can only be 0 if both vectors are empty, already excluded above.
        scores = inter / union

        keep = np.flatnonzero(scores >= threshold)
        if keep.size == 0:
            return []

        order = keep[np.argsort(-scores[keep], kind="stable")]
        if limit is not None:
            order = order[:limit]

        return [
            SimilarityHit(chembl_id=str(self.ids[candidates[i]]), tanimoto=float(scores[i]))
            for i in order
        ]

    def __len__(self) -> int:
        return len(self.ids)

    def __repr__(self) -> str:
        mb = self.matrix.nbytes / 1e6
        return f"<FingerprintIndex n={len(self):,} {FP.signature} {mb:.0f}MB>"


def tanimoto(fp_a: bytes, fp_b: bytes) -> float:
    """Single-pair Tanimoto. Convenience for tests and one-off comparisons."""
    a = np.ascontiguousarray(np.frombuffer(fp_a, dtype=np.uint8)).view(np.uint64)
    b = np.ascontiguousarray(np.frombuffer(fp_b, dtype=np.uint8)).view(np.uint64)
    pa, pb = int(_popcount_u64(a)), int(_popcount_u64(b))
    if pa == 0 and pb == 0:
        return 1.0
    c = int(_popcount_u64(a & b))
    return c / (pa + pb - c)
