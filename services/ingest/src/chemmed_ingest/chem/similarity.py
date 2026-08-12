"""Query-time Tanimoto similarity search over the full compound set.

Similarity is computed per query rather than stored as edges: all-pairs over
2.4M ChEMBL structures is ~2.9e12 comparisons, and the packed matrix is only
~600MB in RAM. This also makes the threshold a live query parameter.

Tanimoto over bit vectors:   Ts(A,B) = c / (a + b - c)
where a = popcount(A), b = popcount(B), c = popcount(A AND B).

An exact popcount bound runs first: Ts(A,B) <= min(a,b)/max(a,b), so anything
outside [t*a, a/t] is skipped untouched. It is free but not what makes this
fast -- measured on real ChEMBL it prunes 51.7% at t=0.85 but only 1.5% at the
0.40 default, since drug-like molecules have similar bit counts.
See docs/decisions/0001.
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


@dataclass(frozen=True)
class SearchResult:
    """Hits plus what the scan actually did, so callers can tell truncation
    and pruning apart from a genuinely small result."""

    hits: list[SimilarityHit]
    searched: int
    candidates: int
    total_matches: int

    @property
    def truncated(self) -> bool:
        return self.total_matches > len(self.hits)


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
        """Compounds scoring >= `threshold`, best first."""
        return self.search_detailed(query_fp, threshold, limit).hits

    def search_detailed(
        self,
        query_fp: bytes,
        threshold: float = 0.40,
        limit: int | None = 500,
    ) -> SearchResult:
        """As `search`, but also reports how many compounds survived the
        popcount bound and how many matched in total before `limit`.

        `limit` truncates the *result*, not the scan, so hits are always the
        globally best matches -- but callers need `total_matches` to know that
        truncation happened at all.
        """
        if not 0.0 < threshold <= 1.0:
            raise ValueError(f"threshold must be in (0, 1], got {threshold}")

        q = np.frombuffer(query_fp, dtype=np.uint8)
        if q.size != FP.n_bytes:
            raise ValueError(f"query is {q.size} bytes, index expects {FP.n_bytes}")
        q = np.ascontiguousarray(q).view(np.uint64)
        q_pop = int(_popcount_u64(q))

        empty = SearchResult(hits=[], searched=len(self), candidates=0, total_matches=0)
        if q_pop == 0:
            return empty  # e.g. a bare atom: no bits set, no meaningful match

        # Exact popcount bound -- prune before touching any bits.
        lo = q_pop * threshold
        hi = q_pop / threshold
        candidates = np.flatnonzero((self.popcounts >= lo) & (self.popcounts <= hi))
        if candidates.size == 0:
            return empty

        inter = _popcount_u64(self.matrix[candidates] & q)
        union = self.popcounts[candidates].astype(np.int32) + q_pop - inter
        # union can only be 0 if both vectors are empty, already excluded above.
        scores = inter / union

        keep = np.flatnonzero(scores >= threshold)
        if keep.size == 0:
            return SearchResult(
                hits=[], searched=len(self), candidates=int(candidates.size), total_matches=0
            )

        order = keep[np.argsort(-scores[keep], kind="stable")]
        if limit is not None:
            order = order[:limit]

        return SearchResult(
            hits=[
                SimilarityHit(chembl_id=str(self.ids[candidates[i]]), tanimoto=float(scores[i]))
                for i in order
            ],
            searched=len(self),
            candidates=int(candidates.size),
            total_matches=int(keep.size),
        )

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
