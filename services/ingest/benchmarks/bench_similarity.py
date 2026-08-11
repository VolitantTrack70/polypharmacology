"""Benchmark the similarity index at full ChEMBL scale.

docs/decisions/0001 claims that dropping precomputed similarity edges is safe
because ~2.4M fingerprints fit in RAM and scan in well under a second. That is
the load-bearing claim of the whole architecture, so it should be measured
rather than asserted.

This builds a synthetic index at ChEMBL scale and times real searches. The
fingerprints are random but realistically sparse: drug-like molecules set
roughly 40-90 of 2048 Morgan bits, and sparsity is what governs both the memory
footprint and how well the popcount bound prunes.

    .venv/Scripts/python benchmarks/bench_similarity.py
    .venv/Scripts/python benchmarks/bench_similarity.py --n 500000
"""

from __future__ import annotations

import argparse
import time

import numpy as np

from chemmed_ingest.chem.similarity import FingerprintIndex, _popcount_u64
from chemmed_ingest.config import FP

CHEMBL_35_COMPOUNDS = 2_400_000
CHUNK = 50_000
# Observed range for drug-like molecules at Morgan r=2 / 2048 bits.
BITS_MIN, BITS_MAX = 40, 90


def synth_fingerprints(n: int, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Generate `n` realistically sparse packed fingerprints.

    Built in chunks: a dense boolean matrix of the full set would need ~5 GB,
    while the packed result needs ~600 MB.
    """
    rng = np.random.default_rng(seed)
    n_bytes = FP.n_bytes
    out = np.empty((n, n_bytes), dtype=np.uint8)

    done = 0
    while done < n:
        size = min(CHUNK, n - done)
        dense = np.zeros((size, FP.n_bits), dtype=np.uint8)
        counts = rng.integers(BITS_MIN, BITS_MAX + 1, size=size)
        for i in range(size):
            idx = rng.choice(FP.n_bits, size=counts[i], replace=False)
            dense[i, idx] = 1
        out[done : done + size] = np.packbits(dense, axis=1)
        done += size
        print(f"\r  generating... {done:,}/{n:,}", end="", flush=True)
    print()

    ids = np.array([f"CHEMBL{i}" for i in range(n)], dtype=object)
    return ids, out


def build(n: int) -> FingerprintIndex:
    print(f"Building synthetic index of {n:,} fingerprints ({FP.signature})")
    t0 = time.perf_counter()
    ids, packed = synth_fingerprints(n)
    matrix = np.ascontiguousarray(packed).view(np.uint64)
    popcounts = _popcount_u64(matrix)
    index = FingerprintIndex(ids, matrix, popcounts)
    print(f"  built in {time.perf_counter() - t0:.1f}s")
    return index


def report_memory(index: FingerprintIndex) -> None:
    mb = index.matrix.nbytes / 1e6
    pc = index.popcounts.nbytes / 1e6
    print("\nMemory")
    print(f"  fingerprint matrix   {mb:8.1f} MB")
    print(f"  popcount vector      {pc:8.1f} MB")
    print(f"  total (excl. ids)    {mb + pc:8.1f} MB")


def bench(index: FingerprintIndex, thresholds: list[float], repeats: int = 5) -> None:
    rng = np.random.default_rng(12345)
    print(f"\nSearch latency over {len(index):,} compounds ({repeats} runs each)")
    print(f"  {'threshold':>10} {'median':>10} {'p95':>10} {'candidates':>12} {'pruned':>8}")

    for t in thresholds:
        times: list[float] = []
        candidates = 0
        for _ in range(repeats):
            # A fresh query each run so no single popcount is favoured.
            q_bits = np.zeros(FP.n_bits, dtype=np.uint8)
            q_bits[rng.choice(FP.n_bits, size=rng.integers(BITS_MIN, BITS_MAX), replace=False)] = 1
            query = np.packbits(q_bits).tobytes()

            t0 = time.perf_counter()
            index.search(query, threshold=t, limit=500)
            times.append((time.perf_counter() - t0) * 1000)

            # Recompute what the popcount bound admitted, to show its effect.
            q_pop = int(_popcount_u64(np.frombuffer(query, dtype=np.uint8).view(np.uint64)))
            lo, hi = q_pop * t, q_pop / t
            candidates = int(((index.popcounts >= lo) & (index.popcounts <= hi)).sum())

        times.sort()
        median = times[len(times) // 2]
        p95 = times[min(len(times) - 1, int(len(times) * 0.95))]
        pruned = 100.0 * (1 - candidates / len(index))
        print(f"  {t:>10.2f} {median:>9.1f}ms {p95:>9.1f}ms {candidates:>12,} {pruned:>7.1f}%")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=CHEMBL_35_COMPOUNDS)
    ap.add_argument("--repeats", type=int, default=5)
    args = ap.parse_args()

    index = build(args.n)
    report_memory(index)
    bench(index, [0.30, 0.40, 0.55, 0.70, 0.85], repeats=args.repeats)

    print(
        "\nADR 0001 claims ~600 MB resident and a sub-second scan at this scale.\n"
        "Compare the figures above against that claim."
    )


if __name__ == "__main__":
    main()
