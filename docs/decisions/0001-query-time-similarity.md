# 0001 — Similarity is computed at query time, not stored as edges

**Status:** Accepted · **Supersedes:** the `IS_STRUCTURALLY_SIMILAR_TO` edge in the original blueprint

## Problem

The blueprint modelled structural similarity as a materialised graph edge:

```
(:Compound)-[:IS_STRUCTURALLY_SIMILAR_TO {tanimoto_coefficient}]->(:Compound)
```

ChEMBL contains ~2.4M distinct structures. All-pairs similarity is

```
(2.4e6 choose 2) ≈ 2.9e12 comparisons
```

That is not a tuning problem, it is an infeasible one. Restricting to compounds
that carry bioactivity data reduces the count by less than an order of magnitude
and leaves it still infeasible.

There is a second, worse problem. A stored edge has to be computed against a
*fixed* threshold. The blueprint's own UI calls for a slider that adjusts the
cutoff interactively — which a materialised edge set at a baked-in threshold
cannot serve. The storage design and the interaction design contradict each other.

## Decision

Do not store similarity edges. Compute similarity at query time against an
in-memory packed fingerprint index.

```
2.4e6 compounds × 2048 bits = ~600 MB as uint64
```

That fits in RAM on any machine that can run the rest of this stack. A vectorised
popcount scan over it completes in well under a second.

The scan is preceded by an **exact** prune. Since `c ≤ min(a,b)` and
`(a + b − c) ≥ max(a,b)`:

```
Ts(A,B) ≤ min(a,b) / max(a,b)
```

so any compound whose popcount falls outside `[t·a, a/t]` cannot reach threshold
`t` and is skipped without its bits ever being touched. The bound is exact, not
heuristic — it never drops a hit that would have qualified, and
`test_popcount_bound_does_not_drop_valid_hits` asserts exactly that against a
brute-force baseline.

Search is exact. No LSH, no ANN. At this scale a linear scan is fast enough that
trading away exactness would buy nothing.

## Measured

`benchmarks/bench_similarity.py --real`, against the loaded ChEMBL 35 index —
**2,474,571 real fingerprints**, 12-core desktop:

| | |
|---|---|
| Fingerprint matrix | **633 MB** |
| Popcount vector | 10 MB |
| Total resident | **643 MB** |

| Threshold | Median | p95 | Pruned by the bound |
|---|---|---|---|
| 0.30 | 481 ms | 499 ms | 2.5% |
| **0.40** (default) | **391 ms** | 417 ms | **1.5%** |
| 0.55 | 252 ms | 375 ms | 41.5% |
| 0.70 | 308 ms | 334 ms | 43.1% |
| 0.85 | 51 ms | 218 ms | 51.7% |

**The memory and latency claims hold** — 643 MB against a predicted ~600 MB, and
sub-second at every threshold with room to spare.

**The pruning claim did not.** An earlier draft asserted the bound "discards
70–90% of the database on a typical query". At the default 0.40 it discards
**1.5%**. Drug-like molecules have similar bit counts, so `min/max` rarely falls
below ~0.44 and nearly everything survives the bound.

The prune stays — free, exact, and it does real work at 0.55+ — but it is not
what makes this fast. A vectorised popcount over 633 MB simply is fast.

Synthetic fingerprints (`--n 2400000`, no `--real`) gave 624 MB and 0.0% pruning
at 0.40. Real data prunes better in the middle of the range (41.5% vs 13.7% at
0.55) because ChEMBL spans fragments to large natural products, but the
conclusion at the default threshold is unchanged.

## Consequences

**Good**

- The threshold becomes a live query parameter. The UI slider genuinely
  re-searches — see [0002](0002-similarity-thresholds.md).
- An entire pipeline stage disappears, along with the storage it needed.
- Adding a compound means appending one row, not recomputing its edges against
  2.4M others.

**Costs**

- The API process must hold ~600 MB resident, and pay a cold-start load.
- Similarity is not expressible in Cypher, so the traversal is two steps:
  fingerprint search produces a compound set, then the graph query expands it.
  The `Compound -> Target -> Pathway` half remains a genuine graph query.

## Note on FPSim2

[FPSim2](https://github.com/chembl/FPSim2) is the ChEMBL team's own library for
exactly this and implements the same popcount-bound strategy against an on-disk
index. It is wired as an optional dependency. The in-repo NumPy backend exists so
the pipeline has no hard dependency on it and so the algorithm stays legible and
directly testable.
