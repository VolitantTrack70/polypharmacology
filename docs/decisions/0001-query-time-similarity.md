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
`t` and is skipped without its bits ever being touched. On a typical drug-like
query this discards most of the database before any real work happens. The bound
is exact, not heuristic — it never drops a hit that would have qualified, and
`test_popcount_bound_does_not_drop_valid_hits` asserts exactly that against a
brute-force baseline.

Search is exact. No LSH, no ANN. At this scale a linear scan is fast enough that
trading away exactness would buy nothing.

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
