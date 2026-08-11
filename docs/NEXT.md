# What's done, what's next

## The whole stack runs

All three services run together and answer a real query end to end:

```
POST /api/offtargets  { smiles: <imatinib>, tanimoto: 0.35 }
  -> 2 similar compounds (imatinib 1.000, nilotinib 0.517)
  -> 3 human targets (ABL1, KIT, PDGFRA)
  -> 16 Reactome pathways (default scope), 22 nodes / 23 edges, ~2 ms
```

`/api/resolve` identifies a pasted structure as `CHEMBL941 / IMATINIB`, and the
non-human target is correctly excluded by the organism filter.

| Component | State |
|---|---|
| Ingest pipeline (Python + RDKit) | Running. 73 tests. |
| PostgreSQL 16 schema + views | Running. Migrations 001–003, 005–006 applied. |
| chemworker (Python gRPC) | Running on `:50051`. 17 tests. |
| Rust API (axum 0.8) | Running on `:8080`. 11 tests. |
| SvelteKit UI | Running on `:5173`, verified against the live API. |
| CI (GitHub Actions) | Green. Three parallel jobs. |
| Reactome data | Loaded: 2,883 pathways, 142,108 protein–pathway links. |
| ChEMBL data | **Fixture only.** Real download in progress. |
| Apache AGE overlay | Skipped — needs Docker. Nothing depends on it. |
| Kafka | Provisioned, unused. Fingerprinting uses a local process pool. |

**101 tests total.** Run them with `pytest` in each Python service and
`cargo test` in `services/api`.

## The finding that mattered

The threshold problem reproduces in the product, not a side script:

```
search <imatinib> --threshold 0.35  ->  imatinib 1.000, nilotinib 0.517
search <imatinib> --threshold 0.85  ->  imatinib only
```

At the blueprint's original 0.85 the textbook polypharmacology pair is
invisible. See [0002](decisions/0002-similarity-thresholds.md). This is now
pinned by tests at three levels: the servicer, the API, and the e2e suite.

## Measured, not assumed

[0001](decisions/0001-query-time-similarity.md) claimed the index would be
~600 MB and scan sub-second. `benchmarks/bench_similarity.py` confirms it —
**624 MB, 214–516 ms** across thresholds at 2.4M compounds.

It also **falsified** part of the same ADR. The claim that the popcount bound
"discards 70–90% of the database" holds only at high thresholds; at the default
0.40 it prunes **0%**, because drug-like molecules have similar bit counts. The
ADR now says so. The prune stays (free, exact, helps where cost is highest) but
it is not what makes the search fast.

## Bugs that only surfaced by running things

Each was invisible to inspection:

1. **`biological_domain` was NULL for all 2,883 pathways.** Migration 003's
   denormalising UPDATE runs once, against a table still empty at migration
   time — and `refresh_derived` rebuilt `target_pathway` (which reads that
   column) *before* `pathway_domain` computed it.
2. **Staging tables silently required every column.** `CREATE TABLE (LIKE ...)`
   copies NOT NULL but *not* defaults.
3. **Re-ingest silently ignored updates** — `ON CONFLICT DO NOTHING` meant a
   newer release would never refresh a corrected structure.
4. **A healthy database reported itself down** — `SELECT 1` is INT4, decoded as
   `i64`.
5. **Relative `DATA_*` paths resolved against cwd**, so `parse` and `load`
   disagreed about where the Parquet lived.
6. **grpcio emits a flat `import chemworker_pb2`** that cannot resolve inside a
   package — the worker could not import at all.
7. **A typo'd SMILES returned 503 "service not responding"** instead of 422.
8. **The fingerprint stage materialised all 2.4M results** via `pool.starmap`,
   and the progress bar would have sat at 0% for an hour then jumped to 100%.
9. **`FPSIM2_INDEX_PATH` defaulted to `.h5`** while every call site rewrote it
   to `.npz` — a config option that quietly ignored what you set.

## Next

1. **Finish the ChEMBL ingest.** Download running. Then `parse --release 35`,
   `fingerprint` (~2.4M structures, 20–60 min), `load`, `index`. The parser now
   preflights the schema, so a release mismatch fails in seconds rather than
   hours in.
2. **Re-run the benchmark against the real index** — real ChEMBL has a wider
   popcount distribution than the synthetic set, so pruning at 0.40 should beat
   0%. Worth recording the true figure in ADR 0001.
3. **Expose `pathway_scope` in the UI.** The API supports it; the UI still uses
   the default.
4. **`project-graph` is unimplemented** — migration 004 creates AGE labels but
   nothing populates them. Only matters if you want Cypher.
5. **No auth, permissive CORS.** Fine on localhost, not beyond it.
6. **UniProt ingestion still stubbed.** ChEMBL already supplies accessions,
   names and organism; UniProt would only add clean gene symbols.
