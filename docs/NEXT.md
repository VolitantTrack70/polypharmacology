# What's done, what's next

## The whole stack runs

All three services run together and answer a real query end to end:

```
POST /api/offtargets  { smiles: <imatinib>, tanimoto: 0.35 }
  -> 2 similar compounds (imatinib 1.000, nilotinib 0.517)
  -> 3 human targets (ABL1, KIT, PDGFRA)
  -> 100 Reactome pathways
  -> 106 nodes / 126 edges in ~2 ms
```

The non-human target is correctly excluded by the organism filter, and
`/api/resolve` identifies a pasted structure as `CHEMBL941 / IMATINIB`.

| Component | State |
|---|---|
| Ingest pipeline (Python + RDKit) | Running. 42 tests. |
| PostgreSQL 16 schema + views | Running. Migrations 001–003, 005 applied. |
| chemworker (Python gRPC) | Running on `:50051`. |
| Rust API (axum 0.8) | Running on `:8080`. 11 tests. |
| SvelteKit UI | Running on `:5173`, confirmed against the live API. |
| Reactome data | Loaded: 2,883 pathways, 142,108 protein–pathway links. |
| ChEMBL data | **Fixture only.** Real download in progress. |
| Apache AGE overlay | Skipped — needs Docker. Nothing depends on it. |
| Kafka | Provisioned, unused. Fingerprinting uses a local process pool. |

## The finding that mattered

The threshold problem now reproduces in the product, not a side script:

```
search <imatinib> --threshold 0.35  ->  imatinib 1.000, nilotinib 0.517
search <imatinib> --threshold 0.85  ->  imatinib only
```

At the blueprint's original 0.85, the textbook polypharmacology pair is
invisible. See [0002](decisions/0002-similarity-thresholds.md).

## Bugs that only surfaced by running things

Worth recording, because each was invisible to inspection:

1. **`biological_domain` was NULL for all 2,883 pathways.** Migration 003's
   denormalising UPDATE runs once, against a table still empty at migration
   time. Worse, `refresh_derived` rebuilt `target_pathway` (which reads that
   column) *before* `pathway_domain` computed it. Refresh is now explicitly
   ordered.
2. **Staging tables silently required every column.** `CREATE TABLE (LIKE ...)`
   copies NOT NULL but *not* defaults, so `withdrawn_flag NOT NULL DEFAULT
   FALSE` became mandatory-with-no-default and any Parquet omitting it failed
   the COPY. Now `CREATE TABLE AS SELECT ... WITH NO DATA`.
3. **Re-ingest silently ignored updates.** `ON CONFLICT DO NOTHING` meant a
   newer ChEMBL release would never refresh a corrected structure. Now upserts.
4. **A healthy database reported itself down.** `SELECT 1` is INT4; it was being
   decoded as `i64`.
5. **Relative `DATA_*` paths resolved against cwd**, so `parse` and `load`
   disagreed about where the Parquet lived depending on where they were run.
6. **grpcio emits a flat `import chemworker_pb2`** that cannot resolve inside a
   package — the worker could not import at all. `gen_proto.py` patches it.

## Next

1. **Finish the ChEMBL ingest.** The download is running. Afterwards:
   `parse --release 35`, then `fingerprint` (~2.4M structures, expect
   20–60 min across a process pool), `load`, `index`.
2. **Pathway noise.** ABL1 alone maps to ~52 pathways, so the cascade graph is
   visually dense. Rolling up to `biological_domain` by default, with
   drill-down, would help more than any styling change.
3. **`project-graph` is still unimplemented** — migration 004 creates AGE labels
   but nothing populates them. Only matters if you want Cypher.
4. **No auth, permissive CORS.** Fine on localhost, not beyond it.
5. **UniProt ingestion still stubbed.** ChEMBL already supplies accessions,
   names and organism; UniProt would only add clean gene symbols.
