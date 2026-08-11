# What's done, what's next

## Verified working, end to end

The pipeline runs for real against a live PostgreSQL 16 instance:

```
parse -> fingerprint -> index -> search        (no database needed)
migrate -> load -> binds_to                    (against Postgres)
```

- **Cheminformatics core** — standardisation (desalting, neutralisation,
  canonicalisation) and Morgan fingerprinting. Imatinib mesylate collapses to
  the free base (Tanimoto 1.000).
- **Similarity search through the CLI, on real structures.** Querying imatinib
  at cutoff 0.35 returns imatinib (1.000) and nilotinib (0.517). The same query
  at the blueprint's original 0.85 returns **only imatinib** — the threshold
  problem from [0002](decisions/0002-similarity-thresholds.md) now reproduces
  in the actual product, not just a side script.
- **ChEMBL parser** — 32 tests against a schema-accurate SQLite fixture
  (`tests/chembl_fixture.py`), pinning down which joins are INNER vs LEFT and
  which rows get dropped.
- **Schema + quality filters.** Migrations 001–003 applied; `binds_to`
  aggregation verified against hand-computed expectations:
  - A compound/target pair with 4 raw activities counts **2** — the `>`-relation
    row (a non-binder) and the low-confidence-assay row are correctly excluded.
  - A pair whose raw maximum pChEMBL is 9.0 shows **7.0** in the view — the
    row flagged with a `data_validity_comment` is correctly excluded.
- **Idempotent loads.** Re-running `load` inserts 0 rows rather than erroring.

## Database setup as it actually stands

PostgreSQL 16 is installed **natively**, not via Docker, so migration 004
(Apache AGE) is **skipped automatically** — the runner detects the missing
extension and reports it rather than failing.

This costs nothing today. Per [0003](decisions/0003-postgres-over-neo4j.md) the
relational tables are the system of record and answer every query the API
makes; `project-graph` was never implemented, so the AGE graph would be empty
regardless. `docker-compose.yml` still provisions the `apache/age` image for
when Cypher is actually wanted.

## Still written but unexecuted

| Component | Blocked on | Likely first friction |
|---|---|---|
| Rust API | `cargo`, `protoc` | axum 0.8 path syntax; sqlx array binding in the `UNNEST` query |
| chemworker gRPC server | `protoc` | Needs stubs generated before it will import |
| Apache AGE overlay | Docker | The property-index syntax in `004` is the least certain part |
| Kafka fan-out | Docker | Provisioned but unused; fingerprinting uses a local process pool |

## Suggested order from here

1. **Real ChEMBL data.** `chemmed-ingest download --source chembl --release 35`
   then `parse --limit 50000`. This is a ~5 GB download expanding to tens of GB,
   so it is a deliberate decision rather than something to run casually.
   Everything upstream of it is already proven on the fixture.
2. **Reactome**, which is small (a few tens of MB) and unlocks the pathway third
   of the cascade — currently the only part of the graph with no data at all.
3. **Rust API**, once `cargo` and `protoc` are installed. Least logic, most
   toolchain setup.

## Known gaps worth tracking

- **Pathway data is entirely absent.** Reactome has not been downloaded, so
  `target_pathway` is empty and the cascade currently stops at targets.
- **`project-graph` is unimplemented.** Migration 004 creates labels; nothing
  populates them.
- **UniProt ingestion is stubbed.** ChEMBL supplies accessions, names and
  organism already; UniProt would only add clean gene symbols.
- **No auth on the API**, and CORS is permissive. Fine on localhost.
- **`known_name` in `/api/resolve` always returns null** — ChEMBL's
  `molecule_dictionary.pref_name` is parsed in the fixture but not yet loaded
  into `chem.compound`.
