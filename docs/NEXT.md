# What's done, what's next

## Verified working

- **Cheminformatics core** — standardisation (desalting, neutralisation,
  canonicalisation) and Morgan fingerprinting. Confirmed against real drugs:
  imatinib mesylate collapses to the free base (Tanimoto 1.000).
- **Similarity engine** — packed-uint64 index with exact popcount-bound pruning.
  15 unit tests pass, including a check that the bound never drops a hit a brute
  force scan would have kept.
- **Threshold calibration** — measured, not assumed. See
  [0002](decisions/0002-similarity-thresholds.md).
- **ChEMBL + Reactome parsers** — stream to Parquet with explicit Arrow schemas.
- **CLI** — all eight stages wired and loading cleanly; `ruff` clean.

## Written but not yet executed

Nothing here is broken as far as I know — it simply hasn't been run, because
the toolchains aren't installed on this machine yet.

| Component | Blocked on | Likely first friction |
|---|---|---|
| Postgres schema (4 migrations) | Docker | AGE image tag; the property-index syntax in `004` is the least certain part |
| Rust API | `cargo`, `protoc` | axum 0.8 path syntax; sqlx array binding in the `UNNEST` query |
| chemworker gRPC server | `protoc` | Needs stubs generated before it will import |
| SvelteKit UI | Node 20+ | `cytoscape-fcose` import shape under Vite SSR |

## Suggested order

1. **Install Docker**, `docker compose up -d`, then `chemmed-ingest migrate`.
   This exercises the schema first, since everything downstream depends on it.
2. **`chemmed-ingest run-all --limit 5000`.** A few thousand compounds end to
   end proves the pipeline in about a minute. Do this before the full ingest.
3. **`chemmed-ingest search "CC(=O)Oc1ccccc1C(=O)O"`.** Similarity search from
   the shell, no API or UI needed — the fastest way to confirm the index is real.
4. **Install Node, run the UI against a stubbed API** if you want to iterate on
   the frontend before the Rust side compiles.
5. **Rust last.** It is the layer with the least logic in it and the most
   toolchain setup.

## Known gaps worth tracking

- **`project-graph` is not implemented.** The AGE migration creates the labels
  but nothing populates them. The relational path answers every current query,
  so this is not blocking — it matters when you want Cypher.
- **UniProt ingestion is stubbed out.** ChEMBL supplies accessions, names, and
  organism already; UniProt would only add clean gene symbols. Listed in the
  blueprint, deliberately deferred as near-redundant.
- **Kafka is provisioned but unused.** The fingerprint stage currently uses a
  local process pool, which is the right tool for a single machine. The topics
  and config exist for when fan-out across hosts is actually wanted.
- **No auth on the API**, and CORS is permissive. Fine on localhost; both need
  fixing before this is reachable from anywhere else.
- **`known_name` in `/api/resolve` always returns null** — ChEMBL's preferred
  compound names aren't loaded yet (`molecule_dictionary.pref_name`).
