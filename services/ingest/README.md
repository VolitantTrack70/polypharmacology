# chemmed-ingest

Ingestion + cheminformatics pipeline. Parses the public source dumps, standardises
and fingerprints structures, and bulk-loads Postgres.

## Setup

```bash
uv venv --python 3.12 .venv
```

```bash
uv pip install --python .venv -e ".[dev,fpsim,kafka]"
```

RDKit wheels do not yet exist for CPython 3.14, so the venv pins 3.12 explicitly.

## Pipeline stages

```
download  ->  parse  ->  fingerprint  ->  load  ->  index  ->  project-graph
```

| Stage | What it does |
|---|---|
| `download` | Fetches ChEMBL SQLite + Reactome TSVs into `data/raw/`. Resumable. |
| `parse` | ChEMBL SQLite / Reactome TSV → normalised Parquet in `data/processed/`. |
| `fingerprint` | Standardise + Morgan FP over every structure, in a process pool (optionally fanned out over Kafka). |
| `load` | `COPY` the Parquet into Postgres. Indexes are built *after*, not before. |
| `index` | Build the packed fingerprint matrix used for query-time similarity search. |
| `project-graph` | Project the relational tables into the Apache AGE property graph. |

Run the whole thing:

```bash
.venv/Scripts/chemmed-ingest run-all --release 35
```

Or start small — this loads a few thousand compounds end-to-end in about a minute
and is the right way to validate the pipeline before committing to the full dump:

```bash
.venv/Scripts/chemmed-ingest run-all --release 35 --limit 5000
```

## Scale expectations

Measured on ChEMBL 35, 12-core desktop:

| | |
|---|---|
| Tarball | 4.99 GB |
| Extracted SQLite | 26 GB |
| Compounds / structures | 2,496,335 / 2,474,590 |
| Activities (all) | 21,123,501 |
| Targets | 16,003 |
| Parse (compounds) | ~50 s |

Use `--limit` to iterate without paying for the full run.

The download is resumable and checksum-verified. Extraction is not — if it is
interrupted you get a truncated database that only reveals itself as
`database disk image is malformed`. Delete `data/raw/chembl_<release>/` and
re-run `download`; the tarball is reused.

## Notes that will save you a debugging session

- **Do not change `FP_RADIUS` / `FP_NBITS` without rebuilding.** The index stamps
  its configuration and refuses to load under a mismatch, but stored BYTEA
  fingerprints in Postgres are not automatically invalidated.
- **Indexes are created after bulk load.** Running `002_indexes_and_views.sql`
  before the `COPY` makes a 21M-row load several times slower.
- **Standardisation is shared between ingest and query.** Both paths call
  `chem.fingerprint.standardize()`. If you add a rule there, the index is stale.

## Tests

```bash
.venv/Scripts/python -m pytest -v
```

The similarity tests deliberately avoid RDKit so the core maths can be checked
without the heavy dependency.
