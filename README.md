# Polypharmacology & Off-Target Binding Knowledge Graph

A research tool for generating off-target hypotheses. Given a compound, it finds
structurally similar molecules, collects the proteins those molecules are known
to bind, and maps those proteins onto biological pathways — surfacing the
`drug → off-target → pathway` cascade.

The underlying premise is the *similarity principle*: structurally similar
molecules tend to bind similar proteins. This is the same reasoning behind the
Similarity Ensemble Approach (Keiser et al., 2007).

> **This tool ranks hypotheses. It does not predict safety.**
> Chemical similarity implies *possible* shared binding, not confirmed binding.
> Absence of a reported interaction is not evidence of absence — ChEMBL is
> heavily biased toward well-studied target families. Every result is a lead to
> investigate, not a finding.

---

## Architecture

```
                    +---------------------------+
                    |   Web UI  (SvelteKit)     |
                    |   cytoscape.js canvas     |
                    +-------------+-------------+
                                  | REST/JSON
                    +-------------v-------------+
                    |   API  (Rust / Axum)      |
                    |   query orchestration     |
                    +--+---------------------+--+
        similarity     |                     |    graph + relational
        search (gRPC)  |                     |    queries
        +--------------v-------+   +---------v------------------+
        |  chemworker (Python) |   |  PostgreSQL 16 + AGE       |
        |  RDKit + FPSim2      |   |  chem.* tables (truth)     |
        |  in-memory FP index  |   |  polypharm graph (proj.)   |
        +----------------------+   +---------^------------------+
                                              | COPY
        +-------------------------------------+------------------+
        |  ingest (Python + RDKit + Polars)                       |
        |  download -> parse -> fingerprint -> load -> index      |
        |  Kafka fans fingerprint work across worker processes    |
        +---------------------------------------------------------+
                    ^              ^               ^
                 ChEMBL         Reactome        UniProt
```

**Where the work actually happens:** the expensive step is scanning ~2.4M
fingerprints per query, which lives in `chemworker` as a vectorised popcount over
a packed in-memory matrix. The Rust API orchestrates and serves; it does not do
chemistry.

## Design decisions worth reading before you touch the code

These three changed materially from the original blueprint, each for a reason
documented in full:

| ADR | Decision | Why it matters |
|---|---|---|
| [0001](docs/decisions/0001-query-time-similarity.md) | Similarity is computed at query time, not stored as edges | All-pairs over 2.4M compounds is ~2.9e12 comparisons — infeasible. It also makes the UI threshold slider genuinely live. |
| [0002](docs/decisions/0002-similarity-thresholds.md) | Default Tanimoto cutoff is **0.40**, not 0.85 | At 0.85, **zero** of 28 real drug pairs match — including imatinib/nilotinib (0.517), the textbook polypharmacology case. |
| [0003](docs/decisions/0003-postgres-over-neo4j.md) | PostgreSQL + Apache AGE, not Neo4j | The workload is half relational (21M rows, numeric filters) and half graph. Postgres covers both. |

## Data sources

| Source | What we take | Notes |
|---|---|---|
| **ChEMBL** (v35+, SQLite) | Compounds, targets, proteins, bioactivity | Ships as a relational DB. MW/logP come from `compound_properties` — we do **not** recompute them. |
| **Reactome** | Pathways, hierarchy, protein→pathway | Use `UniProt2Reactome_All_Levels.txt` and `ReactomePathwaysRelation.txt`. |
| **UniProt** | Gene symbols, richer annotation | *Optional.* ChEMBL already carries UniProt accessions via `target_components → component_sequences.accession`. |

## Getting started

**Prerequisites:** Docker, Rust toolchain, Node 20+, Python 3.12 (via `uv`).

```bash
cp .env.example .env
```

```bash
docker compose up -d
```

```bash
cd services/ingest && uv venv --python 3.12 .venv && uv pip install --python .venv -e ".[dev]"
```

Load a small slice end-to-end first — this takes about a minute and validates the
whole pipeline before you commit to the multi-hour full ingest:

```bash
cd services/ingest && .venv/Scripts/chemmed-ingest run-all --release 35 --limit 5000
```

Then the API and UI:

```bash
cd services/api && cargo run
```

```bash
cd services/web && npm install && npm run dev
```

## Repository layout

```
db/migrations/      Schema. 001 core -> 002 indexes/views -> 003 rollup -> 004 AGE graph.
services/ingest/    Python. Parsing, standardisation, fingerprinting, bulk load.
services/chemworker/Python. gRPC similarity-search service holding the FP index.
services/api/       Rust/Axum. HTTP surface, query orchestration.
services/web/       SvelteKit. Structure input, threshold slider, graph canvas.
docs/decisions/     Architecture decision records.
data/               Source dumps and intermediates. Git-ignored.
```

## Status

Early. Schema, cheminformatics core, and ChEMBL parser are implemented and
tested; the similarity engine has full unit coverage including an exactness
check of the pruning bound against brute force. The API and UI are scaffolded.
See the tracking notes in `docs/` for what's next.
