# 0003 — PostgreSQL + Apache AGE as the storage engine

**Status:** Accepted · **Resolves:** the `Neo4j Enterprise or PostgreSQL + Apache AGE` alternative left open in the blueprint

## Context

The blueprint offered both options without choosing. Both can hold this data.
The question is which fits the actual workload.

## The workload is genuinely two workloads

**Relational.** ~21M activity rows filtered on numeric ranges — pChEMBL bounds,
assay confidence score, activity type, relation operator. Aggregating them into
one edge per (compound, target) pair means grouping millions of rows. This is
textbook relational work and it is precisely where a graph engine is weakest:
Neo4j scanning several million edges to filter on a numeric property is doing
the thing it is worst at.

**Graph.** `Compound → Target → Pathway` expansion, plus the Reactome pathway
hierarchy, which is a genuine variable-depth DAG.

Postgres covers both. Neo4j covers one well and one poorly.

## Scale does not argue for a dedicated graph engine here

After projection the graph holds roughly:

| | Count |
|---|---|
| Compound vertices (with ≥1 qualifying edge) | ~1M |
| Target vertices | ~15k |
| Pathway vertices | ~2.7k |
| `BINDS_TO` edges | ~3M |
| `PARTICIPATES_IN` edges | ~100k |

Neo4j is engineered for billions of edges. Three million is not a scale that
requires index-free adjacency. And critically — per
[0001](0001-query-time-similarity.md) — there are no similarity edges, which is
what would have made this a genuinely large graph.

The traversal is also a *fixed, known* shape: exactly three hops, always the
same labels. Graph databases earn their operational cost on variable-depth,
unknown-shape traversal. Two joins with good indexes is not that.

## Two factors that settled it

**ChEMBL is natively relational.** It ships *as a database*. Loading it
relationally is lossless, which means the full source schema stays available
alongside the derived graph. For a research tool where the interesting questions
are the ones you did not anticipate, being able to drop into unmodified ChEMBL
matters more than it looks. A graph projection is inherently lossy — it encodes
the questions you thought of.

**The RDKit cartridge exists for Postgres.** Native chemical types, substructure
search, GiST-indexed fingerprints, in-database Tanimoto. Not adopted today
(0001 explains why the in-process index wins for full-database search) but it is
a real upgrade path that Neo4j structurally cannot offer.

Also worth noting: Neo4j **Enterprise**, as named in the blueprint, is
commercially licensed. Community would have been the honest comparison.

## Decision

**PostgreSQL 16 + Apache AGE**, with the relational tables as the system of
record and the AGE graph as a rebuildable projection.

## Consequences

- One system to run, back up, and reason about.
- Cypher is still available via AGE where it reads better than SQL.
- AGE is less mature than Neo4j — thinner tooling, rougher docs, and bulk
  loading needs its file loaders rather than `CREATE` statements. Accepted,
  because the graph is a projection that can be dropped and rebuilt at will.
- If this ever grows PPI networks and genuine variable-depth mechanism
  inference, revisit. That is the scenario where Neo4j would start earning its
  keep.
