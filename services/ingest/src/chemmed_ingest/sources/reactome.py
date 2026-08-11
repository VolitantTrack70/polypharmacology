"""Reactome flat files -> normalised Parquet.

File selection matters here. The original blueprint named `Ensemble2Reactome.txt`,
which is (a) spelled `Ensembl2Reactome.txt` and (b) keyed on Ensembl *gene* IDs.
Since targets are identified by UniProt accession everywhere else in this system,
routing through Ensembl would mean an extra lossy ID-mapping hop.

The correct files:

  UniProt2Reactome_All_Levels.txt   accession -> pathway, directly.
                                    "_All_Levels" includes ancestor pathways, not
                                    just the leaf a protein is annotated to.
  ReactomePathwaysRelation.txt      parent -> child hierarchy.
  ReactomePathways.txt              pathway id -> name, species.

All three are headerless TSV.
"""

from __future__ import annotations

import logging
from pathlib import Path

import polars as pl

log = logging.getLogger(__name__)

FILES = {
    "pathways": "ReactomePathways.txt",
    "relations": "ReactomePathwaysRelation.txt",
    "uniprot_map": "UniProt2Reactome_All_Levels.txt",
}

BASE_URL = "https://reactome.org/download/current"

DEFAULT_SPECIES = "Homo sapiens"


def _require(raw_dir: Path, key: str) -> Path:
    path = raw_dir / FILES[key]
    if not path.exists():
        raise FileNotFoundError(
            f"missing {path.name}. Run `chemmed-ingest download --source reactome` "
            f"(fetches from {BASE_URL}/{FILES[key]})."
        )
    return path


def parse_pathways(raw_dir: Path, species: str | None = DEFAULT_SPECIES) -> pl.DataFrame:
    """ReactomePathways.txt -> (reactome_id, pathway_name, species).

    `is_top_level` is derived rather than read: a pathway is top-level exactly
    when it never appears as a child in the relation file.
    """
    pathways = pl.read_csv(
        _require(raw_dir, "pathways"),
        separator="\t",
        has_header=False,
        new_columns=["reactome_id", "pathway_name", "species"],
        quote_char=None,  # pathway names contain unbalanced quotes
    )

    relations = pl.read_csv(
        _require(raw_dir, "relations"),
        separator="\t",
        has_header=False,
        new_columns=["parent_reactome_id", "child_reactome_id"],
        quote_char=None,
    )

    if species:
        pathways = pathways.filter(pl.col("species") == species)

    children = relations.select("child_reactome_id").unique()
    pathways = pathways.with_columns(
        pl.col("reactome_id")
        .is_in(children["child_reactome_id"])
        .not_()
        .alias("is_top_level"),
        pl.lit(None, dtype=pl.Utf8).alias("biological_domain"),  # filled by migration 003
    )

    log.info(
        "reactome: %d pathways (%d top-level)",
        pathways.height,
        int(pathways["is_top_level"].sum()),
    )
    return pathways


def parse_hierarchy(raw_dir: Path, valid_ids: pl.Series | None = None) -> pl.DataFrame:
    """ReactomePathwaysRelation.txt -> (parent_reactome_id, child_reactome_id).

    Filtered to `valid_ids` when given, because the file spans every species and
    the FK to `pathway` would otherwise fail on non-human rows.
    """
    relations = pl.read_csv(
        _require(raw_dir, "relations"),
        separator="\t",
        has_header=False,
        new_columns=["parent_reactome_id", "child_reactome_id"],
        quote_char=None,
    )

    if valid_ids is not None:
        keep = set(valid_ids.to_list())
        relations = relations.filter(
            pl.col("parent_reactome_id").is_in(keep) & pl.col("child_reactome_id").is_in(keep)
        )

    relations = relations.unique()
    log.info("reactome: %d hierarchy edges", relations.height)
    return relations


def parse_protein_pathway(
    raw_dir: Path,
    valid_ids: pl.Series | None = None,
    species: str | None = DEFAULT_SPECIES,
) -> pl.DataFrame:
    """UniProt2Reactome_All_Levels.txt -> (uniprot_accession, reactome_id, evidence_code).

    Columns are positional and headerless:
        0 accession  1 reactome_id  2 url  3 pathway_name  4 evidence_code  5 species

    Isoform suffixes (`P12345-2`) are stripped to the canonical accession so these
    join against ChEMBL's `component_sequences.accession`, which carries none.
    """
    df = pl.read_csv(
        _require(raw_dir, "uniprot_map"),
        separator="\t",
        has_header=False,
        new_columns=[
            "uniprot_accession",
            "reactome_id",
            "url",
            "pathway_name",
            "evidence_code",
            "species",
        ],
        quote_char=None,
    )

    if species:
        df = df.filter(pl.col("species") == species)

    df = df.with_columns(
        pl.col("uniprot_accession").str.split("-").list.first().alias("uniprot_accession")
    )

    if valid_ids is not None:
        df = df.filter(pl.col("reactome_id").is_in(set(valid_ids.to_list())))

    df = df.select("uniprot_accession", "reactome_id", "evidence_code").unique(
        subset=["uniprot_accession", "reactome_id"]
    )

    log.info("reactome: %d protein-pathway links", df.height)
    return df


def export_all(
    raw_dir: Path,
    out_dir: Path,
    species: str | None = DEFAULT_SPECIES,
) -> dict[str, int]:
    """Parse all three files and write Parquet, with referential integrity
    enforced against the pathway set so the Postgres FKs hold."""
    out_dir.mkdir(parents=True, exist_ok=True)

    pathways = parse_pathways(raw_dir, species=species)
    ids = pathways["reactome_id"]

    hierarchy = parse_hierarchy(raw_dir, valid_ids=ids)
    links = parse_protein_pathway(raw_dir, valid_ids=ids, species=species)

    pathways.write_parquet(out_dir / "pathway.parquet", compression="zstd")
    hierarchy.write_parquet(out_dir / "pathway_hierarchy.parquet", compression="zstd")
    links.write_parquet(out_dir / "protein_pathway.parquet", compression="zstd")

    return {
        "pathway": pathways.height,
        "pathway_hierarchy": hierarchy.height,
        "protein_pathway": links.height,
    }
