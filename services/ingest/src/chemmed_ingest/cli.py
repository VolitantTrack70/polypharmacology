"""Command-line entry point for the ingestion pipeline.

    download -> parse -> fingerprint -> load -> index -> project-graph

Every stage is independently runnable and re-runnable. `run-all` chains them.
"""

from __future__ import annotations

import logging
import multiprocessing as mp
import tarfile
from pathlib import Path

import typer
from rich.console import Console
from rich.logging import RichHandler
from rich.table import Table

from chemmed_ingest.config import DATABASE_URL, FP, PATHS, REPO_ROOT, THRESHOLDS

app = typer.Typer(add_completion=False, help="Polypharmacology graph ingestion pipeline.")
console = Console()

CHEMBL_URL = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_{v}/chembl_{v}_sqlite.tar.gz"
REACTOME_BASE = "https://reactome.org/download/current"
REACTOME_FILES = [
    "ReactomePathways.txt",
    "ReactomePathwaysRelation.txt",
    "UniProt2Reactome_All_Levels.txt",
]


def _setup_logging(verbose: bool = True) -> None:
    logging.basicConfig(
        level=logging.INFO if verbose else logging.WARNING,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(console=console, rich_tracebacks=True, show_path=False)],
    )


def _chembl_sqlite_path(release: str) -> Path:
    return PATHS.raw / f"chembl_{release}" / f"chembl_{release}_sqlite" / f"chembl_{release}.db"


# ---------------------------------------------------------------------------


@app.command()
def download(
    source: str = typer.Option("all", help="chembl | reactome | all"),
    release: str = typer.Option("35", help="ChEMBL release number"),
) -> None:
    """Fetch source dumps into data/raw. Resumable; skips files already present."""
    _setup_logging()
    import requests
    from tqdm import tqdm

    PATHS.ensure()

    def fetch(url: str, dest: Path) -> None:
        if dest.exists() and dest.stat().st_size > 0:
            console.print(f"[dim]have[/dim] {dest.name}")
            return
        dest.parent.mkdir(parents=True, exist_ok=True)
        console.print(f"[cyan]GET[/cyan] {url}")
        with requests.get(url, stream=True, timeout=120) as r:
            r.raise_for_status()
            total = int(r.headers.get("content-length", 0))
            tmp = dest.with_suffix(dest.suffix + ".part")
            with open(tmp, "wb") as fh, tqdm(
                total=total, unit="B", unit_scale=True, desc=dest.name
            ) as bar:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    fh.write(chunk)
                    bar.update(len(chunk))
            tmp.rename(dest)

    if source in ("chembl", "all"):
        url = CHEMBL_URL.format(v=release)
        tarball = PATHS.raw / f"chembl_{release}_sqlite.tar.gz"
        fetch(url, tarball)
        target = PATHS.raw / f"chembl_{release}"
        if not target.exists():
            console.print(f"[cyan]extracting[/cyan] {tarball.name} (this takes a while)")
            with tarfile.open(tarball) as tf:
                tf.extractall(target, filter="data")

    if source in ("reactome", "all"):
        for name in REACTOME_FILES:
            fetch(f"{REACTOME_BASE}/{name}", PATHS.raw / name)

    console.print("[green]download complete[/green]")


@app.command()
def parse(
    release: str = typer.Option("35"),
    limit: int | None = typer.Option(None, help="Rows per table. For smoke tests."),
    species: str = typer.Option("Homo sapiens", help="Reactome species filter; '' for all"),
) -> None:
    """Source dumps -> normalised Parquet in data/processed."""
    _setup_logging()
    from chemmed_ingest.sources import chembl, reactome

    PATHS.ensure()
    counts: dict[str, int] = {}

    sqlite_path = _chembl_sqlite_path(release)
    if sqlite_path.exists():
        counts |= chembl.export_all(sqlite_path, PATHS.processed, limit=limit)
    else:
        console.print(f"[yellow]skipping ChEMBL: {sqlite_path} not found[/yellow]")

    if (PATHS.raw / "ReactomePathways.txt").exists():
        counts |= reactome.export_all(PATHS.raw, PATHS.processed, species=species or None)
    else:
        console.print("[yellow]skipping Reactome: files not found[/yellow]")

    _print_counts("Parsed", counts)


@app.command()
def fingerprint(
    workers: int = typer.Option(0, help="Process count; 0 = cpu_count() - 1"),
    limit: int | None = typer.Option(None),
) -> None:
    """Standardise + Morgan-fingerprint every structure, in a process pool."""
    _setup_logging()
    import polars as pl
    import pyarrow as pa
    import pyarrow.parquet as pq
    from tqdm import tqdm

    from chemmed_ingest.chem.fingerprint import smiles_to_record

    src = PATHS.processed / "compound.parquet"
    if not src.exists():
        raise typer.BadParameter(f"{src} not found -- run `parse` first")

    df = pl.read_parquet(src, columns=["chembl_id", "canonical_smiles"])
    if limit:
        df = df.head(limit)
    pairs = list(zip(df["chembl_id"].to_list(), df["canonical_smiles"].to_list(), strict=True))

    n_workers = workers or max(1, (mp.cpu_count() or 2) - 1)
    console.print(f"fingerprinting {len(pairs):,} structures across {n_workers} workers")

    schema = pa.schema([
        ("chembl_id", pa.string()),
        ("radius", pa.int16()),
        ("n_bits", pa.int16()),
        ("fp", pa.binary()),
        ("popcount", pa.int16()),
    ])
    out = PATHS.processed / "compound_fingerprint.parquet"

    failures = 0
    written = 0
    buffer: list[dict] = []

    with mp.Pool(n_workers) as pool, pq.ParquetWriter(out, schema, compression="zstd") as writer:
        results = pool.starmap(smiles_to_record, pairs, chunksize=2000)
        for rec in tqdm(results, total=len(pairs), desc="fingerprint"):
            if rec is None:
                failures += 1
                continue
            buffer.append(rec)
            if len(buffer) >= 100_000:
                writer.write_batch(_fp_batch(buffer, schema))
                written += len(buffer)
                buffer.clear()
        if buffer:
            writer.write_batch(_fp_batch(buffer, schema))
            written += len(buffer)

    console.print(
        f"[green]{written:,} fingerprints[/green] "
        f"({failures:,} structures could not be standardised)"
    )


def _fp_batch(records: list[dict], schema):
    import pyarrow as pa

    return pa.RecordBatch.from_arrays(
        [
            pa.array([r["chembl_id"] for r in records], type=pa.string()),
            pa.array([r["radius"] for r in records], type=pa.int16()),
            pa.array([r["n_bits"] for r in records], type=pa.int16()),
            pa.array([r["fp"] for r in records], type=pa.binary()),
            pa.array([r["popcount"] for r in records], type=pa.int16()),
        ],
        schema=schema,
    )


@app.command()
def migrate() -> None:
    """Apply db/migrations in order. Idempotent."""
    _setup_logging()
    import psycopg

    from chemmed_ingest.load.postgres import apply_migrations

    with psycopg.connect(DATABASE_URL) as conn:
        applied, skipped = apply_migrations(conn, REPO_ROOT / "db" / "migrations")

    console.print(f"[green]applied {len(applied)}[/green]: {', '.join(applied) or '-'}")
    for name, reason in skipped:
        console.print(f"[yellow]skipped[/yellow] {name} ({reason})")


@app.command(name="load")
def load_cmd(release: str = typer.Option("35")) -> None:
    """COPY the Parquet into Postgres, then refresh derived views."""
    _setup_logging()
    import psycopg

    from chemmed_ingest.load.postgres import load_table, record_release, refresh_derived

    # Order matters: parents before children, for the FK guards.
    order = [
        "compound", "target", "protein", "target_component",
        "pathway", "pathway_hierarchy", "protein_pathway",
        "activity", "compound_fingerprint",
    ]

    counts: dict[str, int] = {}
    with psycopg.connect(DATABASE_URL) as conn:
        rid = record_release(conn, "chembl", release, CHEMBL_URL.format(v=release))
        for table in order:
            path = PATHS.processed / f"{table}.parquet"
            counts[table] = load_table(conn, table, path, release_id=rid)
        refresh_derived(conn)

    _print_counts("Loaded", counts)


@app.command()
def index(
    source: str = typer.Option(
        "parquet",
        help="parquet | db. Parquet needs no database and is the default.",
    ),
) -> None:
    """Build the packed fingerprint matrix used for query-time similarity.

    Reads from the Parquet produced by `fingerprint` by default -- the index is
    derived data and there is no reason to route it through Postgres. `--source
    db` exists for rebuilding after a load that came from somewhere else.
    """
    _setup_logging()
    from chemmed_ingest.chem.similarity import FingerprintIndex

    if source == "parquet":
        import polars as pl

        src = PATHS.processed / "compound_fingerprint.parquet"
        if not src.exists():
            raise typer.BadParameter(f"{src} not found -- run `fingerprint` first")
        df = pl.read_parquet(src, columns=["chembl_id", "fp"])
        records = list(zip(df["chembl_id"].to_list(), df["fp"].to_list(), strict=True))
    elif source == "db":
        import psycopg

        with psycopg.connect(DATABASE_URL) as conn, conn.cursor() as cur:
            cur.execute("SELECT chembl_id, fp FROM chem.compound_fingerprint")
            records = [(r[0], bytes(r[1])) for r in cur]
    else:
        raise typer.BadParameter(f"unknown source {source!r}; use 'parquet' or 'db'")

    if not records:
        raise typer.BadParameter(f"no fingerprints found in {source}")

    idx = FingerprintIndex.from_records(records)
    path = PATHS.fpsim_index.with_suffix(".npz")
    idx.save(path)
    console.print(f"[green]{idx}[/green] -> {path}")


@app.command()
def search(
    smiles: str = typer.Argument(..., help="Query structure as SMILES"),
    threshold: float = typer.Option(THRESHOLDS.tanimoto),
    limit: int = typer.Option(20),
) -> None:
    """Similarity search from the shell. Useful for sanity-checking the index."""
    _setup_logging(verbose=False)
    from chemmed_ingest.chem.fingerprint import fingerprint_bytes, standardize
    from chemmed_ingest.chem.similarity import FingerprintIndex

    path = PATHS.fpsim_index.with_suffix(".npz")
    if not path.exists():
        raise typer.BadParameter(f"{path} not found -- run `index` first")

    std = standardize(smiles)
    console.print(f"[dim]canonical:[/dim] {std.canonical_smiles}")

    idx = FingerprintIndex.load(path)
    hits = idx.search(fingerprint_bytes(std.mol), threshold=threshold, limit=limit)

    table = Table(title=f"{len(hits)} hits at Tanimoto >= {threshold}")
    table.add_column("ChEMBL ID")
    table.add_column("Tanimoto", justify="right")
    for hit in hits:
        table.add_row(hit.chembl_id, f"{hit.tanimoto:.3f}")
    console.print(table)


@app.command(name="run-all")
def run_all(
    release: str = typer.Option("35"),
    limit: int | None = typer.Option(None, help="Row cap per table. Use for a fast smoke run."),
    skip_download: bool = typer.Option(False),
) -> None:
    """Run the full pipeline end to end."""
    if not skip_download:
        download(source="all", release=release)
    parse(release=release, limit=limit, species="Homo sapiens")
    fingerprint(workers=0, limit=limit)
    migrate()
    load_cmd(release=release)
    index()
    console.print("\n[bold green]pipeline complete[/bold green]")
    console.print(f"config: {FP.signature}, default cutoff {THRESHOLDS.tanimoto}")


def _print_counts(title: str, counts: dict[str, int]) -> None:
    table = Table(title=title)
    table.add_column("Entity")
    table.add_column("Rows", justify="right")
    for k, v in counts.items():
        table.add_row(k, f"{v:,}")
    console.print(table)


if __name__ == "__main__":
    app()
