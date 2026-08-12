"""Central configuration, loaded from environment / .env.

Fingerprint parameters live here rather than being passed around, because a
mismatch between the radius/n_bits used at index-build time and at query time
produces silently wrong similarity scores rather than an error.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

from dotenv import load_dotenv

# Walk up to the repo root (.../chem-med) and load .env if present.
REPO_ROOT = Path(__file__).resolve().parents[4]
load_dotenv(REPO_ROOT / ".env")


def _env(key: str, default: str) -> str:
    return os.environ.get(key, default)


@dataclass(frozen=True)
class FingerprintConfig:
    """Morgan (ECFP-like) fingerprint parameters.

    radius=2 / n_bits=2048 corresponds to ECFP4 and is the field default.
    Changing either invalidates every stored fingerprint and the search index.
    """

    radius: int = int(_env("FP_RADIUS", "2"))
    n_bits: int = int(_env("FP_NBITS", "2048"))

    @property
    def n_bytes(self) -> int:
        return self.n_bits // 8

    @property
    def signature(self) -> str:
        """Stamped into index files so a mismatched index fails loudly."""
        return f"morgan-r{self.radius}-{self.n_bits}"


@dataclass(frozen=True)
class Thresholds:
    """Default cutoffs.

    tanimoto=0.40 is deliberate and is NOT the 0.85 that intuition suggests.
    Morgan/ECFP4 Tanimoto is a far harsher metric than path- or MACCS-based
    fingerprints: 0.85 means near-identical molecules (a methyl group apart)
    and returns essentially nothing. The range over which compounds
    meaningfully share targets is roughly 0.35-0.55.
    See docs/decisions/0002-similarity-thresholds.md.
    """

    tanimoto: float = float(_env("DEFAULT_TANIMOTO_CUTOFF", "0.40"))
    pchembl: float = float(_env("DEFAULT_PCHEMBL_CUTOFF", "6.0"))
    # ChEMBL assay-to-target confidence, 0-9. The chem.binds_to view applies a
    # hard floor of 7 (see migration 002); this is the queryable floor on top,
    # e.g. 9 for direct single-protein assays only.
    min_confidence: int = int(_env("MIN_CONFIDENCE_SCORE", "7"))


def _path_env(key: str, default: Path) -> Path:
    """Resolve a path from the environment, anchoring relatives to the repo root.

    A bare `./data/processed` in .env would otherwise resolve against whatever
    directory the command happened to be run from -- so `parse` (run from the
    repo root) and `load` (run from services/ingest) would silently disagree
    about where the Parquet lives.
    """
    raw = os.environ.get(key)
    if not raw:
        return default
    path = Path(raw)
    return path if path.is_absolute() else (REPO_ROOT / path).resolve()


@dataclass(frozen=True)
class Paths:
    raw: Path = _path_env("DATA_RAW_DIR", REPO_ROOT / "data" / "raw")
    processed: Path = _path_env("DATA_PROCESSED_DIR", REPO_ROOT / "data" / "processed")
    # Packed fingerprint matrix, written as .npz by FingerprintIndex.save().
    fingerprint_index: Path = _path_env(
        "FINGERPRINT_INDEX_PATH", REPO_ROOT / "data" / "processed" / "chembl_morgan_2048.npz"
    )

    def ensure(self) -> None:
        self.raw.mkdir(parents=True, exist_ok=True)
        self.processed.mkdir(parents=True, exist_ok=True)


@dataclass(frozen=True)
class KafkaConfig:
    """Not yet consumed by any code path.

    docker-compose provisions a broker and these topic names are the intended
    contract, but fingerprinting currently uses a local process pool -- the
    right tool on one machine. Kept so the topics are defined in one place if
    the work is ever distributed across hosts.
    """

    bootstrap_servers: str = _env("KAFKA_BOOTSTRAP_SERVERS", "localhost:9092")
    topic_compounds: str = _env("KAFKA_TOPIC_COMPOUNDS", "chemmed.compounds.raw")
    topic_fingerprints: str = _env("KAFKA_TOPIC_FINGERPRINTS", "chemmed.compounds.fingerprinted")
    group_id: str = _env("KAFKA_INGEST_GROUP", "chemmed-fingerprint-workers")


DATABASE_URL = _env("DATABASE_URL", "postgres://chemmed:change_me_locally@localhost:5432/chemmed")

FP = FingerprintConfig()
THRESHOLDS = Thresholds()
PATHS = Paths()
KAFKA = KafkaConfig()
