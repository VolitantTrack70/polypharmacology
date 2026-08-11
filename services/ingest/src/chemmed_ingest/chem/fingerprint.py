"""Molecule standardisation and Morgan fingerprint generation.

Standardisation matters more than it looks. ChEMBL contains salts, mixtures,
and charge variants of the same parent structure. Fingerprinting a
hydrochloride salt and its free base produces different vectors, so a naive
pipeline reports two "different" compounds that are pharmacologically one.
Every SMILES -- from the dumps and from user input -- goes through the same
`standardize()` path so query and index agree.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize

from chemmed_ingest.config import FP

# RDKit is extremely chatty about malformed input; we count failures ourselves.
RDLogger.DisableLog("rdApp.*")


@lru_cache(maxsize=1)
def _generator() -> rdFingerprintGenerator.FingerprintGenerator64:
    """Morgan generator. Cached because construction is not free and this is
    called per-molecule in tight loops."""
    return rdFingerprintGenerator.GetMorganGenerator(radius=FP.radius, fpSize=FP.n_bits)


@lru_cache(maxsize=1)
def _largest_fragment_chooser() -> rdMolStandardize.LargestFragmentChooser:
    return rdMolStandardize.LargestFragmentChooser()


@lru_cache(maxsize=1)
def _uncharger() -> rdMolStandardize.Uncharger:
    return rdMolStandardize.Uncharger()


class StandardizationError(ValueError):
    """Raised when a SMILES cannot be parsed or standardised."""


@dataclass(frozen=True)
class StandardizedMolecule:
    canonical_smiles: str
    inchikey: str
    mol: Chem.Mol

    @property
    def parent_inchikey(self) -> str:
        """First InChIKey block -- the skeleton hash, ignoring stereo and
        protonation. Useful for grouping stereoisomers."""
        return self.inchikey.split("-")[0]


def standardize(smiles: str) -> StandardizedMolecule:
    """Parse, desalt, neutralise, and canonicalise a SMILES string.

    Raises StandardizationError on unparseable input rather than returning
    None, so callers can't accidentally propagate a null into the index.
    """
    if not smiles or not smiles.strip():
        raise StandardizationError("empty SMILES")

    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        raise StandardizationError(f"RDKit could not parse SMILES: {smiles[:120]!r}")

    try:
        # Order matters: strip counter-ions first, then neutralise what's left.
        mol = _largest_fragment_chooser().choose(mol)
        mol = _uncharger().uncharge(mol)
        Chem.SanitizeMol(mol)
    except Exception as exc:  # RDKit raises a variety of C++-backed types
        raise StandardizationError(f"standardisation failed: {exc}") from exc

    canonical = Chem.MolToSmiles(mol, canonical=True)
    inchikey = Chem.MolToInchiKey(mol)
    if not inchikey:
        raise StandardizationError("InChIKey generation returned empty")

    return StandardizedMolecule(canonical_smiles=canonical, inchikey=inchikey, mol=mol)


def fingerprint(mol: Chem.Mol) -> DataStructs.ExplicitBitVect:
    """Morgan fingerprint as an RDKit bit vector."""
    return _generator().GetFingerprint(mol)


def fingerprint_bytes(mol: Chem.Mol) -> bytes:
    """Fingerprint packed to raw bytes for storage in Postgres BYTEA.

    Bit ordering is whatever RDKit's binary text uses; it is never interpreted
    outside this module, and Tanimoto is order-invariant so long as index and
    query use the same packing.
    """
    return DataStructs.BitVectToBinaryText(fingerprint(mol))


def fingerprint_array(mol: Chem.Mol) -> np.ndarray:
    """Fingerprint as a uint8 array of length n_bits/8."""
    return np.frombuffer(fingerprint_bytes(mol), dtype=np.uint8)


def bytes_to_bitvect(raw: bytes) -> DataStructs.ExplicitBitVect:
    """Inverse of `fingerprint_bytes`, for round-tripping out of Postgres."""
    return DataStructs.CreateFromBinaryText(raw)


def popcount(raw: bytes) -> int:
    """Number of set bits. Stored alongside each fingerprint so the search can
    bound Tanimoto cheaply: Ts(A,B) <= min(|A|,|B|) / max(|A|,|B|)."""
    return int(np.unpackbits(np.frombuffer(raw, dtype=np.uint8)).sum())


def smiles_to_record(chembl_id: str, smiles: str) -> dict[str, object] | None:
    """Full per-compound pipeline used by the parallel fingerprint workers.

    Returns None (rather than raising) on bad input, because a handful of
    unparseable structures in a 2.4M-row dump should not kill the job. The
    caller is responsible for counting and reporting the failures.
    """
    try:
        std = standardize(smiles)
    except StandardizationError:
        return None

    raw = fingerprint_bytes(std.mol)
    return {
        "chembl_id": chembl_id,
        "canonical_smiles": std.canonical_smiles,
        "standard_inchi_key": std.inchikey,
        "fp": raw,
        "popcount": popcount(raw),
        "radius": FP.radius,
        "n_bits": FP.n_bits,
    }
