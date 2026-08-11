"""Tests for checksum manifest parsing and the ChEMBL database finder.

Both guard expensive mistakes: a corrupted 5 GB download that only reveals
itself as a bizarre SQLite error during parse, and a path assumption that
breaks immediately after a 25 GB extraction.
"""

from __future__ import annotations

import hashlib

import pytest

from chemmed_ingest.cli import _sha256, find_checksum

# Verbatim from https://ftp.ebi.ac.uk/.../chembl_35/checksums.txt
REAL_MANIFEST = (
    "381842aae142bf68b035f8d1aba94dfb715ad5d979640f16186b0f9d44fa5320"
    "\tchembl_35_sqlite.tar.gz\n"
    "a1b2c3d4e5f60718293a4b5c6d7e8f90a1b2c3d4e5f60718293a4b5c6d7e8f90"
    "\tchembl_35_mysql.tar.gz\n"
)


class TestFindChecksum:
    def test_parses_the_real_chembl_format(self):
        assert find_checksum(REAL_MANIFEST, "chembl_35_sqlite.tar.gz") == (
            "381842aae142bf68b035f8d1aba94dfb715ad5d979640f16186b0f9d44fa5320"
        )

    def test_picks_the_right_file_among_several(self):
        assert find_checksum(REAL_MANIFEST, "chembl_35_mysql.tar.gz").startswith("a1b2c3")

    def test_missing_file_returns_none(self):
        assert find_checksum(REAL_MANIFEST, "chembl_99_sqlite.tar.gz") is None

    def test_space_separated_coreutils_style(self):
        text = "381842aae142bf68b035f8d1aba94dfb715ad5d979640f16186b0f9d44fa5320  file.tar.gz"
        assert find_checksum(text, "file.tar.gz").startswith("381842")

    def test_binary_mode_asterisk_is_stripped(self):
        text = "381842aae142bf68b035f8d1aba94dfb715ad5d979640f16186b0f9d44fa5320 *file.tar.gz"
        assert find_checksum(text, "file.tar.gz").startswith("381842")

    def test_directory_prefix_is_ignored(self):
        text = (
            "381842aae142bf68b035f8d1aba94dfb715ad5d979640f16186b0f9d44fa5320"
            "  ./dist/file.tar.gz"
        )
        assert find_checksum(text, "file.tar.gz").startswith("381842")

    def test_md5_entry_is_not_mistaken_for_sha256(self):
        """A 32-char MD5 must not be returned and then compared against a
        SHA-256, which would report every download as corrupt."""
        text = "d41d8cd98f00b204e9800998ecf8427e  file.tar.gz"
        assert find_checksum(text, "file.tar.gz") is None

    def test_non_hex_digest_rejected(self):
        text = "z" * 64 + "  file.tar.gz"
        assert find_checksum(text, "file.tar.gz") is None

    def test_partial_name_does_not_match(self):
        """'chembl_35_sqlite.tar.gz' must not match a request for 'sqlite.tar.gz'."""
        assert find_checksum(REAL_MANIFEST, "sqlite.tar.gz") is None

    def test_empty_and_malformed_lines_are_skipped(self):
        text = "\n\n# comment\ngarbage\n" + REAL_MANIFEST
        assert find_checksum(text, "chembl_35_sqlite.tar.gz") is not None


class TestSha256:
    def test_matches_hashlib(self, tmp_path):
        payload = b"polypharmacology" * 5000
        path = tmp_path / "blob.bin"
        path.write_bytes(payload)
        assert _sha256(path) == hashlib.sha256(payload).hexdigest()

    def test_handles_a_file_larger_than_one_chunk(self, tmp_path):
        payload = bytes(range(256)) * 40_000  # ~10 MB
        path = tmp_path / "big.bin"
        path.write_bytes(payload)
        assert _sha256(path, chunk=4096) == hashlib.sha256(payload).hexdigest()

    def test_empty_file(self, tmp_path):
        path = tmp_path / "empty.bin"
        path.write_bytes(b"")
        assert _sha256(path) == hashlib.sha256(b"").hexdigest()


class TestChemblDatabaseFinder:
    """The finder searches rather than assumes, because the tarball supplies
    its own directory level and the layout has varied across releases."""

    @pytest.fixture
    def raw_dir(self, tmp_path, monkeypatch):
        # Paths is a frozen dataclass, so swap the whole object rather than
        # trying to mutate a field on it.
        from dataclasses import replace

        from chemmed_ingest import cli

        monkeypatch.setattr(cli, "PATHS", replace(cli.PATHS, raw=tmp_path))
        return tmp_path

    def _make_db(self, path, size=1000):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"\x00" * size)
        return path

    def test_finds_the_canonical_layout(self, raw_dir):
        from chemmed_ingest.cli import _chembl_sqlite_path

        want = self._make_db(raw_dir / "chembl_35" / "chembl_35_sqlite" / "chembl_35.db")
        assert _chembl_sqlite_path("35") == want

    def test_finds_the_doubled_layout(self, raw_dir):
        """What you get by extracting the archive into a directory you also
        named chembl_35."""
        from chemmed_ingest.cli import _chembl_sqlite_path

        want = self._make_db(
            raw_dir / "chembl_35" / "chembl_35" / "chembl_35_sqlite" / "chembl_35.db"
        )
        assert _chembl_sqlite_path("35") == want

    def test_finds_an_unexpected_nesting_via_glob(self, raw_dir):
        from chemmed_ingest.cli import _chembl_sqlite_path

        want = self._make_db(raw_dir / "a" / "b" / "c" / "chembl_35.db")
        assert _chembl_sqlite_path("35") == want

    def test_prefers_the_largest_match(self, raw_dir):
        """A partial or stray file must not win over the real database."""
        from chemmed_ingest.cli import _chembl_sqlite_path

        self._make_db(raw_dir / "stray" / "chembl_35.db", size=10)
        want = self._make_db(raw_dir / "real" / "chembl_35.db", size=100_000)
        assert _chembl_sqlite_path("35") == want

    def test_missing_release_returns_none(self, raw_dir):
        from chemmed_ingest.cli import _chembl_sqlite_path

        assert _chembl_sqlite_path("99") is None

    def test_does_not_match_a_different_release(self, raw_dir):
        from chemmed_ingest.cli import _chembl_sqlite_path

        self._make_db(raw_dir / "chembl_34" / "chembl_34_sqlite" / "chembl_34.db")
        assert _chembl_sqlite_path("35") is None
