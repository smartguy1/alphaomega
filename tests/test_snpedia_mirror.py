"""
AlphaOmega — SNPedia Mirror Tests

Tests the SNPedia mirror using a small in-memory SQLite fixture.
No network access needed — all tests use local mock data.
"""

import os
import sqlite3
import tempfile

import pytest

from src.snpedia_mirror import SNPediaMirror, SNPediaEntry, SNPediaGenotype


# ── Fixtures ─────────────────────────────────────────────────────────────────

@pytest.fixture
def mock_db_dir(tmp_path):
    """Create a temp directory with a mock SNPedia SQLite database."""
    db_path = tmp_path / "snpedia.sqlite"

    conn = sqlite3.connect(str(db_path))

    # Create tables matching a realistic SNPedia dump schema
    conn.execute("""
        CREATE TABLE snps (
            rsid TEXT PRIMARY KEY,
            chromosome TEXT,
            position INTEGER,
            gene TEXT,
            summary TEXT
        )
    """)

    conn.execute("""
        CREATE TABLE genotypes (
            rsid TEXT,
            genotype TEXT,
            magnitude REAL,
            repute TEXT,
            summary TEXT,
            details TEXT
        )
    """)

    conn.execute("""
        CREATE TABLE categories (
            rsid TEXT,
            category TEXT
        )
    """)

    # Insert test data
    conn.executemany(
        "INSERT INTO snps VALUES (?, ?, ?, ?, ?)",
        [
            ("rs1815739", "11", 66560624, "ACTN3",
             "Muscle performance - sprinting vs endurance"),
            ("rs4988235", "2", 136608646, "MCM6",
             "Lactose intolerance / persistence"),
            ("rs1800497", "11", 113400106, "ANKK1",
             "DRD2/ANKK1 TaqIA - dopamine receptor"),
            ("rs762551", "15", 75041917, "CYP1A2",
             "Caffeine metabolism speed"),
            ("rs53576", "3", 8804371, "OXTR",
             "Oxytocin receptor - empathy and social behavior"),
        ],
    )

    conn.executemany(
        "INSERT INTO genotypes VALUES (?, ?, ?, ?, ?, ?)",
        [
            ("rs1815739", "(C;C)", 3.0, "good",
             "Optimal sprint/power performance", "Both copies of functional ACTN3"),
            ("rs1815739", "(C;T)", 2.0, "neutral",
             "Carrier - mixed muscle type", "One functional copy"),
            ("rs1815739", "(T;T)", 2.5, "bad",
             "Endurance-oriented muscle type", "No functional ACTN3"),
            ("rs4988235", "(G;G)", 0.0, "neutral",
             "Likely lactose intolerant", ""),
            ("rs4988235", "(A;G)", 2.5, "good",
             "Likely lactose tolerant", ""),
            ("rs4988235", "(A;A)", 2.5, "good",
             "Likely lactose tolerant", ""),
            ("rs762551", "(A;A)", 2.0, "good",
             "Fast caffeine metabolizer", ""),
            ("rs762551", "(A;C)", 1.5, "neutral",
             "Moderate caffeine metabolizer", ""),
            ("rs762551", "(C;C)", 2.5, "bad",
             "Slow caffeine metabolizer", "Increased cardiovascular risk with coffee"),
        ],
    )

    conn.executemany(
        "INSERT INTO categories VALUES (?, ?)",
        [
            ("rs1815739", "muscle"),
            ("rs1815739", "athletics"),
            ("rs4988235", "metabolism"),
            ("rs4988235", "diet"),
            ("rs762551", "metabolism"),
            ("rs762551", "caffeine"),
            ("rs53576", "behavior"),
            ("rs1800497", "behavior"),
        ],
    )

    conn.commit()
    conn.close()

    return str(tmp_path)


# ── Tests ────────────────────────────────────────────────────────────────────

class TestSNPediaMirror:

    def test_is_available(self, mock_db_dir):
        mirror = SNPediaMirror(db_dir=mock_db_dir)
        assert mirror.is_available

    def test_not_available(self, tmp_path):
        mirror = SNPediaMirror(db_dir=str(tmp_path / "nonexistent"))
        assert not mirror.is_available

    def test_get_schema(self, mock_db_dir):
        with SNPediaMirror(db_dir=mock_db_dir) as mirror:
            schema = mirror.get_schema()

            assert "snps" in schema
            assert "genotypes" in schema
            assert "categories" in schema
            assert "rsid" in schema["snps"]

    def test_get_stats(self, mock_db_dir):
        with SNPediaMirror(db_dir=mock_db_dir) as mirror:
            stats = mirror.get_stats()

            assert stats["tables"]["snps"] == 5
            assert stats["tables"]["genotypes"] == 9
            assert stats["tables"]["categories"] == 8
            assert stats["db_size_mb"] > 0

    def test_lookup_existing_snp(self, mock_db_dir):
        with SNPediaMirror(db_dir=mock_db_dir) as mirror:
            entry = mirror.lookup("rs1815739")

            assert entry is not None
            assert entry.rsid == "rs1815739"
            assert len(entry.genotypes) == 3

    def test_lookup_case_insensitive(self, mock_db_dir):
        with SNPediaMirror(db_dir=mock_db_dir) as mirror:
            entry = mirror.lookup("RS1815739")
            assert entry is not None

    def test_lookup_nonexistent(self, mock_db_dir):
        with SNPediaMirror(db_dir=mock_db_dir) as mirror:
            entry = mirror.lookup("rs99999999")
            assert entry is None

    def test_lookup_genotype_specific(self, mock_db_dir):
        with SNPediaMirror(db_dir=mock_db_dir) as mirror:
            gt = mirror.lookup_genotype("rs1815739", "CC")

            assert gt is not None
            assert gt.magnitude == 3.0
            assert gt.repute == "good"
            assert "sprint" in gt.summary.lower()

    def test_lookup_genotype_reversed_alleles(self, mock_db_dir):
        """Allele order shouldn't matter: TC == CT."""
        with SNPediaMirror(db_dir=mock_db_dir) as mirror:
            gt = mirror.lookup_genotype("rs1815739", "TC")

            assert gt is not None
            assert gt.genotype == "(C;T)"

    def test_lookup_genotype_nonexistent(self, mock_db_dir):
        with SNPediaMirror(db_dir=mock_db_dir) as mirror:
            gt = mirror.lookup_genotype("rs1815739", "GG")
            assert gt is None

    def test_export_rsid_list(self, mock_db_dir, tmp_path):
        output = str(tmp_path / "rsids.txt")

        with SNPediaMirror(db_dir=mock_db_dir) as mirror:
            count = mirror.export_rsid_list(output)

            assert count == 5
            assert os.path.exists(output)

            with open(output) as f:
                lines = f.read().strip().splitlines()
                assert len(lines) == 5
                assert "rs1815739" in lines

    def test_search_by_category(self, mock_db_dir):
        with SNPediaMirror(db_dir=mock_db_dir) as mirror:
            results = mirror.search_by_category("metabolism")

            assert len(results) == 2  # rs4988235, rs762551

    def test_search_by_category_no_results(self, mock_db_dir):
        with SNPediaMirror(db_dir=mock_db_dir) as mirror:
            results = mirror.search_by_category("nonexistent_category")
            assert len(results) == 0

    def test_context_manager(self, mock_db_dir):
        """Verify context manager properly opens and closes."""
        with SNPediaMirror(db_dir=mock_db_dir) as mirror:
            assert mirror.is_available
            entry = mirror.lookup("rs762551")
            assert entry is not None

        # After exit, connection should be closed
        assert mirror._conn is None


class TestSNPediaGenotype:

    def test_caffeine_fast(self, mock_db_dir):
        with SNPediaMirror(db_dir=mock_db_dir) as mirror:
            gt = mirror.lookup_genotype("rs762551", "AA")

            assert gt is not None
            assert gt.repute == "good"
            assert "fast" in gt.summary.lower()

    def test_caffeine_slow(self, mock_db_dir):
        with SNPediaMirror(db_dir=mock_db_dir) as mirror:
            gt = mirror.lookup_genotype("rs762551", "CC")

            assert gt is not None
            assert gt.repute == "bad"
            assert gt.magnitude == 2.5

    def test_lactose_intolerant(self, mock_db_dir):
        with SNPediaMirror(db_dir=mock_db_dir) as mirror:
            gt = mirror.lookup_genotype("rs4988235", "GG")

            assert gt is not None
            assert "intolerant" in gt.summary.lower()


class TestSNPediaFindColumn:
    """Test the static column-finding helper."""

    def test_exact_match(self):
        cols = ["rsid", "chromosome", "position"]
        assert SNPediaMirror._find_column(cols, ["rsid"]) == "rsid"

    def test_case_insensitive(self):
        cols = ["RSID", "Chromosome", "Position"]
        assert SNPediaMirror._find_column(cols, ["rsid"]) == "RSID"

    def test_fallback_candidates(self):
        cols = ["name", "chr", "pos"]
        assert SNPediaMirror._find_column(cols, ["rsid", "name"]) == "name"

    def test_no_match(self):
        cols = ["foo", "bar"]
        assert SNPediaMirror._find_column(cols, ["rsid", "name"]) is None
