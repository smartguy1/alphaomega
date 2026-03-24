"""
AlphaOmega — ClinVar Mirror Tests

Tests the ClinVar mirror using a mock SQLite database.
No network access needed.
"""

import sqlite3
import pytest

from src.clinvar_mirror import ClinVarMirror, ClinVarEntry, REVIEW_STARS


# ── Fixtures ─────────────────────────────────────────────────────────────────

@pytest.fixture
def mock_clinvar_dir(tmp_path):
    """Create a temp directory with a mock ClinVar SQLite database."""
    db_path = tmp_path / "clinvar.sqlite"
    conn = sqlite3.connect(str(db_path))

    conn.execute("""
        CREATE TABLE variants (
            rsid TEXT,
            chromosome TEXT,
            position INTEGER,
            ref TEXT,
            alt TEXT,
            clinical_significance TEXT,
            review_status TEXT,
            review_stars INTEGER,
            conditions TEXT,
            gene TEXT,
            molecular_consequence TEXT,
            origin TEXT,
            clinvar_id TEXT
        )
    """)
    conn.execute("CREATE INDEX idx_rsid ON variants(rsid)")
    conn.execute("CREATE INDEX idx_gene ON variants(gene)")

    conn.executemany(
        "INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
        [
            ("rs121918506", "17", 7673803, "G", "A",
             "Pathogenic", "reviewed_by_expert_panel", 3,
             "Li-Fraumeni syndrome|Hereditary cancer-predisposing syndrome",
             "TP53", "missense_variant", "germline", "12345"),
            ("rs1800497", "11", 113400106, "G", "A",
             "Benign", "criteria_provided,_single_submitter", 1,
             "not_provided", "ANKK1", "missense_variant", "germline", "67890"),
            ("rs4244285", "10", 96541616, "G", "A",
             "Pathogenic", "practice_guideline", 4,
             "CYP2C19 poor metabolizer",
             "CYP2C19", "splice_donor_variant", "germline", "11111"),
            ("rs1057910", "10", 96702047, "A", "C",
             "Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts", 2,
             "warfarin sensitivity",
             "CYP2C9", "missense_variant", "germline", "22222"),
            ("rs9999999", "1", 100000, "C", "T",
             "Uncertain significance", "criteria_provided,_single_submitter", 1,
             "not_provided", "FAKE", "synonymous_variant", "germline", "99999"),
        ],
    )

    conn.commit()
    conn.close()

    return str(tmp_path)


# ── Tests ────────────────────────────────────────────────────────────────────

class TestClinVarMirror:

    def test_is_available(self, mock_clinvar_dir):
        mirror = ClinVarMirror(db_dir=mock_clinvar_dir)
        assert mirror.is_available

    def test_lookup_existing(self, mock_clinvar_dir):
        with ClinVarMirror(db_dir=mock_clinvar_dir) as mirror:
            entries = mirror.lookup("rs121918506")

            assert len(entries) == 1
            entry = entries[0]
            assert entry.rsid == "rs121918506"
            assert entry.gene == "TP53"
            assert entry.is_pathogenic
            assert entry.review_stars == 3
            assert "Li-Fraumeni" in entry.conditions[0]

    def test_lookup_nonexistent(self, mock_clinvar_dir):
        with ClinVarMirror(db_dir=mock_clinvar_dir) as mirror:
            entries = mirror.lookup("rs00000000")
            assert len(entries) == 0

    def test_lookup_case_insensitive(self, mock_clinvar_dir):
        with ClinVarMirror(db_dir=mock_clinvar_dir) as mirror:
            entries = mirror.lookup("RS121918506")
            assert len(entries) == 1

    def test_lookup_gene(self, mock_clinvar_dir):
        with ClinVarMirror(db_dir=mock_clinvar_dir) as mirror:
            entries = mirror.lookup_gene("CYP2C19")
            assert len(entries) == 1
            assert entries[0].rsid == "rs4244285"

    def test_lookup_pathogenic(self, mock_clinvar_dir):
        with ClinVarMirror(db_dir=mock_clinvar_dir) as mirror:
            entries = mirror.lookup_pathogenic(min_stars=1)
            assert len(entries) == 3  # TP53, CYP2C19, CYP2C9

    def test_lookup_pathogenic_high_confidence(self, mock_clinvar_dir):
        with ClinVarMirror(db_dir=mock_clinvar_dir) as mirror:
            entries = mirror.lookup_pathogenic(min_stars=3)
            assert len(entries) == 2  # only TP53 (3 stars) and CYP2C19 (4 stars)

    def test_get_stats(self, mock_clinvar_dir):
        with ClinVarMirror(db_dir=mock_clinvar_dir) as mirror:
            stats = mirror.get_stats()

            assert stats["total_variants"] == 5
            assert stats["pathogenic"] == 3
            assert stats["benign"] == 1

    def test_context_manager(self, mock_clinvar_dir):
        with ClinVarMirror(db_dir=mock_clinvar_dir) as mirror:
            entries = mirror.lookup("rs1800497")
            assert len(entries) == 1
            assert entries[0].is_benign

        assert mirror._conn is None


class TestClinVarEntry:

    def test_is_pathogenic(self):
        entry = ClinVarEntry(
            rsid="rs1", chromosome="1", position=1, ref="A", alt="T",
            clinical_significance="Pathogenic",
            review_status="", review_stars=3,
            conditions=[], gene="", molecular_consequence="",
            origin="", clinvar_id="",
        )
        assert entry.is_pathogenic
        assert not entry.is_benign
        assert not entry.is_uncertain

    def test_is_benign(self):
        entry = ClinVarEntry(
            rsid="rs1", chromosome="1", position=1, ref="A", alt="T",
            clinical_significance="Likely benign",
            review_status="", review_stars=1,
            conditions=[], gene="", molecular_consequence="",
            origin="", clinvar_id="",
        )
        assert entry.is_benign
        assert not entry.is_pathogenic

    def test_is_uncertain(self):
        entry = ClinVarEntry(
            rsid="rs1", chromosome="1", position=1, ref="A", alt="T",
            clinical_significance="Uncertain significance",
            review_status="", review_stars=1,
            conditions=[], gene="", molecular_consequence="",
            origin="", clinvar_id="",
        )
        assert entry.is_uncertain


class TestReviewStars:

    def test_practice_guideline(self):
        assert REVIEW_STARS["practice_guideline"] == 4

    def test_expert_panel(self):
        assert REVIEW_STARS["reviewed_by_expert_panel"] == 3

    def test_no_assertion(self):
        assert REVIEW_STARS["no_assertion_criteria_provided"] == 0
