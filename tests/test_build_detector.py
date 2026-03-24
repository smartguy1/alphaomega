"""
AlphaOmega — Build Detector Tests

Tests the cascading genome build detection logic using mock BAM headers.
No actual BAM files needed — we mock pysam.AlignmentFile.
"""

import pytest
from unittest.mock import MagicMock, patch

from src.build_detector import (
    BuildDetectionResult,
    ChromosomeStyle,
    GenomeBuild,
    GRCH37_LENGTHS,
    GRCH38_LENGTHS,
    detect_from_dictionary,
    detect_from_header,
    detect_genome_build,
    _normalize_chrom_name,
    _detect_chromosome_style,
)


# ── Helpers ──────────────────────────────────────────────────────────────────

def _mock_bam(references: list[str], lengths: list[int]):
    """Create a mock pysam.AlignmentFile with given references and lengths."""
    mock = MagicMock()
    mock.references = references
    mock.lengths = lengths
    mock.close = MagicMock()
    return mock


def _make_ucsc_refs(length_dict: dict[str, int]) -> tuple[list[str], list[int]]:
    """Convert normalized length dict to UCSC-style (chr1, chr2, ...) lists."""
    refs = [f"chr{k}" for k in length_dict.keys()]
    lens = list(length_dict.values())
    return refs, lens


def _make_ensembl_refs(length_dict: dict[str, int]) -> tuple[list[str], list[int]]:
    """Convert normalized length dict to Ensembl-style (1, 2, ...) lists."""
    refs = list(length_dict.keys())
    lens = list(length_dict.values())
    return refs, lens


# ── Unit tests ───────────────────────────────────────────────────────────────

class TestNormalizeChromName:
    def test_ucsc_style(self):
        assert _normalize_chrom_name("chr1") == "1"
        assert _normalize_chrom_name("chrX") == "X"
        assert _normalize_chrom_name("chr22") == "22"

    def test_ensembl_style(self):
        assert _normalize_chrom_name("1") == "1"
        assert _normalize_chrom_name("X") == "X"

    def test_already_normalized(self):
        assert _normalize_chrom_name("MT") == "MT"


class TestDetectChromosomeStyle:
    def test_ucsc(self):
        refs = ["chr1", "chr2", "chr3", "chrX"]
        assert _detect_chromosome_style(refs) == ChromosomeStyle.UCSC

    def test_ensembl(self):
        refs = ["1", "2", "3", "X"]
        assert _detect_chromosome_style(refs) == ChromosomeStyle.ENSEMBL

    def test_empty(self):
        assert _detect_chromosome_style([]) == ChromosomeStyle.UNKNOWN

    def test_mixed_majority_ucsc(self):
        refs = ["chr1", "chr2", "3", "chrX"]
        assert _detect_chromosome_style(refs) == ChromosomeStyle.UCSC


# ── Header detection tests ───────────────────────────────────────────────────

class TestDetectFromHeader:

    @patch("src.build_detector.pysam")
    def test_grch38_ucsc(self, mock_pysam):
        refs, lens = _make_ucsc_refs(GRCH38_LENGTHS)
        mock_pysam.AlignmentFile.return_value = _mock_bam(refs, lens)

        result = detect_from_header("fake.bam")

        assert result.build == GenomeBuild.GRCH38
        assert result.confidence == 0.95
        assert result.method == "header_chr1"
        assert result.chromosome_style == ChromosomeStyle.UCSC

    @patch("src.build_detector.pysam")
    def test_grch37_ensembl(self, mock_pysam):
        refs, lens = _make_ensembl_refs(GRCH37_LENGTHS)
        mock_pysam.AlignmentFile.return_value = _mock_bam(refs, lens)

        result = detect_from_header("fake.bam")

        assert result.build == GenomeBuild.GRCH37
        assert result.confidence == 0.95
        assert result.chromosome_style == ChromosomeStyle.ENSEMBL

    @patch("src.build_detector.pysam")
    def test_unknown_length(self, mock_pysam):
        mock_pysam.AlignmentFile.return_value = _mock_bam(
            ["chr1"], [999999999]
        )

        result = detect_from_header("fake.bam")

        assert result.build == GenomeBuild.UNKNOWN
        assert result.confidence == 0.0

    @patch("src.build_detector.pysam")
    def test_no_chr1(self, mock_pysam):
        # Partial BAM (e.g., exome with only specific chromosomes)
        mock_pysam.AlignmentFile.return_value = _mock_bam(
            ["chr22"], [50818468]
        )

        result = detect_from_header("fake.bam")

        assert result.build == GenomeBuild.UNKNOWN
        assert "not found" in result.details


# ── Dictionary detection tests ───────────────────────────────────────────────

class TestDetectFromDictionary:

    @patch("src.build_detector.pysam")
    def test_grch38_full_match(self, mock_pysam):
        refs, lens = _make_ucsc_refs(GRCH38_LENGTHS)
        mock_pysam.AlignmentFile.return_value = _mock_bam(refs, lens)

        result = detect_from_dictionary("fake.bam")

        assert result.build == GenomeBuild.GRCH38
        assert result.confidence == 1.0
        assert result.method == "dictionary"

    @patch("src.build_detector.pysam")
    def test_grch37_full_match(self, mock_pysam):
        refs, lens = _make_ucsc_refs(GRCH37_LENGTHS)
        mock_pysam.AlignmentFile.return_value = _mock_bam(refs, lens)

        result = detect_from_dictionary("fake.bam")

        assert result.build == GenomeBuild.GRCH37
        assert result.confidence == 1.0

    @patch("src.build_detector.pysam")
    def test_partial_match(self, mock_pysam):
        # Only chr1 and chr2 from GRCh38
        mock_pysam.AlignmentFile.return_value = _mock_bam(
            ["chr1", "chr2"],
            [GRCH38_LENGTHS["1"], GRCH38_LENGTHS["2"]],
        )

        result = detect_from_dictionary("fake.bam")

        assert result.build == GenomeBuild.GRCH38
        assert result.confidence > 0.6

    @patch("src.build_detector.pysam")
    def test_ambiguous(self, mock_pysam):
        # Mix of lengths from both builds (should not happen in practice)
        mock_pysam.AlignmentFile.return_value = _mock_bam(
            ["chr1", "chr2"],
            [GRCH38_LENGTHS["1"], GRCH37_LENGTHS["2"]],
        )

        result = detect_from_dictionary("fake.bam")

        # Neither build gets majority
        assert result.build == GenomeBuild.UNKNOWN


# ── Master function tests ────────────────────────────────────────────────────

class TestDetectGenomeBuild:

    @patch("src.build_detector.pysam")
    def test_stops_at_header_if_confident(self, mock_pysam):
        refs, lens = _make_ucsc_refs(GRCH38_LENGTHS)
        mock_pysam.AlignmentFile.return_value = _mock_bam(refs, lens)

        result = detect_genome_build("fake.bam")

        assert result.build == GenomeBuild.GRCH38
        assert result.method == "header_chr1"  # Stopped at first step

    @patch("src.build_detector.pysam")
    def test_falls_through_to_dictionary(self, mock_pysam):
        # chr1 has wrong length, but chr2/chr10 match GRCh38
        mock_pysam.AlignmentFile.return_value = _mock_bam(
            ["chr1", "chr2", "chr3", "chr10", "chr22", "chrX"],
            [
                999999999,              # chr1 doesn't match
                GRCH38_LENGTHS["2"],
                GRCH38_LENGTHS["3"],
                GRCH38_LENGTHS["10"],
                GRCH38_LENGTHS["22"],
                GRCH38_LENGTHS["X"],
            ],
        )

        result = detect_genome_build("fake.bam")

        assert result.build == GenomeBuild.GRCH38
        assert result.method == "dictionary"  # Fell through to step 2


# ── Result serialization ─────────────────────────────────────────────────────

class TestBuildDetectionResult:
    def test_to_dict(self):
        result = BuildDetectionResult(
            build=GenomeBuild.GRCH38,
            confidence=0.95,
            method="header_chr1",
            chromosome_style=ChromosomeStyle.UCSC,
            details="test details",
        )
        d = result.to_dict()

        assert d["build"] == "GRCh38"
        assert d["confidence"] == 0.95
        assert d["chromosome_style"] == "ucsc"
        assert d["method"] == "header_chr1"
