"""
AlphaOmega — Genome Build Detector

Auto-detects whether a BAM file was aligned to GRCh37 or GRCh38
using a cascading approach:
  1. Header-based detection (chr1 length)
  2. Multi-chromosome dictionary fingerprinting
  3. SNP position probing fallback
"""

import logging
import subprocess
from dataclasses import dataclass
from enum import Enum
from typing import Optional

try:
    import pysam
except ImportError:
    pysam = None

logger = logging.getLogger(__name__)


class GenomeBuild(Enum):
    GRCH37 = "GRCh37"
    GRCH38 = "GRCh38"
    UNKNOWN = "unknown"


class ChromosomeStyle(Enum):
    UCSC = "ucsc"       # chr1, chr2, chrX
    ENSEMBL = "ensembl" # 1, 2, X
    UNKNOWN = "unknown"


@dataclass
class BuildDetectionResult:
    """Result of genome build detection."""
    build: GenomeBuild
    confidence: float
    method: str
    chromosome_style: ChromosomeStyle
    details: Optional[str] = None

    def to_dict(self) -> dict:
        return {
            "build": self.build.value,
            "confidence": self.confidence,
            "method": self.method,
            "chromosome_style": self.chromosome_style.value,
            "details": self.details,
        }


# Known chromosome lengths for fingerprinting
# Using a subset of chromosomes that have distinct lengths between builds
GRCH37_LENGTHS = {
    "1":  249250621,
    "2":  243199373,
    "3":  198022430,
    "10": 135534747,
    "22": 51304566,
    "X":  155270560,
}

GRCH38_LENGTHS = {
    "1":  248956422,
    "2":  242193529,
    "3":  198295559,
    "10": 133797422,
    "22": 50818468,
    "X":  156040895,
}


def _normalize_chrom_name(name: str) -> str:
    """Strip 'chr' prefix for comparison."""
    return name.replace("chr", "")


def _detect_chromosome_style(references: list[str]) -> ChromosomeStyle:
    """Detect whether BAM uses UCSC (chr1) or Ensembl (1) naming."""
    chr_prefixed = sum(1 for r in references if r.startswith("chr"))
    if chr_prefixed > len(references) / 2:
        return ChromosomeStyle.UCSC
    elif len(references) > 0:
        return ChromosomeStyle.ENSEMBL
    return ChromosomeStyle.UNKNOWN


def _get_ref_lengths(bam_path: str) -> tuple[dict[str, int], ChromosomeStyle]:
    """
    Extract reference sequence lengths from BAM header.
    Returns normalized (no 'chr' prefix) name→length dict and chromosome style.
    """
    if pysam is None:
        raise RuntimeError(
            "pysam is not installed. Install with: pip install pysam"
        )

    bam = pysam.AlignmentFile(bam_path, "rb")
    references = list(bam.references)
    lengths = list(bam.lengths)
    bam.close()

    style = _detect_chromosome_style(references)

    normalized = {}
    for ref, length in zip(references, lengths):
        norm = _normalize_chrom_name(ref)
        normalized[norm] = length

    return normalized, style


def detect_from_header(bam_path: str) -> BuildDetectionResult:
    """
    Step 1: Fast detection using chr1 length alone.
    Works ~90% of the time.
    """
    refs, style = _get_ref_lengths(bam_path)

    chr1_len = refs.get("1")
    if chr1_len is None:
        return BuildDetectionResult(
            build=GenomeBuild.UNKNOWN,
            confidence=0.0,
            method="header_chr1",
            chromosome_style=style,
            details="chr1/1 not found in BAM header",
        )

    if chr1_len == GRCH37_LENGTHS["1"]:
        return BuildDetectionResult(
            build=GenomeBuild.GRCH37,
            confidence=0.95,
            method="header_chr1",
            chromosome_style=style,
            details=f"chr1 length={chr1_len} matches GRCh37",
        )

    if chr1_len == GRCH38_LENGTHS["1"]:
        return BuildDetectionResult(
            build=GenomeBuild.GRCH38,
            confidence=0.95,
            method="header_chr1",
            chromosome_style=style,
            details=f"chr1 length={chr1_len} matches GRCh38",
        )

    return BuildDetectionResult(
        build=GenomeBuild.UNKNOWN,
        confidence=0.0,
        method="header_chr1",
        chromosome_style=style,
        details=f"chr1 length={chr1_len} matches neither build",
    )


def detect_from_dictionary(bam_path: str) -> BuildDetectionResult:
    """
    Step 2: Multi-chromosome fingerprinting for higher confidence.
    Compares lengths of chr1, chr2, chr3, chr10, chr22, chrX.
    """
    refs, style = _get_ref_lengths(bam_path)

    def _score(signature: dict[str, int]) -> tuple[float, int, int]:
        matches = 0
        total = 0
        for chrom, expected_length in signature.items():
            actual = refs.get(chrom)
            if actual is not None:
                total += 1
                if actual == expected_length:
                    matches += 1
        score = matches / total if total > 0 else 0.0
        return score, matches, total

    score37, m37, t37 = _score(GRCH37_LENGTHS)
    score38, m38, t38 = _score(GRCH38_LENGTHS)

    if score37 > score38 and score37 > 0.6:
        return BuildDetectionResult(
            build=GenomeBuild.GRCH37,
            confidence=score37,
            method="dictionary",
            chromosome_style=style,
            details=f"GRCh37: {m37}/{t37} chromosomes matched",
        )

    if score38 > score37 and score38 > 0.6:
        return BuildDetectionResult(
            build=GenomeBuild.GRCH38,
            confidence=score38,
            method="dictionary",
            chromosome_style=style,
            details=f"GRCh38: {m38}/{t38} chromosomes matched",
        )

    return BuildDetectionResult(
        build=GenomeBuild.UNKNOWN,
        confidence=max(score37, score38),
        method="dictionary",
        chromosome_style=style,
        details=f"Ambiguous: GRCh37={m37}/{t37}, GRCh38={m38}/{t38}",
    )


def detect_by_probe(bam_path: str) -> BuildDetectionResult:
    """
    Step 3: Probe known differing SNP positions as a last resort.
    Tests read coverage at positions that differ between builds.
    Requires samtools on PATH.
    """
    _, style = _get_ref_lengths(bam_path)

    chrom_prefix = "chr" if style == ChromosomeStyle.UCSC else ""

    # Position that exists in GRCh37 but shifted in GRCh38
    grch37_region = f"{chrom_prefix}1:752566-752566"
    grch38_region = f"{chrom_prefix}1:752721-752721"

    def _has_coverage(region: str) -> bool:
        try:
            result = subprocess.run(
                ["samtools", "mpileup", "-r", region, bam_path],
                capture_output=True, text=True, timeout=30,
            )
            return len(result.stdout.strip()) > 0
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            logger.warning(f"Probe failed for {region}: {e}")
            return False

    hit37 = _has_coverage(grch37_region)
    hit38 = _has_coverage(grch38_region)

    if hit37 and not hit38:
        return BuildDetectionResult(
            build=GenomeBuild.GRCH37,
            confidence=0.85,
            method="probe",
            chromosome_style=style,
            details="Coverage at GRCh37-specific position",
        )

    if hit38 and not hit37:
        return BuildDetectionResult(
            build=GenomeBuild.GRCH38,
            confidence=0.85,
            method="probe",
            chromosome_style=style,
            details="Coverage at GRCh38-specific position",
        )

    return BuildDetectionResult(
        build=GenomeBuild.UNKNOWN,
        confidence=0.0,
        method="probe",
        chromosome_style=style,
        details=f"Probe inconclusive: hit37={hit37}, hit38={hit38}",
    )


def detect_genome_build(bam_path: str) -> BuildDetectionResult:
    """
    Master detection function. Cascades through methods:
      1. Header-based (fast, ~95% accurate)
      2. Dictionary fingerprinting (robust)
      3. SNP probing (deterministic last resort)

    Returns the first confident result.
    """
    logger.info(f"Detecting genome build for: {bam_path}")

    # Step 1: fast path
    result = detect_from_header(bam_path)
    if result.build != GenomeBuild.UNKNOWN:
        logger.info(f"Build detected via header: {result.build.value} "
                     f"(confidence={result.confidence:.2f})")
        return result

    logger.info("Header detection inconclusive, trying dictionary matching...")

    # Step 2: dictionary matching
    result = detect_from_dictionary(bam_path)
    if result.build != GenomeBuild.UNKNOWN:
        logger.info(f"Build detected via dictionary: {result.build.value} "
                     f"(confidence={result.confidence:.2f})")
        return result

    logger.info("Dictionary matching inconclusive, trying probe...")

    # Step 3: probe fallback
    result = detect_by_probe(bam_path)
    if result.build != GenomeBuild.UNKNOWN:
        logger.info(f"Build detected via probe: {result.build.value} "
                     f"(confidence={result.confidence:.2f})")
    else:
        logger.warning("All detection methods failed — build is unknown")

    return result
