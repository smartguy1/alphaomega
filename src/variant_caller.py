"""
AlphaOmega — Variant Caller

Wraps bcftools to perform BAM → VCF variant calling.
Auto-selects reference genome based on detected build.
Filters output to SNPedia-relevant rsIDs when available.
"""

import hashlib
import logging
import os
import subprocess
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from src.build_detector import (
    BuildDetectionResult,
    ChromosomeStyle,
    GenomeBuild,
    detect_genome_build,
)

logger = logging.getLogger(__name__)


@dataclass
class VariantCallingResult:
    """Result of a variant calling run."""
    vcf_path: str
    build: GenomeBuild
    chromosome_style: ChromosomeStyle
    total_variants: int
    filtered_variants: Optional[int]  # After SNPedia rsID filter
    bam_hash: str
    success: bool
    error: Optional[str] = None

    def to_dict(self) -> dict:
        return {
            "vcf_path": self.vcf_path,
            "build": self.build.value,
            "chromosome_style": self.chromosome_style.value,
            "total_variants": self.total_variants,
            "filtered_variants": self.filtered_variants,
            "bam_hash": self.bam_hash,
            "success": self.success,
            "error": self.error,
        }


def compute_bam_hash(bam_path: str, chunk_size: int = 8192) -> str:
    """
    Compute a fast hash of the BAM file for cache invalidation.
    Uses first + last chunks rather than hashing the entire file (which
    can be 30-100GB).
    """
    file_size = os.path.getsize(bam_path)
    hasher = hashlib.md5()

    with open(bam_path, "rb") as f:
        # Hash first chunk
        hasher.update(f.read(chunk_size))

        # Hash last chunk
        if file_size > chunk_size:
            f.seek(max(0, file_size - chunk_size))
            hasher.update(f.read(chunk_size))

    # Include file size in hash for extra safety
    hasher.update(str(file_size).encode())
    return hasher.hexdigest()


def _get_reference_path(
    build: GenomeBuild,
    reference_dir: str,
) -> Optional[str]:
    """
    Find the reference genome FASTA for the given build.
    Expects files named like: GRCh38.fa or GRCh37.fa
    """
    ref_dir = Path(reference_dir)
    if not ref_dir.exists():
        return None

    patterns = {
        GenomeBuild.GRCH37: [
            "GRCh37.fa", "GRCh37.fasta",
            "hs37d5.fa", "human_g1k_v37.fasta",
            "hg19.fa",
        ],
        GenomeBuild.GRCH38: [
            "GRCh38.fa", "GRCh38.fasta",
            "GRCh38_full_analysis_set_plus_decoy_hla.fa",
            "hg38.fa",
        ],
    }

    for pattern in patterns.get(build, []):
        candidate = ref_dir / pattern
        if candidate.exists():
            return str(candidate)

    # Try any .fa or .fasta file in the directory (fallback)
    for ext in (".fa", ".fasta"):
        for f in ref_dir.glob(f"*{build.value}*{ext}"):
            return str(f)

    return None


def _count_variants(vcf_path: str) -> int:
    """Count the number of variant records in a VCF file."""
    try:
        result = subprocess.run(
            ["bcftools", "view", "-H", vcf_path],
            capture_output=True, text=True, timeout=120,
        )
        return len(result.stdout.strip().splitlines())
    except Exception:
        return -1


def _annotate_rsids(
    vcf_path: str,
    db_dir: str,
    build: str = "GRCh37",
) -> tuple[bool, int]:
    """
    Inject rsIDs into the VCF by matching chromosome + position against
    a build-matched ClinVar VCF. The ClinVar VCF stores rsIDs in
    the INFO/RS field (not the ID column), so we extract rs-prefixed IDs
    into a tab file and use bcftools annotate.

    Returns (success, count_of_annotated_variants).
    """
    import gzip as gzip_mod

    # Use build-matched ClinVar VCF
    clinvar_vcf = os.path.join(db_dir, f"clinvar_{build.lower()}.vcf.gz")

    # Download if not present
    if not os.path.exists(clinvar_vcf):
        url = f"https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_{build}/clinvar.vcf.gz"
        logger.info(f"Downloading {build} ClinVar VCF for rsID annotation...")
        logger.info(f"URL: {url}")
        try:
            import requests
            resp = requests.get(url, stream=True, timeout=300)
            resp.raise_for_status()
            with open(clinvar_vcf, "wb") as f:
                for chunk in resp.iter_content(chunk_size=65536):
                    f.write(chunk)
            logger.info(f"Downloaded ClinVar {build} VCF: {clinvar_vcf}")
        except Exception as e:
            logger.warning(f"Failed to download ClinVar {build} VCF: {e}")
            # Fall back to existing clinvar.vcf.gz if available
            fallback = os.path.join(db_dir, "clinvar.vcf.gz")
            if os.path.exists(fallback):
                clinvar_vcf = fallback
                logger.info("Falling back to existing clinvar.vcf.gz")
            else:
                return False, 0

    # Step 1: Build a tab-delimited rsID annotation file from ClinVar
    # ClinVar VCF has rsID in INFO/RS field, not the ID column
    tab_file = os.path.join(db_dir, f"rsid_map_{build.lower()}.tab.gz")

    if not os.path.exists(tab_file):
        logger.info("Building rsID annotation table from ClinVar...")
        tab_path_raw = tab_file.replace(".gz", "")
        count = 0

        try:
            with gzip_mod.open(clinvar_vcf, "rt", errors="replace") as fin, \
                 open(tab_path_raw, "w") as fout:
                for line in fin:
                    if line.startswith("#"):
                        continue
                    fields = line.strip().split("\t")
                    if len(fields) < 8:
                        continue
                    chrom = fields[0]
                    pos = fields[1]
                    info = fields[7]

                    # Extract RS number from INFO field
                    rs_num = None
                    for item in info.split(";"):
                        if item.startswith("RS="):
                            rs_num = item[3:]
                            break

                    if rs_num and rs_num != ".":
                        fout.write(f"{chrom}\t{pos}\trs{rs_num}\n")
                        count += 1

            logger.info(f"Extracted {count} rsID mappings from ClinVar")

            # Sort, compress with bgzip, and index with tabix
            sorted_file = tab_path_raw + ".sorted"
            subprocess.run(
                f"sort -k1,1V -k2,2n {tab_path_raw} > {sorted_file}",
                shell=True, capture_output=True, timeout=300,
            )
            os.replace(sorted_file, tab_path_raw)

            subprocess.run(
                ["bgzip", "-f", tab_path_raw],
                capture_output=True, timeout=300,
            )
            subprocess.run(
                ["tabix", "-s1", "-b2", "-e2", tab_file],
                capture_output=True, timeout=300,
            )

        except Exception as e:
            logger.warning(f"Failed to build rsID table: {e}")
            # Clean up partial files
            for f in [tab_path_raw, tab_file, tab_file + ".tbi"]:
                if os.path.exists(f):
                    os.remove(f)
            return False, 0

    # Step 2: Use bcftools annotate with the tab file
    tmp_vcf = vcf_path + ".rsid.tmp.vcf.gz"

    try:
        logger.info("Injecting rsIDs into variant calls...")
        result = subprocess.run(
            [
                "bcftools", "annotate",
                "-a", tab_file,
                "-c", "CHROM,POS,ID",
                "-Oz", "-o", tmp_vcf,
                vcf_path,
            ],
            capture_output=True, text=True, timeout=3600,
        )

        if result.returncode != 0:
            logger.warning(f"rsID annotation failed: {result.stderr}")
            if os.path.exists(tmp_vcf):
                os.remove(tmp_vcf)
            return False, 0

        # Replace original with annotated version
        os.replace(tmp_vcf, vcf_path)

        # Re-index
        subprocess.run(
            ["bcftools", "index", "-f", vcf_path],
            capture_output=True, timeout=300,
        )

        # Count how many variants now have rsIDs
        count_result = subprocess.run(
            ["bcftools", "query", "-f", "%ID\n", vcf_path],
            capture_output=True, text=True, timeout=600,
        )
        rsid_count = sum(
            1 for line in count_result.stdout.strip().splitlines()
            if line.startswith("rs")
        )

        logger.info(f"Injected rsIDs: {rsid_count} variants now have rsIDs")
        return True, rsid_count

    except subprocess.TimeoutExpired:
        logger.warning("rsID annotation timed out")
        if os.path.exists(tmp_vcf):
            os.remove(tmp_vcf)
        return False, 0
    except Exception as e:
        logger.warning(f"rsID annotation error: {e}")
        if os.path.exists(tmp_vcf):
            os.remove(tmp_vcf)
        return False, 0


def _run_bcftools(
    bam_path: str,
    reference_path: str,
    output_vcf: str,
    threads: int = 4,
    min_quality: int = 20,
    min_depth: int = 5,
) -> tuple[bool, Optional[str]]:
    """
    Run bcftools mpileup + call pipeline.
    Returns (success, error_message).
    """
    logger.info(f"Starting variant calling: {bam_path}")
    logger.info(f"Reference: {reference_path}")
    logger.info(f"Output: {output_vcf}")

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_vcf), exist_ok=True)

    # Build the pipeline command
    # bcftools mpileup → bcftools call → bcftools filter → output
    mpileup_cmd = [
        "bcftools", "mpileup",
        "-f", reference_path,
        "--threads", str(threads),
        "-a", "FORMAT/AD,FORMAT/DP",    # Allelic depth + total depth
        "-q", str(min_quality),          # Min mapping quality
        "-Q", "20",                      # Min base quality
        bam_path,
    ]

    call_cmd = [
        "bcftools", "call",
        "-mv",                           # Multiallelic + variants only
        "--threads", str(threads),
    ]

    filter_cmd = [
        "bcftools", "filter",
        "-i", f"QUAL>={min_quality} && FORMAT/DP>={min_depth}",
    ]

    view_cmd = [
        "bcftools", "view",
        "-Oz",                           # Output as compressed VCF
        "-o", output_vcf,
    ]

    try:
        # Run as a pipeline: mpileup | call | filter | view
        p_mpileup = subprocess.Popen(
            mpileup_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
        p_call = subprocess.Popen(
            call_cmd,
            stdin=p_mpileup.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
        p_mpileup.stdout.close()

        p_filter = subprocess.Popen(
            filter_cmd,
            stdin=p_call.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
        p_call.stdout.close()

        p_view = subprocess.Popen(
            view_cmd,
            stdin=p_filter.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        p_filter.stdout.close()

        # Background thread to report progress
        stop_event = threading.Event()
        def progress_monitor():
            start_time = time.time()
            while not stop_event.is_set():
                if stop_event.wait(60):  # Update every 60 seconds
                    break
                elapsed = int((time.time() - start_time) / 60)
                size_mb = 0
                if os.path.exists(output_vcf):
                    size_mb = os.path.getsize(output_vcf) / (1024 * 1024)
                logger.info(f"Still calling variants... {elapsed} mins elapsed, {size_mb:.1f} MB written to VCF")

        monitor_thread = threading.Thread(target=progress_monitor, daemon=True)
        monitor_thread.start()

        _, stderr = p_view.communicate(timeout=28800)  # 8 hour timeout
        
        # Stop monitor thread
        stop_event.set()
        monitor_thread.join()

        # Wait for all processes
        p_mpileup.wait()
        p_call.wait()
        p_filter.wait()

        if p_view.returncode != 0:
            return False, f"bcftools pipeline failed: {stderr.decode()}"

        # Index the output VCF
        logger.info("Indexing VCF...")
        idx_result = subprocess.run(
            ["bcftools", "index", output_vcf],
            capture_output=True, text=True, timeout=300,
        )
        if idx_result.returncode != 0:
            return False, f"bcftools index failed: {idx_result.stderr}"

        return True, None

    except subprocess.TimeoutExpired:
        return False, "Variant calling timed out (>8 hours)"
    except FileNotFoundError:
        return False, (
            "bcftools not found. Install via: "
            "apt-get install bcftools (Linux) or "
            "conda install -c bioconda bcftools"
        )


def filter_to_snpedia_rsids(
    input_vcf: str,
    output_vcf: str,
    rsid_file: str,
) -> Optional[int]:
    """
    Filter a VCF to only include variants with rsIDs in the SNPedia database.
    The rsid_file should contain one rsID per line (e.g., rs123).

    Returns the number of filtered variants, or None on failure.
    """
    if not os.path.exists(rsid_file):
        logger.warning(f"rsID filter file not found: {rsid_file}")
        return None

    try:
        result = subprocess.run(
            [
                "bcftools", "view",
                "-i", f"ID=@{rsid_file}",
                "-Oz", "-o", output_vcf,
                input_vcf,
            ],
            capture_output=True, text=True, timeout=600,
        )
        if result.returncode != 0:
            logger.error(f"rsID filtering failed: {result.stderr}")
            return None

        # Index filtered VCF
        subprocess.run(
            ["bcftools", "index", output_vcf],
            capture_output=True, timeout=120,
        )

        return _count_variants(output_vcf)
    except Exception as e:
        logger.error(f"rsID filtering error: {e}")
        return None


def call_variants(
    bam_path: str,
    output_dir: str,
    reference_dir: str,
    profile_name: str = "default",
    threads: int = 4,
    min_quality: int = 20,
    min_depth: int = 5,
    snpedia_rsid_file: Optional[str] = None,
    db_dir: Optional[str] = None,
    force: bool = False,
) -> VariantCallingResult:
    """
    Full variant calling pipeline:
      1. Detect genome build
      2. Find matching reference genome
      3. Run bcftools mpileup | call | filter
      4. Optionally filter to SNPedia rsIDs
      5. Return result with metadata

    Args:
        bam_path: Path to input BAM file
        output_dir: Directory to write VCF output
        reference_dir: Directory containing reference genome FASTA files
        profile_name: Name for output file naming
        threads: Number of threads for bcftools
        min_quality: Minimum variant quality score
        min_depth: Minimum read depth
        snpedia_rsid_file: Optional file with SNPedia rsIDs for filtering
        force: If True, re-run even if cached VCF exists
    """
    # Validate input
    if not os.path.exists(bam_path):
        return VariantCallingResult(
            vcf_path="", build=GenomeBuild.UNKNOWN,
            chromosome_style=ChromosomeStyle.UNKNOWN,
            total_variants=0, filtered_variants=None,
            bam_hash="", success=False,
            error=f"BAM file not found: {bam_path}",
        )

    # Compute BAM hash for caching
    bam_hash = compute_bam_hash(bam_path)
    logger.info(f"BAM hash: {bam_hash}")

    # Output paths
    os.makedirs(output_dir, exist_ok=True)
    raw_vcf = os.path.join(output_dir, f"{profile_name}.raw.vcf.gz")
    filtered_vcf = os.path.join(output_dir, f"{profile_name}.filtered.vcf.gz")

    # Check cache
    if not force and os.path.exists(raw_vcf):
        logger.info(f"Cached VCF exists: {raw_vcf}")
        total = _count_variants(raw_vcf)
        filtered = _count_variants(filtered_vcf) if os.path.exists(filtered_vcf) else None
        # Note: we don't verify bam_hash here — that's the profile manager's job
        build_result = detect_genome_build(bam_path)
        return VariantCallingResult(
            vcf_path=filtered_vcf if os.path.exists(filtered_vcf) else raw_vcf,
            build=build_result.build,
            chromosome_style=build_result.chromosome_style,
            total_variants=total,
            filtered_variants=filtered,
            bam_hash=bam_hash,
            success=True,
        )

    # Step 1: Detect build
    build_result = detect_genome_build(bam_path)
    if build_result.build == GenomeBuild.UNKNOWN:
        return VariantCallingResult(
            vcf_path="", build=GenomeBuild.UNKNOWN,
            chromosome_style=build_result.chromosome_style,
            total_variants=0, filtered_variants=None,
            bam_hash=bam_hash, success=False,
            error="Could not detect genome build. "
                  "Please specify manually or check BAM file integrity.",
        )

    logger.info(f"Detected build: {build_result.build.value} "
                f"(style={build_result.chromosome_style.value})")

    # Step 2: Find reference
    ref_path = _get_reference_path(build_result.build, reference_dir)
    if ref_path is None:
        return VariantCallingResult(
            vcf_path="", build=build_result.build,
            chromosome_style=build_result.chromosome_style,
            total_variants=0, filtered_variants=None,
            bam_hash=bam_hash, success=False,
            error=f"Reference genome for {build_result.build.value} not found "
                  f"in {reference_dir}. Run 'alphaomega setup' to download.",
        )

    # Step 3: Call variants
    success, error = _run_bcftools(
        bam_path=bam_path,
        reference_path=ref_path,
        output_vcf=raw_vcf,
        threads=threads,
        min_quality=min_quality,
        min_depth=min_depth,
    )

    if not success:
        return VariantCallingResult(
            vcf_path="", build=build_result.build,
            chromosome_style=build_result.chromosome_style,
            total_variants=0, filtered_variants=None,
            bam_hash=bam_hash, success=False, error=error,
        )

    total = _count_variants(raw_vcf)
    logger.info(f"Total variants called: {total}")

    # Step 3.5: Annotate with rsIDs from ClinVar
    if db_dir:
        rsid_success, rsid_count = _annotate_rsids(
            raw_vcf, db_dir, build=build_result.build.value
        )
        if rsid_success:
            logger.info(f"rsID annotation complete: {rsid_count} variants tagged")
        else:
            logger.warning("rsID annotation was skipped or failed")

    # Step 4: Filter to SNPedia rsIDs (optional)
    filtered = None
    final_vcf = raw_vcf
    if snpedia_rsid_file:
        filtered = filter_to_snpedia_rsids(raw_vcf, filtered_vcf, snpedia_rsid_file)
        if filtered is not None:
            final_vcf = filtered_vcf
            logger.info(f"Filtered to {filtered} SNPedia variants")

    return VariantCallingResult(
        vcf_path=final_vcf,
        build=build_result.build,
        chromosome_style=build_result.chromosome_style,
        total_variants=total,
        filtered_variants=filtered,
        bam_hash=bam_hash,
        success=True,
    )
