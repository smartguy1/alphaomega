"""
AlphaOmega — PharmCAT Runner

Wraps the PharmCAT (Pharmacogenomics Clinical Annotation Tool) for
extracting pharmacogenomic variants from VCF files and generating
drug-gene interaction reports.

PharmCAT is a Java-based tool from PharmGKB/CPIC that:
  1. Extracts pharmacogenomic variants from VCF
  2. Infers star alleles and diplotypes
  3. Matches to CPIC guideline recommendations
  4. Generates prescribing recommendations

Since PharmCAT requires Java, this module provides:
  - Automatic PharmCAT JAR download
  - VCF preprocessing for PharmCAT compatibility
  - Result parsing into structured Python objects
"""

import json
import logging
import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import requests

logger = logging.getLogger(__name__)

# PharmCAT download
PHARMCAT_VERSION = "2.15.4"
PHARMCAT_JAR_URL = (
    f"https://github.com/PharmGKB/PharmCAT/releases/download/"
    f"v{PHARMCAT_VERSION}/pharmcat-{PHARMCAT_VERSION}-all.jar"
)

DEFAULT_DB_DIR = "db"
DEFAULT_JAR_NAME = "pharmcat.jar"


@dataclass
class DrugRecommendation:
    """A single drug-gene interaction recommendation."""
    drug: str
    gene: str
    diplotype: str
    phenotype: str
    activity_score: Optional[float]
    recommendation: str
    classification: str    # e.g., "Strong", "Moderate"
    guideline_url: str = ""


@dataclass
class GeneResult:
    """PharmCAT result for a single gene."""
    gene: str
    diplotype: str
    phenotype: str
    activity_score: Optional[float]
    allele1: str = ""
    allele2: str = ""
    missing_positions: list[str] = field(default_factory=list)


@dataclass
class PharmCATResult:
    """Full PharmCAT analysis result."""
    gene_results: list[GeneResult]
    drug_recommendations: list[DrugRecommendation]
    missing_genes: list[str]
    report_path: Optional[str]
    success: bool
    error: Optional[str] = None


# Key pharmacogenes tracked by PharmCAT / CPIC
PHARMCAT_GENES = [
    "CYP2D6", "CYP2C19", "CYP2C9", "CYP3A5", "CYP2B6",
    "CYP4F2", "DPYD", "G6PD", "IFNL3", "NUDT15",
    "RYR1", "SLCO1B1", "TPMT", "UGT1A1", "VKORC1",
    "CYP2C8", "NAT2", "ABCG2", "CACNA1S", "HLA-A",
    "HLA-B", "MT-RNR1",
]


class PharmCATRunner:
    """
    Interface to PharmCAT for pharmacogenomic analysis.

    Usage:
        runner = PharmCATRunner(db_dir="db")
        runner.ensure_installed()
        result = runner.analyze("path/to/variants.vcf.gz")
    """

    def __init__(self, db_dir: str = DEFAULT_DB_DIR):
        self.db_dir = Path(db_dir)
        self.jar_path = self.db_dir / DEFAULT_JAR_NAME

    @property
    def is_available(self) -> bool:
        """Check if PharmCAT JAR is downloaded."""
        return self.jar_path.exists()

    @property
    def java_available(self) -> bool:
        """Check if Java is available on the system."""
        try:
            result = subprocess.run(
                ["java", "-version"],
                capture_output=True, text=True, timeout=10,
            )
            return result.returncode == 0
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False

    def ensure_installed(self, force: bool = False) -> bool:
        """
        Download the PharmCAT JAR if not present.
        Returns True if download was performed.
        """
        if self.is_available and not force:
            logger.info(f"PharmCAT already installed: {self.jar_path}")
            return False

        if not self.java_available:
            logger.warning(
                "Java not found. PharmCAT requires Java 17+. "
                "Install with: apt-get install openjdk-17-jre (Linux) "
                "or download from https://adoptium.net/"
            )

        self.db_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Downloading PharmCAT v{PHARMCAT_VERSION}...")
        logger.info(f"URL: {PHARMCAT_JAR_URL}")

        try:
            response = requests.get(
                PHARMCAT_JAR_URL, stream=True, timeout=120,
                headers={"Accept": "application/octet-stream"},
            )
            response.raise_for_status()

            total = int(response.headers.get("content-length", 0))
            logger.info(f"Download size: {total / 1024 / 1024:.1f} MB")

            with open(self.jar_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=65536):
                    f.write(chunk)

            logger.info(f"PharmCAT installed: {self.jar_path}")
            return True

        except requests.RequestException as e:
            logger.error(f"Download failed: {e}")
            raise RuntimeError(f"Failed to download PharmCAT: {e}") from e

    def analyze(
        self,
        vcf_path: str,
        output_dir: Optional[str] = None,
        sample_id: Optional[str] = None,
    ) -> PharmCATResult:
        """
        Run PharmCAT analysis on a VCF file.

        Args:
            vcf_path: Path to the input VCF file
            output_dir: Directory for output files (default: same as VCF)
            sample_id: Sample ID to use (default: auto-detect)

        Returns:
            PharmCATResult with gene results and drug recommendations
        """
        if not self.is_available:
            return PharmCATResult(
                gene_results=[], drug_recommendations=[],
                missing_genes=PHARMCAT_GENES, report_path=None,
                success=False,
                error="PharmCAT not installed. Run ensure_installed() first.",
            )

        if not self.java_available:
            return PharmCATResult(
                gene_results=[], drug_recommendations=[],
                missing_genes=PHARMCAT_GENES, report_path=None,
                success=False,
                error="Java not found. PharmCAT requires Java 17+.",
            )

        if not os.path.exists(vcf_path):
            return PharmCATResult(
                gene_results=[], drug_recommendations=[],
                missing_genes=PHARMCAT_GENES, report_path=None,
                success=False, error=f"VCF file not found: {vcf_path}",
            )

        # Determine output directory
        if output_dir is None:
            output_dir = os.path.dirname(vcf_path)
        os.makedirs(output_dir, exist_ok=True)

        # Build PharmCAT command
        cmd = [
            "java", "-jar", str(self.jar_path),
            "-vcf", vcf_path,
            "-o", output_dir,
        ]

        if sample_id:
            cmd.extend(["-s", sample_id])

        logger.info(f"Running PharmCAT: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True, text=True,
                timeout=600,  # 10 min timeout
            )

            if result.returncode != 0:
                logger.error(f"PharmCAT failed: {result.stderr}")
                return PharmCATResult(
                    gene_results=[], drug_recommendations=[],
                    missing_genes=PHARMCAT_GENES, report_path=None,
                    success=False,
                    error=f"PharmCAT exited with code {result.returncode}: {result.stderr[:500]}",
                )

            # Parse results
            return self._parse_results(output_dir, vcf_path)

        except subprocess.TimeoutExpired:
            return PharmCATResult(
                gene_results=[], drug_recommendations=[],
                missing_genes=PHARMCAT_GENES, report_path=None,
                success=False, error="PharmCAT timed out (>10 min)",
            )

    def _parse_results(
        self, output_dir: str, vcf_path: str
    ) -> PharmCATResult:
        """Parse PharmCAT output files into structured results."""
        base_name = Path(vcf_path).stem
        if base_name.endswith(".vcf"):
            base_name = base_name[:-4]

        # PharmCAT outputs a JSON report
        json_report = None
        html_report = None

        for ext in [".report.json", "_report.json", ".json"]:
            candidate = os.path.join(output_dir, f"{base_name}{ext}")
            if os.path.exists(candidate):
                json_report = candidate
                break

        for ext in [".report.html", "_report.html", ".html"]:
            candidate = os.path.join(output_dir, f"{base_name}{ext}")
            if os.path.exists(candidate):
                html_report = candidate
                break

        # Also search for any JSON report in the output directory
        if json_report is None:
            for f in Path(output_dir).glob("*.json"):
                if "report" in f.name.lower() or "pharmcat" in f.name.lower():
                    json_report = str(f)
                    break

        gene_results = []
        drug_recommendations = []
        missing_genes = list(PHARMCAT_GENES)

        if json_report and os.path.exists(json_report):
            try:
                with open(json_report, "r") as f:
                    report_data = json.load(f)

                gene_results, drug_recommendations, missing_genes = (
                    self._extract_from_json(report_data)
                )
            except (json.JSONDecodeError, KeyError) as e:
                logger.warning(f"Could not parse PharmCAT JSON: {e}")

        return PharmCATResult(
            gene_results=gene_results,
            drug_recommendations=drug_recommendations,
            missing_genes=missing_genes,
            report_path=html_report or json_report,
            success=True,
        )

    def _extract_from_json(
        self, data: dict
    ) -> tuple[list[GeneResult], list[DrugRecommendation], list[str]]:
        """Extract structured results from PharmCAT JSON report."""
        gene_results = []
        drug_recommendations = []
        found_genes = set()

        # Extract gene call results
        gene_calls = data.get("geneResults", data.get("geneCalls", []))
        if isinstance(gene_calls, dict):
            gene_calls = list(gene_calls.values())

        for gc in gene_calls:
            gene = gc.get("gene", gc.get("geneSymbol", ""))
            if not gene:
                continue

            found_genes.add(gene)

            diplotype = gc.get("diplotype", gc.get("printDiplotype", ""))
            phenotype = gc.get("phenotype", gc.get("lookupDiplotype", ""))
            activity = gc.get("activityScore", None)

            if isinstance(activity, str):
                try:
                    activity = float(activity)
                except ValueError:
                    activity = None

            gene_results.append(GeneResult(
                gene=gene,
                diplotype=str(diplotype),
                phenotype=str(phenotype),
                activity_score=activity,
            ))

        # Extract drug recommendations
        drug_recs = data.get("drugRecommendations", data.get("guidelines", []))
        if isinstance(drug_recs, dict):
            drug_recs = list(drug_recs.values())

        for rec in drug_recs:
            drugs = rec.get("drugs", rec.get("drug", []))
            if isinstance(drugs, str):
                drugs = [drugs]
            elif isinstance(drugs, list) and drugs and isinstance(drugs[0], dict):
                drugs = [d.get("name", str(d)) for d in drugs]

            gene = rec.get("gene", rec.get("relatedGene", ""))
            classification = rec.get("classification", rec.get("strength", ""))
            recommendation = rec.get("recommendation", rec.get("implications", ""))
            phenotype = rec.get("phenotype", "")
            diplotype = rec.get("diplotype", "")
            url = rec.get("url", rec.get("guidelineUrl", ""))

            for drug in drugs:
                drug_recommendations.append(DrugRecommendation(
                    drug=str(drug),
                    gene=str(gene),
                    diplotype=str(diplotype),
                    phenotype=str(phenotype),
                    activity_score=None,
                    recommendation=str(recommendation),
                    classification=str(classification),
                    guideline_url=str(url),
                ))

        # Determine missing genes
        missing = [g for g in PHARMCAT_GENES if g not in found_genes]

        return gene_results, drug_recommendations, missing

    def get_supported_genes(self) -> list[str]:
        """Return list of genes tracked by PharmCAT."""
        return list(PHARMCAT_GENES)
