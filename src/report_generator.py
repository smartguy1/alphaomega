"""
AlphaOmega — Report Generator

Renders annotated, clustered genomic data into HTML (and optionally PDF)
reports using Jinja2 templates. No LLM required — this is purely
data-driven, deterministic report generation.
"""

import json
import logging
import os
from datetime import datetime
from pathlib import Path
from typing import Optional

from jinja2 import Environment, FileSystemLoader, select_autoescape

from src.annotator import AnnotatedVariant
from src.avatar import build_avatar
from src.trait_clusterer import Domain, DomainCluster, TraitClusterer

logger = logging.getLogger(__name__)

# Map domain names to emoji icons for the template
DOMAIN_ICONS = {
    "metabolism": "🔥",
    "disease_risk": "🏥",
    "pharmacogenomics": "💊",
    "physical_traits": "🏃",
    "behavior": "🧠",
    "ancestry": "🌍",
    "immunity": "🛡️",
    "uncategorized": "🔬",
}

# Default template directory (relative to project root)
DEFAULT_TEMPLATE_DIR = "templates"
DEFAULT_TEMPLATE_NAME = "report.html"


class ReportGenerator:
    """
    Generates HTML/PDF reports from annotated, clustered genomic data.

    Usage:
        gen = ReportGenerator(template_dir="templates")
        html = gen.render_html(
            profile_name="martin",
            clusters=clusters,
            variants=variants,
        )
        gen.save_html(html, "output/martin_report.html")
    """

    def __init__(
        self,
        template_dir: str = DEFAULT_TEMPLATE_DIR,
        template_name: str = DEFAULT_TEMPLATE_NAME,
    ):
        self.template_dir = Path(template_dir)
        self.template_name = template_name

        # Initialize Jinja2
        if self.template_dir.exists():
            self._env = Environment(
                loader=FileSystemLoader(str(self.template_dir)),
                autoescape=select_autoescape(["html"]),
            )
        else:
            self._env = None

    def render_html(
        self,
        profile_name: str,
        clusters: dict[Domain, DomainCluster],
        variants: Optional[list[AnnotatedVariant]] = None,
        build: str = "",
        drug_recommendations: Optional[list[dict]] = None,
    ) -> str:
        """
        Render the HTML report from clustered data.

        Args:
            profile_name: The profile name for the report header
            clusters: Domain → DomainCluster mapping from TraitClusterer
            variants: Optional full variant list (for stats)
            build: Genome build string (e.g., "GRCh38")
            drug_recommendations: Optional list of drug-gene interactions

        Returns:
            Rendered HTML string
        """
        if self._env is None:
            raise RuntimeError(
                f"Template directory not found: {self.template_dir}. "
                f"Make sure templates/ exists in the project root."
            )

        template = self._env.get_template(self.template_name)

        # Calculate stats
        total_variants = sum(c.count for c in clusters.values())
        domain_count = sum(1 for c in clusters.values() if c.count > 0)

        # Count notable findings (magnitude >= 2 or pathogenic)
        notable_count = 0
        source_counts = {"snpedia": 0, "clinvar": 0, "pharmcat": 0}

        for cluster in clusters.values():
            for v in cluster.variants:
                if v.has_snpedia:
                    source_counts["snpedia"] += 1
                    if v.snpedia.magnitude >= 2.0:
                        notable_count += 1
                if v.has_clinvar:
                    source_counts["clinvar"] += 1
                    if "pathogenic" in (v.clinvar.significance or "").lower():
                        notable_count += 1
                if v.has_pharmcat:
                    source_counts["pharmcat"] += 1

        # Serialize clusters for template (keyed by domain string name)
        template_clusters = {}
        for domain, cluster in clusters.items():
            key = domain.value if isinstance(domain, Domain) else str(domain)
            template_clusters[key] = cluster

        avatar_data = build_avatar(variants) if variants else {}
        logger.warning(f"AVATAR DICT: {avatar_data}")
        for v in variants:
            if v.rsid in ["rs12913832", "rs1815739"]:
                logger.warning(f"FOUND AVATAR VARIANT IN REPORTS LIST: rsid={v.rsid}, gt={v.genotype}")

        html = template.render(
            profile_name=profile_name,
            build=build,
            generated_at=datetime.now().strftime("%Y-%m-%d %H:%M"),
            total_variants=total_variants,
            domain_count=domain_count,
            notable_count=notable_count,
            source_counts=source_counts,
            clusters=template_clusters,
            domain_icons=DOMAIN_ICONS,
            drug_recommendations=drug_recommendations or [],
            avatar=avatar_data,
        )

        return html

    def save_html(self, html: str, output_path: str) -> str:
        """Save rendered HTML to a file. Returns the absolute path."""
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(html)
        logger.info(f"Report saved: {output_path}")
        return os.path.abspath(output_path)

    def generate_report(
        self,
        profile_name: str,
        variants: list[AnnotatedVariant],
        output_path: str,
        build: str = "",
        drug_recommendations: Optional[list[dict]] = None,
    ) -> str:
        """
        Full pipeline: cluster variants → render HTML → save to file.
        Returns the absolute path to the saved report.
        """
        # Cluster variants
        clusterer = TraitClusterer()
        clusters = clusterer.cluster(variants)

        # Render HTML
        html = self.render_html(
            profile_name=profile_name,
            clusters=clusters,
            variants=variants,
            build=build,
            drug_recommendations=drug_recommendations,
        )

        # Save
        return self.save_html(html, output_path)

    def generate_summary_json(
        self,
        profile_name: str,
        clusters: dict[Domain, DomainCluster],
    ) -> dict:
        """
        Generate a structured JSON summary of the report.
        Useful for programmatic consumption or API responses.
        """
        summary = {
            "profile_name": profile_name,
            "generated_at": datetime.now().isoformat(),
            "total_variants": sum(c.count for c in clusters.values()),
            "domains": {},
        }

        for domain, cluster in clusters.items():
            if cluster.count == 0:
                continue

            key = domain.value if isinstance(domain, Domain) else str(domain)
            summary["domains"][key] = {
                "display_name": cluster.display_name,
                "count": cluster.count,
                "top_findings": [],
            }

            for v in cluster.top_variants[:5]:
                finding = {
                    "rsid": v.rsid,
                    "genotype": v.genotype,
                }
                if v.snpedia:
                    finding["magnitude"] = v.snpedia.magnitude
                    finding["repute"] = v.snpedia.repute
                    finding["summary"] = v.snpedia.summary
                if v.clinvar:
                    finding["significance"] = v.clinvar.significance
                    finding["conditions"] = v.clinvar.conditions
                if v.pharmcat:
                    finding["drug_gene"] = v.pharmcat.gene
                    finding["phenotype"] = v.pharmcat.phenotype

                summary["domains"][key]["top_findings"].append(finding)

        return summary

    def save_summary_json(
        self,
        summary: dict,
        output_path: str,
    ) -> str:
        """Save the JSON summary to a file."""
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2, ensure_ascii=False)
        logger.info(f"JSON summary saved: {output_path}")
        return os.path.abspath(output_path)
