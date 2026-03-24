"""
AlphaOmega — FragmentOfUs Export

Thin adapter that packages AlphaOmega genomic analysis into a format
consumable by the FragmentOfUs AI biographer project.

Exports:
  - Domain-clustered genomic highlights
  - Confidence metadata per finding
  - Pre-structured narrative fragments
  - Profile metadata
"""

import json
import logging
import os
from datetime import datetime
from typing import Optional

from src.annotator import AnnotatedVariant
from src.trait_clusterer import Domain, DomainCluster, TraitClusterer
from src.profile_manager import ProfileManager

logger = logging.getLogger(__name__)


def export_for_biographer(
    profile_name: str,
    variants: list[AnnotatedVariant],
    clusters: Optional[dict[Domain, DomainCluster]] = None,
    profile_manager: Optional[ProfileManager] = None,
) -> dict:
    """
    Export genomic analysis for the FragmentOfUs AI biographer.

    Returns a structured dict with:
      - profile metadata
      - domain-clustered highlights
      - confidence-tagged findings
      - suggested narrative hooks

    Args:
        profile_name: Name of the profile
        variants: Annotated variants
        clusters: Pre-computed clusters (will generate if not provided)
        profile_manager: Optional ProfileManager for metadata
    """
    # Cluster if not provided
    if clusters is None:
        clusterer = TraitClusterer()
        clusters = clusterer.cluster_nonempty(variants)

    # Gather profile metadata
    profile_meta = {
        "name": profile_name,
        "exported_at": datetime.now().isoformat(),
        "total_variants_analyzed": len(variants),
    }

    if profile_manager:
        profile = profile_manager.get_profile(profile_name)
        if profile:
            profile_meta["genome_build"] = profile.detected_build
            profile_meta["created_at"] = profile.created_at
            profile_meta["last_analyzed"] = profile.last_analyzed

    # Build highlights per domain
    genomic_highlights = {}

    for domain, cluster in clusters.items():
        if cluster.count == 0:
            continue

        key = domain.value if isinstance(domain, Domain) else str(domain)

        findings = []
        for v in cluster.top_variants[:10]:
            finding = _build_finding(v)
            if finding:
                findings.append(finding)

        genomic_highlights[key] = {
            "display_name": cluster.display_name,
            "description": cluster.description,
            "finding_count": cluster.count,
            "findings": findings,
            "narrative_hooks": _generate_hooks(key, findings),
        }

    return {
        "schema_version": "1.0",
        "source": "AlphaOmega",
        "profile": profile_meta,
        "genomic_highlights": genomic_highlights,
        "disclaimer": (
            "This genomic data is for informational purposes only. "
            "It is not a medical diagnosis. Associations are based on "
            "published research and may not apply to every individual."
        ),
    }


def _build_finding(variant: AnnotatedVariant) -> Optional[dict]:
    """Build a single finding dict from an annotated variant."""
    finding = {
        "rsid": variant.rsid,
        "genotype": variant.genotype,
    }

    # Determine confidence level
    confidence = "low"
    source_count = variant.source_count

    if source_count >= 2:
        confidence = "high"
    elif source_count == 1:
        if variant.snpedia and variant.snpedia.magnitude >= 3.0:
            confidence = "high"
        elif variant.clinvar and variant.clinvar.review_stars >= 3:
            confidence = "high"
        elif variant.snpedia and variant.snpedia.magnitude >= 2.0:
            confidence = "moderate"
        elif variant.clinvar and variant.clinvar.review_stars >= 1:
            confidence = "moderate"
        else:
            confidence = "low"

    finding["confidence"] = confidence

    # Summary
    if variant.snpedia and variant.snpedia.summary:
        finding["summary"] = variant.snpedia.summary
        finding["repute"] = variant.snpedia.repute

    if variant.clinvar:
        finding["clinical_significance"] = variant.clinvar.significance
        if variant.clinvar.conditions:
            finding["conditions"] = variant.clinvar.conditions

    if variant.pharmcat:
        finding["drug_response"] = {
            "gene": variant.pharmcat.gene,
            "phenotype": variant.pharmcat.phenotype,
            "drugs": variant.pharmcat.drugs,
        }

    return finding


# Narrative hook templates per domain — these give the biographer
# starting points for weaving genetics into a life narrative.
NARRATIVE_HOOKS = {
    "metabolism": [
        "Your body's relationship with food and substances is partly written in your DNA.",
        "These metabolic traits influence how you experience everyday substances like caffeine and food.",
    ],
    "disease_risk": [
        "Genetics load the gun, but lifestyle pulls the trigger.",
        "These associations represent predispositions, not certainties.",
    ],
    "pharmacogenomics": [
        "Your DNA can influence how your body responds to medications.",
        "This information may be valuable to share with your healthcare provider.",
    ],
    "physical_traits": [
        "Some of your physical characteristics have roots deep in your genetic code.",
        "These traits connect you to the long chain of ancestors who carried similar variants.",
    ],
    "behavior": [
        "Our behavioral tendencies have both genetic and environmental components.",
        "These genetic influences are just one part of what makes you who you are.",
    ],
    "ancestry": [
        "Your DNA carries echoes of ancient migrations and population histories.",
        "These markers connect you to broader human population patterns.",
    ],
    "immunity": [
        "Your immune system's blueprint is partially encoded in your HLA genes.",
        "These immune system variants reflect your body's unique defense strategy.",
    ],
}


def _generate_hooks(domain_key: str, findings: list[dict]) -> list[str]:
    """Generate narrative hooks for the biographer based on domain and findings."""
    hooks = list(NARRATIVE_HOOKS.get(domain_key, []))

    # Add finding-specific hooks
    for f in findings[:3]:
        summary = f.get("summary", "")
        if summary and f.get("confidence") in ("high", "moderate"):
            hooks.append(f"Notably: {summary}")

    return hooks


def save_export(export_data: dict, output_path: str) -> str:
    """Save the export to a JSON file."""
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(export_data, f, indent=2, ensure_ascii=False)
    logger.info(f"FragmentOfUs export saved: {output_path}")
    return os.path.abspath(output_path)
