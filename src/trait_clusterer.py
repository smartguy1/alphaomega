"""
AlphaOmega — Trait Clusterer

Groups annotated variants into human-readable domains for reporting
and LLM narration. Uses keyword matching on SNPedia summaries, ClinVar
conditions, and known gene-to-domain mappings.
"""

import logging
import re
from collections import defaultdict
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

from src.annotator import AnnotatedVariant

logger = logging.getLogger(__name__)


class Domain(Enum):
    """High-level trait domains for grouping variants."""
    METABOLISM = "metabolism"
    DISEASE_RISK = "disease_risk"
    PHARMACOGENOMICS = "pharmacogenomics"
    PHYSICAL_TRAITS = "physical_traits"
    BEHAVIOR = "behavior"
    ANCESTRY = "ancestry"
    IMMUNITY = "immunity"
    UNCATEGORIZED = "uncategorized"


@dataclass
class DomainCluster:
    """A cluster of variants within a single domain."""
    domain: Domain
    display_name: str
    description: str
    variants: list[AnnotatedVariant] = field(default_factory=list)

    @property
    def count(self) -> int:
        return len(self.variants)

    @property
    def top_variants(self) -> list[AnnotatedVariant]:
        """Return the most significant variants in this cluster (top 10)."""
        return sorted(
            self.variants,
            key=lambda v: v.max_magnitude,
            reverse=True,
        )[:10]

    def to_dict(self) -> dict:
        return {
            "domain": self.domain.value,
            "display_name": self.display_name,
            "description": self.description,
            "count": self.count,
            "variants": [v.to_dict() for v in self.variants],
        }


# Domain display names and descriptions
DOMAIN_INFO = {
    Domain.METABOLISM: (
        "Metabolism",
        "How your body processes food, caffeine, alcohol, and other substances",
    ),
    Domain.DISEASE_RISK: (
        "Health & Disease Risk",
        "Genetic predispositions to health conditions — not predictions, but associations",
    ),
    Domain.PHARMACOGENOMICS: (
        "Drug Response",
        "How your body processes medications — relevant for dosing and drug selection",
    ),
    Domain.PHYSICAL_TRAITS: (
        "Physical Traits",
        "Observable characteristics influenced by genetics",
    ),
    Domain.BEHAVIOR: (
        "Behavioral Tendencies",
        "Genetic influences on sleep, mood, pain, and cognitive patterns",
    ),
    Domain.ANCESTRY: (
        "Ancestry & Population",
        "Markers related to population genetics and ancestral origins",
    ),
    Domain.IMMUNITY: (
        "Immunity & Immune System",
        "Immune system genetics including HLA type and autoimmune tendencies",
    ),
    Domain.UNCATEGORIZED: (
        "Other Findings",
        "Variants that don't fit neatly into the above categories",
    ),
}


# ── Classification rules ─────────────────────────────────────────────────────

# Known gene → domain mappings (pharmacogenomics genes)
PHARMGENES = {
    "CYP2D6", "CYP2C19", "CYP2C9", "CYP3A5", "CYP2B6", "CYP4F2",
    "CYP1A2", "CYP2C8", "CYP3A4",
    "DPYD", "NUDT15", "TPMT", "UGT1A1", "NAT2",
    "SLCO1B1", "ABCG2",
    "VKORC1",  # warfarin
    "RYR1", "CACNA1S",  # malignant hyperthermia
}

# HLA genes → immunity
HLA_GENES = {"HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1", "HLA-DPB1"}

# Keyword patterns for domain classification from summaries/conditions
DOMAIN_KEYWORDS = {
    Domain.METABOLISM: [
        r"caffeine", r"metaboliz", r"lactose", r"alcohol",
        r"folate", r"vitamin", r"iron", r"glucose", r"insulin",
        r"cholesterol", r"triglycerid", r"lipid", r"obesity",
        r"BMI", r"weight", r"fat\b", r"celiac", r"gluten",
        r"taste", r"bitter", r"sweet", r"umami",
        r"diet", r"nutrient", r"amino acid",
    ],
    Domain.DISEASE_RISK: [
        r"cancer", r"tumor", r"carcinoma", r"melanoma", r"leukemia",
        r"cardiovascular", r"heart", r"stroke", r"hypertension",
        r"diabetes", r"alzheimer", r"parkinson", r"dementia",
        r"asthma", r"COPD", r"fibrosis", r"cirrhosis",
        r"arthritis", r"osteoporosis",
        r"pathogenic", r"likely pathogenic",
        r"syndrome", r"disorder",
        r"risk", r"predispos", r"susceptib",
        r"macular degeneration", r"glaucoma",
        r"crohn", r"colitis", r"IBD",
        r"thromb", r"embolism", r"aneurysm",
    ],
    Domain.PHYSICAL_TRAITS: [
        r"eye color", r"hair color", r"skin color", r"freckling",
        r"height", r"tall", r"short",
        r"earwax", r"cerumen",
        r"muscle", r"ACTN3", r"sprint", r"endurance",
        r"baldness", r"alopecia",
        r"blood type", r"ABO",
        r"curly hair", r"straight hair",
        r"pain sensitiv",
        r"circadian", r"morning person", r"night owl",
    ],
    Domain.BEHAVIOR: [
        r"sleep", r"insomnia", r"circadian",
        r"nicotine", r"smoking", r"addiction",
        r"anxiety", r"depress", r"mood",
        r"COMT", r"warrior.{0,5}worrier",
        r"empathy", r"oxytocin", r"OXTR",
        r"dopamine", r"serotonin", r"MAOA",
        r"novelty seeking", r"risk.{0,5}taking",
        r"memory", r"cognitive", r"intelligence",
        r"ADHD", r"attention",
        r"autism", r"ASD",
        r"alcohol.{0,10}depend", r"substance",
    ],
    Domain.ANCESTRY: [
        r"haplogroup", r"mtDNA", r"Y.chromosome",
        r"population", r"ancestry",
        r"african", r"european", r"asian",
        r"neanderthal", r"denisovan",
        r"founder.{0,5}mutation",
    ],
    Domain.IMMUNITY: [
        r"immune", r"autoimmune", r"HLA",
        r"lupus", r"rheumatoid",
        r"psoriasis", r"eczema",
        r"allerg",
        r"inflammat",
        r"interleukin", r"cytokine", r"TNF",
        r"celiac",  # also metabolism
        r"multiple sclerosis",
        r"type.{0,3}1.{0,3}diabetes",
    ],
}


def classify_variant(variant: AnnotatedVariant) -> Domain:
    """
    Classify a single variant into a domain.
    Uses a priority-based approach:
      1. PharmCAT gene → PHARMACOGENOMICS
      2. HLA gene → IMMUNITY
      3. Keyword matching on summaries/conditions
      4. Fallback to UNCATEGORIZED
    """
    # Priority 1: Pharmacogenomics genes
    gene = _get_gene(variant)
    if gene and gene.upper() in PHARMGENES:
        return Domain.PHARMACOGENOMICS

    if variant.has_pharmcat:
        return Domain.PHARMACOGENOMICS

    # Priority 2: HLA/immunity genes
    if gene and gene.upper() in HLA_GENES:
        return Domain.IMMUNITY

    # Priority 3: Keyword matching
    text = _get_searchable_text(variant)
    if not text:
        return Domain.UNCATEGORIZED

    # Score each domain
    scores: dict[Domain, int] = defaultdict(int)
    for domain, patterns in DOMAIN_KEYWORDS.items():
        for pattern in patterns:
            if re.search(pattern, text, re.IGNORECASE):
                scores[domain] += 1

    if scores:
        # Return domain with highest score
        best_domain = max(scores, key=scores.get)
        return best_domain

    return Domain.UNCATEGORIZED


def _get_gene(variant: AnnotatedVariant) -> str:
    """Extract gene name from variant annotations."""
    if variant.clinvar and variant.clinvar.gene:
        return variant.clinvar.gene
    if variant.pharmcat and variant.pharmcat.gene:
        return variant.pharmcat.gene
    return ""


def _get_searchable_text(variant: AnnotatedVariant) -> str:
    """Build a searchable text string from all annotations."""
    parts = []

    if variant.snpedia:
        parts.append(variant.snpedia.summary)
        parts.append(variant.snpedia.details)

    if variant.clinvar:
        parts.append(variant.clinvar.significance)
        parts.extend(variant.clinvar.conditions)
        parts.append(variant.clinvar.gene)

    if variant.pharmcat:
        parts.append(variant.pharmcat.gene)
        parts.append(variant.pharmcat.phenotype)
        parts.extend(variant.pharmcat.drugs)

    return " ".join(filter(None, parts))


class TraitClusterer:
    """
    Groups annotated variants into domain clusters.

    Usage:
        clusterer = TraitClusterer()
        clusters = clusterer.cluster(annotated_variants)
        summary = clusterer.get_summary(clusters)
    """

    def cluster(
        self, variants: list[AnnotatedVariant]
    ) -> dict[Domain, DomainCluster]:
        """
        Classify and group variants into domain clusters.
        Returns a dict of Domain → DomainCluster.
        """
        clusters: dict[Domain, DomainCluster] = {}

        # Initialize all clusters
        for domain in Domain:
            name, desc = DOMAIN_INFO[domain]
            clusters[domain] = DomainCluster(
                domain=domain,
                display_name=name,
                description=desc,
            )

        # Classify each variant
        for variant in variants:
            domain = classify_variant(variant)
            clusters[domain].variants.append(variant)

        logger.info(
            "Clustering results: " +
            ", ".join(f"{d.value}={c.count}" for d, c in clusters.items() if c.count > 0)
        )

        return clusters

    def cluster_nonempty(
        self, variants: list[AnnotatedVariant]
    ) -> dict[Domain, DomainCluster]:
        """Same as cluster() but only returns non-empty domains."""
        all_clusters = self.cluster(variants)
        return {d: c for d, c in all_clusters.items() if c.count > 0}

    def get_summary(
        self, clusters: dict[Domain, DomainCluster]
    ) -> dict:
        """
        Generate a summary of all clusters.
        Returns a dict suitable for JSON serialization.
        """
        summary = {
            "total_variants": sum(c.count for c in clusters.values()),
            "domains": {},
        }

        for domain, cluster in clusters.items():
            if cluster.count == 0:
                continue

            summary["domains"][domain.value] = {
                "display_name": cluster.display_name,
                "description": cluster.description,
                "count": cluster.count,
                "top_findings": [
                    {
                        "rsid": v.rsid,
                        "genotype": v.genotype,
                        "magnitude": v.max_magnitude,
                        "summary": (
                            v.snpedia.summary if v.snpedia else
                            v.clinvar.significance if v.clinvar else ""
                        ),
                    }
                    for v in cluster.top_variants[:5]
                ],
            }

        return summary

    def to_json(
        self, clusters: dict[Domain, DomainCluster]
    ) -> str:
        """Serialize clusters to JSON."""
        import json
        data = {}
        for domain, cluster in clusters.items():
            if cluster.count > 0:
                data[domain.value] = cluster.to_dict()
        return json.dumps(data, indent=2, ensure_ascii=False)
