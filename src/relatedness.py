"""
AlphaOmega — Relatedness Analyzer

Verifies filiation (parent-child, siblings, etc.) between two genomic
profiles using Identity-by-State (IBS) analysis.

Method:
  1. Find overlapping SNPs between two profiles
  2. For each shared SNP, count alleles shared (0, 1, or 2)
  3. Calculate IBS0/IBS1/IBS2 proportions
  4. Compare against known distributions for each relationship type
  5. Apply Mendelian error check for parent-child claims

Expected distributions:
  Parent-child:  IBS0 ≈ 0%, IBS1 ≈ 100%, IBS2 ≈ 0%  (must share ≥1 allele)
  Full siblings:  IBS0 ≈ 25%, IBS1 ≈ 50%, IBS2 ≈ 25%
  Half siblings:  IBS0 ≈ 50%, IBS1 ≈ 50%, IBS2 ≈ 0%
  Unrelated:      IBS0 ≈ depends on allele frequencies, typically higher
"""

import logging
import sqlite3
from dataclasses import dataclass
from enum import Enum
from typing import Optional

logger = logging.getLogger(__name__)


class Relationship(Enum):
    PARENT_CHILD = "parent_child"
    FULL_SIBLINGS = "full_siblings"
    HALF_SIBLINGS = "half_siblings"
    SECOND_DEGREE = "second_degree"
    UNRELATED = "unrelated"
    SELF_OR_TWIN = "self_or_twin"
    INCONCLUSIVE = "inconclusive"


@dataclass
class RelatednessResult:
    """Result of relatedness analysis between two profiles."""
    profile_a: str
    profile_b: str
    relationship: Relationship
    confidence: float
    shared_snps: int
    ibs0: float       # Proportion sharing 0 alleles
    ibs1: float       # Proportion sharing 1 allele
    ibs2: float       # Proportion sharing 2 alleles
    mendelian_errors: int
    mendelian_error_rate: float
    details: str = ""

    def to_dict(self) -> dict:
        return {
            "profile_a": self.profile_a,
            "profile_b": self.profile_b,
            "relationship": self.relationship.value,
            "confidence": round(self.confidence, 3),
            "shared_snps": self.shared_snps,
            "ibs0": round(self.ibs0, 4),
            "ibs1": round(self.ibs1, 4),
            "ibs2": round(self.ibs2, 4),
            "mendelian_errors": self.mendelian_errors,
            "mendelian_error_rate": round(self.mendelian_error_rate, 4),
            "details": self.details,
        }

    @property
    def summary(self) -> str:
        rel_names = {
            Relationship.PARENT_CHILD: "Parent–Child",
            Relationship.FULL_SIBLINGS: "Full Siblings",
            Relationship.HALF_SIBLINGS: "Half Siblings",
            Relationship.SECOND_DEGREE: "Second-Degree Relatives",
            Relationship.UNRELATED: "Unrelated",
            Relationship.SELF_OR_TWIN: "Self/Identical Twin",
            Relationship.INCONCLUSIVE: "Inconclusive",
        }
        name = rel_names.get(self.relationship, self.relationship.value)
        return (
            f"{self.profile_a} ↔ {self.profile_b}: "
            f"{name} (confidence={self.confidence:.1%}, "
            f"shared_snps={self.shared_snps})"
        )


# Minimum SNPs needed for reliable analysis
MIN_SHARED_SNPS = 100

# Thresholds for relationship classification (based on IBS proportions)
# These are approximate and work well with >1000 shared SNPs
RELATIONSHIP_THRESHOLDS = {
    # (ibs0_max, ibs1_min, ibs1_max, ibs2_min)
    Relationship.SELF_OR_TWIN:   (0.02, 0.0, 0.10, 0.88),
    Relationship.PARENT_CHILD:   (0.03, 0.85, 1.00, 0.0),
    Relationship.FULL_SIBLINGS:  (0.35, 0.35, 0.65, 0.15),
    Relationship.HALF_SIBLINGS:  (0.60, 0.35, 0.65, 0.0),
}


def _count_shared_alleles(genotype_a: str, genotype_b: str) -> int:
    """
    Count the number of alleles shared between two genotypes.
    Returns 0, 1, or 2.

    Examples:
      AA vs AA → 2 (both alleles match)
      AA vs AG → 1 (one allele matches)
      AA vs GG → 0 (no alleles match)
      AG vs AG → 2 (both alleles match)
      AG vs AC → 1 (A matches)
    """
    if not genotype_a or not genotype_b:
        return -1  # Missing data

    # Normalize: ensure exactly 2 characters
    a = genotype_a.upper().strip()
    b = genotype_b.upper().strip()

    if len(a) < 2 or len(b) < 2:
        return -1

    a1, a2 = a[0], a[1]
    b1, b2 = b[0], b[1]

    # Check all pairing combinations
    if (a1 == b1 and a2 == b2) or (a1 == b2 and a2 == b1):
        return 2  # Both alleles shared

    if a1 in (b1, b2) or a2 in (b1, b2):
        return 1  # One allele shared

    return 0  # No alleles shared


def _check_mendelian(genotype_parent: str, genotype_child: str) -> bool:
    """
    Check if a parent-child genotype pair is Mendelian-consistent.
    A child MUST inherit one allele from each parent.
    Returns True if consistent, False if Mendelian error.
    """
    if not genotype_parent or not genotype_child:
        return True  # Can't check, assume consistent

    p = genotype_parent.upper().strip()
    c = genotype_child.upper().strip()

    if len(p) < 2 or len(c) < 2:
        return True

    p1, p2 = p[0], p[1]
    c1, c2 = c[0], c[1]

    # Child must have at least one allele from parent
    return c1 in (p1, p2) or c2 in (p1, p2)


def analyze_relatedness(
    profile_a_name: str,
    profile_b_name: str,
    db_path: str,
) -> RelatednessResult:
    """
    Analyze relatedness between two profiles using their stored annotations.

    Args:
        profile_a_name: First profile name
        profile_b_name: Second profile name
        db_path: Path to the profiles SQLite database

    Returns:
        RelatednessResult with relationship classification and statistics
    """
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row

    # Get profile IDs
    row_a = conn.execute(
        "SELECT id FROM profiles WHERE name=?", (profile_a_name,)
    ).fetchone()
    row_b = conn.execute(
        "SELECT id FROM profiles WHERE name=?", (profile_b_name,)
    ).fetchone()

    if row_a is None or row_b is None:
        conn.close()
        missing = profile_a_name if row_a is None else profile_b_name
        return RelatednessResult(
            profile_a=profile_a_name, profile_b=profile_b_name,
            relationship=Relationship.INCONCLUSIVE, confidence=0.0,
            shared_snps=0, ibs0=0, ibs1=0, ibs2=0,
            mendelian_errors=0, mendelian_error_rate=0,
            details=f"Profile '{missing}' not found",
        )

    id_a, id_b = row_a["id"], row_b["id"]

    # Find overlapping SNPs with genotypes in both profiles
    rows = conn.execute(
        """
        SELECT a.rsid, a.genotype AS gt_a, b.genotype AS gt_b
        FROM profile_annotations a
        JOIN profile_annotations b
            ON a.rsid = b.rsid
        WHERE a.profile_id = ? AND b.profile_id = ?
            AND a.genotype IS NOT NULL AND a.genotype != ''
            AND b.genotype IS NOT NULL AND b.genotype != ''
        """,
        (id_a, id_b),
    ).fetchall()

    conn.close()

    shared_snps = len(rows)

    if shared_snps < MIN_SHARED_SNPS:
        return RelatednessResult(
            profile_a=profile_a_name, profile_b=profile_b_name,
            relationship=Relationship.INCONCLUSIVE, confidence=0.0,
            shared_snps=shared_snps, ibs0=0, ibs1=0, ibs2=0,
            mendelian_errors=0, mendelian_error_rate=0,
            details=f"Insufficient shared SNPs ({shared_snps} < {MIN_SHARED_SNPS})",
        )

    # Count IBS states
    ibs_counts = {0: 0, 1: 0, 2: 0}
    mendelian_errors = 0

    for row in rows:
        shared = _count_shared_alleles(row["gt_a"], row["gt_b"])
        if shared < 0:
            continue  # Skip missing data
        ibs_counts[shared] += 1

        if not _check_mendelian(row["gt_a"], row["gt_b"]):
            mendelian_errors += 1

    total_compared = sum(ibs_counts.values())
    if total_compared == 0:
        return RelatednessResult(
            profile_a=profile_a_name, profile_b=profile_b_name,
            relationship=Relationship.INCONCLUSIVE, confidence=0.0,
            shared_snps=shared_snps, ibs0=0, ibs1=0, ibs2=0,
            mendelian_errors=0, mendelian_error_rate=0,
            details="No valid genotype comparisons could be made",
        )

    ibs0 = ibs_counts[0] / total_compared
    ibs1 = ibs_counts[1] / total_compared
    ibs2 = ibs_counts[2] / total_compared
    merr_rate = mendelian_errors / total_compared

    # Classify relationship
    relationship, confidence, details = _classify_relationship(
        ibs0, ibs1, ibs2, merr_rate, total_compared,
    )

    return RelatednessResult(
        profile_a=profile_a_name,
        profile_b=profile_b_name,
        relationship=relationship,
        confidence=confidence,
        shared_snps=shared_snps,
        ibs0=ibs0,
        ibs1=ibs1,
        ibs2=ibs2,
        mendelian_errors=mendelian_errors,
        mendelian_error_rate=merr_rate,
        details=details,
    )


def _classify_relationship(
    ibs0: float, ibs1: float, ibs2: float,
    merr_rate: float, n_snps: int,
) -> tuple[Relationship, float, str]:
    """
    Classify the relationship based on IBS proportions.
    Returns (relationship, confidence, details).
    """

    # Self or identical twin: nearly all IBS2
    if ibs2 > 0.88 and ibs0 < 0.02:
        conf = min(1.0, ibs2)
        return (
            Relationship.SELF_OR_TWIN, conf,
            f"IBS2={ibs2:.1%} indicates same individual or identical twin",
        )

    # Parent-child: nearly all IBS1, very low IBS0
    if ibs0 < 0.03 and ibs1 > 0.70:
        # Mendelian error rate should also be very low for true parent-child
        if merr_rate < 0.02:
            conf = min(1.0, (1.0 - ibs0) * (1.0 - merr_rate))
            return (
                Relationship.PARENT_CHILD, conf,
                f"IBS0={ibs0:.1%}, IBS1={ibs1:.1%}, "
                f"Mendelian errors={merr_rate:.2%} — consistent with parent-child",
            )
        else:
            conf = 0.5
            return (
                Relationship.PARENT_CHILD, conf,
                f"IBS pattern suggests parent-child but Mendelian error rate "
                f"is elevated ({merr_rate:.2%}). May indicate genotyping errors.",
            )

    # Full siblings: roughly equal IBS0/IBS1/IBS2
    if 0.15 < ibs0 < 0.35 and 0.35 < ibs1 < 0.65 and ibs2 > 0.15:
        conf = min(1.0, 1.0 - abs(ibs0 - 0.25) - abs(ibs1 - 0.50) - abs(ibs2 - 0.25))
        return (
            Relationship.FULL_SIBLINGS, max(conf, 0.5),
            f"IBS0={ibs0:.1%}, IBS1={ibs1:.1%}, IBS2={ibs2:.1%} — "
            f"consistent with full siblings",
        )

    # Half siblings: higher IBS0, moderate IBS1, low IBS2
    if 0.10 < ibs0 < 0.60 and 0.35 < ibs1 < 0.65 and ibs2 < 0.15:
        conf = 0.6
        return (
            Relationship.HALF_SIBLINGS, conf,
            f"IBS0={ibs0:.1%}, IBS1={ibs1:.1%}, IBS2={ibs2:.1%} — "
            f"consistent with half siblings or avuncular relationship",
        )

    # Second degree: some sharing but less than siblings
    if 0.15 < ibs0 < 0.50 and ibs1 > 0.30:
        return (
            Relationship.SECOND_DEGREE, 0.4,
            f"IBS pattern suggests second-degree relatives "
            f"(grandparent, uncle/aunt, nephew/niece)",
        )

    # Default: unrelated
    conf = min(1.0, ibs0 * 2)  # Higher IBS0 = more confident unrelated
    return (
        Relationship.UNRELATED, max(conf, 0.3),
        f"IBS0={ibs0:.1%} — no close familial relationship detected",
    )


def analyze_relatedness_from_variants(
    profile_a_name: str,
    profile_b_name: str,
    variants_a: list[dict],
    variants_b: list[dict],
) -> RelatednessResult:
    """
    Analyze relatedness from pre-loaded annotation dicts (in-memory).
    Useful for testing without a database.

    Each variant dict must have 'rsid' and 'genotype' keys.
    """
    # Build lookup
    lookup_b = {v["rsid"]: v["genotype"] for v in variants_b if v.get("genotype")}

    ibs_counts = {0: 0, 1: 0, 2: 0}
    mendelian_errors = 0
    shared_snps = 0

    for va in variants_a:
        rsid = va.get("rsid", "")
        gt_a = va.get("genotype", "")
        if not gt_a or rsid not in lookup_b:
            continue

        gt_b = lookup_b[rsid]
        shared_snps += 1

        shared = _count_shared_alleles(gt_a, gt_b)
        if shared >= 0:
            ibs_counts[shared] += 1
            if not _check_mendelian(gt_a, gt_b):
                mendelian_errors += 1

    total = sum(ibs_counts.values())
    if total < MIN_SHARED_SNPS:
        return RelatednessResult(
            profile_a=profile_a_name, profile_b=profile_b_name,
            relationship=Relationship.INCONCLUSIVE, confidence=0.0,
            shared_snps=shared_snps, ibs0=0, ibs1=0, ibs2=0,
            mendelian_errors=0, mendelian_error_rate=0,
            details=f"Insufficient shared SNPs ({total} < {MIN_SHARED_SNPS})",
        )

    ibs0 = ibs_counts[0] / total
    ibs1 = ibs_counts[1] / total
    ibs2 = ibs_counts[2] / total
    merr_rate = mendelian_errors / total

    relationship, confidence, details = _classify_relationship(
        ibs0, ibs1, ibs2, merr_rate, total,
    )

    return RelatednessResult(
        profile_a=profile_a_name, profile_b=profile_b_name,
        relationship=relationship, confidence=confidence,
        shared_snps=shared_snps, ibs0=ibs0, ibs1=ibs1, ibs2=ibs2,
        mendelian_errors=mendelian_errors, mendelian_error_rate=merr_rate,
        details=details,
    )
