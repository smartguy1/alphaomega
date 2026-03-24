"""
AlphaOmega — Core Annotation Engine

Merges annotations from SNPedia, ClinVar, and PharmCAT into a unified
structure per variant. Reads a VCF file and produces annotated JSON
ready for trait clustering and reporting.
"""

import json
import logging
import os
import subprocess
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

from src.snpedia_mirror import SNPediaMirror
from src.clinvar_mirror import ClinVarMirror
from src.pharmcat_runner import PharmCATRunner, PharmCATResult

logger = logging.getLogger(__name__)


@dataclass
class SNPediaAnnotation:
    """SNPedia annotation for a variant."""
    magnitude: float = 0.0
    repute: str = ""
    summary: str = ""
    details: str = ""


@dataclass
class ClinVarAnnotation:
    """ClinVar annotation for a variant."""
    significance: str = ""
    conditions: list[str] = field(default_factory=list)
    review_stars: int = 0
    gene: str = ""


@dataclass
class PharmCATAnnotation:
    """PharmCAT annotation for a variant."""
    gene: str = ""
    phenotype: str = ""
    drugs: list[str] = field(default_factory=list)
    recommendation: str = ""


@dataclass
class AnnotatedVariant:
    """
    A single variant with merged annotations from all sources.
    This is the core data structure consumed by the trait clusterer
    and report generator.
    """
    rsid: str
    chromosome: str = ""
    position: int = 0
    genotype: str = ""
    ref: str = ""
    alt: str = ""
    snpedia: Optional[SNPediaAnnotation] = None
    clinvar: Optional[ClinVarAnnotation] = None
    pharmcat: Optional[PharmCATAnnotation] = None

    @property
    def has_snpedia(self) -> bool:
        return self.snpedia is not None and bool(self.snpedia.summary or self.snpedia.magnitude > 0)

    @property
    def has_clinvar(self) -> bool:
        return self.clinvar is not None and bool(self.clinvar.significance)

    @property
    def has_pharmcat(self) -> bool:
        return self.pharmcat is not None and bool(self.pharmcat.gene)

    @property
    def source_count(self) -> int:
        return sum([self.has_snpedia, self.has_clinvar, self.has_pharmcat])

    @property
    def max_magnitude(self) -> float:
        mag = 0.0
        if self.snpedia:
            mag = max(mag, self.snpedia.magnitude)
        if self.clinvar and self.clinvar.review_stars:
            mag = max(mag, self.clinvar.review_stars)
        return mag

    def to_dict(self) -> dict:
        d = {
            "rsid": self.rsid,
            "chromosome": self.chromosome,
            "position": self.position,
            "genotype": self.genotype,
            "sources": {},
        }
        if self.snpedia:
            d["sources"]["snpedia"] = asdict(self.snpedia)
        if self.clinvar:
            d["sources"]["clinvar"] = asdict(self.clinvar)
        if self.pharmcat:
            d["sources"]["pharmcat"] = asdict(self.pharmcat)
        return d


def _parse_vcf_genotypes(
    vcf_path: str,
    position_to_rsid: Optional[dict[tuple[str, int], str]] = None,
) -> list[dict]:
    """
    Parse a VCF file and extract rsID + genotype for each variant.
    Returns list of {rsid, chromosome, position, ref, alt, genotype}.
    Uses bcftools query for efficiency, falls back to manual parsing.

    If position_to_rsid is provided, variants without rsIDs will be
    looked up by (chromosome, position) to resolve their rsID.
    """
    variants = []

    try:
        result = subprocess.run(
            [
                "bcftools", "query",
                "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n",
                vcf_path,
            ],
            capture_output=True, text=True, timeout=300,
        )

        if result.returncode == 0:
            for line in result.stdout.strip().splitlines():
                parts = line.split("\t")
                if len(parts) < 6:
                    continue

                chrom, pos, rsid, ref, alt, gt = parts[:6]

                # If rsid is missing, try position-based lookup
                if not rsid.startswith("rs") and position_to_rsid:
                    chrom_norm = chrom.replace("chr", "")
                    try:
                        pos_int = int(pos)
                        rsid = position_to_rsid.get(
                            (chrom_norm, pos_int), rsid
                        )
                    except ValueError:
                        pass

                if not rsid.startswith("rs"):
                    continue

                genotype = _gt_to_alleles(gt, ref, alt)
                variants.append({
                    "rsid": rsid.lower(),
                    "chromosome": chrom.replace("chr", ""),
                    "position": int(pos),
                    "ref": ref,
                    "alt": alt,
                    "genotype": genotype,
                })
            return variants

    except (FileNotFoundError, subprocess.TimeoutExpired):
        logger.info("bcftools not available, falling back to manual VCF parsing")

    # Fallback: manual parsing
    return _parse_vcf_manual(vcf_path, position_to_rsid)


def _parse_vcf_manual(
    vcf_path: str,
    position_to_rsid: Optional[dict[tuple[str, int], str]] = None,
) -> list[dict]:
    """Manual VCF parser fallback when bcftools is not available."""
    import gzip

    variants = []
    open_func = gzip.open if str(vcf_path).endswith(".gz") else open

    try:
        with open_func(vcf_path, "rt", errors="replace") as f:
            for line in f:
                if line.startswith("#"):
                    continue

                parts = line.strip().split("\t")
                if len(parts) < 10:
                    continue

                chrom, pos, rsid, ref, alt = parts[:5]

                # If rsid is missing, try position-based lookup
                if not rsid.startswith("rs") and position_to_rsid:
                    chrom_norm = chrom.replace("chr", "")
                    try:
                        pos_int = int(pos)
                        rsid = position_to_rsid.get(
                            (chrom_norm, pos_int), rsid
                        )
                    except ValueError:
                        pass

                if not rsid.startswith("rs"):
                    continue

                # Parse genotype from sample column
                fmt = parts[8].split(":")
                sample = parts[9].split(":")
                gt_idx = fmt.index("GT") if "GT" in fmt else 0
                gt = sample[gt_idx] if gt_idx < len(sample) else "."

                genotype = _gt_to_alleles(gt, ref, alt)
                variants.append({
                    "rsid": rsid.lower(),
                    "chromosome": chrom.replace("chr", ""),
                    "position": int(pos),
                    "ref": ref,
                    "alt": alt,
                    "genotype": genotype,
                })

    except Exception as e:
        logger.error(f"Error parsing VCF: {e}")

    return variants


def _gt_to_alleles(gt_str: str, ref: str, alt: str) -> str:
    """
    Convert VCF GT field to allele string.
    E.g., "0/1" with ref=A, alt=G → "AG"
    """
    alts = alt.split(",")
    alleles = [ref] + alts

    gt_clean = gt_str.replace("|", "/").replace(".", "0")
    indices = gt_clean.split("/")

    result = ""
    for idx in indices:
        try:
            i = int(idx)
            if i < len(alleles):
                result += alleles[i]
        except ValueError:
            result += "?"

    return result


class Annotator:
    """
    Core annotation engine that merges data from multiple sources.

    Usage:
        annotator = Annotator(
            snpedia=SNPediaMirror("db"),
            clinvar=ClinVarMirror("db"),
        )
        variants = annotator.annotate_vcf("path/to/variants.vcf.gz")
    """

    def __init__(
        self,
        snpedia: Optional[SNPediaMirror] = None,
        clinvar: Optional[ClinVarMirror] = None,
        pharmcat_result: Optional[PharmCATResult] = None,
        min_magnitude: float = 0.0,
        min_review_stars: int = 0,
        include_unannotated: bool = False,
        force_include_rsids: Optional[set[str]] = None,
    ):
        self.snpedia = snpedia
        self.clinvar = clinvar
        self.pharmcat_result = pharmcat_result
        self.min_magnitude = min_magnitude
        self.min_review_stars = min_review_stars
        self.include_unannotated = include_unannotated
        self.force_include_rsids = force_include_rsids or set()

        # Build position → rsid index from SNPedia for resolving
        # variants that bcftools outputs without rsIDs
        self._position_index: dict[tuple[str, int], str] = {}
        if snpedia:
            try:
                self._position_index = snpedia.build_position_index()
            except Exception as e:
                logger.warning(f"Failed to build position index: {e}")

        # Build pharmcat lookup index by gene
        self._pharmcat_by_gene: dict[str, dict] = {}
        if pharmcat_result and pharmcat_result.success:
            for gr in pharmcat_result.gene_results:
                self._pharmcat_by_gene[gr.gene.upper()] = {
                    "gene": gr.gene,
                    "diplotype": gr.diplotype,
                    "phenotype": gr.phenotype,
                    "activity_score": gr.activity_score,
                }

            for dr in pharmcat_result.drug_recommendations:
                gene_key = dr.gene.upper()
                if gene_key in self._pharmcat_by_gene:
                    self._pharmcat_by_gene[gene_key].setdefault("drugs", []).append({
                        "drug": dr.drug,
                        "recommendation": dr.recommendation,
                        "classification": dr.classification,
                    })

    def annotate_vcf(self, vcf_path: str) -> list[AnnotatedVariant]:
        """
        Annotate all variants in a VCF file.
        Returns list of AnnotatedVariant sorted by source count (desc) then magnitude (desc).
        """
        logger.info(f"Parsing VCF: {vcf_path}")
        raw_variants = _parse_vcf_genotypes(
            vcf_path, position_to_rsid=self._position_index or None
        )
        logger.info(f"Found {len(raw_variants)} variants with rsIDs")

        annotated = []
        for var in raw_variants:
            av = self.annotate_single(
                rsid=var["rsid"],
                genotype=var["genotype"],
                chromosome=var["chromosome"],
                position=var["position"],
                ref=var["ref"],
                alt=var["alt"],
            )

            if av is not None:
                annotated.append(av)

        # Sort: most-annotated first, then by magnitude
        annotated.sort(
            key=lambda v: (v.source_count, v.max_magnitude),
            reverse=True,
        )

        logger.info(f"Annotated {len(annotated)} variants "
                     f"(from {len(raw_variants)} total)")
        return annotated

    def annotate_single(
        self,
        rsid: str,
        genotype: str = "",
        chromosome: str = "",
        position: int = 0,
        ref: str = "",
        alt: str = "",
    ) -> Optional[AnnotatedVariant]:
        """
        Annotate a single variant by looking up in all sources.
        Returns None if no annotations found and include_unannotated is False.
        """
        rsid = rsid.lower().strip()
        av = AnnotatedVariant(
            rsid=rsid,
            chromosome=chromosome,
            position=position,
            genotype=genotype,
            ref=ref,
            alt=alt,
        )

        # SNPedia lookup
        if self.snpedia:
            try:
                if genotype:
                    gt_result = self.snpedia.lookup_genotype(rsid, genotype)
                    if gt_result:
                        av.snpedia = SNPediaAnnotation(
                            magnitude=gt_result.magnitude,
                            repute=gt_result.repute,
                            summary=gt_result.summary,
                            details=gt_result.details if hasattr(gt_result, 'details') else "",
                        )
                else:
                    entry = self.snpedia.lookup(rsid)
                    if entry and entry.summary:
                        av.snpedia = SNPediaAnnotation(
                            summary=entry.summary,
                        )
            except Exception as e:
                logger.debug(f"SNPedia lookup failed for {rsid}: {e}")

        # ClinVar lookup
        if self.clinvar:
            try:
                entries = self.clinvar.lookup(rsid)
                if entries:
                    # Use the entry with highest review stars
                    best = max(entries, key=lambda e: e.review_stars)
                    av.clinvar = ClinVarAnnotation(
                        significance=best.clinical_significance,
                        conditions=best.conditions,
                        review_stars=best.review_stars,
                        gene=best.gene,
                    )
            except Exception as e:
                logger.debug(f"ClinVar lookup failed for {rsid}: {e}")

        # PharmCAT lookup (by gene from ClinVar annotation)
        gene = ""
        if av.clinvar and av.clinvar.gene:
            gene = av.clinvar.gene.upper()
        if gene and gene in self._pharmcat_by_gene:
            pg = self._pharmcat_by_gene[gene]
            drugs = [d["drug"] for d in pg.get("drugs", [])]
            av.pharmcat = PharmCATAnnotation(
                gene=pg["gene"],
                phenotype=pg.get("phenotype", ""),
                drugs=drugs,
                recommendation=pg.get("drugs", [{}])[0].get("recommendation", "") if pg.get("drugs") else "",
            )

        # Apply filters
        # If it's a forced inclusion rsID, never filter it out
        if self.force_include_rsids and av.rsid in self.force_include_rsids:
            return av
            
        if not self.include_unannotated and av.source_count == 0:
            return None

        if av.snpedia and av.snpedia.magnitude < self.min_magnitude and not av.has_clinvar:
            if not self.include_unannotated:
                return None

        if av.clinvar and av.clinvar.review_stars < self.min_review_stars and not av.has_snpedia:
            if not self.include_unannotated:
                return None

        return av

    def annotate_rsids(self, rsids: list[str]) -> list[AnnotatedVariant]:
        """
        Annotate a list of rsIDs without VCF context (no genotype info).
        Useful for quick lookups.
        """
        results = []
        for rsid in rsids:
            av = self.annotate_single(rsid=rsid)
            if av is not None:
                results.append(av)
        return results

    def to_json(self, variants: list[AnnotatedVariant]) -> str:
        """Serialize annotated variants to JSON string."""
        return json.dumps(
            [v.to_dict() for v in variants],
            indent=2,
            ensure_ascii=False,
        )

    def save_json(
        self, variants: list[AnnotatedVariant], output_path: str
    ):
        """Save annotated variants to a JSON file."""
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        with open(output_path, "w") as f:
            json.dump(
                [v.to_dict() for v in variants],
                f, indent=2, ensure_ascii=False,
            )
        logger.info(f"Saved {len(variants)} annotations to {output_path}")
