"""
AlphaOmega — Phase 5 Tests

Tests for relatedness analysis, FragmentOfUs export, and CLI commands.
"""

import json
import os
import sqlite3
import pytest
from click.testing import CliRunner

from src.annotator import AnnotatedVariant, SNPediaAnnotation, ClinVarAnnotation, PharmCATAnnotation
from src.relatedness import (
    _count_shared_alleles,
    _check_mendelian,
    _classify_relationship,
    analyze_relatedness_from_variants,
    Relationship,
    MIN_SHARED_SNPS,
)
from src.export import export_for_biographer, _build_finding
from src.cli import cli


# ── Helpers ──────────────────────────────────────────────────────────────────

def _make_variant(rsid, genotype="", snpedia_summary="", snpedia_mag=0.0,
                  snpedia_repute="", clinvar_sig="", clinvar_gene="",
                  clinvar_conditions=None, clinvar_stars=0,
                  pharmcat_gene="", pharmcat_phenotype=""):
    v = AnnotatedVariant(rsid=rsid, genotype=genotype)
    if snpedia_summary or snpedia_mag > 0:
        v.snpedia = SNPediaAnnotation(
            magnitude=snpedia_mag, repute=snpedia_repute,
            summary=snpedia_summary,
        )
    if clinvar_sig:
        v.clinvar = ClinVarAnnotation(
            significance=clinvar_sig, gene=clinvar_gene,
            conditions=clinvar_conditions or [], review_stars=clinvar_stars,
        )
    if pharmcat_gene:
        v.pharmcat = PharmCATAnnotation(
            gene=pharmcat_gene, phenotype=pharmcat_phenotype,
        )
    return v


# ── Relatedness: Allele Sharing ──────────────────────────────────────────────

class TestCountSharedAlleles:

    def test_identical_homozygous(self):
        assert _count_shared_alleles("AA", "AA") == 2

    def test_identical_heterozygous(self):
        assert _count_shared_alleles("AG", "AG") == 2

    def test_reversed_heterozygous(self):
        assert _count_shared_alleles("AG", "GA") == 2

    def test_one_shared(self):
        assert _count_shared_alleles("AA", "AG") == 1

    def test_none_shared(self):
        assert _count_shared_alleles("AA", "GG") == 0

    def test_missing_data(self):
        assert _count_shared_alleles("", "AA") == -1
        assert _count_shared_alleles("A", "AA") == -1

    def test_different_bases(self):
        assert _count_shared_alleles("CT", "AG") == 0

    def test_partial_overlap(self):
        assert _count_shared_alleles("AC", "AG") == 1


class TestCheckMendelian:

    def test_consistent_het_to_homo(self):
        # Parent AG, Child AA → child got A from parent ✓
        assert _check_mendelian("AG", "AA") is True

    def test_consistent_homo_to_het(self):
        # Parent AA, Child AG → child got A from parent ✓
        assert _check_mendelian("AA", "AG") is True

    def test_violation(self):
        # Parent AA, Child GG → child has no A from parent ✗
        assert _check_mendelian("AA", "GG") is False

    def test_identical(self):
        assert _check_mendelian("AA", "AA") is True

    def test_missing(self):
        assert _check_mendelian("", "AA") is True  # Can't check, assume OK


class TestClassifyRelationship:

    def test_self_or_twin(self):
        rel, conf, _ = _classify_relationship(0.01, 0.05, 0.94, 0.0, 10000)
        assert rel == Relationship.SELF_OR_TWIN

    def test_parent_child(self):
        rel, conf, _ = _classify_relationship(0.01, 0.95, 0.04, 0.005, 10000)
        assert rel == Relationship.PARENT_CHILD
        assert conf > 0.8

    def test_full_siblings(self):
        rel, conf, _ = _classify_relationship(0.25, 0.50, 0.25, 0.01, 10000)
        assert rel == Relationship.FULL_SIBLINGS

    def test_unrelated(self):
        rel, conf, _ = _classify_relationship(0.70, 0.25, 0.05, 0.10, 10000)
        assert rel == Relationship.UNRELATED


class TestAnalyzeRelatednessFromVariants:

    def _make_parent_child_data(self, n=200):
        """Generate synthetic parent-child genotype data."""
        import random
        random.seed(42)
        bases = ["A", "C", "G", "T"]

        parent_variants = []
        child_variants = []

        for i in range(n):
            rsid = f"rs{i}"
            # Parent heterozygous
            b1, b2 = random.sample(bases, 2)
            parent_gt = b1 + b2

            # Child inherits one allele from parent + random from other parent
            inherited = random.choice([b1, b2])
            other = random.choice(bases)
            child_gt = inherited + other

            parent_variants.append({"rsid": rsid, "genotype": parent_gt})
            child_variants.append({"rsid": rsid, "genotype": child_gt})

        return parent_variants, child_variants

    def _make_unrelated_data(self, n=200):
        """Generate synthetic unrelated genotype data."""
        import random
        random.seed(99)
        bases = ["A", "C", "G", "T"]

        va, vb = [], []
        for i in range(n):
            rsid = f"rs{i}"
            va.append({"rsid": rsid, "genotype": random.choice(bases) + random.choice(bases)})
            vb.append({"rsid": rsid, "genotype": random.choice(bases) + random.choice(bases)})

        return va, vb

    def _make_identical_data(self, n=200):
        """Generate identical genotype data (self/twin)."""
        import random
        random.seed(7)
        bases = ["A", "C", "G", "T"]

        variants = []
        for i in range(n):
            rsid = f"rs{i}"
            gt = random.choice(bases) + random.choice(bases)
            variants.append({"rsid": rsid, "genotype": gt})

        return variants, [dict(v) for v in variants]

    def test_parent_child_detected(self):
        parent, child = self._make_parent_child_data(500)
        result = analyze_relatedness_from_variants("parent", "child", parent, child)

        assert result.relationship == Relationship.PARENT_CHILD
        assert result.mendelian_error_rate < 0.05

    def test_unrelated_detected(self):
        va, vb = self._make_unrelated_data(500)
        result = analyze_relatedness_from_variants("alice", "stranger", va, vb)

        assert result.relationship in (
            Relationship.UNRELATED, Relationship.SECOND_DEGREE,
            Relationship.HALF_SIBLINGS,
        )

    def test_self_detected(self):
        va, vb = self._make_identical_data(500)
        result = analyze_relatedness_from_variants("alice", "alice_copy", va, vb)

        assert result.relationship == Relationship.SELF_OR_TWIN

    def test_insufficient_snps(self):
        va = [{"rsid": "rs1", "genotype": "AA"}]
        vb = [{"rsid": "rs1", "genotype": "AA"}]
        result = analyze_relatedness_from_variants("a", "b", va, vb)

        assert result.relationship == Relationship.INCONCLUSIVE

    def test_result_to_dict(self):
        parent, child = self._make_parent_child_data(300)
        result = analyze_relatedness_from_variants("p", "c", parent, child)

        d = result.to_dict()
        assert "profile_a" in d
        assert "relationship" in d
        assert "ibs0" in d

    def test_result_summary(self):
        parent, child = self._make_parent_child_data(300)
        result = analyze_relatedness_from_variants("alice", "bob", parent, child)

        s = result.summary
        assert "alice" in s
        assert "bob" in s


# ── Export ───────────────────────────────────────────────────────────────────

class TestExport:

    def test_export_structure(self):
        variants = [
            _make_variant("rs1", "AA", snpedia_summary="Fast caffeine", snpedia_mag=2.5),
            _make_variant("rs2", "GG", clinvar_sig="Pathogenic", clinvar_gene="TP53"),
        ]
        data = export_for_biographer("martin", variants)

        assert data["schema_version"] == "1.0"
        assert data["source"] == "AlphaOmega"
        assert data["profile"]["name"] == "martin"
        assert "disclaimer" in data
        assert len(data["genomic_highlights"]) > 0

    def test_export_has_narrative_hooks(self):
        variants = [
            _make_variant("rs1", "AA", snpedia_summary="Fast caffeine metabolizer", snpedia_mag=2.5),
        ]
        data = export_for_biographer("test", variants)

        # Metabolism domain should have hooks
        for key, domain in data["genomic_highlights"].items():
            assert "narrative_hooks" in domain
            assert len(domain["narrative_hooks"]) > 0

    def test_export_confidence_tagging(self):
        high_conf = _make_variant(
            "rs1", "AA", snpedia_summary="Well-known", snpedia_mag=4.0,
            clinvar_sig="Pathogenic", clinvar_gene="TP53",
        )
        finding = _build_finding(high_conf)
        assert finding["confidence"] == "high"

    def test_export_low_confidence(self):
        low_conf = _make_variant("rs1", snpedia_summary="Unknown", snpedia_mag=0.5)
        finding = _build_finding(low_conf)
        assert finding["confidence"] == "low"

    def test_save_export(self, tmp_path):
        from src.export import save_export

        variants = [_make_variant("rs1", "AA", snpedia_summary="Test", snpedia_mag=1.0)]
        data = export_for_biographer("test", variants)

        output = str(tmp_path / "export.json")
        path = save_export(data, output)

        assert os.path.exists(path)
        with open(path) as f:
            loaded = json.load(f)
            assert loaded["profile"]["name"] == "test"

    def test_export_empty_variants(self):
        data = export_for_biographer("empty", [])
        assert data["profile"]["name"] == "empty"
        assert len(data["genomic_highlights"]) == 0


# ── CLI ──────────────────────────────────────────────────────────────────────

class TestCLI:

    def test_help(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "AlphaOmega" in result.output

    def test_status_command(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["status"])
        assert result.exit_code == 0
        assert "AlphaOmega Status" in result.output

    def test_setup_help(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["setup", "--help"])
        assert result.exit_code == 0
        assert "reference" in result.output.lower()

    def test_analyze_help(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["analyze", "--help"])
        assert result.exit_code == 0
        assert "--bam" in result.output

    def test_compare_help(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["compare", "--help"])
        assert result.exit_code == 0
        assert "--profiles" in result.output

    def test_export_help(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["export", "--help"])
        assert result.exit_code == 0
        assert "FragmentOfUs" in result.output

    def test_report_missing_profile(self, tmp_path):
        runner = CliRunner()
        result = runner.invoke(cli, [
            "-c", str(tmp_path / "empty.yaml"),
            "report", "-p", "nonexistent",
        ])
        assert result.exit_code != 0
