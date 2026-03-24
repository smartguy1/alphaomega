"""
AlphaOmega — Report Generator Tests

Tests HTML report rendering using mock annotated data.
Verifies template rendering, stats calculation, and output structure.
"""

import json
import os
import pytest

from src.annotator import (
    AnnotatedVariant,
    SNPediaAnnotation,
    ClinVarAnnotation,
    PharmCATAnnotation,
)
from src.trait_clusterer import Domain, DomainCluster, TraitClusterer
from src.report_generator import ReportGenerator, DOMAIN_ICONS


# ── Fixtures ─────────────────────────────────────────────────────────────────

def _make_variant(
    rsid, genotype="", snpedia_mag=0.0, snpedia_repute="",
    snpedia_summary="", clinvar_sig="", clinvar_gene="",
    clinvar_conditions=None, clinvar_stars=0,
    pharmcat_gene="", pharmcat_phenotype="",
):
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


@pytest.fixture
def sample_variants():
    """Set of realistic annotated variants for testing."""
    return [
        _make_variant(
            "rs1815739", "CC",
            snpedia_mag=3.0, snpedia_repute="good",
            snpedia_summary="Sprint/power muscle performance",
        ),
        _make_variant(
            "rs762551", "AA",
            snpedia_mag=2.0, snpedia_repute="good",
            snpedia_summary="Fast caffeine metabolizer",
        ),
        _make_variant(
            "rs4988235", "GG",
            snpedia_mag=0.0, snpedia_repute="neutral",
            snpedia_summary="Likely lactose intolerant",
        ),
        _make_variant(
            "rs121918506", "",
            clinvar_sig="Pathogenic", clinvar_gene="TP53",
            clinvar_conditions=["Li-Fraumeni syndrome"],
            clinvar_stars=3,
        ),
        _make_variant(
            "rs4244285", "",
            clinvar_sig="Pathogenic", clinvar_gene="CYP2C19",
            clinvar_conditions=["CYP2C19 poor metabolizer"],
            clinvar_stars=4,
            pharmcat_gene="CYP2C19", pharmcat_phenotype="Poor metabolizer",
        ),
    ]


@pytest.fixture
def sample_clusters(sample_variants):
    """Pre-clustered variants."""
    clusterer = TraitClusterer()
    return clusterer.cluster(sample_variants)


@pytest.fixture
def template_dir():
    """Return the real template directory path."""
    # Look for templates/ relative to project root
    candidates = [
        os.path.join(os.path.dirname(os.path.dirname(__file__)), "templates"),
        "templates",
    ]
    for d in candidates:
        if os.path.exists(d) and os.path.isfile(os.path.join(d, "report.html")):
            return d
    pytest.skip("templates/ directory not found")


# ── Tests ────────────────────────────────────────────────────────────────────

class TestReportGenerator:

    def test_render_html_basic(self, template_dir, sample_clusters):
        gen = ReportGenerator(template_dir=template_dir)
        html = gen.render_html(
            profile_name="martin",
            clusters=sample_clusters,
            build="GRCh38",
        )

        assert isinstance(html, str)
        assert len(html) > 500
        assert "martin" in html
        assert "GRCh38" in html
        assert "AlphaOmega" in html

    def test_contains_disclaimer(self, template_dir, sample_clusters):
        gen = ReportGenerator(template_dir=template_dir)
        html = gen.render_html("martin", sample_clusters)

        assert "not a medical diagnostic tool" in html.lower()
        assert "disclaimer" in html.lower()

    def test_contains_domain_sections(self, template_dir, sample_clusters):
        gen = ReportGenerator(template_dir=template_dir)
        html = gen.render_html("test", sample_clusters)

        # Should contain at least some domain headers
        assert "variant" in html.lower()

    def test_contains_rsid_links(self, template_dir, sample_clusters):
        gen = ReportGenerator(template_dir=template_dir)
        html = gen.render_html("test", sample_clusters)

        # rsIDs should be rendered as SNPedia links
        assert "snpedia.com" in html
        assert "rs1815739" in html

    def test_pathogenic_highlighted(self, template_dir, sample_clusters):
        gen = ReportGenerator(template_dir=template_dir)
        html = gen.render_html("test", sample_clusters)

        assert "Pathogenic" in html

    def test_save_html(self, template_dir, sample_clusters, tmp_path):
        gen = ReportGenerator(template_dir=template_dir)
        html = gen.render_html("test", sample_clusters)

        output = str(tmp_path / "report.html")
        path = gen.save_html(html, output)

        assert os.path.exists(path)

        with open(path, "r", encoding="utf-8") as f:
            content = f.read()
            assert "test" in content
            assert len(content) > 500

    def test_generate_report_full_pipeline(
        self, template_dir, sample_variants, tmp_path
    ):
        gen = ReportGenerator(template_dir=template_dir)
        output = str(tmp_path / "full_report.html")

        path = gen.generate_report(
            profile_name="martin",
            variants=sample_variants,
            output_path=output,
            build="GRCh38",
        )

        assert os.path.exists(path)

        with open(path, "r", encoding="utf-8") as f:
            content = f.read()
            assert "martin" in content
            assert "rs1815739" in content

    def test_drug_recommendations_table(self, template_dir, sample_clusters):
        gen = ReportGenerator(template_dir=template_dir)

        drugs = [
            {
                "drug": "Clopidogrel",
                "gene": "CYP2C19",
                "phenotype": "Poor metabolizer",
                "recommendation": "Consider alternative antiplatelet therapy",
            },
        ]

        html = gen.render_html(
            "test", sample_clusters, drug_recommendations=drugs,
        )

        assert "Clopidogrel" in html
        assert "CYP2C19" in html
        assert "alternative" in html

    def test_empty_clusters(self, template_dir):
        gen = ReportGenerator(template_dir=template_dir)

        empty = {
            Domain.METABOLISM: DomainCluster(
                domain=Domain.METABOLISM,
                display_name="Metabolism",
                description="",
            ),
        }

        html = gen.render_html("empty_test", empty)
        assert "empty_test" in html
        assert isinstance(html, str)

    def test_missing_template_dir(self):
        gen = ReportGenerator(template_dir="/nonexistent/path")
        assert gen._env is None

        with pytest.raises(RuntimeError, match="Template directory"):
            gen.render_html("test", {})


class TestReportSummaryJSON:

    def test_summary_structure(self, template_dir, sample_clusters):
        gen = ReportGenerator(template_dir=template_dir)
        summary = gen.generate_summary_json("martin", sample_clusters)

        assert summary["profile_name"] == "martin"
        assert "generated_at" in summary
        assert summary["total_variants"] == 5
        assert len(summary["domains"]) > 0

    def test_summary_top_findings(self, template_dir, sample_clusters):
        gen = ReportGenerator(template_dir=template_dir)
        summary = gen.generate_summary_json("test", sample_clusters)

        # Each domain should have top findings
        for domain_key, domain_data in summary["domains"].items():
            assert "top_findings" in domain_data
            assert "display_name" in domain_data
            assert "count" in domain_data

    def test_save_json(self, template_dir, sample_clusters, tmp_path):
        gen = ReportGenerator(template_dir=template_dir)
        summary = gen.generate_summary_json("test", sample_clusters)

        output = str(tmp_path / "summary.json")
        path = gen.save_summary_json(summary, output)

        assert os.path.exists(path)

        with open(path) as f:
            loaded = json.load(f)
            assert loaded["profile_name"] == "test"


class TestDomainIcons:

    def test_all_domains_have_icons(self):
        for domain in Domain:
            assert domain.value in DOMAIN_ICONS, f"Missing icon for {domain.value}"
