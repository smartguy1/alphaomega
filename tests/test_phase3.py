"""
AlphaOmega — Annotator & Trait Clusterer Tests

Tests annotation merging and domain classification using mock data.
No network access or real VCF files needed.
"""

import gzip
import json
import os
import shutil
import sqlite3
import subprocess
import pytest

from src.annotator import (
    Annotator,
    AnnotatedVariant,
    SNPediaAnnotation,
    ClinVarAnnotation,
    PharmCATAnnotation,
    _gt_to_alleles,
)
from src.trait_clusterer import (
    TraitClusterer,
    Domain,
    classify_variant,
)
from src.snpedia_mirror import SNPediaMirror
from src.clinvar_mirror import ClinVarMirror

# Check if bcftools is available for integration tests
_bcftools_available = shutil.which("bcftools") is not None


# ── Fixtures ─────────────────────────────────────────────────────────────────

@pytest.fixture
def mock_snpedia_dir(tmp_path):
    """Mock SNPedia database with test data."""
    db_path = tmp_path / "snpedia.sqlite"
    conn = sqlite3.connect(str(db_path))
    conn.execute("""
        CREATE TABLE snps (
            rsid TEXT PRIMARY KEY, chromosome TEXT,
            position INTEGER, gene TEXT, summary TEXT
        )
    """)
    conn.execute("""
        CREATE TABLE genotypes (
            rsid TEXT, genotype TEXT, magnitude REAL,
            repute TEXT, summary TEXT, details TEXT
        )
    """)
    conn.executemany("INSERT INTO snps VALUES (?,?,?,?,?)", [
        ("rs1815739", "11", 66560624, "ACTN3", "Muscle performance"),
        ("rs762551", "15", 75041917, "CYP1A2", "Caffeine metabolism"),
        ("rs4988235", "2", 136608646, "MCM6", "Lactose intolerance"),
    ])
    conn.executemany("INSERT INTO genotypes VALUES (?,?,?,?,?,?)", [
        ("rs1815739", "(C;C)", 3.0, "good", "Sprint performance", ""),
        ("rs1815739", "(T;T)", 2.5, "bad", "Endurance type", ""),
        ("rs762551", "(A;A)", 2.0, "good", "Fast caffeine metabolizer", ""),
        ("rs762551", "(C;C)", 2.5, "bad", "Slow caffeine metabolizer", ""),
        ("rs4988235", "(A;A)", 2.5, "good", "Lactose tolerant", ""),
    ])
    conn.commit()
    conn.close()
    return str(tmp_path)


@pytest.fixture
def mock_clinvar_dir(tmp_path):
    """Mock ClinVar database with test data."""
    db_path = tmp_path / "clinvar.sqlite"
    conn = sqlite3.connect(str(db_path))
    conn.execute("""
        CREATE TABLE variants (
            rsid TEXT, chromosome TEXT, position INTEGER,
            ref TEXT, alt TEXT, clinical_significance TEXT,
            review_status TEXT, review_stars INTEGER,
            conditions TEXT, gene TEXT,
            molecular_consequence TEXT, origin TEXT, clinvar_id TEXT
        )
    """)
    conn.execute("CREATE INDEX idx_rsid ON variants(rsid)")
    conn.executemany("INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)", [
        ("rs121918506", "17", 7673803, "G", "A", "Pathogenic",
         "reviewed_by_expert_panel", 3, "Li-Fraumeni syndrome",
         "TP53", "missense_variant", "germline", "12345"),
        ("rs4244285", "10", 96541616, "G", "A", "Pathogenic",
         "practice_guideline", 4, "CYP2C19 poor metabolizer",
         "CYP2C19", "splice_donor_variant", "germline", "11111"),
    ])
    conn.commit()
    conn.close()
    return str(tmp_path)


def _make_variant(
    rsid: str,
    snpedia_summary: str = "",
    snpedia_mag: float = 0.0,
    snpedia_repute: str = "",
    clinvar_sig: str = "",
    clinvar_gene: str = "",
    clinvar_conditions: list = None,
    clinvar_stars: int = 0,
    pharmcat_gene: str = "",
    pharmcat_phenotype: str = "",
) -> AnnotatedVariant:
    """Helper to build AnnotatedVariant for testing."""
    v = AnnotatedVariant(rsid=rsid)
    if snpedia_summary or snpedia_mag > 0:
        v.snpedia = SNPediaAnnotation(
            magnitude=snpedia_mag, repute=snpedia_repute, summary=snpedia_summary,
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


# ── Annotator tests ──────────────────────────────────────────────────────────

class TestAnnotatedVariant:

    def test_source_count_all(self):
        v = _make_variant(
            "rs1", snpedia_summary="test", snpedia_mag=2.0,
            clinvar_sig="Pathogenic", pharmcat_gene="CYP2D6",
        )
        assert v.source_count == 3

    def test_source_count_snpedia_only(self):
        v = _make_variant("rs1", snpedia_summary="test", snpedia_mag=1.0)
        assert v.source_count == 1
        assert v.has_snpedia
        assert not v.has_clinvar
        assert not v.has_pharmcat

    def test_source_count_none(self):
        v = AnnotatedVariant(rsid="rs1")
        assert v.source_count == 0

    def test_max_magnitude(self):
        v = _make_variant("rs1", snpedia_mag=3.5, clinvar_sig="Pathogenic", clinvar_stars=2)
        assert v.max_magnitude == 3.5

    def test_to_dict(self):
        v = _make_variant("rs1", snpedia_summary="test", snpedia_mag=2.0)
        d = v.to_dict()
        assert d["rsid"] == "rs1"
        assert "snpedia" in d["sources"]
        assert d["sources"]["snpedia"]["magnitude"] == 2.0


class TestGtToAlleles:

    def test_homozygous_ref(self):
        assert _gt_to_alleles("0/0", "A", "G") == "AA"

    def test_heterozygous(self):
        assert _gt_to_alleles("0/1", "A", "G") == "AG"

    def test_homozygous_alt(self):
        assert _gt_to_alleles("1/1", "A", "G") == "GG"

    def test_phased(self):
        assert _gt_to_alleles("0|1", "C", "T") == "CT"

    def test_multi_allelic(self):
        assert _gt_to_alleles("1/2", "A", "G,T") == "GT"


class TestAnnotatorWithMocks:

    def test_annotate_single_snpedia(self, mock_snpedia_dir):
        with SNPediaMirror(db_dir=mock_snpedia_dir) as snpedia:
            annotator = Annotator(snpedia=snpedia)
            v = annotator.annotate_single("rs1815739", genotype="CC")

            assert v is not None
            assert v.has_snpedia
            assert v.snpedia.magnitude == 3.0
            assert "sprint" in v.snpedia.summary.lower()

    def test_annotate_single_clinvar(self, mock_clinvar_dir):
        with ClinVarMirror(db_dir=mock_clinvar_dir) as clinvar:
            annotator = Annotator(clinvar=clinvar)
            v = annotator.annotate_single("rs121918506")

            assert v is not None
            assert v.has_clinvar
            assert "pathogenic" in v.clinvar.significance.lower()
            assert v.clinvar.gene == "TP53"

    def test_annotate_single_both_sources(self, mock_snpedia_dir, mock_clinvar_dir):
        with SNPediaMirror(db_dir=mock_snpedia_dir) as snpedia:
            with ClinVarMirror(db_dir=mock_clinvar_dir) as clinvar:
                annotator = Annotator(snpedia=snpedia, clinvar=clinvar)
                # This rsid exists in ClinVar but not SNPedia mock
                v = annotator.annotate_single("rs121918506")
                assert v is not None
                assert v.has_clinvar

    def test_annotate_nonexistent_returns_none(self, mock_snpedia_dir):
        with SNPediaMirror(db_dir=mock_snpedia_dir) as snpedia:
            annotator = Annotator(snpedia=snpedia)
            v = annotator.annotate_single("rs99999999")
            assert v is None

    def test_annotate_nonexistent_with_unannotated(self, mock_snpedia_dir):
        with SNPediaMirror(db_dir=mock_snpedia_dir) as snpedia:
            annotator = Annotator(snpedia=snpedia, include_unannotated=True)
            v = annotator.annotate_single("rs99999999")
            assert v is not None
            assert v.source_count == 0

    def test_annotate_rsids_batch(self, mock_snpedia_dir):
        with SNPediaMirror(db_dir=mock_snpedia_dir) as snpedia:
            annotator = Annotator(snpedia=snpedia)
            results = annotator.annotate_rsids(
                ["rs1815739", "rs762551", "rs99999999"]
            )
            # rs99999999 has no annotation and should be excluded
            assert len(results) == 2

    def test_json_serialization(self, mock_snpedia_dir):
        with SNPediaMirror(db_dir=mock_snpedia_dir) as snpedia:
            annotator = Annotator(snpedia=snpedia)
            results = annotator.annotate_rsids(["rs1815739"])
            json_str = annotator.to_json(results)

            parsed = json.loads(json_str)
            assert len(parsed) == 1
            assert parsed[0]["rsid"] == "rs1815739"


# ── Trait Clusterer tests ────────────────────────────────────────────────────

class TestClassifyVariant:

    def test_pharma_gene(self):
        v = _make_variant("rs1", clinvar_gene="CYP2D6", clinvar_sig="Pathogenic")
        assert classify_variant(v) == Domain.PHARMACOGENOMICS

    def test_pharmcat_annotation(self):
        v = _make_variant("rs1", pharmcat_gene="CYP2C19", pharmcat_phenotype="poor metabolizer")
        assert classify_variant(v) == Domain.PHARMACOGENOMICS

    def test_hla_gene(self):
        v = _make_variant("rs1", clinvar_gene="HLA-B", clinvar_sig="risk factor")
        assert classify_variant(v) == Domain.IMMUNITY

    def test_metabolism_keyword(self):
        v = _make_variant("rs1", snpedia_summary="Fast caffeine metabolizer")
        assert classify_variant(v) == Domain.METABOLISM

    def test_disease_risk_keyword(self):
        v = _make_variant("rs1", snpedia_summary="Increased cancer risk", snpedia_mag=3.0)
        assert classify_variant(v) == Domain.DISEASE_RISK

    def test_physical_trait_keyword(self):
        v = _make_variant("rs1", snpedia_summary="Determines eye color")
        assert classify_variant(v) == Domain.PHYSICAL_TRAITS

    def test_behavior_keyword(self):
        v = _make_variant("rs1", snpedia_summary="Oxytocin receptor and empathy")
        assert classify_variant(v) == Domain.BEHAVIOR

    def test_muscle_is_physical(self):
        v = _make_variant("rs1", snpedia_summary="Muscle performance - ACTN3 sprint type")
        assert classify_variant(v) == Domain.PHYSICAL_TRAITS

    def test_clinvar_condition_disease(self):
        v = _make_variant(
            "rs1", clinvar_sig="Pathogenic",
            clinvar_conditions=["Li-Fraumeni syndrome"],
        )
        assert classify_variant(v) == Domain.DISEASE_RISK

    def test_uncategorized_fallback(self):
        v = _make_variant("rs1", snpedia_summary="Unknown significance", snpedia_mag=1.0)
        assert classify_variant(v) == Domain.UNCATEGORIZED

    def test_empty_variant(self):
        v = AnnotatedVariant(rsid="rs1")
        assert classify_variant(v) == Domain.UNCATEGORIZED


class TestTraitClusterer:

    def test_cluster_basic(self):
        variants = [
            _make_variant("rs1", snpedia_summary="Fast caffeine metabolizer", snpedia_mag=2.0),
            _make_variant("rs2", snpedia_summary="Sprint muscle type ACTN3", snpedia_mag=3.0),
            _make_variant("rs3", clinvar_gene="CYP2D6", clinvar_sig="Pathogenic"),
            _make_variant("rs4", snpedia_summary="Increased cancer risk", snpedia_mag=4.0),
        ]

        clusterer = TraitClusterer()
        clusters = clusterer.cluster_nonempty(variants)

        assert Domain.METABOLISM in clusters
        assert Domain.PHYSICAL_TRAITS in clusters
        assert Domain.PHARMACOGENOMICS in clusters
        assert Domain.DISEASE_RISK in clusters

        assert clusters[Domain.METABOLISM].count == 1
        assert clusters[Domain.PHYSICAL_TRAITS].count == 1

    def test_summary(self):
        variants = [
            _make_variant("rs1", snpedia_summary="Fast caffeine metabolizer", snpedia_mag=2.0),
            _make_variant("rs2", snpedia_summary="Increased cancer risk", snpedia_mag=4.0),
        ]

        clusterer = TraitClusterer()
        clusters = clusterer.cluster(variants)
        summary = clusterer.get_summary(clusters)

        assert summary["total_variants"] == 2
        assert "metabolism" in summary["domains"]
        assert "disease_risk" in summary["domains"]

    def test_empty_input(self):
        clusterer = TraitClusterer()
        clusters = clusterer.cluster_nonempty([])
        assert len(clusters) == 0

    def test_json_output(self):
        variants = [
            _make_variant("rs1", snpedia_summary="Caffeine metabolizer", snpedia_mag=2.0),
        ]
        clusterer = TraitClusterer()
        clusters = clusterer.cluster(variants)
        json_str = clusterer.to_json(clusters)

        parsed = json.loads(json_str)
        assert "metabolism" in parsed


# ── Profile Manager tests ────────────────────────────────────────────────────

class TestProfileManager:

    def test_create_and_get(self, tmp_path):
        from src.profile_manager import ProfileManager
        with ProfileManager(
            db_dir=str(tmp_path / "db"),
            profiles_dir=str(tmp_path / "profiles"),
        ) as pm:
            pid = pm.create_profile("martin", bam_path="/data/martin.bam")
            assert pid > 0

            profile = pm.get_profile("martin")
            assert profile is not None
            assert profile.name == "martin"
            assert profile.bam_path == "/data/martin.bam"

    def test_list_profiles(self, tmp_path):
        from src.profile_manager import ProfileManager
        with ProfileManager(
            db_dir=str(tmp_path / "db"),
            profiles_dir=str(tmp_path / "profiles"),
        ) as pm:
            pm.create_profile("alice")
            pm.create_profile("bob")

            profiles = pm.list_profiles()
            names = [p.name for p in profiles]
            assert "alice" in names
            assert "bob" in names

    def test_store_and_retrieve_annotations(self, tmp_path):
        from src.profile_manager import ProfileManager
        with ProfileManager(
            db_dir=str(tmp_path / "db"),
            profiles_dir=str(tmp_path / "profiles"),
        ) as pm:
            pm.create_profile("martin")

            variants = [
                _make_variant("rs1", snpedia_summary="Caffeine", snpedia_mag=2.5),
                _make_variant("rs2", clinvar_sig="Pathogenic", clinvar_gene="TP53"),
            ]
            pm.store_annotations("martin", variants, domain_map={
                "rs1": "metabolism", "rs2": "disease_risk",
            })

            annots = pm.get_annotations("martin")
            assert len(annots) == 2

    def test_get_annotations_by_domain(self, tmp_path):
        from src.profile_manager import ProfileManager
        with ProfileManager(
            db_dir=str(tmp_path / "db"),
            profiles_dir=str(tmp_path / "profiles"),
        ) as pm:
            pm.create_profile("martin")
            variants = [
                _make_variant("rs1", snpedia_summary="Caffeine", snpedia_mag=2.5),
                _make_variant("rs2", clinvar_sig="Pathogenic", clinvar_gene="TP53"),
            ]
            pm.store_annotations("martin", variants, domain_map={
                "rs1": "metabolism", "rs2": "disease_risk",
            })

            metab = pm.get_annotations("martin", domain="metabolism")
            assert len(metab) == 1
            assert metab[0]["rsid"] == "rs1"

    def test_shared_variants(self, tmp_path):
        from src.profile_manager import ProfileManager
        with ProfileManager(
            db_dir=str(tmp_path / "db"),
            profiles_dir=str(tmp_path / "profiles"),
        ) as pm:
            pm.create_profile("alice")
            pm.create_profile("bob")

            shared_v = _make_variant("rs100", snpedia_summary="Shared trait", snpedia_mag=2.0)
            unique_v = _make_variant("rs200", snpedia_summary="Unique trait", snpedia_mag=1.0)

            pm.store_annotations("alice", [shared_v, unique_v])
            pm.store_annotations("bob", [shared_v])

            shared = pm.find_shared_variants(["alice", "bob"])
            assert len(shared) == 1
            assert shared[0]["rsid"] == "rs100"

    def test_unique_variants(self, tmp_path):
        from src.profile_manager import ProfileManager
        with ProfileManager(
            db_dir=str(tmp_path / "db"),
            profiles_dir=str(tmp_path / "profiles"),
        ) as pm:
            pm.create_profile("alice")
            pm.create_profile("bob")

            pm.store_annotations("alice", [
                _make_variant("rs100", snpedia_summary="Shared", snpedia_mag=2.0),
                _make_variant("rs200", snpedia_summary="Only Alice", snpedia_mag=3.0),
            ])
            pm.store_annotations("bob", [
                _make_variant("rs100", snpedia_summary="Shared", snpedia_mag=2.0),
            ])

            unique = pm.find_unique_variants("alice", ["bob"])
            assert len(unique) == 1
            assert unique[0]["rsid"] == "rs200"

    def test_variant_across_profiles(self, tmp_path):
        from src.profile_manager import ProfileManager
        with ProfileManager(
            db_dir=str(tmp_path / "db"),
            profiles_dir=str(tmp_path / "profiles"),
        ) as pm:
            pm.create_profile("alice")
            pm.create_profile("bob")

            pm.store_annotations("alice", [
                _make_variant("rs100", snpedia_summary="Test", snpedia_mag=2.0),
            ])
            pm.store_annotations("bob", [
                _make_variant("rs100", snpedia_summary="Test", snpedia_mag=2.0),
            ])

            results = pm.get_variant_across_profiles("rs100")
            assert len(results) == 2
            names = {r["profile_name"] for r in results}
            assert names == {"alice", "bob"}

    def test_delete_profile(self, tmp_path):
        from src.profile_manager import ProfileManager
        with ProfileManager(
            db_dir=str(tmp_path / "db"),
            profiles_dir=str(tmp_path / "profiles"),
        ) as pm:
            pm.create_profile("temp")
            pm.store_annotations("temp", [
                _make_variant("rs1", snpedia_summary="Test", snpedia_mag=1.0),
            ])

            assert pm.delete_profile("temp")
            assert pm.get_profile("temp") is None
            assert pm.get_annotation_count("temp") == 0

    def test_profile_dir_created(self, tmp_path):
        from src.profile_manager import ProfileManager
        with ProfileManager(
            db_dir=str(tmp_path / "db"),
            profiles_dir=str(tmp_path / "profiles"),
        ) as pm:
            pm.create_profile("martin")
            assert pm.get_profile_dir("martin").exists()


# ── rsID Annotation tests ────────────────────────────────────────────────────

class TestAnnotateRsids:

    @pytest.mark.skipif(not _bcftools_available, reason="bcftools not installed")
    def test_annotate_rsids_injects_ids(self, tmp_path):
        """End-to-end test: create a raw VCF and a ClinVar VCF, verify rsIDs are injected."""
        from src.variant_caller import _annotate_rsids

        db_dir = str(tmp_path / "db")
        os.makedirs(db_dir)

        # Create a minimal raw VCF (no rsIDs — just ".")
        raw_vcf_content = (
            "##fileformat=VCFv4.2\n"
            '##contig=<ID=1,length=249250621>\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            "1\t752566\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
            "1\t752721\t.\tC\tT\t60\tPASS\t.\tGT\t1/1\n"
        )
        raw_vcf_path = str(tmp_path / "test.vcf.gz")
        with gzip.open(raw_vcf_path, "wt") as f:
            f.write(raw_vcf_content)

        # Index the raw VCF
        subprocess.run(["bcftools", "index", "-f", raw_vcf_path],
                      capture_output=True, timeout=30)

        # Create a minimal ClinVar-style VCF with rsIDs at matching positions
        clinvar_content = (
            "##fileformat=VCFv4.2\n"
            '##contig=<ID=1,length=249250621>\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t752566\trs3094315\tA\tG\t.\t.\t.\n"
        )
        clinvar_path = os.path.join(db_dir, "clinvar.vcf.gz")
        with gzip.open(clinvar_path, "wt") as f:
            f.write(clinvar_content)

        # Index ClinVar VCF
        subprocess.run(["bcftools", "index", "-f", "-t", clinvar_path],
                      capture_output=True, timeout=30)

        # Run _annotate_rsids
        success, count = _annotate_rsids(raw_vcf_path, db_dir)

        assert success is True
        assert count >= 1  # At least one variant should get an rsID

        # Verify the VCF now contains the rsID
        result = subprocess.run(
            ["bcftools", "query", "-f", "%ID\\n", raw_vcf_path],
            capture_output=True, text=True, timeout=30,
        )
        ids = result.stdout.strip().splitlines()
        assert "rs3094315" in ids

    def test_annotate_rsids_missing_clinvar(self, tmp_path):
        """Verify graceful fallback when ClinVar VCF is not present."""
        from src.variant_caller import _annotate_rsids

        success, count = _annotate_rsids("fake.vcf.gz", str(tmp_path))
        assert success is False
        assert count == 0
