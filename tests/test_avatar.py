import pytest
from src.annotator import AnnotatedVariant, SNPediaAnnotation
from src.avatar import build_avatar

def test_avatar_eye_color():
    v1 = AnnotatedVariant(rsid="rs12913832", genotype="A/A")
    avatar = build_avatar([v1])
    assert "eye_color" in avatar
    assert avatar["eye_color"]["label"] == "Brown Eyes"

    v2 = AnnotatedVariant(rsid="rs12913832", genotype="G/G")
    avatar = build_avatar([v2])
    assert avatar["eye_color"]["label"] == "Blue/Green Eyes"

def test_avatar_muscle_type():
    v = AnnotatedVariant(rsid="rs1815739", genotype="C/T")
    avatar = build_avatar([v])
    assert "muscle" in avatar
    assert avatar["muscle"]["label"] == "Sprinter/Power"

    v2 = AnnotatedVariant(rsid="rs1815739", genotype="T/T")
    avatar = build_avatar([v2])
    assert avatar["muscle"]["label"] == "Endurance"

def test_avatar_lactose_tolerance_by_genotype():
    v = AnnotatedVariant(rsid="rs4988235", genotype="C/T")
    avatar = build_avatar([v])
    assert "lactose" in avatar
    assert avatar["lactose"]["label"] == "Lactose Tolerant"

def test_avatar_lactose_tolerance_by_summary():
    v = AnnotatedVariant(rsid="rs4988235", genotype="G/G")
    v.snpedia = SNPediaAnnotation(magnitude=2.0, repute="good", summary="Lactose tolerant")
    avatar = build_avatar([v])
    assert "lactose" in avatar
    assert avatar["lactose"]["label"] == "Lactose Tolerant"

def test_avatar_flush():
    v = AnnotatedVariant(rsid="rs671", genotype="A/G")
    avatar = build_avatar([v])
    assert "flush" in avatar
    assert avatar["flush"]["label"] == "Alcohol Flush"

def test_empty_variants():
    avatar = build_avatar([])
    assert avatar == {}

def test_irrelevant_variants():
    v = AnnotatedVariant(rsid="rs999999", genotype="A/A")
    avatar = build_avatar([v])
    assert avatar == {}
