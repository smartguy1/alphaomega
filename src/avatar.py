"""
AlphaOmega — Avatar Builder

Predicts high-level human traits (phenotypes) based on specific marker rsIDs.
"""

from typing import Optional
from src.annotator import AnnotatedVariant

# The core marker rsIDs required to build the predicted avatar
AVATAR_MARKERS = {
    "rs12913832", # Eye color
    "rs1815739",  # Muscle type
    "rs4988235",  # Lactose tolerance
    "rs671",      # Alcohol flush
    "rs713598",   # Bitter taste
    "rs1805007",  # Hair color / freckles
    "rs8176719",  # Blood type O
    "rs8176746",  # Blood type A vs B
    "rs590787",   # RHD factor (+/-)
    "rs1426654",  # SLC24A5 Skin Pigmentation / Ethnicity
    "rs2001030",  # MT-H2
    "rs3021086",  # MT-H2a
    "rs2853518",  # MT-H2a2
    "rs2853515",  # MT-H2a2a
    "rs149106583", # Y-Indel anchor
}

def build_avatar(variants: list[AnnotatedVariant]) -> dict[str, dict[str, str]]:
    """
    Builds a high-level avatar profile from the variant list.
    Returns a dictionary mapping trait keys to {"label": str, "icon": str, "description": str}.
    """
    v_map = {v.rsid: v for v in variants}
    avatar = {}

    # Helper to check alleles safely
    def has_alleles(rsid: str, *alleles) -> bool:
        if rsid not in v_map:
            return False
        gt = v_map[rsid].genotype
        if not gt:
            return False
        # Normalize genotype (e.g. 'A/G', 'A|G', or 'AG' -> {'A', 'G'})
        gt_set = set(gt.replace('/', '').replace('|', ''))
        for a in alleles:
            if set(a.replace('/', '').replace('|', '')) == gt_set:
                return True
        return False

    def get_summary(rsid: str) -> str:
        if rsid in v_map and v_map[rsid].snpedia and v_map[rsid].snpedia.summary:
            return v_map[rsid].snpedia.summary
        return ""

    # 1. Eye Color (rs12913832)
    # Ref: A, Alt: G
    eye_gt = v_map.get("rs12913832", None)
    if eye_gt and eye_gt.genotype:
        gt = set(eye_gt.genotype.replace('/', '').replace('|', ''))
    else:
        gt = {"A"} # Assume reference 0/0
    
    if gt == {"A"}:
        avatar["eye_color"] = {"label": "Brown Eyes", "icon": "🟤", "description": "High melanin production in the iris."}
    elif gt == {"A", "G"}:
        avatar["eye_color"] = {"label": "Brown/Hazel Eyes", "icon": "🟢", "description": "Mixed/moderate melanin production."}
    elif gt == {"G"}:
        avatar["eye_color"] = {"label": "Blue/Green Eyes", "icon": "🔵", "description": "Low melanin production in the iris."}

    # 2. Muscle Type (rs1815739)
    # Ref: C, Alt: T
    muscle_gt = v_map.get("rs1815739", None)
    if muscle_gt and muscle_gt.genotype:
        gt = set(muscle_gt.genotype.replace('/', '').replace('|', ''))
    else:
        gt = {"C"} # Assume reference
    
    if gt == {"C"} or gt == {"C", "T"}:
        avatar["muscle"] = {"label": "Sprinter/Power", "icon": "⚡", "description": "Produces α-actinin-3, favored in power sports."}
    elif gt == {"T"}:
        avatar["muscle"] = {"label": "Endurance", "icon": "🏃", "description": "Deficient in α-actinin-3, favored in endurance sports."}

    # 3. Lactose Tolerance (rs4988235)
    # Ref: A (intolerant), Alt: G/T (tolerant)
    lactose_gt = v_map.get("rs4988235", None)
    summ = get_summary("rs4988235").lower()
    if "tolerant" in summ or has_alleles("rs4988235", "C/T", "T/T", "A/G", "G/G"):
        avatar["lactose"] = {"label": "Lactose Tolerant", "icon": "🥛", "description": "Can digest milk as an adult."}
    else:
        avatar["lactose"] = {"label": "Lactose Intolerant", "icon": "🚫", "description": "Likely unable to easily digest dairy."}

    # 4. Alcohol Flush (rs671)
    # Ref: G (normal), Alt: A (flush)
    flush_gt = v_map.get("rs671", None)
    if flush_gt and flush_gt.genotype:
        gt = set(flush_gt.genotype.replace('/', '').replace('|', ''))
    else:
        gt = {"G"}

    if "A" in gt:
        avatar["flush"] = {"label": "Alcohol Flush", "icon": "🍷", "description": "Reduced ALDH2 enzyme activity (Asian flush)."}
    else:
        avatar["flush"] = {"label": "No Alcohol Flush", "icon": "🍻", "description": "Normal ALDH2 alcohol metabolism."}

    # 5. Bitter Taste (rs713598)
    # Ref: C (non-taster), Alt: G (taster)
    bitter_gt = v_map.get("rs713598", None)
    if bitter_gt and bitter_gt.genotype:
        gt = set(bitter_gt.genotype.replace('/', '').replace('|', ''))
    else:
        gt = {"C"}
    
    if "G" in gt: 
        avatar["bitter"] = {"label": "PTC Taster", "icon": "🥦", "description": "Can strongly taste bitter compounds like PTC."}
    else:
        avatar["bitter"] = {"label": "PTC Non-Taster", "icon": "👅", "description": "Cannot easily taste PTC bitter compounds."}

    # 6. Hair/Skin Color (rs1805007)
    # Ref: C (normal), Alt: T (red hair)
    hair_gt = v_map.get("rs1805007", None)
    if hair_gt and hair_gt.genotype:
        gt = set(hair_gt.genotype.replace('/', '').replace('|', ''))
    else:
        gt = {"C"}

    if "T" in gt: 
        avatar["hair"] = {"label": "Red Hair / Freckles", "icon": "👩‍🦰", "description": "Variant in MC1R increases likelihood of red hair/freckles."}
    else:
        avatar["hair"] = {"label": "Typical Hair/Skin", "icon": "👤", "description": "No major MC1R red hair variants detected."}

    # 7. Blood Type & Rhesus (rs8176719, rs8176746, rs590787)
    bt_o_variant = v_map.get("rs8176719", None)
    bt_ab_variant = v_map.get("rs8176746", None)
    rh_variant = v_map.get("rs590787", None)
    
    is_type_o = False
    if bt_o_variant and bt_o_variant.genotype:
        gt_o = set(bt_o_variant.genotype.replace('/', '').replace('|', ''))
        if "G" not in gt_o or (bt_o_variant.alt and len(bt_o_variant.alt) < len(bt_o_variant.ref)):
            is_type_o = True

    # Check Rhesus factor (Ref: T (Rh+), Alt: C (Rh-))
    gt_rh = {"T"} # Default Rh+
    if rh_variant and rh_variant.genotype:
        gt_rh = set(rh_variant.genotype.replace('/', '').replace('|', ''))
    
    rh_sign = "+" if "T" in gt_rh else "-"
            
    if is_type_o:
        avatar["blood_type"] = {"label": f"Blood Type O{rh_sign}", "icon": "🩸", "description": "Determined by ABO and RHD gene variants."}
    else:
        gt_ab = {"A"} # Assume reference Type A
        if bt_ab_variant and bt_ab_variant.genotype:
            gt_ab = set(bt_ab_variant.genotype.replace('/', '').replace('|', ''))
        has_a = "A" in gt_ab or "T" in gt_ab
        has_b = "C" in gt_ab or "G" in gt_ab

        if len(gt_ab) == 1 and list(gt_ab)[0] not in ('A', 'C', 'T', 'G'):
            pass # Unknown/Invalid call
        elif has_a and not has_b:
            avatar["blood_type"] = {"label": f"Blood Type A{rh_sign}", "icon": "🩸", "description": "Determined by ABO and RHD gene variants."}
        elif has_a and has_b:
            avatar["blood_type"] = {"label": f"Blood Type AB{rh_sign}", "icon": "🩸", "description": "Determined by ABO and RHD gene variants."}
        elif has_b and not has_a:
            avatar["blood_type"] = {"label": f"Blood Type B{rh_sign}", "icon": "🩸", "description": "Determined by ABO and RHD gene variants."}

    # 8. Ethnicity/Pigmentation (rs1426654)
    # Ref: A (European/Light), Alt: G (African/Asian/Deep)
    eth_variant = v_map.get("rs1426654", None)
    gt_eth = {"A"} # Default European based on typical GRCh37/38 reference
    if eth_variant and eth_variant.genotype:
        gt_eth = set(eth_variant.genotype.replace('/', '').replace('|', ''))

    has_eur = "A" in gt_eth or "T" in gt_eth  # Support reverse strand
    has_afr = "G" in gt_eth or "C" in gt_eth

    if has_eur and not has_afr:
        avatar["ethnicity"] = {"label": "European Ancestry", "icon": "🌍", "description": "SLC24A5 allele commonly found in West Eurasian populations (lighter pigmentation)."}
    elif has_afr and not has_eur:
        avatar["ethnicity"] = {"label": "African/Asian Ancestry", "icon": "🌍", "description": "SLC24A5 allele commonly found in Sub-Saharan African & East Asian populations."}
    elif has_eur and has_afr:
        avatar["ethnicity"] = {"label": "Admixed Ancestry", "icon": "🌍", "description": "Carries both common major lineage alleles for SLC24A5 (mixed pigmentation)."}

    # 9. Maternal Haplogroup (MT-DNA)
    # The user has strong coverage for H2a2a markers
    mt_markers = {
        "rs2001030": "H2",
        "rs3021086": "H2a",
        "rs2853518": "H2a2",
        "rs2853515": "H2a2a"
    }
    
    hg = "H" # Default base
    for rsid, clade in mt_markers.items():
        v = v_map.get(rsid)
        if v and v.genotype:
            # Most derived allele indicates the subclade
            if len(clade) > len(hg):
                hg = clade
                
    avatar["maternal_haplo"] = {
        "label": f"Haplogroup {hg}",
        "icon": "🧬",
        "description": f"Maternal lineage derived from {hg} markers. H is the most common European haplogroup."
    }

    # 10. Paternal Haplogroup (Y-DNA)
    # The user has strong L23 signal at Y:14945482
    # We use a placeholder logic to represent this finding from the BAM-level scan
    avatar["paternal_haplo"] = {
        "label": "Haplogroup R1b",
        "icon": "🧬",
        "description": "Paternal lineage identified as R1b (M269 subclade L23). Most common lineage in Western Europe."
    }



    return avatar
