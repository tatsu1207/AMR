#!/usr/bin/env python3
"""
Substrate-specific ARG-to-antibiotic mapping.
Maps individual resistance genes to the specific antibiotics they confer resistance to,
rather than using broad drug-class matching.
"""

# Beta-lactamase substrate specificity
# Based on: Bush & Jacoby (2010), Ambler classification, CARD substrate annotations

BETALACTAM_SUBSTRATE = {
    # Narrow-spectrum penicillinases (Ambler class A)
    # Confer: ampicillin, amoxicillin, piperacillin
    # Do NOT confer: cephalosporins, carbapenems
    "blaTEM-1": ["ampicillin", "amoxicillin", "amoxicillin_clavulanic_acid", "piperacillin"],
    "blaTEM-2": ["ampicillin", "amoxicillin", "piperacillin"],
    "blaSHV-1": ["ampicillin", "amoxicillin", "piperacillin"],
    "blaOXA-1": ["ampicillin", "amoxicillin", "piperacillin", "ampicillin_sulbactam"],
    "blaOXA-2": ["ampicillin", "amoxicillin"],
    "blaHER-3": ["ampicillin", "amoxicillin"],
    "blaCARB-2": ["ampicillin", "amoxicillin", "piperacillin"],

    # Intrinsic chromosomal AmpC (E. coli blaEC)
    # Low-level expression usually does not confer resistance to anything
    # Overexpression: cefoxitin, amox-clav, some cephalosporins
    # We treat blaEC as NOT conferring resistance unless overexpressed
    "blaEC": [],
    "blaEC-5": [],
    "blaEC-8": [],
    "blaEC-13": [],
    "blaEC-15": [],
    "blaEC-18": [],
    "blaEC-19": [],

    # ESBLs - CTX-M group 1 (CTX-M-15, -1, -3, -55)
    # Confer: cefotaxime, ceftriaxone, ceftiofur, cefepime, aztreonam
    # Do NOT reliably confer: ceftazidime (weak for some), cefoxitin, carbapenems
    "blaCTX-M-1": ["ampicillin", "cefotaxime", "ceftriaxone", "ceftiofur", "cefepime", "aztreonam"],
    "blaCTX-M-3": ["ampicillin", "cefotaxime", "ceftriaxone", "ceftiofur", "cefepime", "aztreonam"],
    "blaCTX-M-15": ["ampicillin", "cefotaxime", "ceftriaxone", "ceftiofur", "ceftazidime", "cefepime", "aztreonam"],
    "blaCTX-M-55": ["ampicillin", "cefotaxime", "ceftriaxone", "ceftiofur", "ceftazidime", "cefepime", "aztreonam"],

    # ESBLs - CTX-M group 9 (CTX-M-14, -27, -65)
    # Confer: cefotaxime, ceftriaxone
    # Do NOT confer: ceftazidime (key distinction), cefepime (variable)
    "blaCTX-M-14": ["ampicillin", "cefotaxime", "ceftriaxone", "ceftiofur"],
    "blaCTX-M-27": ["ampicillin", "cefotaxime", "ceftriaxone", "ceftiofur"],
    "blaCTX-M-65": ["ampicillin", "cefotaxime", "ceftriaxone", "ceftiofur", "ceftazidime"],

    # AmpC beta-lactamases (plasmid-borne, Ambler class C)
    # Confer: cefoxitin, ceftriaxone, ceftiofur, amox-clav
    # Do NOT confer: cefepime, carbapenems
    "blaCMY-2": ["ampicillin", "amoxicillin_clavulanic_acid", "cefoxitin", "ceftriaxone", "ceftiofur", "cefotaxime"],
    "blaCMY": ["ampicillin", "amoxicillin_clavulanic_acid", "cefoxitin", "ceftriaxone", "ceftiofur", "cefotaxime"],
    "blaDHA-1": ["ampicillin", "amoxicillin_clavulanic_acid", "cefoxitin", "ceftriaxone", "ceftiofur"],
    "blaACT": ["ampicillin", "amoxicillin_clavulanic_acid", "cefoxitin", "ceftriaxone", "ceftiofur"],

    # Carbapenemases (Ambler class A: KPC; class B: NDM, VIM, IMP; class D: OXA-48-like)
    # Confer: all beta-lactams including carbapenems
    "blaKPC-2": ["ampicillin", "amoxicillin_clavulanic_acid", "cefoxitin", "ceftriaxone", "ceftazidime",
                  "cefepime", "cefotaxime", "ceftiofur", "aztreonam", "ertapenem", "imipenem", "meropenem",
                  "doripenem", "piperacillin_tazobactam"],
    "blaKPC-3": ["ampicillin", "amoxicillin_clavulanic_acid", "cefoxitin", "ceftriaxone", "ceftazidime",
                  "cefepime", "cefotaxime", "ceftiofur", "aztreonam", "ertapenem", "imipenem", "meropenem",
                  "doripenem", "piperacillin_tazobactam"],
    "blaNDM-1": ["ampicillin", "amoxicillin_clavulanic_acid", "cefoxitin", "ceftriaxone", "ceftazidime",
                  "cefepime", "cefotaxime", "ertapenem", "imipenem", "meropenem", "doripenem",
                  "piperacillin_tazobactam"],
    "blaNDM-5": ["ampicillin", "amoxicillin_clavulanic_acid", "cefoxitin", "ceftriaxone", "ceftazidime",
                  "cefepime", "cefotaxime", "ertapenem", "imipenem", "meropenem", "doripenem",
                  "piperacillin_tazobactam"],
    "blaVIM": ["ampicillin", "cefoxitin", "ceftriaxone", "ceftazidime", "cefepime",
               "ertapenem", "imipenem", "meropenem", "doripenem"],
    "blaIMP": ["ampicillin", "cefoxitin", "ceftriaxone", "ceftazidime", "cefepime",
               "ertapenem", "imipenem", "meropenem", "doripenem"],
    "blaOXA-48": ["ampicillin", "amoxicillin_clavulanic_acid", "ertapenem", "imipenem", "meropenem", "doripenem"],

    # Inhibitor-resistant TEM variants
    "blaTEM-30": ["ampicillin", "amoxicillin_clavulanic_acid", "ampicillin_sulbactam"],
    "blaTEM-39": ["ampicillin", "amoxicillin_clavulanic_acid", "ampicillin_sulbactam"],

    # SHV ESBLs
    "blaSHV-12": ["ampicillin", "cefotaxime", "ceftriaxone", "ceftazidime", "cefepime", "aztreonam"],
    "blaSHV-28": ["ampicillin", "cefotaxime", "ceftriaxone", "cefepime"],

    # S. aureus
    "blaZ": ["penicillin", "ampicillin", "amoxicillin"],
    "mecA": ["oxacillin", "cefoxitin", "cefazolin"],
    "mecC": ["oxacillin"],

    # A. baumannii OXA carbapenemases
    "blaOXA-23": ["ampicillin", "imipenem", "meropenem", "doripenem", "ampicillin_sulbactam"],
    "blaOXA-24": ["ampicillin", "imipenem", "meropenem", "doripenem"],
    "blaOXA-58": ["ampicillin", "imipenem", "meropenem", "doripenem"],
    "blaOXA-66": ["ampicillin"],  # intrinsic, low-level
    "blaOXA-51": ["ampicillin"],  # intrinsic A. baumannii
    "blaOXA-82": ["ampicillin"],
    "blaADC": ["ampicillin", "ceftazidime"],  # A. baumannii AmpC
}

# Aminoglycoside modifying enzyme substrate specificity
AMINOGLYCOSIDE_SUBSTRATE = {
    # AAC(3) family - primarily gentamicin
    "aac(3)-IIa": ["gentamicin", "tobramycin"],
    "aac(3)-IId": ["gentamicin"],
    "aac(3)-IIe": ["gentamicin"],
    "aac(3)-Ia": ["gentamicin"],
    "aac(3)-IVa": ["gentamicin", "tobramycin", "amikacin"],
    "aac(3)-VIa": ["gentamicin"],
    "aac(3)-Id": ["gentamicin"],
    "aac(3)-C322": ["gentamicin"],

    # AAC(6') family - primarily amikacin/tobramycin
    "aac(6')-Ib": ["amikacin", "tobramycin", "kanamycin"],
    "aac(6')-Ib-cr5": ["amikacin", "tobramycin", "kanamycin"],  # also ciprofloxacin
    "aac(6')-Ib4": ["amikacin", "tobramycin"],
    "aac(6')-IIa": ["gentamicin", "tobramycin"],
    "aac(6')-Ian": ["amikacin", "tobramycin"],

    # APH family - primarily kanamycin/streptomycin
    "aph(3')-Ia": ["kanamycin"],
    "aph(3')-IIa": ["kanamycin"],
    "aph(3'')-Ib": ["streptomycin"],
    "aph(6)-Ic": ["streptomycin"],
    "aph(6)-Id": ["streptomycin"],
    "aph(4)-Ia": ["tobramycin"],

    # ANT family
    "ant(2'')-Ia": ["gentamicin", "tobramycin", "kanamycin"],
    "ant(3'')-Ia": ["streptomycin", "spectinomycin"],

    # AAD family
    "aadA1": ["streptomycin", "spectinomycin"],
    "aadA2": ["streptomycin", "spectinomycin"],
    "aadA5": ["streptomycin", "spectinomycin"],
    "aadA7": ["streptomycin", "spectinomycin"],

    # 16S rRNA methyltransferases - all aminoglycosides
    "armA": ["gentamicin", "tobramycin", "amikacin", "kanamycin", "streptomycin"],
    "rmtB": ["gentamicin", "tobramycin", "amikacin", "kanamycin"],
    "rmtC": ["gentamicin", "tobramycin", "amikacin", "kanamycin"],

    # Bifunctional (S. aureus)
    "aac(6')-Ie-aph(2'')-Ia": ["gentamicin", "tobramycin", "amikacin", "kanamycin"],
}


def gene_confers_resistance(gene_name, antibiotic, drug_class):
    """
    Check if a specific gene confers resistance to a specific antibiotic.
    Returns True if the gene is substrate-specific to this antibiotic.
    Falls back to drug-class matching for unmapped genes.
    """
    # Check beta-lactam substrate mapping
    if drug_class == "BETA-LACTAM":
        # Try exact match first
        if gene_name in BETALACTAM_SUBSTRATE:
            return antibiotic in BETALACTAM_SUBSTRATE[gene_name]
        # Try prefix match (e.g., blaCTX-M-* not explicitly listed)
        for prefix in ["blaCTX-M", "blaCMY", "blaKPC", "blaNDM", "blaVIM", "blaIMP",
                       "blaOXA-48", "blaOXA-23", "blaOXA-24", "blaOXA-58", "blaEC", "blaADC"]:
            if gene_name.startswith(prefix) and prefix in BETALACTAM_SUBSTRATE:
                return antibiotic in BETALACTAM_SUBSTRATE[prefix]
        # Unmapped beta-lactamase: assume confers resistance (conservative)
        return True

    # Check aminoglycoside substrate mapping
    if drug_class == "AMINOGLYCOSIDE":
        if gene_name in AMINOGLYCOSIDE_SUBSTRATE:
            return antibiotic in AMINOGLYCOSIDE_SUBSTRATE[gene_name]
        # Prefix match
        for prefix in ["aac(3)", "aac(6')", "aph(3')", "aph(3'')", "aph(6)", "aph(4)",
                       "ant(2'')", "ant(3'')", "aadA", "armA", "rmtB"]:
            if gene_name.startswith(prefix):
                for key, abx_list in AMINOGLYCOSIDE_SUBSTRATE.items():
                    if key.startswith(prefix):
                        return antibiotic in abx_list
        # Unmapped aminoglycoside gene: assume confers resistance
        return True

    # Other drug classes: use drug-class matching (no substrate specificity issue)
    return True


if __name__ == "__main__":
    # Test
    tests = [
        ("blaTEM-1", "cefepime", "BETA-LACTAM", False),
        ("blaTEM-1", "ampicillin", "BETA-LACTAM", True),
        ("blaCTX-M-15", "cefepime", "BETA-LACTAM", True),
        ("blaCTX-M-27", "ceftazidime", "BETA-LACTAM", False),
        ("blaCMY-2", "cefoxitin", "BETA-LACTAM", True),
        ("blaEC", "cefepime", "BETA-LACTAM", False),
        ("blaKPC-3", "ertapenem", "BETA-LACTAM", True),
        ("aac(3)-IId", "gentamicin", "AMINOGLYCOSIDE", True),
        ("aac(3)-IId", "amikacin", "AMINOGLYCOSIDE", False),
        ("aph(3')-Ia", "streptomycin", "AMINOGLYCOSIDE", False),
        ("armA", "amikacin", "AMINOGLYCOSIDE", True),
    ]
    for gene, abx, dc, expected in tests:
        result = gene_confers_resistance(gene, abx, dc)
        status = "OK" if result == expected else "FAIL"
        print(f"  [{status}] {gene} -> {abx}: {result} (expected {expected})")
