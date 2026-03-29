#!/usr/bin/env python3
"""
Build per-species feature matrices for AMR phenotype prediction from WGS.

Reads annotation outputs (from annotate_one.sh) and phenotype labels,
then produces a fixed-width feature matrix per species with one row per
isolate × antibiotic.

Usage:
    python3 build_feature_matrix_v2.py [--results-dir DIR] [--phenotype CSV]
                                       [--output-dir DIR] [--species NAME]

Output per species:
    {output_dir}/{species}/feature_matrix.tsv
    {output_dir}/{species}/antibiotic_summary.tsv
"""

import argparse
import csv
import glob
import os
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

# ═══════════════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════════════

# Default paths (override via CLI args)
DEFAULT_RESULTS_DIR = "/home/unnot/WGS/phenotype/results"
DEFAULT_PHENOTYPE_CSV = "/home/unnot/WGS/phenotype/vault/phenotype/data/total_phenotype.csv"
DEFAULT_OUTPUT_DIR = "/home/unnot/WGS/phenotype/features"

# Species mapping: keyword in skani Ref_name → canonical species name
TARGET_SPECIES = {
    "salmonella": "Salmonella_enterica",
    "escherichia": "Escherichia_coli",
    # "campylobacter": "Campylobacter_jejuni",  # excluded: only 3 antibiotics qualify (1 ARG-only, 2 point-mutation-only)
    "klebsiella": "Klebsiella_pneumoniae",
    "mycobacterium": "Mycobacterium_tuberculosis",
    "shigella": "Shigella_sonnei",
    "staphylococcus": "Staphylococcus_aureus",
    "acinetobacter": "Acinetobacter_baumannii",
}

# Antibiotic → AMRFinder drug class keyword
DRUG_CLASS_MAP = {
    # Aminoglycosides
    "amikacin": "AMINOGLYCOSIDE", "gentamicin": "AMINOGLYCOSIDE",
    "kanamycin": "AMINOGLYCOSIDE", "tobramycin": "AMINOGLYCOSIDE",
    "streptomycin": "AMINOGLYCOSIDE", "spectinomycin": "AMINOGLYCOSIDE",
    "neomycin": "AMINOGLYCOSIDE", "apramycin": "AMINOGLYCOSIDE",
    "netilmicin": "AMINOGLYCOSIDE", "plazomicin": "AMINOGLYCOSIDE",
    "capreomycin": "AMINOGLYCOSIDE",
    # Beta-lactams
    "ampicillin": "BETA-LACTAM", "amoxicillin": "BETA-LACTAM",
    "amoxicillin_clavulanic_acid": "BETA-LACTAM",
    "ampicillin_sulbactam": "BETA-LACTAM",
    "piperacillin": "BETA-LACTAM", "piperacillin_tazobactam": "BETA-LACTAM",
    "ticarcillin": "BETA-LACTAM", "ticarcillin_clavulanic_acid": "BETA-LACTAM",
    "cefazolin": "BETA-LACTAM", "cefalexin": "BETA-LACTAM",
    "cephalexin": "BETA-LACTAM", "cefalotin": "BETA-LACTAM",
    "cephalothin": "BETA-LACTAM", "cefuroxime": "BETA-LACTAM",
    "cefoxitin": "BETA-LACTAM", "cefotetan": "BETA-LACTAM",
    "cefotaxime": "BETA-LACTAM", "cefotaxime_clavulanic_acid": "BETA-LACTAM",
    "ceftriaxone": "BETA-LACTAM", "ceftazidime": "BETA-LACTAM",
    "ceftazidime_avibactam": "BETA-LACTAM",
    "ceftazidime_clavulanic_acid": "BETA-LACTAM",
    "cefepime": "BETA-LACTAM", "ceftiofur": "BETA-LACTAM",
    "ceftolozane_tazobactam": "BETA-LACTAM", "ceftaroline": "BETA-LACTAM",
    "cefiderocol": "BETA-LACTAM", "cefixime": "BETA-LACTAM",
    "cefoperazone": "BETA-LACTAM", "cefpirome": "BETA-LACTAM",
    "aztreonam": "BETA-LACTAM",
    "imipenem": "BETA-LACTAM", "imipenem_relebactam": "BETA-LACTAM",
    "meropenem": "BETA-LACTAM", "meropenem_vaborbactam": "BETA-LACTAM",
    "ertapenem": "BETA-LACTAM", "doripenem": "BETA-LACTAM",
    "penicillin": "BETA-LACTAM", "benzylpenicillin": "BETA-LACTAM",
    "oxacillin": "BETA-LACTAM", "temocillin": "BETA-LACTAM",
    "sulbactam": "BETA-LACTAM",
    # Quinolones
    "ciprofloxacin": "QUINOLONE", "levofloxacin": "QUINOLONE",
    "moxifloxacin": "QUINOLONE", "nalidixic_acid": "QUINOLONE",
    "norfloxacin": "QUINOLONE", "ofloxacin": "QUINOLONE",
    "enrofloxacin": "QUINOLONE", "delafloxacin": "QUINOLONE",
    "pefloxacin": "QUINOLONE",
    # Tetracyclines
    "tetracycline": "TETRACYCLINE", "doxycycline": "TETRACYCLINE",
    "minocycline": "TETRACYCLINE", "tigecycline": "TETRACYCLINE",
    "eravacycline": "TETRACYCLINE", "omadacycline": "TETRACYCLINE",
    "chlortetracycline": "TETRACYCLINE", "oxytetracycline": "TETRACYCLINE",
    # Phenicols
    "chloramphenicol": "PHENICOL", "florfenicol": "PHENICOL",
    # Macrolides
    "erythromycin": "MACROLIDE", "azithromycin": "MACROLIDE",
    "clarithromycin": "MACROLIDE", "telithromycin": "MACROLIDE",
    "spiramycin": "MACROLIDE", "tilmicosin": "MACROLIDE", "tylosin": "MACROLIDE",
    # Lincosamides
    "clindamycin": "LINCOSAMIDE", "lincomycin": "LINCOSAMIDE",
    # Glycopeptides
    "vancomycin": "GLYCOPEPTIDE", "teicoplanin": "GLYCOPEPTIDE",
    # Oxazolidinones
    "linezolid": "OXAZOLIDINONE", "tedizolid": "OXAZOLIDINONE",
    # Folate pathway
    "trimethoprim": "TRIMETHOPRIM",
    "trimethoprim_sulfamethoxazole": "TRIMETHOPRIM",
    "sulfamethoxazole_trimethoprim": "TRIMETHOPRIM",
    "sulfamethoxazole": "SULFONAMIDE", "sulfisoxazole": "SULFONAMIDE",
    "sulfonamide": "SULFONAMIDE", "sulfonamides": "SULFONAMIDE",
    "sulphadimethoxine": "SULFONAMIDE",
    # Polymyxins
    "colistin": "COLISTIN", "polymyxin_b": "COLISTIN", "colomycin": "COLISTIN",
    # Others
    "fosfomycin": "FOSFOMYCIN", "phosphomycin": "FOSFOMYCIN",
    "nitrofurantoin": "NITROFURANTOIN", "strofurantoin": "NITROFURANTOIN",
    "furazolidone": "NITROFURANTOIN",
    "rifampin": "RIFAMYCIN", "daptomycin": "LIPOPEPTIDE",
    "mupirocin": "MUPIROCIN", "fusidic_acid": "FUSIDIC ACID",
    "metronidazole": "NITROIMIDAZOLE",
    # TB drugs
    "isoniazid": "ISONIAZID", "ethambutol": "ETHAMBUTOL",
    "pyrazinamide": "PYRAZINAMIDE", "ethionamide": "ETHIONAMIDE",
    "cycloserine": "CYCLOSERINE",
    # Streptogramins / Pleuromutilins
    "quinupristin_dalfopristin": "STREPTOGRAMIN",
    "synercid": "STREPTOGRAMIN", "pristimycin": "STREPTOGRAMIN",
    "tiamulin": "PLEUROMUTILIN",
    # Antifungals
    "fluconazole": "ANTIFUNGAL", "itraconazole": "ANTIFUNGAL",
    "voriconazole": "ANTIFUNGAL", "posaconazole": "ANTIFUNGAL",
    "isavuconazole": "ANTIFUNGAL", "amphotericin_b": "ANTIFUNGAL",
    "caspofungin": "ANTIFUNGAL", "anidulafungin": "ANTIFUNGAL",
    "micafungin": "ANTIFUNGAL", "flucytosine": "ANTIFUNGAL",
    "fidaxomicin": "ANTIFUNGAL",
}

# Filtering thresholds
MIN_MINORITY = 30       # Minimum samples in minority class
MAX_IMBALANCE = 20      # Maximum majority/minority ratio

# Numeric columns from arg_context_summary.tsv to aggregate per ARG
NUMERIC_COLS = [
    "identity", "coverage",
    "promoter_ldf", "promoter_tf_sites", "promoter_distance", "promoter_up_at_ratio",
    "rbs_expression", "rbs_dg_total", "rbs_dg_mrna",
    "cai", "rare_codon_pct", "rare_codon_clusters",
    "gene_gc", "genome_gc", "gc_deviation",
    "gene_copies",
    "operon_size", "operon_position",
    "nearest_IS_distance_bp",
    "nearest_sRNA_distance_bp",
    "nearest_integron_distance_bp",
    "synteny_n_transposase", "synteny_n_amr", "synteny_n_stress",
    "synteny_n_hypothetical", "synteny_n_virulence",
]

# Binary/categorical columns from arg_context_summary.tsv
BINARY_COLS = ["on_plasmid", "on_prophage", "in_integron"]


# ═══════════════════════════════════════════════════════════════════════════════
# Helper functions
# ═══════════════════════════════════════════════════════════════════════════════

def safe_float(val):
    """Convert to float, returning NaN for empty/invalid values."""
    if val is None or val == "" or val == "NA":
        return np.nan
    try:
        return float(val)
    except (ValueError, TypeError):
        return np.nan


def detect_species(results_dir, sample_id):
    """Detect species from skani results, with MLST fallback."""
    # Primary: skani
    skani_path = os.path.join(results_dir, sample_id, "species", "skani_results.tsv")
    if os.path.exists(skani_path) and os.path.getsize(skani_path) > 0:
        with open(skani_path) as f:
            for line in f:
                line_lower = line.lower()
                for key, species in TARGET_SPECIES.items():
                    if key in line_lower:
                        return species

    # Fallback: MLST scheme
    mlst_path = os.path.join(results_dir, sample_id, "mlst", "mlst_results.tsv")
    if os.path.exists(mlst_path) and os.path.getsize(mlst_path) > 0:
        with open(mlst_path) as f:
            parts = f.readline().strip().split("\t")
            if len(parts) >= 2:
                scheme = parts[1].lower()
                for key, species in TARGET_SPECIES.items():
                    if key in scheme:
                        return species

    return "unknown"


def parse_summary(summary_path):
    """Parse arg_context_summary.tsv into a list of row dicts."""
    rows = []
    if not os.path.exists(summary_path):
        return rows
    with open(summary_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rows.append(row)
    return rows


def parse_pointfinder(results_dir, sample_id):
    """Parse PointFinder results into per-sample mutation features."""
    pf_path = os.path.join(results_dir, sample_id, "pointfinder", "PointFinder_results.txt")
    features = {
        "pf_has_gyrA": 0, "pf_has_gyrB": 0,
        "pf_has_parC": 0, "pf_has_parE": 0,
        "pf_has_ampC_promoter": 0,
        "pf_n_mutations": 0,
        "pf_has_quinolone_mutation": 0,
        "pf_n_quinolone_mutations": 0,
    }
    if not os.path.exists(pf_path) or os.path.getsize(pf_path) == 0:
        return features
    try:
        with open(pf_path) as f:
            f.readline()  # skip header
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue
                mutation, resistance = parts[0], parts[3].lower()
                features["pf_n_mutations"] += 1
                if mutation.startswith("gyrA"):
                    features["pf_has_gyrA"] = 1
                elif mutation.startswith("gyrB"):
                    features["pf_has_gyrB"] = 1
                elif mutation.startswith("parC"):
                    features["pf_has_parC"] = 1
                elif mutation.startswith("parE"):
                    features["pf_has_parE"] = 1
                elif "ampc" in mutation.lower():
                    features["pf_has_ampC_promoter"] = 1
                if "nalidixic" in resistance or "ciprofloxacin" in resistance:
                    features["pf_has_quinolone_mutation"] = 1
                    features["pf_n_quinolone_mutations"] += 1
    except Exception:
        pass
    return features


def aggregate_features(arg_rows, target_drug_class=None):
    """
    Aggregate per-ARG annotation features into a fixed-width feature vector.

    For each isolate × antibiotic pair, this computes:
    - Count features (n_amr_genes, n_target_args, has_target_arg)
    - Per-ARG numeric features aggregated as max/mean/min for target-class ARGs
    - Genome-wide aggregate features (total gene copies, synteny counts)
    - Binary features (plasmid, prophage, integron)
    - Contig type distribution (n_plasmid, n_prophage, n_chromosome)
    - sRNA proximity features
    - Interaction features (arg_functional_score, arg_truncated, etc.)
    - Isolate-level features (genome_gc, MLST, gene list)
    """
    features = {}

    # ── Count features ──
    amr_args = [r for r in arg_rows if r.get("element_type") == "AMR"]
    features["n_total_elements"] = len(arg_rows)
    features["n_amr_genes"] = len(amr_args)

    # Target drug class matching
    if target_drug_class:
        target_args = [r for r in amr_args
                       if target_drug_class.upper() in r.get("drug_class", "").upper()]
    else:
        target_args = amr_args

    features["n_target_args"] = len(target_args)
    features["has_target_arg"] = 1 if target_args else 0

    # ── Numeric features for target ARGs (max/mean, some also min) ──
    for col in NUMERIC_COLS:
        vals = [v for v in (safe_float(r.get(col, "")) for r in target_args)
                if not np.isnan(v)]
        if vals:
            features[f"target_{col}_max"] = max(vals)
            features[f"target_{col}_mean"] = np.mean(vals)
            if col in ("promoter_ldf", "rbs_expression", "cai"):
                features[f"target_{col}_min"] = min(vals)
        else:
            features[f"target_{col}_max"] = np.nan
            features[f"target_{col}_mean"] = np.nan
            if col in ("promoter_ldf", "rbs_expression", "cai"):
                features[f"target_{col}_min"] = np.nan

    # ── Genome-wide aggregate features (across ALL ARGs) ──
    for col in ["gene_copies", "synteny_n_amr", "synteny_n_transposase"]:
        vals = [v for v in (safe_float(r.get(col, "")) for r in amr_args)
                if not np.isnan(v)]
        features[f"genome_{col}_sum"] = sum(vals) if vals else 0
        features[f"genome_{col}_max"] = max(vals) if vals else 0

    # ── Binary features for target ARGs ──
    for col in BINARY_COLS:
        vals = [r.get(col, "False") for r in target_args]
        features[f"target_any_{col}"] = 1 if any(v == "True" for v in vals) else 0

    # ── Contig type distribution ──
    contig_types = [r.get("contig_type", "chromosome") for r in target_args]
    features["target_n_plasmid"] = sum(1 for c in contig_types if c == "plasmid")
    features["target_n_prophage"] = sum(1 for c in contig_types if c == "prophage")
    features["target_n_chromosome"] = sum(1 for c in contig_types if c == "chromosome")

    # ── IS element orientation ──
    is_orients = [r.get("nearest_IS_orientation", "") for r in target_args]
    features["target_has_nearby_IS"] = 1 if any(r.get("nearest_IS_distance_bp", "") for r in target_args) else 0
    features["target_IS_same_strand"] = sum(1 for o in is_orients if o == "same")
    features["target_IS_opposite_strand"] = sum(1 for o in is_orients if o == "opposite")

    # ── sRNA proximity ──
    srna_dists = [v for v in (safe_float(r.get("nearest_sRNA_distance_bp", ""))
                              for r in target_args) if not np.isnan(v)]
    features["target_has_nearby_srna"] = 1 if any(d <= 5000 for d in srna_dists) else 0
    features["target_min_srna_distance"] = min(srna_dists) if srna_dists else np.nan

    # ── Isolate-level features ──
    if arg_rows:
        features["mlst_scheme"] = arg_rows[0].get("mlst_scheme", "")
        features["mlst_st"] = arg_rows[0].get("mlst_st", "")
        features["genome_gc"] = safe_float(arg_rows[0].get("genome_gc", ""))
    else:
        features["mlst_scheme"] = ""
        features["mlst_st"] = ""
        features["genome_gc"] = np.nan

    # ── Gene name list (for potential downstream use) ──
    target_gene_names = sorted(set(r.get("gene", "") for r in target_args))
    features["target_gene_list"] = ";".join(target_gene_names)

    # ── Interaction features: ARG functionality indicators ──
    has_arg = features["has_target_arg"]
    if target_args:
        rbs_vals = [v for v in (safe_float(r.get("rbs_expression", "")) for r in target_args)
                    if not np.isnan(v)]
        cov_vals = [v for v in (safe_float(r.get("coverage", "")) for r in target_args)
                    if not np.isnan(v)]
        gc_dev_vals = [v for v in (safe_float(r.get("gc_deviation", "")) for r in target_args)
                       if not np.isnan(v)]
        pldf_vals = [v for v in (safe_float(r.get("promoter_ldf", "")) for r in target_args)
                     if not np.isnan(v)]

        # Binary: ARG present but likely non-functional
        features["arg_low_rbs"] = 1 if (has_arg and rbs_vals and max(rbs_vals) < 50) else 0
        features["arg_low_promoter"] = 1 if (has_arg and pldf_vals and max(pldf_vals) < 3.0) else 0
        features["arg_truncated"] = 1 if (has_arg and cov_vals and max(cov_vals) < 98.0) else 0
        features["arg_poor_adaptation"] = 1 if (has_arg and gc_dev_vals and max(gc_dev_vals) > 0.08) else 0

        # Composite: high = likely functional ARG
        rbs_norm = min(max(rbs_vals) / 500.0, 1.0) if rbs_vals else 0
        cov_norm = max(cov_vals) / 100.0 if cov_vals else 0
        gc_dev_pen = 1.0 - min(max(gc_dev_vals) / 0.15, 1.0) if gc_dev_vals else 0.5
        pldf_norm = min(max(pldf_vals) / 5.0, 1.0) if pldf_vals else 0
        features["arg_functional_score"] = has_arg * rbs_norm * cov_norm * gc_dev_pen * pldf_norm
    else:
        features["arg_low_rbs"] = 0
        features["arg_low_promoter"] = 0
        features["arg_truncated"] = 0
        features["arg_poor_adaptation"] = 0
        features["arg_functional_score"] = 0.0

    return features


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def process_species(species_name, samples, pheno, results_dir, output_dir):
    """Build feature matrix for one species."""
    print(f"\n{'='*60}")
    print(f"Building feature matrix for {species_name} (n={len(samples)})")
    print(f"{'='*60}")

    species_dir = os.path.join(output_dir, species_name)
    os.makedirs(species_dir, exist_ok=True)

    antibiotics = list(pheno.columns)

    # Load annotation data
    print("  Loading annotation summaries...")
    all_summaries = {}
    all_pointfinder = {}
    for sample_id in samples:
        summary_path = os.path.join(results_dir, sample_id, "arg_context_summary.tsv")
        all_summaries[sample_id] = parse_summary(summary_path)
        all_pointfinder[sample_id] = parse_pointfinder(results_dir, sample_id)

    # Select antibiotics with sufficient data
    print("  Selecting antibiotics...")
    valid_antibiotics = []
    for abx in antibiotics:
        if abx not in DRUG_CLASS_MAP:
            continue
        vals = pheno.loc[samples, abx].dropna()
        if len(vals) == 0:
            continue
        labels = vals.map({0: 0, 1: 1, 2: 1}).dropna()
        if len(labels) == 0:
            continue
        n_r, n_s = int(labels.sum()), int(len(labels) - labels.sum())
        minority = min(n_r, n_s)
        if minority < MIN_MINORITY:
            continue
        imbalance = max(n_r, n_s) / minority if minority > 0 else float("inf")
        if imbalance > MAX_IMBALANCE:
            continue
        valid_antibiotics.append(abx)
        print(f"    {abx}: n={len(labels)}, R={n_r}, S={n_s}, ratio={imbalance:.1f}")

    if not valid_antibiotics:
        print("  No valid antibiotics — skipping")
        return

    # Build feature matrix
    print(f"  Building features for {len(valid_antibiotics)} antibiotics...")
    rows = []
    for abx in valid_antibiotics:
        drug_class = DRUG_CLASS_MAP.get(abx, "")
        for sample_id in samples:
            val = pheno.loc[sample_id, abx]
            if pd.isna(val):
                continue
            # S (0) = susceptible, I (1) and R (2) = non-susceptible
            label = 0 if val == 0 else 1
            feats = aggregate_features(all_summaries[sample_id], target_drug_class=drug_class)
            feats.update(all_pointfinder[sample_id])
            feats["sample_id"] = sample_id
            feats["antibiotic"] = abx
            feats["drug_class"] = drug_class
            feats["label"] = label
            rows.append(feats)

    df = pd.DataFrame(rows)

    # Encode MLST ST as one-hot (top 30 + other)
    if "mlst_st" in df.columns:
        top_sts = df["mlst_st"].value_counts().head(30).index.tolist()
        df["mlst_st_encoded"] = df["mlst_st"].apply(
            lambda x: x if x in top_sts else "other"
        )
        st_dummies = pd.get_dummies(df["mlst_st_encoded"], prefix="st")
        df = pd.concat([df, st_dummies], axis=1)

    # Encode per-gene identity as one-hot (top N genes by frequency)
    if "target_gene_list" in df.columns:
        from collections import Counter
        gene_counts = Counter()
        for gene_list in df["target_gene_list"].dropna():
            for g in str(gene_list).split(";"):
                g = g.strip()
                if g:
                    gene_counts[g] += 1
        # Keep genes present in at least 20 samples
        frequent_genes = [g for g, n in gene_counts.items() if n >= 20]
        for gene in sorted(frequent_genes):
            col_name = f"gene_{gene}"
            df[col_name] = df["target_gene_list"].apply(
                lambda x: 1 if isinstance(x, str) and gene in x.split(";") else 0
            )
        print(f"  Gene identity features: {len(frequent_genes)} genes (present in >=20 samples)")

    # Save feature matrix
    out_path = os.path.join(species_dir, "feature_matrix.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"  Saved: {out_path}")
    print(f"  Shape: {df.shape[0]} rows × {df.shape[1]} columns")
    print(f"  Antibiotics: {len(valid_antibiotics)}")

    # Save antibiotic summary
    summary_path = os.path.join(species_dir, "antibiotic_summary.tsv")
    with open(summary_path, "w") as f:
        f.write("antibiotic\tdrug_class\tn_total\tn_resistant\tn_susceptible\n")
        for abx in valid_antibiotics:
            sub = df[df["antibiotic"] == abx]
            n_r = int(sub["label"].sum())
            n_s = int(len(sub) - n_r)
            f.write(f"{abx}\t{DRUG_CLASS_MAP.get(abx,'')}\t{len(sub)}\t{n_r}\t{n_s}\n")
    print(f"  Saved: {summary_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Build per-species feature matrices for AMR phenotype prediction."
    )
    parser.add_argument("--results-dir", default=DEFAULT_RESULTS_DIR,
                        help="Directory containing per-sample annotation results")
    parser.add_argument("--phenotype", default=DEFAULT_PHENOTYPE_CSV,
                        help="Phenotype CSV file (Run × antibiotic, 0=S, 2=R)")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR,
                        help="Output directory for feature matrices")
    parser.add_argument("--species", default=None,
                        help="Process only this species (e.g. Salmonella_enterica)")
    parser.add_argument("--min-samples", type=int, default=50,
                        help="Minimum samples to process a species (default: 50)")
    args = parser.parse_args()

    print("Loading phenotype data...")
    pheno = pd.read_csv(args.phenotype).set_index("Run")
    print(f"  {len(pheno)} isolates, {len(pheno.columns)} antibiotics")

    # Find completed samples
    completed = []
    for sample_dir in glob.glob(os.path.join(args.results_dir, "*/arg_context_summary.tsv")):
        sample_id = os.path.basename(os.path.dirname(sample_dir))
        if sample_id in pheno.index:
            completed.append(sample_id)
    print(f"  {len(completed)} samples with completed annotations and phenotype data")

    # Detect species
    print("Detecting species...")
    species_map = {}
    for i, sample_id in enumerate(completed):
        species_map[sample_id] = detect_species(args.results_dir, sample_id)
        if (i + 1) % 1000 == 0:
            print(f"  {i+1}/{len(completed)}...")

    species_counts = defaultdict(int)
    for sp in species_map.values():
        species_counts[sp] += 1
    print("  Species distribution:")
    for sp, n in sorted(species_counts.items(), key=lambda x: -x[1]):
        print(f"    {sp}: {n}")

    # Process species
    target_species = sorted(set(TARGET_SPECIES.values()))
    if args.species:
        target_species = [args.species]

    for species_name in target_species:
        samples = [s for s, sp in species_map.items() if sp == species_name]
        if len(samples) < args.min_samples:
            print(f"\nSkipping {species_name} — only {len(samples)} samples")
            continue
        process_species(species_name, samples, pheno, args.results_dir, args.output_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
