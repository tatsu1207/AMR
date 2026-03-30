#!/usr/bin/env python3
"""Export pre-trained RF models for web deployment."""
import os
import sys
import json
import numpy as np
import pandas as pd
import joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score

RANDOM_STATE = 42
FEATURES_DIR = "/home/unnot/WGS/phenotype/features"
OUTPUT_DIR = "/home/unnot/github/AMR/models"

SPECIES = [
    "Salmonella_enterica", "Escherichia_coli", "Klebsiella_pneumoniae",
    "Staphylococcus_aureus", "Acinetobacter_baumannii",
]

# Minimal feature groups (no mobility/regulation)
MINIMAL_PREFIXES = [
    # ARG_only
    "has_target_arg", "n_target_args", "n_amr_genes",
    # promoter
    "target_promoter_",
    # RBS
    "target_rbs_",
    # codon
    "target_cai_", "target_rare_codon_",
    # genomic_context
    "target_gene_gc_", "target_gc_deviation_", "genome_gc",
    "target_gene_copies_", "genome_gene_copies_",
    # point_mutations
    "target_has_point_mutation", "genome_has_point_mutation",
    # population (identity/coverage)
    "target_identity_", "target_coverage_",
    # functionality
    "arg_low_rbs", "arg_low_promoter", "arg_truncated",
    "arg_poor_adaptation", "arg_functional_score",
    # point_mutations_detail
    "pf_has_gyrA", "pf_has_gyrB", "pf_has_parC", "pf_has_parE",
    "pf_has_ampC_promoter", "pf_n_mutations",
    "pf_has_quinolone_mutation", "pf_n_quinolone_mutations",
    # A. baumannii custom
    "pf_has_gyrA_Ser83Leu",
    # MLST
    "st_",
]

META_COLS = {"sample_id", "antibiotic", "label", "species", "drug_class", "target_gene_list"}


def is_minimal_feature(col):
    for prefix in MINIMAL_PREFIXES:
        if col.startswith(prefix) or col == prefix:
            return True
    return False


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    manifest = {"species": {}, "total_models": 0}

    for sp in SPECIES:
        fm_path = os.path.join(FEATURES_DIR, sp, "feature_matrix.tsv")
        if not os.path.exists(fm_path):
            print(f"Skipping {sp}: no feature matrix")
            continue

        print(f"\n{'='*60}")
        print(f"  {sp}")
        print(f"{'='*60}")

        df = pd.read_csv(fm_path, sep='\t', low_memory=False)

        # Get minimal feature columns
        all_cols = [c for c in df.columns if c not in META_COLS]
        feat_cols = [c for c in all_cols if is_minimal_feature(c)
                     and df[c].dtype in ['float64', 'int64', 'float32', 'int32', 'uint8']]

        print(f"  Features: {len(feat_cols)} minimal features")

        sp_dir = os.path.join(OUTPUT_DIR, sp)
        os.makedirs(sp_dir, exist_ok=True)

        sp_manifest = {"antibiotics": {}, "n_features": len(feat_cols)}

        for abx in sorted(df['antibiotic'].unique()):
            abx_df = df[df['antibiotic'] == abx]
            n = len(abx_df)
            if n < 30:
                continue

            y = abx_df['label'].values
            n_r = int(y.sum())
            n_s = n - n_r
            if n_r < 5 or n_s < 5:
                continue

            # Get available features for this antibiotic
            avail_feats = [c for c in feat_cols if c in abx_df.columns
                          and abx_df[c].notna().any()]

            X = abx_df[avail_feats].values.copy()
            medians = {}
            for i, col in enumerate(avail_feats):
                med = np.nanmedian(X[:, i])
                if np.isnan(med):
                    med = 0.0
                medians[col] = float(med)
                X[np.isnan(X[:, i]), i] = med

            # Get drug class
            drug_class = "Unknown"
            if 'drug_class' in abx_df.columns:
                dc_vals = abx_df['drug_class'].dropna().unique()
                if len(dc_vals) > 0:
                    drug_class = str(dc_vals[0])

            # 5-fold CV for F1 estimate
            clf = RandomForestClassifier(
                n_estimators=300, max_depth=None, class_weight='balanced',
                random_state=RANDOM_STATE, n_jobs=4
            )
            try:
                cv_scores = cross_val_score(clf, X, y, cv=StratifiedKFold(5, shuffle=True, random_state=RANDOM_STATE),
                                           scoring='f1')
                cv_f1 = float(np.mean(cv_scores))
            except Exception:
                cv_f1 = 0.0

            # Train final model on ALL data
            final_model = RandomForestClassifier(
                n_estimators=300, max_depth=None, class_weight='balanced',
                random_state=RANDOM_STATE, n_jobs=4
            )
            final_model.fit(X, y)

            # Save model
            model_path = os.path.join(sp_dir, f"{abx}_rf.joblib")
            joblib.dump(final_model, model_path, compress=3)

            # Save metadata
            meta = {
                "feature_columns": avail_feats,
                "median_values": medians,
                "drug_class": drug_class,
                "n_samples": n,
                "n_resistant": n_r,
                "n_susceptible": n_s,
                "cv_f1": cv_f1,
            }
            meta_path = os.path.join(sp_dir, f"{abx}_meta.json")
            with open(meta_path, 'w') as f:
                json.dump(meta, f, indent=2)

            sp_manifest["antibiotics"][abx] = {
                "cv_f1": round(cv_f1, 3),
                "n_samples": n,
                "drug_class": drug_class,
            }
            manifest["total_models"] += 1

            print(f"  {abx:35s} F1={cv_f1:.3f}  n={n:5d} ({n_r} R / {n_s} S)  {len(avail_feats)} feats")

        manifest["species"][sp] = sp_manifest

    # Save manifest
    with open(os.path.join(OUTPUT_DIR, "model_manifest.json"), 'w') as f:
        json.dump(manifest, f, indent=2)

    print(f"\n{'='*60}")
    print(f"  Exported {manifest['total_models']} models to {OUTPUT_DIR}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
