#!/usr/bin/env python3
"""
Train per-species, per-antibiotic ML models for AMR phenotype prediction.

Includes:
1. Full model training (LR, RF, XGBoost, MLP)
2. Ablation study (incremental feature groups)
3. Learning curve analysis (sample size sensitivity)
4. SHAP feature importance

Usage:
    python3 train_models.py

Requires: build_feature_matrix.py to have been run first.
"""

import os
import sys
import json
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from collections import defaultdict
from pathlib import Path

from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import (
    balanced_accuracy_score, f1_score, matthews_corrcoef,
    roc_auc_score, average_precision_score, brier_score_loss,
    precision_score, recall_score
)

import xgboost as xgb

# ── Configuration ──
FEATURES_DIR = "/home/unnot/WGS/phenotype/features"
OUTPUT_DIR = "/home/unnot/WGS/phenotype/ml_results"
RANDOM_STATE = 42
N_FOLDS = 5

# Feature groups for ablation study
FEATURE_GROUPS = {
    "ARG_only": [
        "has_target_arg", "n_target_args", "n_amr_genes",
    ],
    "promoter": [
        "target_promoter_ldf_max", "target_promoter_ldf_mean", "target_promoter_ldf_min",
        "target_promoter_tf_sites_max", "target_promoter_tf_sites_mean",
        "target_promoter_distance_max", "target_promoter_distance_mean",
        "target_promoter_up_at_ratio_max", "target_promoter_up_at_ratio_mean",
    ],
    "RBS": [
        "target_rbs_expression_max", "target_rbs_expression_mean", "target_rbs_expression_min",
        "target_rbs_dg_total_max", "target_rbs_dg_total_mean",
        "target_rbs_dg_mrna_max", "target_rbs_dg_mrna_mean",
    ],
    "codon": [
        "target_cai_max", "target_cai_mean", "target_cai_min",
        "target_rare_codon_pct_max", "target_rare_codon_pct_mean",
        "target_rare_codon_clusters_max", "target_rare_codon_clusters_mean",
    ],
    "mobility": [
        "target_nearest_IS_distance_bp_max", "target_nearest_IS_distance_bp_mean",
        "target_has_nearby_IS", "target_IS_same_strand", "target_IS_opposite_strand",
        "target_any_on_plasmid", "target_any_on_prophage",
        "target_n_plasmid", "target_n_prophage", "target_n_chromosome",
        "genome_synteny_n_transposase_sum", "genome_synteny_n_transposase_max",
    ],
    "regulation": [
        "target_has_nearby_srna", "target_min_srna_distance",
        "target_nearest_sRNA_distance_bp_max", "target_nearest_sRNA_distance_bp_mean",
        "target_operon_size_max", "target_operon_size_mean",
        "target_operon_position_max", "target_operon_position_mean",
    ],
    "integron": [
        "target_any_in_integron",
        "target_nearest_integron_distance_bp_max", "target_nearest_integron_distance_bp_mean",
    ],
    "synteny": [
        "target_synteny_n_amr_max", "target_synteny_n_amr_mean",
        "target_synteny_n_stress_max", "target_synteny_n_stress_mean",
        "target_synteny_n_hypothetical_max", "target_synteny_n_hypothetical_mean",
        "target_synteny_n_transposase_max", "target_synteny_n_transposase_mean",
        "target_synteny_n_virulence_max", "target_synteny_n_virulence_mean",
    ],
    "genomic_context": [
        "target_gene_gc_max", "target_gene_gc_mean",
        "target_gc_deviation_max", "target_gc_deviation_mean",
        "genome_gc",
        "target_gene_copies_max", "target_gene_copies_mean",
        "genome_gene_copies_sum", "genome_gene_copies_max",
    ],
    "point_mutations": [
        "target_has_point_mutation", "genome_has_point_mutation",
    ],
    "population": [
        "target_identity_max", "target_identity_mean",
        "target_coverage_max", "target_coverage_mean",
    ],
    "functionality": [
        "arg_low_rbs", "arg_low_promoter", "arg_truncated",
        "arg_poor_adaptation", "arg_functional_score",
    ],
    "point_mutations_detail": [
        "pf_has_gyrA", "pf_has_gyrB", "pf_has_parC", "pf_has_parE",
        "pf_has_ampC_promoter", "pf_n_mutations",
        "pf_has_quinolone_mutation", "pf_n_quinolone_mutations",
    ],
}

# Ablation configurations: cumulative feature groups
ABLATION_CONFIGS = {
    "1_ARG_only": ["ARG_only"],
    "2_ARG+expression": ["ARG_only", "promoter", "RBS", "codon"],
    "3_ARG+expr+mobility": ["ARG_only", "promoter", "RBS", "codon", "mobility", "integron"],
    "4_ARG+expr+mob+regulation": ["ARG_only", "promoter", "RBS", "codon", "mobility", "integron", "regulation", "synteny"],
    "5_full_model": list(FEATURE_GROUPS.keys()),
    "6_full+gene_identity": list(FEATURE_GROUPS.keys()) + ["gene_identity"],
}

# Learning curve sample sizes
LEARNING_CURVE_SIZES = [50, 100, 200, 300, 500, 750, 1000, 2000, 3000, 5000, 7500, 10000]


def get_models():
    """Return dict of model name -> model instance."""
    return {
        "LR": LogisticRegression(
            max_iter=1000, class_weight="balanced", random_state=RANDOM_STATE, n_jobs=4
        ),
        "RF": RandomForestClassifier(
            n_estimators=300, max_depth=None, class_weight="balanced",
            random_state=RANDOM_STATE, n_jobs=4
        ),
        "XGB": xgb.XGBClassifier(
            n_estimators=100, max_depth=6, learning_rate=0.1,
            use_label_encoder=False, eval_metric="logloss",
            random_state=RANDOM_STATE, n_jobs=4
        ),
        "MLP": MLPClassifier(
            hidden_layer_sizes=(64, 32), max_iter=200, early_stopping=True,
            random_state=RANDOM_STATE
        ),
    }


def evaluate(y_true, y_pred, y_prob):
    """Compute evaluation metrics."""
    metrics = {
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "precision": precision_score(y_true, y_pred, zero_division=0),
        "recall": recall_score(y_true, y_pred, zero_division=0),
        "f1": f1_score(y_true, y_pred, zero_division=0),
        "mcc": matthews_corrcoef(y_true, y_pred),
        "roc_auc": roc_auc_score(y_true, y_prob) if len(set(y_true)) > 1 else np.nan,
        "pr_auc": average_precision_score(y_true, y_prob) if len(set(y_true)) > 1 else np.nan,
        "brier": brier_score_loss(y_true, y_prob),
        "n_samples": len(y_true),
        "n_positive": int(sum(y_true)),
        "n_negative": int(len(y_true) - sum(y_true)),
    }
    return metrics


def get_feature_columns(df, groups):
    """Get column names for given feature groups, filtering to columns that exist."""
    cols = []
    for g in groups:
        if g in FEATURE_GROUPS:
            cols.extend(FEATURE_GROUPS[g])
    # Also add MLST ST dummies if they exist
    if "population" in groups:
        st_cols = [c for c in df.columns if c.startswith("st_")]
        cols.extend(st_cols)
    # Add per-gene identity features if requested
    if "gene_identity" in groups:
        gene_cols = [c for c in df.columns if c.startswith("gene_")]
        cols.extend(gene_cols)
    # Filter to columns that actually exist in the dataframe
    return [c for c in cols if c in df.columns]


def prepare_xy(df, feature_cols):
    """Prepare X, y arrays with NaN handling."""
    X = df[feature_cols].copy()
    y = df["label"].values
    # Fill NaN with median
    for col in X.columns:
        median_val = X[col].median()
        X[col] = X[col].fillna(median_val if not np.isnan(median_val) else 0)
    return X.values, y


def train_and_evaluate(X, y, model_name, model_factory):
    """Train with stratified 5-fold CV and return mean metrics."""
    skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    fold_metrics = []

    for train_idx, test_idx in skf.split(X, y):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Scale
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        # Handle XGBoost class weight
        model = model_factory()
        if model_name == "XGB":
            n_neg = sum(y_train == 0)
            n_pos = sum(y_train == 1)
            if n_pos > 0:
                model.set_params(scale_pos_weight=n_neg / n_pos)

        try:
            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)
            y_prob = model.predict_proba(X_test)[:, 1]
            fold_metrics.append(evaluate(y_test, y_pred, y_prob))
        except Exception as e:
            print(f"      Warning: {model_name} failed on fold: {e}")
            continue

    if not fold_metrics:
        return None

    # Average across folds
    avg = {}
    for key in fold_metrics[0]:
        vals = [m[key] for m in fold_metrics if not np.isnan(m[key])]
        avg[key] = np.mean(vals) if vals else np.nan
        avg[f"{key}_std"] = np.std(vals) if vals else np.nan
    return avg


def run_learning_curve(df, feature_cols, antibiotic, species_dir):
    """Run learning curve analysis with subsampling."""
    print(f"    Learning curve for {antibiotic}...")
    results = []

    X, y = prepare_xy(df, feature_cols)
    n_total = len(y)

    for size in LEARNING_CURVE_SIZES:
        if size >= n_total:
            size = n_total

        # Subsample (stratified)
        if size < n_total:
            try:
                _, X_sub, _, y_sub = train_test_split(
                    X, y, test_size=size, stratify=y, random_state=RANDOM_STATE
                )
            except ValueError:
                continue
        else:
            X_sub, y_sub = X, y

        # Check minimum minority class
        minority = min(sum(y_sub == 0), sum(y_sub == 1))
        if minority < 10:
            continue

        # Train RF for learning curve
        metrics = train_and_evaluate(
            X_sub, y_sub, "RF",
            lambda: RandomForestClassifier(
                n_estimators=300, max_depth=None, class_weight='balanced',
                random_state=RANDOM_STATE, n_jobs=4
            )
        )
        if metrics:
            metrics["sample_size"] = size
            metrics["antibiotic"] = antibiotic
            results.append(metrics)

        if size == n_total:
            break

    return results


def run_shap_analysis(X_train, X_test, y_train, feature_names, output_path):
    """Run SHAP analysis on RF model."""
    try:
        import shap
        model = RandomForestClassifier(
            n_estimators=300, max_depth=None, class_weight='balanced',
            random_state=RANDOM_STATE, n_jobs=-1
        )
        model.fit(X_train, y_train)

        # Subsample for speed if dataset is large
        if len(X_test) > 1000:
            rng = np.random.RandomState(RANDOM_STATE)
            idx = rng.choice(len(X_test), 1000, replace=False)
            X_shap = X_test[idx] if hasattr(X_test, 'iloc') else X_test[idx]
        else:
            X_shap = X_test
        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(X_shap)

        # Handle multi-dimensional SHAP output from RF
        # Can be list [class0, class1] or 3D array (samples, features, classes)
        if isinstance(shap_values, list):
            shap_values = shap_values[1]
        elif shap_values.ndim == 3:
            shap_values = shap_values[:, :, 1]

        # Save mean absolute SHAP values
        mean_shap = np.abs(shap_values).mean(axis=0)
        shap_df = pd.DataFrame({
            "feature": feature_names,
            "mean_abs_shap": mean_shap
        }).sort_values("mean_abs_shap", ascending=False)
        shap_df.to_csv(output_path, sep="\t", index=False)
        return shap_df
    except ImportError:
        print("      SHAP not installed — skipping")
        return None
    except Exception as e:
        print(f"      SHAP failed: {e}")
        return None


def process_species(species_name, feature_matrix_path):
    """Run full analysis for one species."""
    print(f"\n{'='*70}")
    print(f"Processing {species_name}")
    print(f"{'='*70}")

    df = pd.read_csv(feature_matrix_path, sep="\t")
    species_dir = os.path.join(OUTPUT_DIR, species_name)
    os.makedirs(species_dir, exist_ok=True)

    antibiotics = df["antibiotic"].unique()
    print(f"  Antibiotics: {len(antibiotics)}")

    all_results = []
    all_ablation = []
    all_learning = []
    all_shap = []

    for abx in sorted(antibiotics):
        abx_df = df[df["antibiotic"] == abx].copy()
        n_total = len(abx_df)
        n_r = int(abx_df["label"].sum())
        n_s = n_total - n_r
        print(f"\n  {abx} (n={n_total}, R={n_r}, S={n_s})")

        # Get full feature set
        full_features = get_feature_columns(abx_df, list(FEATURE_GROUPS.keys()))
        if not full_features:
            print("    No features available — skipping")
            continue

        # ── 1. Full model comparison ──
        print("    Training full models...")
        X, y = prepare_xy(abx_df, full_features)

        models = get_models()
        best_model = None
        best_f1 = -1

        for model_name, model in models.items():
            metrics = train_and_evaluate(
                X, y, model_name, lambda mn=model_name: get_models()[mn]
            )
            if metrics:
                metrics["model"] = model_name
                metrics["antibiotic"] = abx
                metrics["config"] = "full"
                all_results.append(metrics)
                print(f"      {model_name}: F1={metrics['f1']:.3f}, AUC={metrics.get('roc_auc', 0):.3f}")
                if metrics["f1"] > best_f1:
                    best_f1 = metrics["f1"]
                    best_model = model_name

        # ── 2. Ablation study ──
        print("    Ablation study...")
        for config_name, groups in ABLATION_CONFIGS.items():
            config_features = get_feature_columns(abx_df, groups)
            if not config_features:
                continue
            X_abl, y_abl = prepare_xy(abx_df, config_features)
            metrics = train_and_evaluate(
                X_abl, y_abl, "RF",
                lambda: RandomForestClassifier(
                    n_estimators=300, max_depth=None, class_weight='balanced',
                    random_state=RANDOM_STATE, n_jobs=4
                )
            )
            if metrics:
                metrics["model"] = "RF"
                metrics["antibiotic"] = abx
                metrics["config"] = config_name
                metrics["n_features"] = len(config_features)
                all_ablation.append(metrics)
                print(f"      {config_name}: F1={metrics['f1']:.3f} ({len(config_features)} features)")

        # ── 3. Learning curve ──
        if n_total >= 100:
            lc_results = run_learning_curve(abx_df, full_features, abx, species_dir)
            all_learning.extend(lc_results)

        # ── 4. SHAP analysis (on best model) ──
        print("    SHAP analysis...")
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, stratify=y, random_state=RANDOM_STATE
        )
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)

        shap_path = os.path.join(species_dir, f"shap_{abx}.tsv")
        shap_df = run_shap_analysis(
            X_train_scaled, X_test_scaled, y_train,
            full_features, shap_path
        )
        if shap_df is not None:
            shap_df["antibiotic"] = abx
            all_shap.append(shap_df)

    # ── Save results ──
    if all_results:
        pd.DataFrame(all_results).to_csv(
            os.path.join(species_dir, "model_comparison.tsv"), sep="\t", index=False
        )
    if all_ablation:
        pd.DataFrame(all_ablation).to_csv(
            os.path.join(species_dir, "ablation_study.tsv"), sep="\t", index=False
        )
    if all_learning:
        pd.DataFrame(all_learning).to_csv(
            os.path.join(species_dir, "learning_curve.tsv"), sep="\t", index=False
        )
    if all_shap:
        pd.concat(all_shap).to_csv(
            os.path.join(species_dir, "shap_summary.tsv"), sep="\t", index=False
        )

    print(f"\n  Results saved to {species_dir}/")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Find all species with feature matrices
    species_dirs = sorted(Path(FEATURES_DIR).glob("*/feature_matrix.tsv"))
    if not species_dirs:
        print("No feature matrices found. Run build_feature_matrix.py first.")
        sys.exit(1)

    print(f"Found {len(species_dirs)} species with feature matrices")

    for fm_path in species_dirs:
        species_name = fm_path.parent.name
        process_species(species_name, str(fm_path))

    # ── Cross-species learning curve summary ──
    print("\n" + "=" * 70)
    print("Cross-species learning curve summary")
    print("=" * 70)
    all_lc = []
    for species_dir in Path(OUTPUT_DIR).iterdir():
        lc_path = species_dir / "learning_curve.tsv"
        if lc_path.exists():
            lc_df = pd.read_csv(str(lc_path), sep="\t")
            lc_df["species"] = species_dir.name
            all_lc.append(lc_df)
    if all_lc:
        combined = pd.concat(all_lc)
        combined.to_csv(
            os.path.join(OUTPUT_DIR, "learning_curve_all_species.tsv"),
            sep="\t", index=False
        )
        # Summary: mean F1 at each sample size across antibiotics
        summary = combined.groupby("sample_size").agg(
            mean_f1=("f1", "mean"),
            std_f1=("f1", "std"),
            mean_auc=("roc_auc", "mean"),
            n_models=("f1", "count"),
        ).reset_index()
        print(summary.to_string(index=False))

    print("\nDone. All results in:", OUTPUT_DIR)


if __name__ == "__main__":
    main()
