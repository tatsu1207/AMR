"""Load pre-trained models and predict phenotypes."""
import os
import json
import joblib
import numpy as np
from app.config import MODELS_DIR, SPECIES_DISPLAY
from app.models.schemas import AntibioticPrediction, GeneResult, PredictionResults

# Global model cache
_models = {}  # {species: {antibiotic: {'model': rf, 'meta': dict}}}
_manifest = {}


def load_all_models():
    """Load all pre-trained models at startup."""
    global _models, _manifest

    manifest_path = os.path.join(MODELS_DIR, "model_manifest.json")
    if os.path.exists(manifest_path):
        with open(manifest_path) as f:
            _manifest = json.load(f)

    for species_dir in os.listdir(MODELS_DIR):
        sp_path = os.path.join(MODELS_DIR, species_dir)
        if not os.path.isdir(sp_path) or species_dir == '__pycache__':
            continue

        _models[species_dir] = {}
        for fname in os.listdir(sp_path):
            if fname.endswith('_rf.joblib'):
                abx = fname.replace('_rf.joblib', '')
                model_path = os.path.join(sp_path, fname)
                meta_path = os.path.join(sp_path, f"{abx}_meta.json")

                try:
                    model = joblib.load(model_path)
                    meta = {}
                    if os.path.exists(meta_path):
                        with open(meta_path) as f:
                            meta = json.load(f)
                    _models[species_dir][abx] = {'model': model, 'meta': meta}
                except Exception as e:
                    print(f"Failed to load model {species_dir}/{abx}: {e}")

    total = sum(len(v) for v in _models.values())
    print(f"Loaded {total} models across {len(_models)} species")


def predict_phenotypes(species: str, feature_data: dict, mlst_st: str,
                       amr_results: dict = None) -> dict:
    """Run predictions for all antibiotics of a species."""
    if species not in _models:
        raise ValueError(f"No models available for {species}")

    species_models = _models[species]
    predictions = []

    for abx, model_data in sorted(species_models.items()):
        model = model_data['model']
        meta = model_data['meta']
        feature_cols = meta.get('feature_columns', [])
        medians = meta.get('median_values', {})
        drug_class = meta.get('drug_class', 'Unknown')

        # Build feature vector
        X = np.array([feature_data.get(col, medians.get(col, 0.0)) for col in feature_cols])
        X = np.nan_to_num(X, nan=0.0).reshape(1, -1)

        # Predict
        pred_class = model.predict(X)[0]
        pred_proba = model.predict_proba(X)[0]
        # Ensure we use the correct class index for Resistant (class 1)
        resistant_idx = list(model.classes_).index(1) if 1 in model.classes_ else -1
        if resistant_idx >= 0:
            prob_resistant = pred_proba[resistant_idx]
        else:
            prob_resistant = pred_proba[1] if len(pred_proba) > 1 else pred_proba[0]

        is_resistant = prob_resistant >= 0.5

        # Confidence based on distance from decision boundary
        # Show as confidence in the predicted outcome
        pred_confidence = prob_resistant if is_resistant else (1.0 - prob_resistant)
        if pred_confidence > 0.85:
            confidence = "High"
        elif pred_confidence > 0.7:
            confidence = "Moderate"
        else:
            confidence = "Low"

        # Identify key determinants from feature importance
        # Only show genes that are relevant to this specific antibiotic
        key_genes = []
        key_mutations = []
        if amr_results:
            # Get feature importances to find which genes actually matter
            importances = dict(zip(feature_cols, model.feature_importances_))
            class_genes = amr_results.get('genes_by_class', {}).get(drug_class, [])

            # Only include genes whose corresponding features have nonzero importance
            for gene in class_genes:
                gene_lower = gene.lower()
                # Check if any feature referencing this gene has importance
                relevant = any(
                    gene_lower in col.lower() and importances.get(col, 0) > 0.01
                    for col in feature_cols
                )
                if relevant:
                    key_genes.append(gene)
                if len(key_genes) >= 3:
                    break

            # If no gene-specific features found, fall back to showing
            # class-level genes only when predicted resistant
            if not key_genes and is_resistant:
                key_genes = class_genes[:3]

            key_mutations = amr_results.get('mutations', [])
            if 'QUINOLONE' in drug_class.upper():
                key_mutations = [m for m in key_mutations if 'gyr' in m.lower() or 'par' in m.lower()]
            else:
                key_mutations = []

        predictions.append(AntibioticPrediction(
            antibiotic=abx,
            drug_class=drug_class,
            prediction="Resistant" if is_resistant else "Susceptible",
            probability=round(float(prob_resistant), 3),
            confidence=confidence,
            key_genes=key_genes,
            key_mutations=key_mutations,
        ))

    # Build detected genes list
    detected_genes = []
    if amr_results and 'all_genes' in amr_results:
        for g in amr_results['all_genes']:
            detected_genes.append(GeneResult(
                gene=g.get('gene', ''),
                drug_class=g.get('drug_class', ''),
                identity=g.get('identity', 0),
                coverage=g.get('coverage', 0),
                on_plasmid=g.get('on_plasmid', False),
            ))

    n_r = sum(1 for p in predictions if p.prediction == "Resistant")
    n_s = sum(1 for p in predictions if p.prediction == "Susceptible")

    return PredictionResults(
        job_id="",
        species=species,
        species_display=SPECIES_DISPLAY.get(species, species),
        mlst_st=mlst_st or "unknown",
        n_antibiotics=len(predictions),
        n_resistant=n_r,
        n_susceptible=n_s,
        predictions=predictions,
        detected_genes=detected_genes,
    ).model_dump()
