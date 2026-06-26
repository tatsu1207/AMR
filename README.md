# AMR: Machine Learning-based AMR Phenotype Prediction from WGS

Uncovering species- and drug-class-specific antimicrobial resistance mechanisms from large-scale whole-genome sequencing data using discordance analysis and machine learning.

## Overview

A 16-step genomic annotation pipeline that characterizes resistance gene expression context (promoter strength, ribosome binding site efficiency, codon adaptation), mobile element context (plasmid, IS elements, integrons), and chromosomal point mutations to predict antimicrobial resistance phenotype from assembled bacterial genomes. The pipeline was applied to 20,296 genomes across five WHO priority pathogens, producing 107 species-antibiotic prediction models.

## Supported Species

| Species | Genomes | Antibiotics | Mean F1 (RF) |
|---------|---------|-------------|--------------|
| *Salmonella enterica* | 10,070 | 20 | 0.924 |
| *Escherichia coli* | 6,668 | 35 | 0.857 |
| *Klebsiella pneumoniae* | 1,824 | 26 | 0.873 |
| *Staphylococcus aureus* | 977 | 13 | 0.804 |
| *Acinetobacter baumannii* | 757 | 13 | 0.932 |

## Repository Structure

```
AMR/
├── annotate_one.sh          # 16-step annotation pipeline for a single assembly
├── build_feature_matrix.py  # Aggregate per-ARG annotations into ML feature matrix
├── train_models.py          # ML training: ablation, comparison, SHAP, learning curves
├── substrate_mapping.py     # Gene-variant-specific ARG-to-antibiotic mapping
├── regen_summary_fast.py    # Regenerate annotation summary from existing results
├── models/                  # Pre-trained RF models (.joblib) + metadata per species
│   ├── Salmonella_enterica/
│   ├── Escherichia_coli/
│   ├── Klebsiella_pneumoniae/
│   ├── Staphylococcus_aureus/
│   ├── Acinetobacter_baumannii/
│   └── model_manifest.json
├── databases/               # Reference databases (AMRFinderPlus, BPROM, etc.)
├── scripts/                 # Utility scripts (model export, etc.)
├── Dockerfile               # Docker container for streamlined prediction
├── docker-compose.yml
└── manuscript/              # Supplementary data tables
```

## Quick Start

### 1. Annotate a single genome assembly

```bash
# Usage: bash annotate_one.sh <input.fasta> <output_dir> <threads>
bash annotate_one.sh input.fasta results_dir 4
```

This runs the 16-step pipeline:
1. AMRFinderPlus (ARG detection + point mutations)
2. BPROM (promoter prediction)
3. OSTIR (RBS / translation initiation rate)
4. MobileElementFinder (IS elements)
5. MOB-recon (plasmid identification)
6. geNomad (mobile element classification)
7. skani (species identification)
8. Prodigal (gene prediction for operon/CAI)
9. Infernal/Rfam (sRNA detection in ARG flanks)
10. Operon structure analysis
11. Gene dosage calculation
12. Codon adaptation index (CAI)
13. GC content deviation
14. MLST (sequence typing)
15. *(Folded into Step 1 — PointFinder via AMRFinderPlus --organism)*
16. IntegronFinder
17. Gene synteny analysis

> **Note:** Step 15 (PointFinder) is no longer executed separately; point mutation detection is handled by AMRFinderPlus `--organism` in Step 1. The pipeline has 16 active steps (17 numbered for historical compatibility).

### 2. Build feature matrix

```bash
# Aggregate per-ARG annotations into a fixed-width feature vector per isolate-antibiotic pair
python3 build_feature_matrix.py --results-dir results --output-dir features
```

Output: `features/<species>/feature_matrix.tsv.gz` (~111 features per isolate-antibiotic pair)

### 3. Train and evaluate models

```bash
# Trains LR, RF, XGBoost, MLP with 5-fold stratified CV
# Runs ablation study, SHAP analysis, learning curves
python3 train_models.py
```

Output in `ml_results/<species>/`:
- `model_comparison.tsv` — all 4 classifiers, all metrics
- `ablation_study.tsv` — cumulative feature ablation
- `shap_<antibiotic>.tsv` — SHAP feature importance
- `learning_curve.tsv` — sample size vs performance

### 4. Predict with pre-trained models (Docker)

```bash
# Build and run the Docker container
docker-compose up --build

# The container runs a streamlined 6-step annotation pipeline
# and applies pre-trained RF models for all 107 prediction tasks
```

Pre-trained models are in `models/<species>/<antibiotic>_rf.joblib` with metadata in `<antibiotic>_meta.json`.

## Feature Groups

The ~111 features per isolate-antibiotic pair are organized into 12 biologically motivated groups:

| Group | Features | Description |
|-------|----------|-------------|
| ARG presence | 3 | has_target_arg, n_target_args, n_amr_genes |
| Promoter | 3 | BPROM-predicted promoter strength (LDF score) |
| RBS | 3 | OSTIR-predicted translation initiation rate |
| Codon adaptation | 3 | Codon adaptation index of target ARG |
| Mobility | 9 | IS elements, transposons near ARG |
| Regulation | 3 | Regulatory elements in ARG flanks |
| Integron | 3 | Integron cassette context |
| Synteny | 3 | ARG neighborhood conservation |
| Genomic context | 6 | Contig size, position, genome normalization |
| Point mutations | 2 | Binary quinolone mutation indicators |
| Point mutations (detail) | 6 | Locus-specific (gyrA, gyrB, parC, parE) |
| Population structure | 35 | MLST sequence types (top 30 + other) + gene identity |
| Functionality | 3 | Gene coverage, truncation indicators |

## Hyperparameters

| Classifier | Key Parameters |
|------------|---------------|
| Logistic Regression | L2 penalty, max_iter=1000, balanced class weights |
| Random Forest | 300 trees, unlimited depth, balanced class weights |
| XGBoost | 100 rounds, max_depth=6, learning_rate=0.1 |
| MLP | Hidden layers (64, 32), early stopping, max_iter=200 |

All models use `class_weight='balanced'` and `random_state=42`. No hyperparameter tuning was performed.

## Data Availability

- **Genome assemblies**: NCBI Sequence Read Archive (SRA) — accession list in `manuscript/table_S7_sra_accessions.tsv`
- **Pre-trained models**: `models/` directory (107 RF models as joblib files)
- **Feature matrices**: Available upon request or via the processed data release

## Requirements

- Python 3.12+
- scikit-learn 1.4, XGBoost 2.0, SHAP 0.45
- BBDuk, SPAdes, AMRFinderPlus, BPROM, OSTIR, ISEScan, MOB-suite, geNomad, IntegronFinder, Prodigal, Infernal, skani, MLST
- See `environment.yml` for full conda environment specification

## Hardware

All model training was performed on commodity CPU hardware (224 cores, 755 GB RAM). No GPU required. Training all 107 species-antibiotic models with 5-fold CV completes in approximately 4 hours.

## Citation

Cheon NJ, Nguyen XC, Unno T. Uncovering species- and drug-class-specific antimicrobial resistance mechanisms from large-scale whole-genome sequencing data using discordance analysis and machine learning. *Briefings in Bioinformatics* (under revision).

## License

[MIT License](LICENSE)
