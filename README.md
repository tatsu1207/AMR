# AMR: Machine Learning-based AMR Phenotype Prediction from WGS

Machine learning uncovers species-specific and drug-class-specific mechanisms of antimicrobial resistance from whole-genome sequencing.

## Overview

A 17-step genomic annotation pipeline that characterizes resistance gene expression context (promoter strength, ribosome binding site, codon adaptation), mobile element context (plasmid, IS elements, integrons), and chromosomal point mutations to predict antimicrobial resistance phenotype from assembled bacterial genomes.

## Supported Species

- *Salmonella enterica* (20 antibiotics)
- *Escherichia coli* (35 antibiotics)
- *Klebsiella pneumoniae* (26 antibiotics)
- *Staphylococcus aureus* (13 antibiotics)
- *Acinetobacter baumannii* (13 antibiotics)

## Scripts

| Script | Description |
|--------|-------------|
| `annotate_one.sh` | 17-step annotation pipeline for a single assembly |
| `build_feature_matrix.py` | Aggregate per-ARG annotations into ML feature matrix |
| `train_models.py` | ML training: ablation, algorithm comparison, SHAP, learning curves |
| `substrate_mapping.py` | Gene-variant-specific ARG-to-antibiotic mapping |

## Usage

```bash
# Annotate a single assembly
bash annotate_one.sh input.fasta results_dir 4

# Build feature matrix
python3 build_feature_matrix.py --results-dir results --output-dir features

# Train models
python3 train_models.py
```

## Citation

[Manuscript in preparation]
