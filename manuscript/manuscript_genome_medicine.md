---
title: "Machine learning uncovers species-specific and drug-class-specific mechanisms of antimicrobial resistance from whole-genome sequencing"
---

## Table of contents

- [Abstract](#abstract)
- [Background](#background)
- [Methods](#methods)
  - [Data collection](#data-collection)
  - [Species and antibiotic selection](#species-and-antibiotic-selection)
  - [Genomic annotation pipeline](#genomic-annotation-pipeline)
  - [Feature matrix construction and machine learning](#feature-matrix-construction-and-machine-learning)
  - [Discordance analysis](#discordance-analysis)
  - [Software and reproducibility](#software-and-reproducibility)
- [Results](#results)
  - [Dataset composition: source and geographic distribution](#dataset-composition-source-and-geographic-distribution)
  - [ARG presence alone fails to predict phenotype](#arg-presence-alone-fails-to-predict-phenotype)
  - [Expression features resolve silent gene discordance](#expression-features-resolve-silent-gene-discordance)
  - [PointFinder mutations resolve missing mechanism discordance](#pointfinder-mutations-resolve-missing-mechanism-discordance)
  - [Full model performance](#full-model-performance)
  - [What remains unexplained](#what-remains-unexplained)
  - [SHAP feature importance reflects mechanism types](#shap-feature-importance-reflects-mechanism-types)
  - [Learning curves reveal species-specific data requirements](#learning-curves-reveal-species-specific-data-requirements)
- [Discussion](#discussion)
  - [A causal framework for WGS-based AMR prediction](#a-causal-framework-for-wgs-based-amr-prediction)
  - [Expression features address a fundamental gap in gene catalogs](#expression-features-address-a-fundamental-gap-in-gene-catalogs)
  - [PointFinder integration is non-negotiable for quinolone prediction](#pointfinder-integration-is-non-negotiable-for-quinolone-prediction)
  - [High resistance prevalence inflates prediction performance for A. baumannii](#high-resistance-prevalence-inflates-prediction-performance-for-a-baumannii)
  - [Remaining challenges and future directions](#remaining-challenges-and-future-directions)
- [Conclusions](#conclusions)
- [Declarations](#declarations)

---

# Abstract

**Background:** The presence of an antimicrobial resistance gene does not guarantee phenotypic resistance — silent genes, chromosomal mutations, and expression context create widespread genotype-phenotype discordance. However, the molecular mechanisms underlying this discordance have not been systematically quantified across diverse bacterial species and drug classes.

**Methods:** We developed a 17-step genomic annotation pipeline capturing promoter strength, ribosome binding site (RBS) expression, codon adaptation, mobile element context, gene synteny, and chromosomal point mutations for each resistance gene. The pipeline was applied to 27,856 bacterial genomes across five WHO priority species (*S. enterica*, *E. coli*, *K. pneumoniae*, *S. aureus*, *A. baumannii*), covering 107 species-antibiotic combinations. Machine learning models were trained to predict resistance phenotype, and genotype-phenotype discordance was analyzed by matching each resistance gene variant to the specific antibiotics it is known to hydrolyze or modify, rather than treating all genes within a drug class as equivalent.

**Results:** Gene-variant-specific discordance analysis revealed that 3-43% of genotype-phenotype predictions are discordant when using gene presence alone. Silent resistance genes — carrying the correct gene variant but remaining susceptible — showed 30-85% lower predicted translation initiation rates, 4-14% lower gene coverage, and weaker promoter activity than functional copies, indicating that insufficient protein production through poor ribosome binding, gene truncation, and weak transcription is a major source of false-positive genotype-based resistance predictions. Quinolone resistance was 85-100% mutation-mediated across all five species, with chromosomal gyrA/parC mutations explaining resistance in the absence of acquired ARGs. Prediction performance varied by species, with *A. baumannii* achieving the highest F1 (0.941) — driven partly by high baseline resistance prevalence — while mechanistically diverse *E. coli* required the full expression feature set to achieve F1=0.871. Overall, genomic context features improved prediction by 0.062-0.168 F1 over gene presence alone.

**Conclusions:** Resistance mechanisms fall into three categories: antibiotics predictable by gene presence alone (tetracyclines, sulfonamides, ampicillin), those requiring expression context to distinguish functional from silent genes (cephalosporins, aminoglycosides, carbapenems), and those dependent on chromosomal mutations invisible to gene catalogs (quinolones, colistin). The annotation pipeline and pre-trained models are released as a Docker image for clinical application.

# Background

The relationship between antimicrobial resistance (AMR) genotype and phenotype is imperfect. While whole-genome sequencing (WGS) has enabled rapid identification of resistance genes across bacterial species, the presence of an acquired resistance gene does not guarantee phenotypic resistance. This genotype-phenotype discordance — where bacteria carry resistance genes yet remain susceptible, or lack known genes yet are resistant — represents a fundamental gap in our understanding of how resistance is expressed and regulated at the molecular level.

Three biological mechanisms contribute to this discordance. First, resistance genes may be "silent" — present in the genome but not producing functional protein at levels sufficient to confer resistance. Gene silencing can arise from weak promoter activity, inefficient ribosome binding, poor codon adaptation to the host, gene truncation, or insertion sequence disruption (Deekshit & Srikumar, 2022). For aminoglycosides in particular, different modifying enzyme variants have overlapping but distinct substrate specificities, meaning gene detection alone cannot predict which antibiotics are affected (Ramirez & Tolmasky, 2010). Second, for several drug classes — most notably fluoroquinolones — resistance arises primarily through chromosomal point mutations (gyrA, parC) rather than acquired genes (Redgrave et al., 2014), making gene catalog approaches fundamentally blind to the dominant resistance mechanism. Third, the genomic context of a resistance gene modulates its phenotypic impact: location on plasmid versus chromosome, proximity to insertion sequences that may provide outward-facing promoters, and integration within mobile genetic elements all influence gene expression (Partridge et al., 2018).

Current WGS-based AMR tools — ResFinder (Bortolaia et al., 2020), AMRFinderPlus (Feldgarden et al., 2021), CARD RGI (Alcock et al., 2020) — operate as gene catalogs that cannot distinguish functional from silent copies, and machine learning approaches using k-mer or pangenome features (Drouin et al., 2019; Nguyen et al., 2019) achieve high accuracy but cannot explain mechanistically why a prediction was made. Individual genomic context features have been explored in isolation — IS element-driven promoters in *A. baumannii* (Mugnier et al., 2009), codon usage as a horizontal transfer indicator (Lawrence & Ochman, 1997) — but no study has systematically quantified the contribution of gene silencing, chromosomal mutations, and expression context to genotype-phenotype discordance across diverse bacterial species and drug classes.

In this study, we developed a 17-step genomic annotation pipeline that characterizes each detected resistance gene with features spanning transcription (promoter strength), translation (RBS expression, codon adaptation), mobility (plasmid location, IS element proximity, integron context), genomic neighborhood (gene synteny), and chromosomal mutations (PointFinder). By applying this pipeline to 27,856 genomes from five WHO priority pathogens and systematically analyzing genotype-phenotype discordance using gene-variant-specific matching, we identify the molecular mechanisms that explain why resistance genes fail to predict phenotype and demonstrate that these mechanisms can be captured computationally from assembled genome sequences alone.

# Methods

## Data collection

A total of 27,856 bacterial whole-genome sequencing datasets with associated antimicrobial susceptibility testing (AST) phenotypes were retrieved from the NCBI Sequence Read Archive (SRA) and BioSample databases. These represent clinical and surveillance isolates that had been tested against at least one antimicrobial agent, with phenotypes coded as susceptible (S), intermediate (I), or resistant (R). For modeling, intermediate and resistant phenotypes were grouped as non-susceptible, consistent with clinical interpretation where standard-dose treatment is expected to be ineffective for both categories.

Raw paired-end reads (FASTQ format) were downloaded from the SRA and quality-trimmed using BBDuk (BBTools suite) to remove adapter sequences and low-quality bases. Trimmed reads were then assembled *de novo* using SPAdes, producing contiguous assemblies for downstream annotation. Species identification was performed on all assemblies using skani v0.2+ against the Genome Taxonomy Database (GTDB) release 220, based on average nucleotide identity (ANI). The dataset spans 30 bacterial species across 18 genera (Additional file 1: Table S1), dominated by Enterobacterales (*Salmonella*, *E. coli*, *Klebsiella*) from surveillance programs (NARMS, GLASS), and including all three WHO critical-priority pathogens (*A. baumannii*, *P. aeruginosa*, *K. pneumoniae*).

Antimicrobial susceptibility phenotypes were available for 126 distinct antimicrobial agents (Additional file 2: Table S2), with resistance rates ranging from 5.4% (azithromycin) to 50.0% (ceftazidime). The number of isolates tested against each agent varied substantially, reflecting differences in clinical testing panels across surveillance programs and geographic regions.

## Species and antibiotic selection

For each species, an antibiotic was included for analysis if at least 30 isolates were present in the minority class (susceptible or non-susceptible) and the majority-to-minority ratio did not exceed 20:1. Species were then selected if they had at least 10 qualifying antibiotics, sufficient isolates (>50), and clinical relevance as WHO priority pathogens (Additional file 3: Table S3).

Three species with sufficient isolates were excluded: *C. jejuni* (n=3,967; only 3 antibiotics qualified, with quinolone resistance >99% mutation-mediated), *M. tuberculosis* (n=1,282; only 1 antibiotic passed filtering), and *N. gonorrhoeae* (n=765; incomplete PointFinder database coverage). In total, 107 species-antibiotic prediction tasks were defined. Beta-lactams constituted the largest group for Gram-negative species (penicillins, cephalosporins, carbapenems, and inhibitor combinations), while *S. aureus* had a distinct panel of Gram-positive-specific agents (oxacillin, vancomycin, clindamycin, linezolid, erythromycin). The full list of antibiotics per species is provided in Additional file 4: Table S4.

## Genomic annotation pipeline

Each assembled genome was processed through a 17-step annotation pipeline to characterize ARGs and their genomic context. ARGs were identified using AMRFinderPlus (Feldgarden et al., 2021) and chromosomal point mutations using PointFinder via ResFinder (Zankari et al., 2017). Promoter strength was predicted using BPROM (Solovyev & Salamov, 2011) from the 500 bp region upstream of each ARG start codon (strand-aware extraction). Translation initiation rates were predicted using OSTIR (Sadler et al., 2023), which models the thermodynamic free energy of mRNA-16S rRNA interaction at the ribosome binding site. Mobile genetic element context was determined using geNomad (Camargo et al., 2024) for plasmid/prophage contig classification and transposase detection, and IntegronFinder (Cury et al., 2016) for integron identification. Small regulatory RNAs within 5 kb of ARGs were detected using Infernal/Rfam (Nawrocki & Eddy, 2013; Kalvari et al., 2021). Gene prediction was performed with Prodigal (Hyatt et al., 2010), and multi-locus sequence typing with mlst (Seemann T).

Additional per-ARG features were computed: codon adaptation index and rare codon frequency relative to the host genome (Sharp & Li, 1987), GC content deviation from the genome average as a horizontal transfer indicator (Lawrence & Ochman, 1997), gene copy number, operon size and position, and gene synteny (classification of the five flanking genes on each side as AMR, transposase, stress, virulence, hypothetical, or other). The complete pipeline output per ARG comprises 46 annotation columns, enabling downstream aggregation into feature matrices.

## Feature matrix construction and machine learning

For each isolate-antibiotic pair, per-ARG annotation features were aggregated into a fixed-width numeric vector. Each antibiotic was mapped to its corresponding AMRFinderPlus drug class (e.g., ciprofloxacin to QUINOLONE), and "target" features were computed only from ARGs matching that drug class, while "genome-wide" features (total gene copies, total neighboring transposase count) were computed across all ARGs in the isolate. For numeric features, max and mean aggregations were used (with min also included for promoter LDF, RBS expression, and CAI). When no target ARG was detected, target-level features were set to NaN and imputed with per-fold median during training. Features were organized into 12 biologically motivated groups (Table 1) spanning ARG presence, expression (promoter, RBS, codon usage), mobility (plasmid, IS elements, integrons), regulation (sRNA, operon, synteny), genomic context (GC deviation, gene copies), point mutations (PointFinder), and population structure (MLST). MLST sequence types were encoded as binary presence/absence columns for the 30 most frequent STs per species (e.g., a ST131 isolate receives 1 in the ST131 column and 0 in all others). The final feature matrix contains one row per isolate-antibiotic pair with ~116 features per row.

**Table 1. Feature groups and their components.**

| Group | n | Features | Biological rationale |
|-------|---|----------|---------------------|
| ARG presence | 3 | has_target_arg, n_target_args, n_amr_genes | Does the isolate carry a resistance gene for this drug class? |
| Promoter | 9 | promoter_ldf, tf_sites, distance, AT ratio (max/mean/min) | Transcription initiation strength |
| RBS | 7 | rbs_expression, dG_total, dG_mRNA (max/mean/min) | Translation initiation efficiency |
| Codon usage | 7 | CAI, rare_codon_pct, rare_codon_clusters (max/mean/min) | Translational elongation efficiency |
| Mobility | 12 | IS_distance, IS_strand, on_plasmid, on_prophage, n_plasmid/prophage/chromosome, genome_transposase_count | Mobile element context |
| Regulation | 8 | sRNA_proximity, operon_size, operon_position (max/mean) | Post-transcriptional regulation |
| Integron | 3 | in_integron, integron_distance (max/mean) | Gene cassette context |
| Synteny | 10 | n_transposase, n_AMR, n_stress, n_hypothetical, n_virulence (max/mean) | Genomic neighborhood |
| Genomic context | 9 | gene_gc, genome_gc, gc_deviation, gene_copies (max/mean) | Horizontal transfer indicators |
| Point mutations | 8 | pf_gyrA/gyrB/parC/parE, pf_ampC, pf_n_quinolone_mutations | Chromosomal resistance mutations |
| Functionality | 5 | arg_low_rbs, arg_truncated, arg_poor_adaptation, arg_functional_score | Non-functional ARG indicators |
| Population | 4+ | identity, coverage (max/mean), MLST ST dummies | Population structure proxies |

Four classification algorithms were evaluated: Logistic Regression (LR; L2 penalty, balanced class weights) as a linear baseline, Random Forest (RF; 300 trees, balanced class weights) as a non-linear ensemble robust to feature noise, XGBoost (XGB; 100 trees, max_depth=6, scale_pos_weight for class imbalance) as a gradient-boosted approach, and Multilayer Perceptron (MLP; two hidden layers of 64 and 32 neurons, early stopping) to test non-linear feature interactions.

All algorithms used a fixed random state (seed=42) for reproducibility. Feature scaling (zero mean, unit variance) was applied within each cross-validation fold to prevent data leakage. All models were evaluated using 5-fold stratified cross-validation, with models trained independently per species and per antibiotic. F1 score was used as the primary metric, as it balances sensitivity and precision. Balanced accuracy, MCC, ROC-AUC, PR-AUC, and Brier score were also computed.

To quantify the contribution of each annotation layer, a cumulative feature ablation study was performed using Random Forest. Five configurations were evaluated, each adding feature groups incrementally: (1) ARG presence only (3 features), (2) plus expression features including promoter, RBS, and codon usage (26 features), (3) plus mobility features including IS elements, plasmid/prophage, and integrons (38 features), (4) plus regulation features including sRNA, operon, and synteny (56 features), and (5) the full model with all features including genomic context, point mutations, and population structure (~116 features). Learning curve analyses were performed by training Random Forest on stratified subsamples of increasing size (n = 50 to 10,000) to determine minimum sample size requirements. Feature importance was assessed using SHAP (SHapley Additive exPlanations; Lundberg & Lee, 2017) TreeExplainer on Random Forest, with the mean absolute SHAP value per feature computed per antibiotic from the non-susceptible class.

## Discordance analysis

To validate the biological relevance of genomic context features, a systematic analysis of phenotype-genotype discordance was performed. For each species-antibiotic pair, isolates were categorized into four groups:

Isolates were classified as concordant (ARG present and resistant, or ARG absent and susceptible) or discordant: either the gene is present but the isolate is susceptible (silent gene), or the isolate is resistant without a detected gene (missing mechanism). ARG presence was assessed using gene-variant-specific matching rather than broad drug-class matching — each gene was mapped to the specific antibiotics it is known to affect (e.g., TEM-1 confers ampicillin but not cefepime resistance; aac(3)-IId confers gentamicin but not amikacin resistance; chromosomal blaEC was excluded as it does not confer resistance at basal expression). For discordant cases, expression-level features were compared between concordant and discordant groups to identify mechanisms of gene silencing, and PointFinder mutation data was examined to determine whether chromosomal point mutations explain resistance in the absence of acquired ARGs.

## Software and reproducibility

All analyses were implemented in Python 3.12 using scikit-learn 1.4, XGBoost 2.0, and SHAP 0.45. The genomic annotation pipeline was implemented as a Bash/Python script (`annotate_one.sh`) using GNU Parallel for batch processing. Feature matrix construction was performed with `build_feature_matrix_v2.py`. Model training, evaluation, and SHAP analysis were performed with `train_models.py`.

All scripts, feature matrices, and model outputs are available at https://github.com/tatsu1207/AMR.

# Results

## Dataset composition: source and geographic distribution

Sample metadata retrieved from NCBI BioSample revealed distinct source compositions across the five species, reflecting the surveillance programs and clinical settings from which isolates originated (Table 2).

**Table 2. Source distribution by species.**

| Source | Salmonella (n=8,937) | E. coli (n=5,236) | Klebsiella (n=1,337) | S. aureus (n=911) | A. baumannii (n=720) |
|--------|---------------------|-------------------|---------------------|-------------------|---------------------|
| Human clinical | 1.2% | 21.6% | 81.3% | 88.9% | 79.0% |
| Poultry | 34.0% | 33.2% | -- | -- | -- |
| Turkey | 18.0% | 17.8% | -- | -- | -- |
| Swine | 7.3% | 11.1% | 6.1% | -- | -- |
| Cattle | 3.9% | 11.0% | -- | -- | -- |
| Unknown | 34.6% | 4.5% | 5.8% | 10.9% | 5.4% |

*Salmonella* and *E. coli* isolates originated predominantly from food-animal surveillance programs (poultry, turkey, swine, cattle), consistent with the NARMS (National Antimicrobial Resistance Monitoring System) origin of many of these isolates. In contrast, *Klebsiella*, *S. aureus*, and *A. baumannii* were overwhelmingly from human clinical specimens (79-89%). Geographically, 68.2% of all isolates were from the USA, with Canada (9.2%) and Australia (2.4%) as the next most common origins.

### Resistance rates vary by source

The host source of isolation was strongly associated with resistance prevalence, particularly for *E. coli* (Table 3).

**Table 3. Resistance rates by source for selected antibiotics.**

| Species | Antibiotic | Human clinical | Poultry | Turkey | Swine |
|---------|-----------|---------------|---------|--------|-------|
| *E. coli* | Ciprofloxacin | 61% (n=984) | 5% (n=1,706) | 1% (n=910) | -- |
| *E. coli* | Ceftriaxone | 84% (n=627) | 26% (n=1,711) | 1% (n=914) | -- |
| *E. coli* | Ampicillin | 83% (n=960) | 38% (n=1,711) | 45% (n=914) | -- |
| *Salmonella* | Nalidixic acid | -- | 18% (n=3,037) | 5% (n=1,608) | 2% (n=652) |
| *Salmonella* | Ampicillin | -- | 16% (n=3,038) | 29% (n=1,609) | 16% (n=652) |
| *Salmonella* | Tetracycline | -- | 52% (n=3,038) | 47% (n=1,608) | 55% (n=652) |
| *S. aureus* | Oxacillin (MRSA) | 14% (n=609) | -- | -- | -- |
| *S. aureus* | Erythromycin | 38% (n=602) | -- | -- | -- |

Human clinical *E. coli* had dramatically higher fluoroquinolone (61% vs 1-5%) and cephalosporin (84% vs 1-26%) resistance than food-animal isolates, reflecting the selective pressure from clinical antibiotic use. Conversely, *Salmonella* tetracycline resistance was similar across sources (47-55%), consistent with widespread tetracycline use in both human medicine and animal production.

### ARG profiles differ by source

ARG distribution also varied by source, with clinical and food-animal isolates carrying distinct resistance gene repertoires (Table 4).

**Table 4. ARG prevalence by source (selected genes, *E. coli*).**

| Gene | Function | Poultry | Human clinical | Turkey | Max Diff |
|------|----------|---------|---------------|--------|----------|
| blaCTX-M-15 | ESBL (3rd/4th gen cephalosporins) | 0.0% | 9.4% | 0.0% | 9.4% |
| blaCMY-2 | AmpC (cephamycins) | 23.0% | 7.0% | 1.2% | 22.9% |
| aac(6')-Ib-cr5 | Aminoglycoside + ciprofloxacin | 0.0% | 12.6% | 0.0% | 12.6% |
| aadA5 | Streptomycin/spectinomycin | 0.0% | 23.0% | 0.0% | 23.0% |
| blaHER-3 | Narrow-spectrum beta-lactamase | 0.3% | 0.0% | 14.7% | 14.7% |

CTX-M-15 (the globally dominant ESBL) was found almost exclusively in human clinical isolates (9.4% vs 0% in food animals), while CMY-2 (plasmid-borne AmpC) was concentrated in poultry (23.0% vs 7.0% clinical). This source-dependent ARG distribution has implications for model generalizability: a model trained predominantly on food-animal *Salmonella* may underperform on human clinical isolates, and vice versa.

### Population structure is strongly associated with source

MLST analysis revealed that different sources harbor phylogenetically distinct lineages, with minimal overlap between food-animal and human clinical populations (Table 5).

**Table 5. Dominant MLST sequence types by source.**

| Species | Human clinical | Poultry | Turkey | Swine |
|---------|---------------|---------|--------|-------|
| *E. coli* | ST131 (383), ST1193 (74), ST69 (52), ST95 (50) | ST10 (172), ST117 (95), ST115 (92) | ST3580 (65), ST58 (64), ST10 (52) | -- |
| *Salmonella* | -- | ST152/Kentucky (616), ST32/Infantis (384) | ST412/Reading (124), ST33 (57) | ST471 (6), ST34 (5) |

For *E. coli*, the separation was near-complete: ST131 — the pandemic fluoroquinolone-resistant, ESBL-producing clone — comprised 34% of human clinical isolates but was absent from poultry and turkey. Conversely, ST117 and ST115 (avian-adapted pathogenic lineages) were common in poultry but rare in clinical specimens. Only ST10, a ubiquitous generalist lineage, was shared across sources.

For *Salmonella*, source-serovar associations were strong: ST152 (serovar Kentucky) dominated poultry (20%), ST412 (serovar Reading) dominated turkey (8%), and ST19 (serovar Typhimurium) was the most common in the "Unknown" category, likely representing human surveillance isolates.

This phylogenetic stratification by source was further confounded with study origin and geography. Canada contributed 1,132 *E. coli* isolates that were 99.6% poultry from a single surveillance program. A cluster of 992 *E. coli* ST152 isolates with unknown metadata (non-SRR sample identifiers) likely originated from a single genomic epidemiology study. For *Salmonella*, 97% of isolates were from USA-based NARMS surveillance.

These confoundings have three implications for AMR prediction:

1. Lineage-resistance confounding: The model may learn lineage-specific associations (e.g., ST131 = fluoroquinolone resistant) rather than mechanistic features. While this improves prediction within the training population, it may reduce generalizability to novel lineages not represented in training data. MLST features were intentionally included to capture these associations explicitly.

2. Source bias: A model trained predominantly on food-animal *Salmonella* (70% of our dataset) may underperform on human clinical isolates, which carry different resistance gene repertoires and lineages. Conversely, *E. coli* models may be biased toward the high-resistance ST131 clone that dominates the clinical fraction.

3. Study batch effects: Large single-study clusters (e.g., 992 ST152 *E. coli*) may introduce batch-specific patterns that the model learns as generalizable features. Cross-validation within such clusters overestimates performance because training and test folds share the same study population.

These limitations are inherent to publicly available genomic surveillance data and should be considered when interpreting prediction performance and deploying the model across different contexts. Future work should evaluate model performance using geographically and temporally independent validation sets.

## ARG presence alone fails to predict phenotype

We first examined whether resistance gene presence alone could predict antimicrobial resistance phenotype. To avoid inflating discordance with gene-variant mismatches, we used gene-variant-specific matching: each detected gene was evaluated against a curated mapping of beta-lactamase and aminoglycoside substrate spectra (e.g., TEM-1 confers ampicillin but not cefepime resistance; aac(3)-IId confers gentamicin but not amikacin resistance; chromosomal blaEC was excluded as it does not confer resistance at basal expression). For other drug classes, standard drug-class matching was used.

With gene-variant-specific matching, isolate-antibiotic observations were classified into four categories: (1) ARG-explained resistance (gene-variant-matched ARG present in resistant isolate), (2) mutation-explained resistance (no gene-variant-matched ARG but PointFinder chromosomal mutation detected), (3) silent genes (gene-variant-matched ARG present but isolate susceptible), and (4) truly unexplained resistance (no ARG, no detected mutation, but resistant). Categories 3 and 4 constitute genuine genotype-phenotype discordance (Table 6).

**Table 6. Genotype-phenotype discordance by species (gene-variant-specific matching + PointFinder mutations).**

| Species | Observations | ARG-explained R | Mutation-explained R | Silent genes (ARG+S) | Unexplained R | True disc. % |
|---------|-------------|----------------|---------------------|----------------------|--------------|-------------|
| *S. enterica* | 117,327 | 21,240 | 1,820 | 1,744 (1.5%) | 2,111 (1.8%) | 3.3% |
| *E. coli* | 83,172 | 18,248 | 1,483 | 2,451 (2.9%) | 3,800 (4.6%) | 7.5% |
| *S. aureus* | 7,052 | 1,279 | 217 | 1,354 (19.2%) | 259 (3.7%) | 22.9% |
| *K. pneumoniae* | 15,640 | 7,875 | 19 | 3,745 (23.9%) | 1,193 (7.6%) | 31.6% |
| *A. baumannii*\* | 8,531 | 3,321 | 0 | 1,063 (12.5%) | 2,584 (30.3%) | 42.7% |

*\*No PointFinder database exists for A. baumannii; quinolone resistance (87% of ciprofloxacin-resistant isolates) is classified as "unexplained" despite being likely gyrA-mediated.*

Two levels of correction were critical for accurate discordance estimation:

Gene-variant-specific matching eliminated false discordance from drug-class-level misattribution. For example, *E. coli* discordance dropped from 33.3% (drug-class matching) to 9.3% (gene-variant-specific) because narrow-spectrum penicillinases (TEM-1) and intrinsic chromosomal AmpC (blaEC) were no longer counted as resistance determinants for cephalosporins and carbapenems.

PointFinder mutation accounting further reduced discordance by reclassifying quinolone-resistant isolates lacking ARGs. In *Salmonella*, 1,820 quinolone-resistant isolates were explained by gyrA/parC mutations, reducing discordance from 4.8% to 3.3%. In *E. coli*, 1,483 were mutation-explained, reducing from 9.3% to 7.5%.

After both corrections, the remaining true discordance comprised two mechanistically distinct types: silent genes (gene-variant-matched ARG present but not conferring resistance) and truly unexplained resistance (no genetic determinant detected). Silent genes dominated in *Klebsiella* (23.9%) and *S. aureus* (19.2%), while unexplained resistance was highest in *A. baumannii* (30.3%, inflated by missing PointFinder coverage; see Table 15 footnote) and *E. coli* (4.6%). These two types of discordance required different feature layers to resolve: expression features addressed silent genes, while PointFinder mutations addressed unexplained quinolone resistance.

### Discordance varies by drug class and mechanism

Discordance rates varied by drug class, driven by the dominant resistance mechanism (Table 7).

**Table 7. High-discordance antibiotics by species (gene-variant-specific matching, >15% discordance, n>200).**

| Species | Antibiotic | Disc. % | Silent | Missing | Mechanism |
|---------|-----------|---------|--------|---------|-----------|
| *S. aureus* | Tetracycline | 87.7% | 670 | 1 | tet(K)/tet(M) expression variability |
| *A. baumannii* | Ciprofloxacin | 87.4% | 0 | 659 | 100% gyrA point mutations |
| *A. baumannii* | Levofloxacin | 82.6% | 0 | 589 | 100% gyrA point mutations |
| *E. coli* | Cefazolin | 82.5% | 4 | 433 | Unknown mechanism (no ESBL detected) |
| *A. baumannii* | TMP-SMX | 80.8% | 1 | 590 | Unknown (no ARG detected) |
| *Klebsiella* | Colistin | 78.1% | 0 | 368 | mgrB inactivation |
| *Salmonella* | Colistin | 75.3% | 0 | 980 | mcr-independent |
| *A. baumannii* | Ampicillin-sulbactam | 73.3% | 326 | 182 | OXA expression + unknown |
| *S. aureus* | Oxacillin | 62.6% | 448 | 3 | mecA expression/SCCmec variability |
| *Klebsiella* | Ceftazidime-avibactam | 60.5% | 127 | 17 | Novel resistance mechanisms |
| *Klebsiella* | Cefepime | 44.7% | 304 | 1 | ESBL substrate specificity |
| *A. baumannii* | Imipenem | 45.3% | 313 | 8 | OXA expression variability |
| *E. coli* | Levofloxacin | 41.8% | 28 | 410 | gyrA/parC mutations |
| *Klebsiella* | Piperacillin-tazobactam | 42.5% | 338 | 118 | Inhibitor resistance complexity |
| *Klebsiella* | Ertapenem | 40.4% | 270 | 6 | Carbapenemase expression variability |
| *E. coli* | Amikacin | 33.9% | 161 | 248 | AAC substrate specificity + unknown |
| *Klebsiella* | Meropenem | 37.9% | 304 | 6 | Carbapenemase expression variability |

Three distinct discordance patterns emerged:

Missing mechanism (no gene-variant-matched ARG): Dominant for quinolones (ciprofloxacin, levofloxacin), colistin, and glycopeptides. Resistance is mediated by chromosomal point mutations (gyrA/parC), regulatory changes (mgrB inactivation), or structural modifications (cell wall thickening) that are not captured by ARG databases.

True silent genes (gene-variant-matched ARG present but not conferring resistance): Dominant for *S. aureus* oxacillin (mecA present but SCCmec type affects expression), *Klebsiella* carbapenem/cephalosporin resistance (carbapenemase/ESBL present but expression levels vary), and *S. aureus* tetracycline (tet(K)/tet(M) present but variably expressed).

Mixed: Some antibiotics show both patterns, such as *E. coli* amikacin (161 silent + 248 missing) and *A. baumannii* ampicillin-sulbactam (326 silent + 182 missing).

Quinolones had the most missing mechanisms (4,861 cases). Across all species, 69% of quinolone-resistant isolates had no acquired quinolone resistance gene — resistance was mediated by chromosomal gyrA/parC point mutations invisible to gene catalog approaches.

Colistin showed the highest overall discordance (68.6%), with 1,437 resistant isolates lacking mcr genes. Colistin resistance in *Klebsiella* and *Salmonella* is predominantly mediated by mgrB inactivation and other chromosomal modifications not captured by ARG databases.

### High-discordance antibiotics

After gene-variant-specific and mutation corrections, the highest true discordance concentrated in specific species-antibiotic pairs (Table 8).

**Table 8. Antibiotics with highest true discordance (gene-variant-specific matching + PointFinder, >25%, n>200).**

| Species | Antibiotic | True Disc. % | True Silent | Unexplained R | Root cause |
|---------|-----------|-------------|-------------|--------------|-----------|
| *S. aureus* | Tetracycline | 87.7% | 670 | 1 | tet(K)/tet(M) expression variability |
| *E. coli* | Cefazolin | 82.5% | 4 | 433 | Missing ESBL not in database |
| *A. baumannii* | TMP-SMX | 80.8% | 1 | 590 | Unknown mechanism (no PointFinder DB) |
| *A. baumannii* | Ampicillin-sulbactam | 73.3% | 326 | 182 | OXA expression variability |
| *S. aureus* | Oxacillin | 62.6% | 448 | 3 | mecA expression / SCCmec type |
| *Klebsiella* | Ceftazidime-avibactam | 60.5% | 127 | 17 | Novel resistance mechanism |
| *Klebsiella* | Cefepime | 44.7% | 304 | 1 | ESBL present but wrong spectrum |
| *A. baumannii* | Imipenem | 45.3% | 313 | 8 | OXA carbapenemase expression |
| *Klebsiella* | Piperacillin-tazobactam | 42.5% | 338 | 118 | Inhibitor resistance complexity |
| *Klebsiella* | Ertapenem | 40.4% | 270 | 6 | Carbapenemase expression variability |
| *Klebsiella* | Meropenem | 37.9% | 304 | 6 | Carbapenemase expression variability |
| *E. coli* | Amikacin | 33.9% | 161 | 248 | AAC substrate specificity |
| *S. aureus* | Levofloxacin | 32.8% | 0 | 145 | gyrA mutations (explained by PointFinder) |
| *Klebsiella* | Ceftazidime | 32.2% | 382 | 41 | ESBL present but wrong spectrum |

Two distinct patterns of true discordance remained:

True silent genes (gene-variant-matched ARG present but not conferring resistance): Dominated in *S. aureus* oxacillin (448 cases — mecA present but SCCmec type or regulatory context affects expression), *Klebsiella* cephalosporins/carbapenems (270-382 cases — ESBL or carbapenemase present but expression insufficient), and *S. aureus* tetracycline (670 cases — tet(K)/tet(M) variably expressed). These are the cases where expression features provide the most value.

Truly unexplained resistance: Highest in *A. baumannii* (TMP-SMX: 590, ampicillin-sulbactam: 182) due to missing PointFinder database, and *E. coli* cefazolin (433) likely representing novel or uncharacterized resistance mechanisms.

## Expression features resolve silent gene discordance

Having identified silent genes as the dominant form of discordance, we examined whether genomic context features could distinguish functional from non-functional ARG copies.

### Expression differences between concordant and discordant cases

For antibiotics with true silent genes (gene-variant-matched ARG present but isolate susceptible), we compared expression-level features between concordant (ARG+ R) and discordant (ARG+ S) isolates. Only substrate-matched ARGs were included — e.g., for gentamicin, only isolates carrying aac(3)/ant(2'') family genes (which confer gentamicin resistance) were compared, excluding isolates carrying aph(3') genes (which confer kanamycin but not gentamicin resistance) (Table 9).

**Table 9. Expression features in concordant vs discordant cases with gene-variant-specific matching (*S. enterica*).**

| Antibiotic | Feature | ARG+ R | ARG+ S | Difference | Interpretation |
|-----------|---------|--------|--------|------------|----------------|
| Chloramphenicol | RBS expression | 319.7 | 46.7 | -85.4% | Near-silent translation |
| Gentamicin | RBS expression | 1381.7 | 873.7 | -36.8% | Weaker translation initiation |
| Ampicillin | RBS expression | 641.2 | 445.8 | -30.5% | Reduced translation |
| Ceftriaxone | RBS expression | 660.5 | 463.2 | -29.9% | Reduced translation |
| Ceftriaxone | Gene coverage | 99.7% | 85.5% | -14.3% | Truncated / fragmented gene |
| Ampicillin | Gene coverage | 99.6% | 94.5% | -5.1% | Partial gene copy |
| Tetracycline | Gene coverage | 99.9% | 95.6% | -4.3% | Partial gene copy |
| Ampicillin | Promoter LDF | 4.05 | 3.37 | -16.9% | Weaker transcription |
| Tetracycline | Promoter LDF | 4.52 | 3.83 | -15.4% | Weaker transcription |
| Chloramphenicol | Promoter LDF | 3.45 | 3.11 | -9.7% | Weaker transcription |
| Gentamicin | Promoter LDF | 4.60 | 4.21 | -8.6% | Weaker transcription |

Gene-variant-specific matching ensured that these comparisons reflect genuine differences in gene functionality rather than drug-class-level misattribution. The differences were consistent across drug classes: silent ARGs had lower predicted translation initiation rates (11-85% reduction, reflecting poor ribosome binding), lower gene coverage (indicating truncation or fragmentation), and weaker promoter activity (reflecting reduced transcription). The most dramatic difference was chloramphenicol, where silent copies had 85% lower predicted translation initiation — indicating near-complete translational silencing of the resistance gene.

These differences correspond to known molecular mechanisms of gene silencing and provide biological validation that the expression features capture real functional variation, not statistical noise.

### Ablation confirms expression features as the key solution

The feature ablation study showed that expression features (promoter, RBS, codon adaptation) provided the largest single improvement over ARG-only prediction for drug classes dominated by silent gene discordance (Table 10).

**Table 10. Ablation gains by feature layer and drug class.**

| Drug Class | ARG only | + Expression | + Mobility | + Regulation | Full model | Expression gain |
|-----------|----------|-------------|-----------|-------------|-----------|----------------|
| Beta-lactam | 0.727 | 0.829 | 0.833 | 0.835 | 0.892 | +0.102 |
| Aminoglycoside | 0.740 | 0.852 | 0.858 | 0.860 | 0.900 | +0.112 |
| Quinolone | 0.608 | 0.614 | 0.614 | 0.613 | 0.945 | +0.006 (expression irrelevant) |
| Tetracycline | 0.920 | 0.924 | 0.924 | 0.925 | 0.930 | +0.004 (ARG sufficient) |
| Sulfonamide | 0.968 | 0.966 | 0.967 | 0.967 | 0.967 | -0.002 (ARG sufficient) |

For beta-lactams and aminoglycosides — the drug classes with 27-36% silent gene discordance — expression features improved mean F1 by +0.10 to +0.11. For quinolones, expression features were irrelevant (+0.006) because the discordance is missing-mechanism type, not silent-gene type. For tetracyclines and sulfonamides, ARG presence was already sufficient (discordance <8%).

### Species-level ablation patterns

The overall ablation gain varied by species, driven by the drug class composition of each species' antibiotic panel (Table 11).

**Table 11. Mean F1 by feature configuration across species (RF, 5-fold CV).**

| Configuration | Salmonella | E. coli | Klebsiella | S. aureus | A. baumannii |
|--------------|-----------|---------|------------|-----------|-------------|
| 1. ARG only | 0.756 | 0.734 | 0.774 | 0.691 | 0.870 |
| 2. + Expression | 0.861 | 0.804 | 0.849 | 0.734 | 0.919 |
| 3. + Mobility | 0.868 | 0.807 | 0.850 | 0.736 | 0.919 |
| 4. + Regulation | 0.869 | 0.810 | 0.851 | 0.736 | 0.920 |
| 5. Full model | 0.924 | 0.857 | 0.873 | 0.804 | 0.932 |
| Expression gain | +0.105 | +0.070 | +0.074 | +0.043 | +0.049 |
| Full model gain | +0.168 | +0.123 | +0.099 | +0.113 | +0.062 |

*Salmonella* showed the largest expression gain (+0.103) because its panel includes aminoglycosides and cephalosporins with high silent gene rates. *A. baumannii* showed the smallest full model gain (+0.062) because its high baseline resistance prevalence and strong ARG-phenotype concordance already provide high ARG-only prediction (F1=0.870).

## PointFinder mutations resolve missing mechanism discordance

### Quinolone resistance is almost entirely mutation-mediated

Across all five species, 69% of quinolone-resistant isolates had no acquired resistance gene — resistance was mediated by chromosomal gyrA/parC/parE mutations (Table 12).

**Table 12. Quinolone resistance mechanisms by species.**

| Species | Antibiotic | Resistant | ARG+ | ARG- (mutation) | ARG- rate | Best F1 |
|---------|-----------|-----------|------|-----------------|-----------|---------|
| *S. aureus* | Ciprofloxacin | 76 | 0 | 76 | 100% | 0.961 |
| *S. aureus* | Levofloxacin | 145 | 0 | 145 | 100% | 0.983 |
| *A. baumannii* | Ciprofloxacin | 659 | 1 | 658 | 99.8% | 0.972 |
| *A. baumannii* | Levofloxacin | 485 | 1 | 484 | 99.8% | 0.958 |
| *S. enterica* | Nalidixic acid | 796 | 113 | 683 | 85.8% | 0.973 |
| *S. enterica* | Ciprofloxacin | 712 | 160 | 552 | 77.5% | 0.960 |
| *E. coli* | Nalidixic acid | 547 | 109 | 438 | 80.1% | 0.926 |
| *E. coli* | Levofloxacin | 570 | 160 | 410 | 71.9% | 0.942 |
| *E. coli* | Ciprofloxacin | 771 | 204 | 567 | 73.5% | 0.921 |
| *K. pneumoniae* | Ciprofloxacin | 601 | 354 | 247 | 41.1% | 0.955 |
| *K. pneumoniae* | Levofloxacin | 334 | 113 | 221 | 66.2% | 0.925 |

Despite the near-complete absence of acquired ARGs, all quinolone predictions achieved F1 > 0.92 with the full model. The ablation showed that this was entirely attributable to PointFinder features — the jump from configuration 4 (no mutations) to configuration 5 (full model including PointFinder) was +0.07 to +0.48 for quinolones.

### PointFinder mutation landscape

PointFinder identified chromosomal mutations in the majority of isolates across all species:

| Species | parC mutations | gyrA mutations | parE mutations | Samples |
|---------|---------------|---------------|---------------|---------|
| *S. enterica* | 6,565 | 806 | -- | 9,111 |
| *E. coli* | 660 | 885 | 411 | 5,236 |
| *S. aureus* | -- | detectable via pf_has_gyrA | -- | 977 |
| *A. baumannii* | detectable | detectable | -- | 757 |
| *K. pneumoniae* | detectable | detectable | detectable | 1,808 |

SHAP analysis confirmed that PointFinder features were among the most important globally: `pf_n_quinolone_mutations` ranked #1 in *Klebsiella* and *S. aureus*, and `pf_has_gyrA` ranked in the top 5 for *E. coli*.

## Full model performance

With both expression features (for silent genes) and PointFinder mutations (for missing mechanisms) integrated, the full model achieved mean best F1 scores (the highest F1 among the four tested algorithms, averaged across antibiotics) of 0.828 to 0.941 across the five species (Table 13).

**Table 13. Full model performance summary.**

| Species | Antibiotics | Mean Best F1 | F1 >= 0.95 | F1 0.80-0.95 | F1 < 0.80 | Best Algorithm |
|---------|-------------|-------------|-----------|-------------|----------|----------------|
| *A. baumannii* | 13 | 0.941 | 6 | 7 | 0 | RF(6), LR(3), XGB(3), MLP(1) |
| *S. enterica* | 20 | 0.930 | 10 | 8 | 2 | RF(7), XGB(6), MLP(5), LR(2) |
| *K. pneumoniae* | 26 | 0.877 | 4 | 18 | 4 | RF(17), XGB(6), LR(3) |
| *E. coli* | 35 | 0.871 | 12 | 18 | 5 | RF(18), XGB(9), LR(6), MLP(2) |
| *S. aureus* | 13 | 0.828 | 3 | 5 | 5 | XGB(6), RF(6), LR(1) |

Performance differences across species were explained by three factors:

Resistance gene diversity: *A. baumannii* (23 unique ARGs per 100 isolates) achieved the best prediction despite highest diversity, largely due to high baseline resistance prevalence for most antibiotics and strong ARG-phenotype concordance (ARG-only F1=0.870). *Salmonella* (1.1 per 100) had simplest gene-phenotype relationships.

Antibiotic panel breadth: *E. coli* (35 antibiotics) and *Klebsiella* (26) included difficult targets (carbapenems, novel inhibitor combinations, tigecycline) that lowered their mean F1, while *A. baumannii* (13) and *S. aureus* (13) had narrower panels.

Mechanism coverage: *S. aureus* had the lowest F1 (0.828) partly because it includes antibiotics with mechanisms not covered by current databases — vancomycin intermediate resistance (cell wall thickening, F1=0.755), linezolid (23S rRNA mutations, F1=0.570).

## What remains unexplained

Despite the improvements from the full model (Table 14), eight species-antibiotic pairs achieved F1 < 0.70, representing the current limits of WGS-based prediction (Table 14).

**Table 14. Poorly predicted antibiotics and root causes (F1 < 0.70).**

| Cause | Species | Antibiotic | F1 | n_R | Specific mechanism |
|-------|---------|-----------|-----|-----|-------------------|
| No database coverage | *E. coli* | Nitrofurantoin | 0.217 | 42 | nfsA/nfsB inactivation; no PointFinder DB |
| | *K. pneumoniae* | Tigecycline | 0.564 | 121 | ramA/ramR efflux regulation; not in databases |
| Complex mechanism | *E. coli* | Piperacillin-tazo | 0.522 | 127 | Requires specific inhibitor-resistant beta-lactamase variants |
| | *E. coli* | Cefepime | 0.622 | 275 | Multiple ESBL variants with different spectra |
| Small sample | *S. aureus* | Linezolid | 0.570 | 31 | Rare; 23S rRNA mutations + cfr |
| | *E. coli* | Doripenem | 0.591 | 31 | Too few resistant isolates for 5-fold CV |
| | *S. aureus* | Doxycycline | 0.630 | 38 | Insufficient training data |
| | *E. coli* | Imipenem | 0.669 | 39 | Small R class + diverse mechanisms |

These failures cluster into three addressable categories: (1) expanding mutation databases to cover regulatory regions (ramA/ramR, nfsA/nfsB), (2) improving resolution of enzyme functional classification within drug classes, and (3) collecting larger datasets for rare resistance phenotypes.

## SHAP feature importance reflects mechanism types

SHAP analysis confirmed that the most important features per species directly reflect the dominant resistance mechanism type in that species' antibiotic panel (Table 15). Per-antibiotic SHAP values are provided in Additional file 5: Table S5.

**Table 15. Top 5 global SHAP features per species.**

| Rank | Salmonella | E. coli | Klebsiella | S. aureus | A. baumannii |
|------|-----------|---------|------------|-----------|-------------|
| 1 | has_target_arg | n_target_args | pf_n_quinolone_mut | pf_n_quinolone_mut | n_amr_genes |
| 2 | n_target_args | target_gc_dev | pf_has_quinolone_mut | pf_n_mutations | genome_gene_copies |
| 3 | genome_gc | has_target_arg | n_amr_genes | genome_gene_copies | genome_gc |
| 4 | target_n_plasmid | target_n_plasmid | genome_gene_copies | pf_has_gyrA | n_target_args |
| 5 | target_any_on_plasmid | arg_functional_score | pf_has_parC | n_amr_genes | target_gc_dev |

### Species-level interpretation

- *E. coli*: Expression features (GC deviation, CAI, gene GC) dominated globally because 17 of 35 prediction tasks involved beta-lactams, where distinguishing functional from silent beta-lactamases is the key challenge.
- *Klebsiella* and *S. aureus*: PointFinder features dominated because quinolone resistance prediction — where mutation detection is decisive — constitutes a large fraction of the overall SHAP signal.
- *A. baumannii* and *Salmonella*: Population features (genome_gc, n_amr_genes) dominated, reflecting the high baseline resistance prevalence where total AMR gene burden correlates with resistance status.
- Universal: Gene count features (n_target_args, n_amr_genes) ranked in the top 5 for all species, confirming that gene dosage contributes to phenotype beyond simple presence/absence.

### Per-antibiotic SHAP reveals drug-class-specific drivers

Examination of per-antibiotic SHAP values revealed that the top predictive feature was determined by the resistance mechanism, not the species (Table 16).

**Table 16. Top SHAP feature by antibiotic class (selected examples).**

| Drug class | Antibiotic | Species | Top SHAP feature | Value | Mechanism |
|-----------|-----------|---------|------------------|-------|-----------|
| Quinolone | Nalidixic acid | *Salmonella* | pf_has_gyrA | 0.164 | gyrA point mutations |
| Quinolone | Nalidixic acid | *E. coli* | pf_has_gyrA | 0.163 | Same mechanism, same feature |
| Quinolone | Ciprofloxacin | *Klebsiella* | pf_has_quinolone_mutation | 0.082 | Same mechanism |
| Quinolone | Levofloxacin | *S. aureus* | pf_n_quinolone_mutations | 0.133 | Same mechanism |
| Beta-lactam | Ceftriaxone | *E. coli* | target_gc_deviation_max | 0.084 | ESBL GC signature |
| Beta-lactam | Ceftiofur | *E. coli* | target_cai_max | 0.079 | ESBL codon adaptation |
| Beta-lactam | Ceftiofur | *Salmonella* | target_gene_gc_max | 0.054 | Same — foreign gene GC |
| Aminoglycoside | Gentamicin | *E. coli* | target_gene_gc_max | 0.075 | AME gene GC signature |
| Aminoglycoside | Gentamicin | *Klebsiella* | target_gc_deviation_max | 0.074 | Same mechanism |
| Beta-lactam | Oxacillin | *S. aureus* | target_rbs_dg_total_max | 0.066 | mecA RBS strength |
| Beta-lactam | Doripenem | *A. baumannii* | target_rbs_dg_total_max | 0.044 | OXA expression level |
| Colistin | Colistin | *Salmonella* | st_17 | 0.104 | Lineage-specific (ST17) |
| Colistin | Colistin | *Klebsiella* | genome_gc | 0.066 | Lineage proxy |
| Any | Various | *A. baumannii* | genome_gene_copies_sum | 0.083-0.125 | MDR clone signature |

Three patterns were consistent across species:

Quinolones: The top SHAP feature was always a PointFinder mutation feature (pf_has_gyrA, pf_n_quinolone_mutations, or pf_has_quinolone_mutation) regardless of species, with SHAP values of 0.07-0.16. This confirms that the same biological mechanism (target modification) is captured by the same computational feature across all organisms.

Beta-lactams: GC deviation and codon adaptation features dominated in Enterobacterales (*E. coli*, *Salmonella*, *Klebsiella*), reflecting the horizontally acquired nature of most beta-lactamases. In *S. aureus* and *A. baumannii*, RBS thermodynamic features were more important, reflecting the role of mecA and OXA expression regulation.

Colistin: MLST and population features dominated in all three species with colistin data (*Salmonella*, *E. coli*, *Klebsiella*). Because colistin resistance is predominantly chromosomal (mgrB inactivation, lipid A modification) without a well-characterized ARG, the model relies on lineage-specific signatures rather than gene-level features.

## Learning curves reveal species-specific data requirements

Learning curve analysis showed that sample size requirements depend on resistance mechanism complexity and resistance prevalence (Table 17). Per-antibiotic learning curve data are provided in Additional file 6: Table S6.

**Table 17. Learning curve summary by species.**

| Species | Plateau at n<=200 | Plateau at n<=1000 | Still improving | Data-hungry antibiotics |
|---------|------------------|-------------------|----------------|------------------------|
| *Salmonella* (n=9,111) | 13/20 (65%) | 18/20 (90%) | 2/20 | ceftazidime, ceftiofur |
| *E. coli* (n=5,236) | 14/35 (40%) | 21/35 (60%) | 14/35 | ertapenem, imipenem, pip-tazo, nitrofurantoin |
| *Klebsiella* (n=1,808) | 16/26 (62%) | 22/26 (85%) | 4/26 | tigecycline, fosfomycin |
| *S. aureus* (n=977) | 5/13 (38%) | 5/13 (38%) | 8/13 | linezolid, TMP-SMX, vancomycin |
| *A. baumannii* (n=757) | 8/13 (62%) | 10/13 (77%) | 3/13 | doripenem, imipenem |

### Species-specific patterns

*Salmonella*: The most "data-efficient" species. Thirteen of 20 antibiotics reached 95% of maximum F1 within 200 samples. Tetracycline, sulfisoxazole, and colistin plateaued at n=50 — the smallest tested size. This efficiency reflects both the large dataset available and the limited resistance gene diversity (1.1 unique ARGs per 100 isolates), meaning a small sample is sufficient to represent the gene-phenotype landscape.

*E. coli*: The most data-hungry species. Only 40% of antibiotics plateaued by n=200, and 14 of 35 were still improving at maximum sample size. Carbapenems (ertapenem: n=1,000 to reach 95% max; imipenem: never reached) and beta-lactam/inhibitor combinations (piperacillin-tazobactam: never reached) were the most data-demanding, consistent with the mechanistic diversity of resistance in this species.

*Klebsiella*: Despite a moderate dataset (1,808 isolates), 62% of antibiotics plateaued by n=200 — nearly as data-efficient as *Salmonella*. This reflects the clonal population structure of *K. pneumoniae*, where dominant clones (ST258, ST512) carry characteristic resistance profiles that are captured with relatively few representatives.

*S. aureus*: The most data-limited species. Only 38% of antibiotics plateaued by n=200, and 8 of 13 were still improving at maximum sample size. Linezolid (F1 improved from 0.067 at n=200 to 0.462 at n=531) and TMP-SMX (0.217 to 0.621) showed the steepest learning curves, indicating that these predictions would benefit substantially from larger datasets. This reflects both the small available sample (977 isolates) and the diversity of resistance mechanisms in Gram-positive organisms.

*A. baumannii*: Highly data-efficient despite the smallest dataset (757 isolates). Eight of 13 antibiotics plateaued by n=200, and most quinolone and aminoglycoside predictions were stable by n=100. The exception was carbapenem resistance (doripenem: did not reach 95% max; imipenem: required n=500), consistent with the greater mechanistic complexity of carbapenem resistance.

### Data efficiency varies by resistance mechanism complexity

Comparing across species, data efficiency correlated with the dominance of single-gene resistance mechanisms. *A. baumannii* (62% plateau by n=200) and *Klebsiella* (62%) were as data-efficient as *Salmonella* (65%), while *E. coli* (40%) and *S. aureus* (38%) required more samples. However, for *A. baumannii*, the high data efficiency may also reflect the high baseline resistance prevalence — when most isolates are resistant, fewer training examples are needed to learn the majority class.

This has practical implications for implementing WGS-based AMR prediction in surveillance: for species and antibiotics with well-characterized single-gene resistance mechanisms, prediction systems can be deployed with relatively small training datasets. For mechanistically diverse species like *E. coli*, training data must be pooled across institutions and geographies to capture the full spectrum of resistance mechanisms.

# Discussion

## A causal framework for WGS-based AMR prediction

Our systematic evaluation across five species and 107 antibiotic prediction tasks reveals that prediction performance is governed by three principal factors: the dominant resistance mechanism (ARG-mediated vs mutation-mediated vs complex), the resistance gene diversity of the species, and the available sample size. This framework explains both the successes and failures of genomic context features and provides actionable guidance for future tool development.

### Resistance mechanism type determines which features matter

The single most important determinant of both baseline performance and feature utility is whether resistance is primarily mediated by acquired genes, chromosomal mutations, or complex multi-factor mechanisms. For ARG-dominated resistance (tetracycline, ampicillin, sulfonamides), gene presence alone achieves F1 > 0.95, and expression features provide marginal additional value. For mutation-dominated resistance (quinolones, colistin), gene presence is uninformative, and PointFinder features are indispensable. For complex resistance (carbapenems in *E. coli*, tigecycline in *Klebsiella*), neither gene presence nor currently available mutation features are sufficient.

This explains why a "one-size-fits-all" feature set cannot optimally predict all antibiotics. The 17-step annotation pipeline provides a superset of features from which the ML models implicitly select the relevant subset per antibiotic. The ablation study confirms this: the gain from expression features is largest for aminoglycosides and beta-lactams (where silent genes are the prediction bottleneck), while the gain from PointFinder features is largest for quinolones (where chromosomal mutations are the bottleneck).

### Species gene diversity modulates difficulty, not mechanism

The five species span a wide range of resistance gene diversity, from 1.1 unique ARGs per 100 isolates (*Salmonella*) to 23.0 (*A. baumannii*). Naively, one might expect higher diversity to impair prediction. Instead, prediction difficulty is modulated by whether diversity translates into mechanistic ambiguity.

*A. baumannii* had the highest gene diversity yet best prediction (F1=0.941). However, this likely reflects high baseline resistance prevalence (ciprofloxacin 99% R, ceftriaxone 100% R) and strong ARG-phenotype concordance rather than the model learning clonal associations — MLST features contributed only +0.012 to mean F1. The dominant clone (ST2, 34% of isolates) showed substantial within-lineage phenotypic variation (e.g., imipenem 22% R, amikacin 49% R), confirming that lineage alone does not determine resistance.

Conversely, *E. coli* has moderate gene diversity (2.4 per 100) but includes the most mechanistically diverse resistance — ESBL spectrum variants, multiple carbapenemase families, porin/efflux-mediated resistance — resulting in the widest F1 range (0.217-0.990). The prediction difficulty in *E. coli* stems not from the number of genes but from the functional heterogeneity within drug classes.

### Sample size interacts with mechanism complexity

Our learning curves show that sample size requirements scale with mechanism complexity: 100-200 samples suffice for single-gene ARG-mediated resistance, 500-1,000 for multi-gene contexts, and larger datasets are needed for complex phenotypes. This has practical implications for surveillance programs — WGS-based prediction for common resistance phenotypes can be implemented with existing datasets, while rare resistance phenotypes (linezolid, tigecycline) require deliberate enrichment of resistant isolates in training sets.

## Expression features address a fundamental gap in gene catalogs

The improvement from expression features suggests that gene catalogs systematically overestimate resistance — they detect gene presence but cannot distinguish translated, functional protein from truncated or poorly expressed copies. This distinction is not trivial: in our data, 19-24% of *Klebsiella* and *S. aureus* observations involved gene-variant-matched silent genes, where the correct resistance enzyme was detected but did not confer phenotypic resistance.

Our discordance analysis provides direct biological evidence for the mechanisms underlying these silent genes. Gene-variant-matched ARGs in susceptible isolates showed 30-85% lower predicted translation initiation rates (reflecting poor ribosome binding), 4-14% lower gene coverage (indicating truncation), and weaker promoter activity (indicating reduced transcription) than identical genes in resistant isolates. Critically, these features can be computed from assembled genome sequences without transcriptomic data, using RBS thermodynamic modeling (OSTIR), promoter prediction (BPROM), and codon usage analysis — making them suitable for routine clinical application.

## PointFinder integration is non-negotiable for quinolone prediction

The near-universal absence of acquired quinolone resistance genes in resistant isolates (85-100% across five species) means that gene catalog approaches are fundamentally blind to the most common resistance mechanism for this drug class. Our ablation study quantifies the cost of this blindness: quinolone F1 scores ranged from 0.41-0.90 without mutation features, rising to 0.92-0.98 with PointFinder integration.

While the role of gyrA/parC mutations in quinolone resistance is well established, the quantitative impact on ML prediction performance had not been measured systematically across multiple species. Our results provide the strongest evidence to date that mutation detection must be integrated into — not separated from — the ML feature set. Notably, the PointFinder contribution was smallest in *A. baumannii* (ciprofloxacin: +0.033), where high resistance prevalence (87% ciprofloxacin-R) provides sufficient predictive signal even without explicit mutation encoding, suggesting that PointFinder is most critical for species with intermediate resistance prevalence.

## High resistance prevalence inflates prediction performance for *A. baumannii*

*A. baumannii* achieved the highest mean F1 (0.941) despite the smallest sample size, but this should be interpreted with caution. Analysis of the dominant clone (ST2, 34% of isolates) revealed substantial phenotypic variation within the lineage — amikacin resistance ranged from 49%, imipenem 22%, to ceftriaxone 100% — indicating that lineage alone does not determine resistance. Accordingly, MLST population features contributed only +0.012 to mean F1 over the model without them.

The high prediction performance is more likely driven by two factors: (1) the high baseline resistance prevalence for most antibiotics (ciprofloxacin 99% R, ceftriaxone 100% R, tetracycline 97% R), which makes classification inherently easier, and (2) the strong ARG-phenotype concordance for common resistance mechanisms (ARG-only baseline F1=0.870). Furthermore, 87% of isolates originated from US hospitals, raising the possibility that study-specific resistance patterns contribute to the apparent predictability.

These observations caution against interpreting high F1 as evidence of biological understanding. For antibiotics with near-universal resistance, prediction accuracy is dominated by class prevalence rather than genomic features. The more informative test of the model's value lies in species and antibiotics with intermediate resistance prevalence, where distinguishing resistant from susceptible isolates requires genuine feature discrimination.

## Remaining challenges and future directions

Three categories of resistance remain poorly predicted and represent targets for future work:

Regulatory mutations: Tigecycline resistance (*Klebsiella* F1=0.554) is mediated by ramA/ramR mutations that upregulate AcrAB-TolC efflux. These regulatory regions are not covered by current mutation databases. Expanding PointFinder or similar tools to include regulatory mutations — particularly in ramA/ramR, marA/marR, soxR/soxS, and ompR/envZ — would directly address this gap.

Chromosomal resistance without database coverage: Nitrofurantoin resistance (*E. coli* F1=0.217) involves nfsA/nfsB mutations that inactivate the prodrug activation pathway. No PointFinder database exists for these targets. Creating species-specific mutation databases for non-standard drug targets is a feasible but labor-intensive task.

Permeability-based resistance: Carbapenem resistance through porin loss (OmpC/OmpF) or efflux pump overexpression in *E. coli* represents a resistance category that is intrinsically difficult to predict from genome sequence alone, as it may involve regulatory rather than structural mutations, and the phenotypic effect depends on the combination of multiple genetic changes.

# Conclusions

This study reveals that antimicrobial resistance mechanisms fall into three distinct categories that determine whether gene presence alone can predict phenotype:

For tetracyclines, sulfonamides, and ampicillin, resistance gene presence is both necessary and sufficient for phenotype prediction (F1 > 0.95 from gene detection alone). These drug classes are mediated by dominant single genes (tet(B), sul1/sul2, blaTEM-1) whose presence reliably confers resistance.

For beta-lactam/inhibitor combinations, cephalosporins, aminoglycosides, and carbapenems, gene presence is necessary but not sufficient — expression context determines whether a detected gene confers resistance. Silent genes with poor ribosome binding (30-85% lower predicted translation initiation), truncated copies (4-14% lower coverage), and weak promoter activity are widespread in these drug classes. This gene silencing is the primary source of false-positive genotype-based predictions, and expression-level features computed from genome sequence alone can resolve it.

For quinolones, colistin, and vancomycin, resistance is primarily mediated by chromosomal point mutations or regulatory changes rather than acquired genes, making gene catalog approaches fundamentally blind to the dominant mechanism. PointFinder mutation detection is essential for these drug classes, accounting for >85% of quinolone resistance across all five species tested.

Antibiotics that remain poorly predicted — including nitrofurantoin, tigecycline, and linezolid — involve regulatory mutations and permeability changes not yet covered by existing databases, representing priorities for future database expansion.

The 17-step annotation pipeline, pre-trained models for 107 species-antibiotic combinations, and a Docker image for local deployment are released at https://github.com/tatsu1207/AMR.

# Declarations

## Availability of data and materials

All genome assemblies used in this study are publicly available from the NCBI Sequence Read Archive under the accession numbers listed in Additional file 7: Table S7. Antimicrobial susceptibility phenotypes were retrieved from NCBI BioSample metadata. Constructed feature matrices for all five species (107 species-antibiotic combinations, ~116 features per row), pre-trained models, and per-ARG annotation summaries are deposited at https://github.com/tatsu1207/AMR.

## Code availability

The annotation pipeline (`annotate_one.sh`), feature matrix builder (`build_feature_matrix_v2.py`), model training script (`train_models.py`), and pre-trained models are available at https://github.com/tatsu1207/AMR. A Docker image containing all dependencies and databases is available at https://github.com/tatsu1207/AMR.

## Authors' contributions

NJ and CV contributed equally to this work. NJ developed the genomic annotation pipeline, performed the ML analysis, and drafted the manuscript. CV developed per-gene feature extraction and contributed to model training. TU conceived the study, supervised the analysis, and revised the manuscript. All authors read and approved the final manuscript.

## Competing interests

The authors declare that they have no competing interests.

## Funding

[To be added]

## Acknowledgements

[To be added]
