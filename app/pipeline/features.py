"""Extract features from annotation results for a single sample."""
import os
import json
import csv


def extract_features(job_dir: str, species: str, mlst_st: str) -> dict:
    """Build feature dictionary from all annotation outputs."""
    features = {}

    # 1. Load AMR results from amrfinderplus.tsv
    amr_tsv = os.path.join(job_dir, "amrfinderplus.tsv")
    n_amr = 0
    if os.path.exists(amr_tsv):
        with open(amr_tsv) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                etype = row.get('Type', '')
                if etype in ('AMR', 'STRESS', 'VIRULENCE') or 'AMR' in str(etype):
                    n_amr += 1
    features['n_amr_genes'] = n_amr

    # 2. Load per-ARG features from arg_context_summary.tsv
    summary_path = os.path.join(job_dir, "arg_context_summary.tsv")
    if os.path.exists(summary_path):
        summary_features = parse_arg_summary(summary_path)
        features.update(summary_features)

    # 3. Load genomic features (fallback if not in summary)
    gf_path = os.path.join(job_dir, "genomic_features.json")
    if os.path.exists(gf_path):
        with open(gf_path) as f:
            gf = json.load(f)
        for k, v in gf.items():
            if k not in features:
                features[k] = v

    # 4. Load PointFinder results
    pf_path = os.path.join(job_dir, "pointfinder", "PointFinder_results.txt")
    if os.path.exists(pf_path):
        pf_features = parse_pointfinder_results(pf_path)
        features.update(pf_features)
    else:
        features['target_has_point_mutation'] = 0
        features['genome_has_point_mutation'] = 0
        features['pf_has_gyrA'] = 0
        features['pf_has_parC'] = 0
        features['pf_has_parE'] = 0
        features['pf_has_gyrB'] = 0
        features['pf_has_ampC_promoter'] = 0
        features['pf_n_mutations'] = 0
        features['pf_has_quinolone_mutation'] = 0
        features['pf_n_quinolone_mutations'] = 0

    # 5. A. baumannii custom gyrA
    gyra_path = os.path.join(job_dir, "abaumannii_gyra.json")
    if os.path.exists(gyra_path):
        with open(gyra_path) as f:
            gyra = json.load(f)
        features['pf_has_gyrA_Ser83Leu'] = gyra.get('has_gyrA_Ser83Leu', 0)
        if gyra.get('has_gyrA_Ser83Leu', 0):
            features['pf_has_gyrA'] = 1
            features['pf_has_quinolone_mutation'] = 1
            features['pf_n_quinolone_mutations'] = max(features.get('pf_n_quinolone_mutations', 0), 1)
            features['genome_has_point_mutation'] = 1

    # 6. MLST features
    features['mlst_st'] = mlst_st
    # ST one-hot features will be added by the prediction step using model metadata

    return features


def parse_arg_summary(summary_path: str) -> dict:
    """Parse arg_context_summary.tsv to extract per-ARG features.

    The model uses features like target_promoter_ldf, target_rbs_expression,
    target_cai, target_gene_gc, etc. These are aggregated across all ARGs
    (e.g., max promoter LDF, mean CAI).
    """
    features = {}
    rows = []
    try:
        with open(summary_path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                rows.append(row)
    except Exception:
        return features

    if not rows:
        return features

    def safe_float(val, default=None):
        try:
            return float(val) if val not in ('', None) else default
        except (ValueError, TypeError):
            return default

    # Aggregate per-ARG features
    cai_vals = [safe_float(r['cai']) for r in rows if safe_float(r['cai']) is not None]
    rare_vals = [safe_float(r['rare_codon_pct']) for r in rows if safe_float(r['rare_codon_pct']) is not None]
    gc_vals = [safe_float(r['gene_gc']) for r in rows if safe_float(r['gene_gc']) is not None]
    gc_dev_vals = [safe_float(r['gc_deviation']) for r in rows if safe_float(r['gc_deviation']) is not None]
    ldf_vals = [safe_float(r['promoter_ldf']) for r in rows if safe_float(r['promoter_ldf']) is not None]
    rbs_vals = [safe_float(r['rbs_expression']) for r in rows if safe_float(r['rbs_expression']) is not None]
    identity_vals = [safe_float(r['identity']) for r in rows if safe_float(r['identity']) is not None]
    coverage_vals = [safe_float(r['coverage']) for r in rows if safe_float(r['coverage']) is not None]
    copies_vals = [safe_float(r['gene_copies']) for r in rows if safe_float(r['gene_copies']) is not None]

    # Target-level features (use first/representative ARG or aggregate)
    if cai_vals:
        features['target_cai_mean'] = sum(cai_vals) / len(cai_vals)
        features['target_cai_max'] = max(cai_vals)
    if rare_vals:
        features['target_rare_codon_pct'] = sum(rare_vals) / len(rare_vals)
    if gc_vals:
        features['target_gene_gc_mean'] = sum(gc_vals) / len(gc_vals)
    if gc_dev_vals:
        features['target_gc_deviation_mean'] = sum(gc_dev_vals) / len(gc_dev_vals)
        features['target_gc_deviation_max'] = max(gc_dev_vals, key=abs)
    if ldf_vals:
        features['target_promoter_ldf'] = max(ldf_vals)
    if rbs_vals:
        features['target_rbs_expression'] = max(rbs_vals)
    if identity_vals:
        features['target_identity_mean'] = sum(identity_vals) / len(identity_vals)
        features['target_identity_min'] = min(identity_vals)
    if coverage_vals:
        features['target_coverage_mean'] = sum(coverage_vals) / len(coverage_vals)
        features['target_coverage_min'] = min(coverage_vals)
    if copies_vals:
        features['target_gene_copies_max'] = max(copies_vals)
        features['genome_gene_copies_total'] = sum(copies_vals)

    # Genome-level from first row
    genome_gc = safe_float(rows[0].get('genome_gc'))
    if genome_gc is not None:
        features['genome_gc'] = genome_gc

    # has_target_arg
    features['has_target_arg'] = 1
    features['n_target_args'] = len(rows)

    # Functionality flags
    low_rbs = sum(1 for r in rows if safe_float(r.get('rbs_expression'), 999) < 1.0)
    low_prom = sum(1 for r in rows if safe_float(r.get('promoter_ldf')) is None)
    features['arg_low_rbs'] = 1 if low_rbs > 0 else 0
    features['arg_low_promoter'] = 1 if low_prom > 0 else 0

    return features


def parse_pointfinder_results(pf_path: str) -> dict:
    """Parse PointFinder output into feature dict."""
    features = {
        'target_has_point_mutation': 0,
        'genome_has_point_mutation': 0,
        'pf_has_gyrA': 0,
        'pf_has_parC': 0,
        'pf_has_parE': 0,
        'pf_has_gyrB': 0,
        'pf_has_ampC_promoter': 0,
        'pf_n_mutations': 0,
        'pf_has_quinolone_mutation': 0,
        'pf_n_quinolone_mutations': 0,
    }

    try:
        with open(pf_path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene = row.get('Mutation', row.get('Gene_ID', '')).lower()
                features['pf_n_mutations'] += 1
                features['genome_has_point_mutation'] = 1
                features['target_has_point_mutation'] = 1

                if 'gyra' in gene:
                    features['pf_has_gyrA'] = 1
                    features['pf_has_quinolone_mutation'] = 1
                    features['pf_n_quinolone_mutations'] += 1
                elif 'parc' in gene:
                    features['pf_has_parC'] = 1
                    features['pf_has_quinolone_mutation'] = 1
                    features['pf_n_quinolone_mutations'] += 1
                elif 'pare' in gene:
                    features['pf_has_parE'] = 1
                elif 'gyrb' in gene:
                    features['pf_has_gyrB'] = 1
                elif 'ampc' in gene or 'promoter' in gene:
                    features['pf_has_ampC_promoter'] = 1
    except Exception:
        pass

    return features
