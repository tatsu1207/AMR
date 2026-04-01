#!/bin/bash
# ============================================================================
# Minimal annotation pipeline for AMR phenotype prediction
# Produces only the features needed by the pre-trained RF models (~52 features)
#
# Skipped (not needed): geNomad, IntegronFinder, Infernal/sRNA, synteny,
#                       operon, IS orientation, IS distance
#
# Usage: annotate_minimal.sh <fasta_path> <output_dir> [threads]
# Output: <output_dir>/arg_context_summary.tsv
# ============================================================================
set -euo pipefail

FASTA_PATH="$1"
OUTPUT_DIR="$2"
THREADS="${3:-4}"

# ── Tool paths (use conda env PATH by default) ──
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# Detect conda env bin directory — prefer the 'amr' env explicitly
if [[ -d "/opt/conda/envs/amr/bin" ]]; then
    CONDA_BIN="/opt/conda/envs/amr/bin"
elif [[ -n "${CONDA_PREFIX:-}" && "$CONDA_PREFIX" == */envs/* ]]; then
    CONDA_BIN="$CONDA_PREFIX/bin"
else
    CONDA_BIN="$(dirname "$(command -v python3 2>/dev/null || echo /usr/bin/python3)")"
fi
AMRFINDER="${AMRFINDER:-${CONDA_BIN}/amrfinder}"
PRODIGAL="${PRODIGAL:-${CONDA_BIN}/prodigal}"
MLST_BIN="${MLST_BIN:-${CONDA_BIN}/mlst}"
# Ensure Perl uses conda env's modules (override any inherited PERL5LIB)
export PERL5LIB="${CONDA_BIN}/../lib/perl5/site_perl:${CONDA_BIN}/../lib/perl5/vendor_perl"
unset PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT 2>/dev/null || true
export PATH="${CONDA_BIN}:${PATH}"
export CONDA_PREFIX="${CONDA_BIN%/bin}"
OSTIR_BIN="${OSTIR_BIN:-ostir}"
BPROM="${BPROM:-${SCRIPT_DIR}/bin/bprom}"
BPROM_DATA="${BPROM_DATA:-${SCRIPT_DIR}/databases/bprom_data}"
POINTFINDER_DB="${POINTFINDER_DB:-${SCRIPT_DIR}/databases/pointfinder_db}"
RESFINDER_DB="${RESFINDER_DB:-${SCRIPT_DIR}/databases/resfinder_db}"

export BPROM
export OSTIR_BIN
export BPROM_DATA
export TSS_DATA="$BPROM_DATA"

BASENAME=$(basename "$FASTA_PATH")
SAMPLE_ID="${BASENAME%.fasta.gz}"
SAMPLE_ID="${SAMPLE_ID%.fasta}"
SAMPLE_ID="${SAMPLE_ID%.fa}"
SAMPLE_ID="${SAMPLE_ID%.fna}"

mkdir -p "$OUTPUT_DIR"

# Use original fasta or copy
if [[ "$FASTA_PATH" == *.gz ]]; then
    ASSEMBLY="${OUTPUT_DIR}/assembly.fasta"
    gunzip -c "$FASTA_PATH" > "$ASSEMBLY"
else
    ASSEMBLY="$FASTA_PATH"
fi

echo "[$(date +%H:%M:%S)] Starting minimal annotation for $SAMPLE_ID"

# ════════════════════════════════════════════════════════════════════════
# Step 1: MLST (species identification + sequence type)
# ════════════════════════════════════════════════════════════════════════
MLST_TSV="${OUTPUT_DIR}/mlst_results.tsv"
if [[ ! -s "$MLST_TSV" ]]; then
    echo "[$(date +%H:%M:%S)] Step 1/6: MLST"
    $MLST_BIN "$ASSEMBLY" > "$MLST_TSV" 2>/dev/null || true
fi

# Parse species and ST
MLST_SCHEME=$(awk -F'\t' '{print $2}' "$MLST_TSV" 2>/dev/null || echo "unknown")
MLST_ST=$(awk -F'\t' '{print $3}' "$MLST_TSV" 2>/dev/null || echo "unknown")
echo "  Species scheme: $MLST_SCHEME, ST: $MLST_ST"

# Map MLST scheme to AMRFinderPlus organism name (for point mutation detection)
AMR_ORGANISM=""
case "$MLST_SCHEME" in
    ecoli*) AMR_ORGANISM="Escherichia" ;;
    senterica*) AMR_ORGANISM="Salmonella" ;;
    saureus*) AMR_ORGANISM="Staphylococcus_aureus" ;;
    klebsiella*|kpneumoniae*) AMR_ORGANISM="Klebsiella_pneumoniae" ;;
    abaumannii*) AMR_ORGANISM="Acinetobacter_baumannii" ;;
esac

# ════════════════════════════════════════════════════════════════════════
# Step 2: AMRFinderPlus (ARG detection + point mutations)
# ════════════════════════════════════════════════════════════════════════
AMR_TSV="${OUTPUT_DIR}/amrfinderplus.tsv"
AMR_MUTATIONS="${OUTPUT_DIR}/amrfinder_mutations.tsv"
if [[ ! -s "$AMR_TSV" ]]; then
    echo "[$(date +%H:%M:%S)] Step 2/6: AMRFinderPlus"
    AMR_ARGS=(-n "$ASSEMBLY" -o "$AMR_TSV" --plus --threads "$THREADS")
    if [[ -n "$AMR_ORGANISM" ]]; then
        AMR_ARGS+=(--organism "$AMR_ORGANISM" --mutation_all "$AMR_MUTATIONS")
    fi
    $AMRFINDER "${AMR_ARGS[@]}" >> "${OUTPUT_DIR}/amrfinder.log" 2>&1 || true
fi

# ════════════════════════════════════════════════════════════════════════
# Step 3: Prodigal (gene prediction — needed for CAI, GC)
# ════════════════════════════════════════════════════════════════════════
GENES_GFF="${OUTPUT_DIR}/genes.gff"
GENES_FNA="${OUTPUT_DIR}/genes.fna"
if [[ ! -s "$GENES_FNA" ]]; then
    echo "[$(date +%H:%M:%S)] Step 3/6: Prodigal"
    $PRODIGAL -i "$ASSEMBLY" -o "$GENES_GFF" -f gff \
        -d "$GENES_FNA" -p meta -q >> "${OUTPUT_DIR}/prodigal.log" 2>&1 || true
fi

# ════════════════════════════════════════════════════════════════════════
# Steps 4-6: Promoter (BPROM), RBS (OSTIR), and genomic features
# All computed per-ARG in a single Python pass
# ════════════════════════════════════════════════════════════════════════
SUMMARY="${OUTPUT_DIR}/arg_context_summary.tsv"
if [[ ! -f "$SUMMARY" ]]; then
    echo "[$(date +%H:%M:%S)] Steps 4-6/6: Promoter + RBS + Genomic features"
    python3 - "$ASSEMBLY" "$AMR_TSV" "$GENES_FNA" "$GENES_GFF" "$AMR_MUTATIONS" "$MLST_SCHEME" "$MLST_ST" "$SUMMARY" <<'PYEOF'
import sys, os, csv, json, subprocess, re, math
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

assembly_path = sys.argv[1]
amr_tsv = sys.argv[2]
genes_fna = sys.argv[3]
genes_gff = sys.argv[4]
amr_mutations_tsv = sys.argv[5]
mlst_scheme = sys.argv[6]
mlst_st = sys.argv[7]
output_tsv = sys.argv[8]

BPROM = os.environ.get("BPROM", "bprom")
OSTIR_BIN = os.environ.get("OSTIR_BIN", "ostir")
BPROM_DATA = os.environ.get("TSS_DATA", "/tmp/bprom/data")

# ── Load assembly ──
contigs = {}
genome_gc_count = 0
genome_len = 0
for rec in SeqIO.parse(assembly_path, "fasta"):
    seq = str(rec.seq).upper()
    contigs[rec.id.split()[0]] = seq
    genome_len += len(seq)
    genome_gc_count += seq.count('G') + seq.count('C')
genome_gc = genome_gc_count / genome_len if genome_len > 0 else 0.5

# ── Load Prodigal genes for CAI reference ──
all_cds_seqs = []
if os.path.exists(genes_fna):
    for rec in SeqIO.parse(genes_fna, "fasta"):
        s = str(rec.seq).upper()
        if len(s) >= 300 and len(s) % 3 == 0:
            all_cds_seqs.append(s)

# Compute genome codon frequencies for CAI
codon_count = defaultdict(int)
aa_codons = defaultdict(list)
CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
for codon, aa in CODON_TABLE.items():
    if aa != '*':
        aa_codons[aa].append(codon)
for seq in all_cds_seqs[:2000]:
    for i in range(0, len(seq)-2, 3):
        c = seq[i:i+3]
        if c in CODON_TABLE:
            codon_count[c] += 1

# Relative adaptiveness (w values) for CAI
w_values = {}
for aa, codons in aa_codons.items():
    if aa == 'M' or aa == 'W':
        for c in codons:
            w_values[c] = 1.0
        continue
    counts = [codon_count.get(c, 0) for c in codons]
    max_count = max(counts) if counts else 1
    for c, cnt in zip(codons, counts):
        w_values[c] = cnt / max_count if max_count > 0 else 0

def compute_cai(seq):
    """Compute Codon Adaptation Index."""
    if len(seq) < 6:
        return 0.0
    log_w_sum = 0
    n = 0
    for i in range(0, len(seq) - 2, 3):
        c = seq[i:i+3]
        if c in w_values and w_values[c] > 0:
            log_w_sum += math.log(w_values[c])
            n += 1
    return math.exp(log_w_sum / n) if n > 0 else 0.0

def compute_rare_codons(seq):
    """Count rare codons (w < 0.1) and clusters."""
    if len(seq) < 6:
        return 0.0, 0
    rare = 0
    total = 0
    clusters = 0
    in_cluster = False
    for i in range(0, len(seq) - 2, 3):
        c = seq[i:i+3]
        if c in w_values:
            total += 1
            if w_values[c] < 0.1:
                rare += 1
                if not in_cluster:
                    clusters += 1
                    in_cluster = True
            else:
                in_cluster = False
    return rare / total if total > 0 else 0.0, clusters

def run_bprom(upstream_seq):
    """Run BPROM on upstream sequence, return LDF score and TF sites."""
    if not os.path.exists(BPROM) or len(upstream_seq) < 50:
        return None, None, None, None
    try:
        proc = subprocess.run([BPROM], input=f">query\n{upstream_seq}\n",
                            capture_output=True, text=True, timeout=30,
                            env={**os.environ, 'TSS_DATA': BPROM_DATA})
        ldf = None
        tf_count = 0
        distance = None
        for line in proc.stdout.split('\n'):
            if 'Linear Discriminant Function' in line or 'LDF' in line:
                m = re.search(r'LDF[=: ]+([0-9.]+)', line)
                if m:
                    ldf = float(m.group(1))
            if 'Promoter Pos' in line:
                m = re.search(r'Pos:\s+(\d+)', line)
                if m:
                    distance = int(m.group(1))
            if re.match(r'\s+[A-Za-z]+\s+\d+', line):
                tf_count += 1
        # AT ratio of upstream region
        at_ratio = (upstream_seq.count('A') + upstream_seq.count('T')) / len(upstream_seq) if upstream_seq else 0.5
        return ldf, tf_count, distance, at_ratio
    except Exception:
        return None, None, None, None

def run_ostir(rbs_region):
    """Run OSTIR on RBS region, return expression and dG values."""
    if len(rbs_region) < 20:
        return None, None, None
    try:
        proc = subprocess.run([OSTIR_BIN, "-i", rbs_region, "-t", "string", "-p"],
                            capture_output=True, text=True, timeout=30)
        for line in proc.stdout.strip().split('\n'):
            if line.startswith('start') or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 4:
                expression = float(parts[1]) if parts[1] != 'NA' else None
                dg_total = float(parts[2]) if parts[2] != 'NA' else None
                dg_mrna = float(parts[3]) if parts[3] != 'NA' else None
                return expression, dg_total, dg_mrna
        return None, None, None
    except Exception:
        return None, None, None

# ── Load AMRFinderPlus results ──
amr_hits = []
if os.path.exists(amr_tsv):
    with open(amr_tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            contig = row.get('Contig id', '').split()[0]
            start = int(row.get('Start', 0))
            stop = int(row.get('Stop', 0))
            strand = row.get('Strand', '+')
            gene = row.get('Element symbol', row.get('Gene symbol', ''))
            etype = row.get('Type', '')
            drug_class = row.get('Class', '')
            subclass = row.get('Subclass', '')
            identity = float(row.get('% Identity to reference', 0))
            coverage = float(row.get('% Coverage of reference', 0))
            amr_hits.append({
                'gene': gene, 'element_type': etype, 'drug_class': drug_class,
                'subclass': subclass, 'contig': contig,
                'start': min(start, stop), 'end': max(start, stop),
                'strand': strand, 'identity': identity, 'coverage': coverage,
            })

# ── Load point mutations from AMRFinderPlus --mutation_all output ──
point_mutations = []
if os.path.exists(amr_mutations_tsv):
    with open(amr_mutations_tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            subtype = row.get('Subtype', '')
            element = row.get('Element symbol', '')
            # Skip wildtype entries, only keep actual mutations
            if subtype == 'POINT' and 'WILDTYPE' not in row.get('Element name', ''):
                point_mutations.append(element)

# ── Process each ARG: promoter, RBS, CAI, GC ──
rows = []
gene_copies = defaultdict(int)
for hit in amr_hits:
    gene_copies[hit['gene']] += 1

for hit in amr_hits:
    contig_seq = contigs.get(hit['contig'], '')
    if not contig_seq:
        continue

    start, end, strand = hit['start'], hit['end'], hit['strand']
    gene_seq = contig_seq[start:end]

    # Strand-aware upstream extraction for promoter/RBS
    if strand == '-':
        # Gene is on minus strand: upstream is AFTER the end position
        upstream_start = end
        upstream_end = min(end + 500, len(contig_seq))
        upstream_seq = str(Seq(contig_seq[upstream_start:upstream_end]).reverse_complement())
        rbs_start = end
        rbs_end = min(end + 70, len(contig_seq))
        rbs_seq = str(Seq(contig_seq[rbs_start:rbs_end]).reverse_complement())
        gene_seq = str(Seq(gene_seq).reverse_complement())
    else:
        upstream_end = start
        upstream_start = max(0, start - 500)
        upstream_seq = contig_seq[upstream_start:upstream_end]
        rbs_start = max(0, start - 70)
        rbs_end = start
        rbs_seq = contig_seq[rbs_start:rbs_end]

    # Promoter (BPROM)
    ldf, tf_sites, prom_distance, at_ratio = run_bprom(upstream_seq)

    # RBS (OSTIR)
    expression, dg_total, dg_mrna = run_ostir(rbs_seq)

    # CAI and rare codons
    cai = compute_cai(gene_seq)
    rare_pct, rare_clusters = compute_rare_codons(gene_seq)

    # GC content
    if len(gene_seq) > 0:
        gene_gc = (gene_seq.count('G') + gene_seq.count('C')) / len(gene_seq)
        gc_deviation = gene_gc - genome_gc
    else:
        gene_gc = genome_gc
        gc_deviation = 0

    rows.append({
        'gene': hit['gene'],
        'element_type': hit['element_type'],
        'drug_class': hit['drug_class'],
        'contig': hit['contig'],
        'start': hit['start'],
        'end': hit['end'],
        'identity': hit['identity'],
        'coverage': hit['coverage'],
        'promoter_ldf': ldf if ldf is not None else '',
        'promoter_tf_sites': tf_sites if tf_sites is not None else '',
        'promoter_distance': prom_distance if prom_distance is not None else '',
        'promoter_up_at_ratio': f"{at_ratio:.4f}" if at_ratio is not None else '',
        'rbs_expression': expression if expression is not None else '',
        'rbs_dg_total': dg_total if dg_total is not None else '',
        'rbs_dg_mrna': dg_mrna if dg_mrna is not None else '',
        'cai': f"{cai:.4f}" if cai else '',
        'rare_codon_pct': f"{rare_pct:.4f}" if rare_pct is not None else '',
        'rare_codon_clusters': rare_clusters,
        'gene_gc': f"{gene_gc:.4f}",
        'genome_gc': f"{genome_gc:.4f}",
        'gc_deviation': f"{gc_deviation:.4f}",
        'gene_copies': gene_copies[hit['gene']],
        'mlst_scheme': mlst_scheme,
        'mlst_st': mlst_st,
        'point_mutations': ';'.join(point_mutations) if point_mutations else '',
    })

# Write summary
COLUMNS = [
    'gene', 'element_type', 'drug_class', 'contig', 'start', 'end',
    'identity', 'coverage',
    'promoter_ldf', 'promoter_tf_sites', 'promoter_distance', 'promoter_up_at_ratio',
    'rbs_expression', 'rbs_dg_total', 'rbs_dg_mrna',
    'cai', 'rare_codon_pct', 'rare_codon_clusters',
    'gene_gc', 'genome_gc', 'gc_deviation',
    'gene_copies', 'mlst_scheme', 'mlst_st', 'point_mutations',
]

with open(output_tsv, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=COLUMNS, delimiter='\t')
    writer.writeheader()
    for row in rows:
        writer.writerow(row)

print(f"  Summary: {len(rows)} ARG annotations written to {output_tsv}")
PYEOF
fi

echo "[$(date +%H:%M:%S)] Done. Output: $SUMMARY"
