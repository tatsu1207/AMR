#!/bin/bash
# annotate_one_v2.sh — Annotate a single WGS assembly FASTA file (17 steps + summary)
# Usage: annotate_one_v2.sh <fasta_path> <results_dir> <threads>
set -euo pipefail

FASTA_PATH="$1"
RESULTS_DIR="$2"
THREADS="${3:-2}"

# ── Tool paths ──
CONDA_BASE="/home/unnot/miniforge3"
AMRFINDER="${CONDA_BASE}/envs/radar-amr/bin/amrfinder"
MEFINDER="${CONDA_BASE}/envs/radar/bin/mefinder"
MOB_RECON="${CONDA_BASE}/envs/radar-plasmid/bin/mob_recon"
GENOMAD="${CONDA_BASE}/envs/radar-genomad/bin/genomad"
SKANI="${CONDA_BASE}/envs/radar-skani/bin/skani"
BPROM="${CONDA_BASE}/envs/radar/bin/bprom"
OSTIR="${CONDA_BASE}/envs/radar/bin/ostir"
PRODIGAL="/usr/bin/prodigal"
CMSCAN="/usr/bin/cmscan"
MLST="${CONDA_BASE}/envs/radar-mlst/bin/mlst"
INTEGRON_FINDER="${CONDA_BASE}/envs/radar-integron/bin/integron_finder"

# ── Database paths ──
GENOMAD_DB="/home/unnot/github/radar/databases/genomad_db"
SKANI_DB="/home/unnot/WGS/phenotype/skani_db/skani-gtdb-r220-sketch"
RFAM_CM="/home/unnot/github/radar/databases/Rfam.cm"
POINTFINDER_DB="/home/unnot/github/radar/databases/pointfinder_db"
RESFINDER_DB="/home/unnot/github/radar/databases/resfinder_db"
export TSS_DATA="/tmp/bprom/data"

UPSTREAM_BP=500
PROXIMITY_BP=5000

# ── Sample setup ──
BASENAME=$(basename "$FASTA_PATH")
SAMPLE_ID="${BASENAME%.fasta.gz}"
SAMPLE_DIR="${RESULTS_DIR}/${SAMPLE_ID}"
mkdir -p "$SAMPLE_DIR"

if [[ -f "${SAMPLE_DIR}/arg_context_summary.tsv" ]]; then
    echo "[SKIP] $SAMPLE_ID — already complete"
    exit 0
fi

echo "[START] $SAMPLE_ID"

# Decompress if gzipped
if [[ "$FASTA_PATH" == *.gz ]]; then
    ASSEMBLY="${SAMPLE_DIR}/assembly.fasta"
    if [[ ! -f "$ASSEMBLY" ]]; then
        gunzip -c "$FASTA_PATH" > "$ASSEMBLY"
    fi
else
    ASSEMBLY="$FASTA_PATH"
fi

ASSEMBLY_NAME=$(basename "$ASSEMBLY" .fasta)

# ════════════════════════════════════════════════════════════════════════════
# Step 1: AMRFinderPlus
# ════════════════════════════════════════════════════════════════════════════
AMR_DIR="${SAMPLE_DIR}/amr"
AMR_TSV="${AMR_DIR}/amrfinderplus.tsv"
mkdir -p "$AMR_DIR"

if [[ ! -f "$AMR_TSV" ]]; then
    if ! $AMRFINDER -n "$ASSEMBLY" -o "$AMR_TSV" --plus --threads "$THREADS" \
        >> "${AMR_DIR}/amrfinder.log" 2>&1; then
        echo "[FAIL] $SAMPLE_ID — AMRFinderPlus failed"
        exit 1
    fi
fi

# ════════════════════════════════════════════════════════════════════════════
# Step 2: BPROM promoter prediction
# ════════════════════════════════════════════════════════════════════════════
PROM_DIR="${SAMPLE_DIR}/promoter"
PROM_TSV="${PROM_DIR}/promoter_results.tsv"
mkdir -p "$PROM_DIR"

if [[ ! -f "$PROM_TSV" ]] && [[ -f "$AMR_TSV" ]]; then
    python3 - "$ASSEMBLY" "$AMR_TSV" "$PROM_DIR" "$PROM_TSV" "$BPROM" "$UPSTREAM_BP" <<'PYEOF' 2>/dev/null || true
import sys, os, re, subprocess, csv
from Bio import SeqIO
from Bio.Seq import Seq

assembly_path, amr_tsv, prom_dir, out_tsv, bprom_bin, upstream_bp = sys.argv[1:7]
upstream_bp = int(upstream_bp)
contigs = {r.id: str(r.seq) for r in SeqIO.parse(assembly_path, "fasta")}
results = []
with open(amr_tsv) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        etype = row.get("Element type", row.get("Type", ""))
        if etype == "VIRULENCE":
            continue
        gene = row.get("Gene symbol", row.get("Element symbol", "unknown"))
        contig_id = row.get("Contig id", "")
        strand = row.get("Strand", "+")
        try:
            start = int(row.get("Start", ""))
            stop = int(row.get("Stop", ""))
        except (ValueError, TypeError):
            continue
        seq = contigs.get(contig_id, "")
        if not seq:
            continue
        # Extract upstream region based on strand
        if strand == "-":
            # Minus strand: promoter is after the stop coordinate (reverse complement)
            up_start = stop
            up_end = min(len(seq), stop + upstream_bp)
            up_seq = str(Seq(seq[up_start:up_end]).reverse_complement())
        else:
            # Plus strand: promoter is before the start coordinate
            up_start = max(0, start - upstream_bp)
            up_seq = seq[up_start:start]
        if len(up_seq) < 50:
            continue
        tmp_fa = os.path.join(prom_dir, f"tmp_{os.getpid()}_{gene}.fasta")
        tmp_out = os.path.join(prom_dir, f"tmp_{os.getpid()}_{gene}.txt")
        with open(tmp_fa, "w") as fout:
            fout.write(f">upstream_{gene}\n")
            for i in range(0, len(up_seq), 60):
                fout.write(up_seq[i:i+60] + "\n")
        try:
            subprocess.run([bprom_bin, tmp_fa, tmp_out], capture_output=True, text=True, timeout=60)
            output = ""
            if os.path.exists(tmp_out):
                with open(tmp_out) as fh:
                    output = fh.read()
            ldf_score = None; tf_count = 0; prom_dist = None; up_at = None
            for line in output.splitlines():
                m = re.search(r'Promoter Pos:\s+(\d+)\s+LDF-?\s+([\d.]+)', line.strip())
                if m and ldf_score is None:
                    pos = int(m.group(1))
                    ldf_score = float(m.group(2))
                    prom_dist = len(up_seq) - pos
                if re.match(r'\s*\w+:\s+[ACGT]+ at position', line.strip()):
                    tf_count += 1
            if prom_dist is not None:
                pp = len(up_seq) - prom_dist
                us, ue = max(0, pp - 20), pp
                if ue > us:
                    region = up_seq[us:ue]
                    up_at = sum(1 for c in region.upper() if c in "AT") / len(region)
            results.append([gene, contig_id, start, ldf_score or "", tf_count,
                            prom_dist or "", f"{up_at:.3f}" if up_at else ""])
        except Exception:
            pass
        finally:
            for tmp in (tmp_fa, tmp_out):
                if os.path.exists(tmp):
                    os.remove(tmp)
with open(out_tsv, "w") as f:
    f.write("gene\tcontig\tstart\tldf_score\ttf_binding_sites\tpromoter_distance\tup_element_at_ratio\n")
    for r in results:
        f.write("\t".join(str(x) for x in r) + "\n")
PYEOF
fi

# ════════════════════════════════════════════════════════════════════════════
# Step 3: OSTIR RBS prediction
# ════════════════════════════════════════════════════════════════════════════
RBS_DIR="${SAMPLE_DIR}/rbs"
RBS_TSV="${RBS_DIR}/rbs_results.tsv"
mkdir -p "$RBS_DIR"

if [[ ! -f "$RBS_TSV" ]] && [[ -f "$AMR_TSV" ]]; then
    python3 - "$ASSEMBLY" "$AMR_TSV" "$RBS_TSV" "$OSTIR" <<'PYEOF' 2>/dev/null || true
import sys, os, subprocess, csv
from Bio import SeqIO
from Bio.Seq import Seq

assembly_path, amr_tsv, out_tsv, ostir_bin = sys.argv[1:5]
contigs = {r.id: str(r.seq) for r in SeqIO.parse(assembly_path, "fasta")}
results = []
with open(amr_tsv) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        etype = row.get("Element type", row.get("Type", ""))
        if etype == "VIRULENCE":
            continue
        gene = row.get("Gene symbol", row.get("Element symbol", "unknown"))
        contig_id = row.get("Contig id", "")
        strand = row.get("Strand", "+")
        try:
            start = int(row.get("Start", ""))
            stop = int(row.get("Stop", ""))
        except (ValueError, TypeError):
            continue
        seq = contigs.get(contig_id, "")
        if not seq:
            continue
        # Extract RBS region based on strand
        # RBS is ~20 bp upstream of start codon + first ~19 bp of gene
        if strand == "-":
            # Minus strand: start codon is at stop coordinate
            rbs_start = stop
            rbs_end = min(len(seq), stop + 51)
            rbs_seq = str(Seq(seq[rbs_start:rbs_end]).reverse_complement())
        else:
            # Plus strand: start codon is at start coordinate
            rbs_start = max(0, start - 51)
            rbs_end = min(len(seq), start + 19)
            rbs_seq = seq[rbs_start:rbs_end]
        if len(rbs_seq) < 30:
            continue
        try:
            result = subprocess.run(["conda", "run", "-n", "radar", ostir_bin, "-i", rbs_seq, "-t", "string", "-p"],
                                    capture_output=True, text=True, timeout=60)
            if result.returncode != 0:
                continue
            lines = result.stdout.strip().splitlines()
            header = None; data_line = None
            for line in lines:
                fields = line.split()
                if len(fields) >= 5 and "expression" in line.lower():
                    header = fields
                elif header and len(fields) >= len(header) and fields[0] in ("ATG", "GTG", "TTG"):
                    data_line = fields
                    break
            if header and data_line:
                col_map = {h.lower(): i for i, h in enumerate(header)}
                expression = data_line[col_map.get("expression", 2)]
                dg_total = data_line[col_map.get("dg_total", 4)]
                dg_mrna = data_line[col_map.get("dg_mrna", 6)]
                results.append([gene, contig_id, start, expression, dg_total, dg_mrna])
        except Exception:
            pass
with open(out_tsv, "w") as f:
    f.write("gene\tcontig\tstart\texpression\tdg_total\tdg_mrna\n")
    for r in results:
        f.write("\t".join(str(x) for x in r) + "\n")
PYEOF
fi

# ════════════════════════════════════════════════════════════════════════════
# (MEFinder and MOB-recon removed — geNomad provides plasmid + IS/transposase detection)
MOB_DIR="${SAMPLE_DIR}/mobility"
MEF_CSV="${MOB_DIR}/mefinder_output.csv"
MOB_REPORT=""

# ════════════════════════════════════════════════════════════════════════════
# Steps 4-8: Assembly-level annotations (run concurrently where possible)
#   geNomad is the bottleneck (286s); skani/Prodigal/MLST run in parallel with it
# ════════════════════════════════════════════════════════════════════════════

# ── geNomad (prophage, plasmid classification, transposase detection) ──
GENOMAD_DIR="${SAMPLE_DIR}/genomad"
VIRUS_SUMMARY="${GENOMAD_DIR}/${ASSEMBLY_NAME}_summary/${ASSEMBLY_NAME}_virus_summary.tsv"
mkdir -p "$GENOMAD_DIR"

if [[ ! -f "$VIRUS_SUMMARY" ]]; then
    $GENOMAD end-to-end "$ASSEMBLY" "$GENOMAD_DIR" "$GENOMAD_DB" \
        --threads "$THREADS" --cleanup \
        >> "${GENOMAD_DIR}/genomad.log" 2>&1 || true
fi &
PID_GENOMAD=$!

# ── skani species identification (runs in parallel with geNomad) ──
SPECIES_DIR="${SAMPLE_DIR}/species"
SKANI_TSV="${SPECIES_DIR}/skani_results.tsv"
mkdir -p "$SPECIES_DIR"

if [[ ! -f "$SKANI_TSV" ]] && [[ -d "$SKANI_DB" ]]; then
    $SKANI search "$ASSEMBLY" -d "$SKANI_DB" -o "$SKANI_TSV" -t 1 -n 5 \
        >> "${SPECIES_DIR}/skani.log" 2>&1 || true
fi &
PID_SKANI=$!

# ════════════════════════════════════════════════════════════════════════════
# Step 8: Prodigal gene prediction
# ════════════════════════════════════════════════════════════════════════════
PRODIGAL_DIR="${SAMPLE_DIR}/prodigal"
PRODIGAL_GFF="${PRODIGAL_DIR}/genes.gff"
PRODIGAL_FNA="${PRODIGAL_DIR}/genes.fna"
mkdir -p "$PRODIGAL_DIR"

if [[ ! -f "$PRODIGAL_GFF" ]]; then
    $PRODIGAL -i "$ASSEMBLY" -o "$PRODIGAL_GFF" -f gff -d "$PRODIGAL_FNA" -p meta \
        >> "${PRODIGAL_DIR}/prodigal.log" 2>&1 || true
fi &
PID_PRODIGAL=$!

# Wait for skani (needed for PointFinder species detection)
wait $PID_SKANI 2>/dev/null || true

# ════════════════════════════════════════════════════════════════════════════
# Step 9: sRNA detection (Infernal/Rfam) — scan ARG flanks only
# ════════════════════════════════════════════════════════════════════════════
SRNA_DIR="${SAMPLE_DIR}/srna"
SRNA_TBLOUT="${SRNA_DIR}/cmscan_hits.tblout"
mkdir -p "$SRNA_DIR"

if [[ ! -f "$SRNA_TBLOUT" ]] && [[ -f "$RFAM_CM" ]] && [[ -f "$AMR_TSV" ]]; then
    SRNA_FLANKS="${SRNA_DIR}/arg_flanks.fasta"
    if [[ ! -f "$SRNA_FLANKS" ]]; then
        python3 - "$ASSEMBLY" "$AMR_TSV" "$SRNA_FLANKS" "$PROXIMITY_BP" <<'PYEOF' 2>/dev/null || true
import sys, csv
from Bio import SeqIO
from collections import defaultdict

assembly_path, amr_tsv, out_fa, flank_bp = sys.argv[1:5]
flank_bp = int(flank_bp)
contigs = {r.id: r for r in SeqIO.parse(assembly_path, "fasta")}

regions = []
with open(amr_tsv) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        contig = row.get("Contig id", "")
        try:
            start = int(row.get("Start", ""))
            end = int(row.get("Stop", ""))
        except (ValueError, TypeError):
            continue
        if contig not in contigs:
            continue
        clen = len(contigs[contig].seq)
        regions.append((contig, max(0, start - flank_bp), min(clen, end + flank_bp)))

by_contig = defaultdict(list)
for c, s, e in regions:
    by_contig[c].append((s, e))

with open(out_fa, "w") as fout:
    for contig, intervals in by_contig.items():
        intervals.sort()
        merged = [intervals[0]]
        for s, e in intervals[1:]:
            if s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        for s, e in merged:
            seq = str(contigs[contig].seq[s:e])
            fout.write(f">{contig}:{s}-{e}\n{seq}\n")
PYEOF
    fi
    if [[ -s "$SRNA_FLANKS" ]]; then
        $CMSCAN --rfam --cut_ga --noali --tblout "$SRNA_TBLOUT" --cpu "$THREADS" \
            "$RFAM_CM" "$SRNA_FLANKS" > /dev/null 2>> "${SRNA_DIR}/cmscan.log" || true
    fi
    rm -f "$SRNA_FLANKS"
fi

# Wait for Prodigal (needed by operon, CAI, synteny steps)
wait $PID_PRODIGAL 2>/dev/null || true

# ════════════════════════════════════════════════════════════════════════════
# Step 10: Operon structure
# ════════════════════════════════════════════════════════════════════════════
OPERON_DIR="${SAMPLE_DIR}/operon"
OPERON_TSV="${OPERON_DIR}/operon_results.tsv"
mkdir -p "$OPERON_DIR"

if [[ ! -f "$OPERON_TSV" ]] && [[ -f "$PRODIGAL_GFF" ]] && [[ -f "$AMR_TSV" ]]; then
    python3 - "$PRODIGAL_GFF" "$AMR_TSV" "$OPERON_TSV" <<'PYEOF' 2>/dev/null || true
import sys, csv, re
from collections import defaultdict

gff_path, amr_tsv, out_tsv = sys.argv[1:4]

genes = []
with open(gff_path) as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9 or parts[2] != 'CDS':
            continue
        genes.append((parts[0], int(parts[3]), int(parts[4]), parts[6]))

contig_genes = defaultdict(list)
for g in genes:
    contig_genes[g[0]].append(g)
for c in contig_genes:
    contig_genes[c].sort(key=lambda x: x[1])

MAX_GAP = 100
operons = []
gene_to_operon = {}

for contig, cgenes in contig_genes.items():
    if not cgenes:
        continue
    current_operon = [cgenes[0]]
    for i in range(1, len(cgenes)):
        prev = cgenes[i-1]
        curr = cgenes[i]
        gap = curr[1] - prev[2]
        if curr[3] == prev[3] and 0 <= gap <= MAX_GAP:
            current_operon.append(curr)
        else:
            operon_idx = len(operons)
            operons.append(current_operon)
            for pos, g in enumerate(current_operon):
                gene_to_operon[(g[0], g[1], g[2])] = (operon_idx, pos, len(current_operon))
            current_operon = [curr]
    operon_idx = len(operons)
    operons.append(current_operon)
    for pos, g in enumerate(current_operon):
        gene_to_operon[(g[0], g[1], g[2])] = (operon_idx, pos, len(current_operon))

results = []
with open(amr_tsv) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        gene = row.get('Gene symbol', row.get('Element symbol', 'unknown'))
        contig = row.get('Contig id', '')
        try:
            start = int(row.get('Start', ''))
            end = int(row.get('Stop', ''))
        except (ValueError, TypeError):
            results.append([gene, contig, '', '1', '1'])
            continue
        best_overlap = 0
        best_key = None
        for g in contig_genes.get(contig, []):
            ov = max(0, min(end, g[2]) - max(start, g[1]))
            if ov > best_overlap:
                best_overlap = ov
                best_key = (g[0], g[1], g[2])
        if best_key and best_key in gene_to_operon:
            _, pos, size = gene_to_operon[best_key]
            results.append([gene, contig, str(start), str(size), str(pos + 1)])
        else:
            results.append([gene, contig, str(start), '1', '1'])

with open(out_tsv, 'w') as f:
    f.write("gene\tcontig\tstart\toperon_size\toperon_position\n")
    for r in results:
        f.write('\t'.join(r) + '\n')
PYEOF
fi

# ════════════════════════════════════════════════════════════════════════════
# Step 11: Gene dosage
# ════════════════════════════════════════════════════════════════════════════
DOSAGE_DIR="${SAMPLE_DIR}/dosage"
DOSAGE_TSV="${DOSAGE_DIR}/gene_dosage.tsv"
mkdir -p "$DOSAGE_DIR"

if [[ ! -f "$DOSAGE_TSV" ]] && [[ -f "$AMR_TSV" ]]; then
    python3 - "$AMR_TSV" "$DOSAGE_TSV" <<'PYEOF' 2>/dev/null || true
import sys, csv
from collections import Counter

amr_tsv, out_tsv = sys.argv[1:3]
gene_counts = Counter()
rows = []
with open(amr_tsv) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        gene = row.get('Gene symbol', row.get('Element symbol', 'unknown'))
        gene_counts[gene] += 1
        contig = row.get('Contig id', '')
        start = row.get('Start', '')
        rows.append((gene, contig, start))

with open(out_tsv, 'w') as f:
    f.write("gene\tcontig\tstart\tgene_copies\n")
    for gene, contig, start in rows:
        f.write(f"{gene}\t{contig}\t{start}\t{gene_counts[gene]}\n")
PYEOF
fi

# ════════════════════════════════════════════════════════════════════════════
# Step 12: Codon adaptation index + rare codons
# ════════════════════════════════════════════════════════════════════════════
CAI_DIR="${SAMPLE_DIR}/cai"
CAI_TSV="${CAI_DIR}/cai_results.tsv"
mkdir -p "$CAI_DIR"

if [[ ! -f "$CAI_TSV" ]] && [[ -f "$PRODIGAL_FNA" ]] && [[ -f "$AMR_TSV" ]]; then
    python3 - "$ASSEMBLY" "$PRODIGAL_FNA" "$AMR_TSV" "$CAI_TSV" <<'PYEOF' 2>/dev/null || true
import sys, csv, math
from collections import Counter, defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

assembly_path, prodigal_fna, amr_tsv, out_tsv = sys.argv[1:5]

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

def count_codons(seq):
    seq = str(seq).upper().replace('U', 'T')
    counts = Counter()
    for i in range(0, len(seq) - 2, 3):
        c = seq[i:i+3]
        if len(c) == 3 and all(b in 'ACGT' for b in c):
            counts[c] += 1
    return counts

genome_codons = Counter()
for rec in SeqIO.parse(prodigal_fna, 'fasta'):
    genome_codons += count_codons(str(rec.seq))

aa_codons = defaultdict(list)
for codon, aa in CODON_TABLE.items():
    if aa != '*':
        aa_codons[aa].append(codon)

rscu = {}
for aa, codons in aa_codons.items():
    total = sum(genome_codons.get(c, 0) for c in codons)
    n = len(codons)
    for c in codons:
        rscu[c] = (genome_codons.get(c, 0) / total * n) if total > 0 else 1.0

w = {}
for aa, codons in aa_codons.items():
    max_rscu = max(rscu.get(c, 0) for c in codons)
    for c in codons:
        w[c] = rscu[c] / max_rscu if max_rscu > 0 else 1.0

def calc_cai(seq):
    codons = count_codons(seq)
    log_sum = 0.0
    count = 0
    for c, n in codons.items():
        aa = CODON_TABLE.get(c, '*')
        if aa in ('*', 'M', 'W'):
            continue
        wc = w.get(c, 0)
        if wc > 0:
            log_sum += n * math.log(wc)
            count += n
    return math.exp(log_sum / count) if count > 0 else None

# Rare codon threshold: bottom 10% of w values (excluding stops, M, W)
w_values_sorted = sorted([v for c, v in w.items()
                          if CODON_TABLE.get(c,'*') not in ('*','M','W') and v > 0])
rare_threshold = w_values_sorted[max(0, len(w_values_sorted) // 10)] if w_values_sorted else 0.1

def count_rare_clusters(seq, window=10, min_rare=3):
    """Count windows with >= min_rare rare codons in a sliding window."""
    seq = str(seq).upper().replace('U', 'T')
    codons_list = []
    for i in range(0, len(seq) - 2, 3):
        c = seq[i:i+3]
        if len(c) == 3 and all(b in 'ACGT' for b in c):
            aa = CODON_TABLE.get(c, '*')
            if aa not in ('*', 'M', 'W'):
                is_rare = 1 if w.get(c, 0) <= rare_threshold else 0
            else:
                is_rare = 0
            codons_list.append(is_rare)
        else:
            codons_list.append(0)
    clusters = 0
    for i in range(len(codons_list) - window + 1):
        if sum(codons_list[i:i+window]) >= min_rare:
            clusters += 1
    return clusters

def calc_rare_codon_pct(seq):
    """Percentage of codons that are rare."""
    seq = str(seq).upper().replace('U', 'T')
    total = 0
    rare = 0
    for i in range(0, len(seq) - 2, 3):
        c = seq[i:i+3]
        if len(c) == 3 and all(b in 'ACGT' for b in c):
            aa = CODON_TABLE.get(c, '*')
            if aa not in ('*', 'M', 'W'):
                total += 1
                if w.get(c, 0) <= rare_threshold:
                    rare += 1
    return rare / total if total > 0 else None

contigs = {r.id: str(r.seq) for r in SeqIO.parse(assembly_path, 'fasta')}

results = []
with open(amr_tsv) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        gene = row.get('Gene symbol', row.get('Element symbol', 'unknown'))
        contig = row.get('Contig id', '')
        strand = row.get('Strand', '+')
        try:
            start = int(row.get('Start', ''))
            end = int(row.get('Stop', ''))
        except (ValueError, TypeError):
            results.append([gene, contig, '', '', '', ''])
            continue
        seq = contigs.get(contig, '')
        if not seq:
            results.append([gene, contig, str(start), '', '', ''])
            continue
        gene_seq = seq[start-1:end]
        if strand == '-':
            gene_seq = str(Seq(gene_seq).reverse_complement())
        cai = calc_cai(gene_seq)
        rare_pct = calc_rare_codon_pct(gene_seq)
        rare_clusters = count_rare_clusters(gene_seq)
        results.append([gene, contig, str(start),
                        f"{cai:.4f}" if cai else '',
                        f"{rare_pct:.4f}" if rare_pct is not None else '',
                        str(rare_clusters)])

with open(out_tsv, 'w') as f:
    f.write("gene\tcontig\tstart\tcai\trare_codon_pct\trare_codon_clusters\n")
    for r in results:
        f.write('\t'.join(r) + '\n')
PYEOF
fi

# ════════════════════════════════════════════════════════════════════════════
# Step 13: GC content deviation
# ════════════════════════════════════════════════════════════════════════════
GC_DIR="${SAMPLE_DIR}/gc"
GC_TSV="${GC_DIR}/gc_deviation.tsv"
mkdir -p "$GC_DIR"

if [[ ! -f "$GC_TSV" ]] && [[ -f "$AMR_TSV" ]]; then
    python3 - "$ASSEMBLY" "$AMR_TSV" "$GC_TSV" <<'PYEOF' 2>/dev/null || true
import sys, csv
from Bio import SeqIO
from Bio.Seq import Seq

assembly_path, amr_tsv, out_tsv = sys.argv[1:4]

contigs = {r.id: str(r.seq) for r in SeqIO.parse(assembly_path, 'fasta')}

def gc_content(seq):
    seq = seq.upper()
    total = sum(1 for c in seq if c in 'ACGT')
    if total == 0:
        return None
    gc = sum(1 for c in seq if c in 'GC')
    return gc / total

genome_seq = ''.join(contigs.values())
genome_gc = gc_content(genome_seq)

results = []
with open(amr_tsv) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        gene = row.get('Gene symbol', row.get('Element symbol', 'unknown'))
        contig = row.get('Contig id', '')
        try:
            start = int(row.get('Start', ''))
            end = int(row.get('Stop', ''))
        except (ValueError, TypeError):
            results.append([gene, contig, '', '', '', ''])
            continue
        seq = contigs.get(contig, '')
        if not seq or start < 1:
            results.append([gene, contig, str(start), '', '', ''])
            continue
        gene_seq = seq[start-1:end]
        gene_gc = gc_content(gene_seq)
        if gene_gc is not None and genome_gc is not None:
            deviation = gene_gc - genome_gc
            results.append([gene, contig, str(start),
                            f"{gene_gc:.4f}", f"{genome_gc:.4f}", f"{deviation:.4f}"])
        else:
            results.append([gene, contig, str(start), '', '', ''])

with open(out_tsv, 'w') as f:
    f.write("gene\tcontig\tstart\tgene_gc\tgenome_gc\tgc_deviation\n")
    for r in results:
        f.write('\t'.join(r) + '\n')
PYEOF
fi

# ════════════════════════════════════════════════════════════════════════════
# Step 14: MLST
# ════════════════════════════════════════════════════════════════════════════
MLST_DIR="${SAMPLE_DIR}/mlst"
MLST_TSV="${MLST_DIR}/mlst_results.tsv"
mkdir -p "$MLST_DIR"

if [[ ! -f "$MLST_TSV" ]]; then
    conda run -n radar-mlst mlst "$ASSEMBLY" > "$MLST_TSV" 2>> "${MLST_DIR}/mlst.log" || true
fi

# ════════════════════════════════════════════════════════════════════════════
# Step 15: PointFinder (point mutations)
#   CRITICAL: skani Ref_name is column index 5 (not 7)
#   CRITICAL: Must pass -db_res to resfinder
#   CRITICAL: Campylobacter maps to 'campylobacter' (not 'campylobacter_jejuni')
# ════════════════════════════════════════════════════════════════════════════
POINT_DIR="${SAMPLE_DIR}/pointfinder"
POINT_TSV="${POINT_DIR}/PointFinder_results.txt"
mkdir -p "$POINT_DIR"

if [[ ! -f "$POINT_TSV" ]] && [[ -d "$POINTFINDER_DB" ]]; then
    PF_SPECIES=""
    if [[ -f "$SKANI_TSV" ]]; then
        PF_SPECIES=$(python3 - "$SKANI_TSV" <<'PYEOF' 2>/dev/null || true
import sys, csv
skani_tsv = sys.argv[1]
species_map = {
    'escherichia': 'escherichia_coli',
    'salmonella': 'salmonella',
    'campylobacter jejuni': 'campylobacter',
    'campylobacter coli': 'campylobacter',
    'staphylococcus aureus': 'staphylococcus_aureus',
    'mycobacterium tuberculosis': 'mycobacterium_tuberculosis',
    'enterococcus faecalis': 'enterococcus_faecalis',
    'enterococcus faecium': 'enterococcus_faecium',
    'klebsiella': 'klebsiella',
    'neisseria gonorrhoeae': 'neisseria_gonorrhoeae',
    'helicobacter pylori': 'helicobacter_pylori',
    'plasmodium falciparum': 'plasmodium_falciparum',
}
try:
    with open(skani_tsv) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) < 6:
                continue
            ref_name = row[5].lower()
            for key, val in species_map.items():
                if key in ref_name:
                    print(val)
                    sys.exit(0)
except Exception:
    pass
PYEOF
        )
    fi
    if [[ -n "$PF_SPECIES" ]]; then
        python3 -m resfinder -ifa "$ASSEMBLY" -o "$POINT_DIR" \
            -c -db_point "$POINTFINDER_DB" -s "$PF_SPECIES" \
            -db_res "$RESFINDER_DB" \
            --ignore_missing_species \
            >> "${POINT_DIR}/pointfinder.log" 2>&1 || true
    fi
fi

# ════════════════════════════════════════════════════════════════════════════
# Step 16: IntegronFinder
# ════════════════════════════════════════════════════════════════════════════
INTEGRON_DIR="${SAMPLE_DIR}/integron"
INTEGRON_SUMMARY="${INTEGRON_DIR}/Results_Integron_Finder_assembly/assembly.summary"
mkdir -p "$INTEGRON_DIR"

if [[ ! -f "$INTEGRON_SUMMARY" ]]; then
    $INTEGRON_FINDER --local-max --cpu "$THREADS" \
        --outdir "${INTEGRON_DIR}/Results_Integron_Finder_assembly" \
        "$ASSEMBLY" >> "${INTEGRON_DIR}/integron_finder.log" 2>&1 || true
fi

# ════════════════════════════════════════════════════════════════════════════
# Wait for geNomad (needed by synteny and summary for plasmid/IS/transposase)
wait $PID_GENOMAD 2>/dev/null || true

# Step 17: Gene synteny
#   CRITICAL: Use geNomad gene annotations for IS/transposase detection (not MEFinder)
#   geNomad genes at: genomad/${ASSEMBLY_NAME}_annotate/${ASSEMBLY_NAME}_genes.tsv
#   Fall back to MEFinder CSV only if geNomad unavailable
# ════════════════════════════════════════════════════════════════════════════
SYNTENY_DIR="${SAMPLE_DIR}/synteny"
SYNTENY_TSV="${SYNTENY_DIR}/synteny_results.tsv"
mkdir -p "$SYNTENY_DIR"

GENOMAD_GENES="${GENOMAD_DIR}/${ASSEMBLY_NAME}_annotate/${ASSEMBLY_NAME}_genes.tsv"

if [[ ! -f "$SYNTENY_TSV" ]] && [[ -f "$PRODIGAL_GFF" ]] && [[ -f "$AMR_TSV" ]]; then
    python3 - "$PRODIGAL_GFF" "$AMR_TSV" "$MEF_CSV" "$SYNTENY_TSV" "$GENOMAD_GENES" <<'PYEOF' 2>/dev/null || true
import sys, csv, re, os
from collections import defaultdict

gff_path, amr_tsv, mef_csv, out_tsv = sys.argv[1:5]
genomad_genes_path = sys.argv[5] if len(sys.argv) > 5 else ""
FLANK_N = 5

# Parse prodigal GFF for CDS positions
genes = []
with open(gff_path) as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9 or parts[2] != 'CDS':
            continue
        genes.append((parts[0], int(parts[3]), int(parts[4]), parts[6]))

contig_genes = defaultdict(list)
for g in genes:
    contig_genes[g[0]].append(g)
for c in contig_genes:
    contig_genes[c].sort(key=lambda x: x[1])

# Build annotation lookup from AMRFinder — classify by coordinate overlap
amr_regions = []
with open(amr_tsv) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        contig = row.get('Contig id', '')
        etype = row.get('Element type', row.get('Type', ''))
        gene_name = row.get('Gene symbol', row.get('Element symbol', ''))
        try:
            s = int(row.get('Start', ''))
            e = int(row.get('Stop', ''))
        except (ValueError, TypeError):
            continue
        if etype == 'AMR':
            cat = 'AMR'
        elif etype == 'STRESS':
            cat = 'stress'
        elif etype == 'VIRULENCE':
            cat = 'virulence'
        else:
            cat = 'AMR'
        amr_regions.append((contig, s, e, cat, gene_name))

# IS elements / transposases from geNomad gene annotations (PRIMARY)
is_regions = []
if genomad_genes_path and os.path.exists(genomad_genes_path):
    with open(genomad_genes_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            desc = row.get("annotation_description", "").lower()
            if "transposase" in desc or "insertion element" in desc or "insertion sequence" in desc:
                gene_id = row.get("gene", "")
                cid = gene_id.rsplit("_", 1)[0] if "_" in gene_id else gene_id
                try:
                    is_regions.append((cid, int(row.get("start", 0)),
                                       int(row.get("end", 0)), desc[:50]))
                except (ValueError, TypeError):
                    pass
# Fallback to MEFinder only if geNomad unavailable
if not is_regions and os.path.exists(mef_csv):
    try:
        with open(mef_csv) as f:
            for row in csv.DictReader(f):
                is_regions.append((row.get('contig',''), int(row.get('start',0)),
                                   int(row.get('end',0)), row.get('name', '')))
    except Exception:
        pass

def classify_gene(contig, gstart, gend):
    """Classify a prodigal CDS by overlap with known annotations."""
    # Check AMRFinder hits
    for ac, as_, ae, acat, aname in amr_regions:
        if ac == contig and max(gstart, as_) < min(gend, ae):
            return acat
    # Check IS elements
    for ic, is_, ie, iname in is_regions:
        if ic == contig and max(gstart, is_) < min(gend, ie):
            return 'transposase'
    # Short genes near mobile elements are often hypothetical
    gene_len = gend - gstart
    if gene_len < 200:
        return 'hypothetical'
    return 'other'

# For each ARG, find its position in prodigal genes and classify neighbors
results = []
with open(amr_tsv) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        gene = row.get('Gene symbol', row.get('Element symbol', 'unknown'))
        contig = row.get('Contig id', '')
        try:
            start = int(row.get('Start', ''))
            end = int(row.get('Stop', ''))
        except (ValueError, TypeError):
            results.append([gene, contig, '', '', '', '', '', '0', '0', '0', '0'])
            continue

        cgenes = contig_genes.get(contig, [])
        best_idx = -1
        best_overlap = 0
        for i, g in enumerate(cgenes):
            ov = max(0, min(end, g[2]) - max(start, g[1]))
            if ov > best_overlap:
                best_overlap = ov
                best_idx = i

        if best_idx < 0:
            results.append([gene, contig, str(start), '', '', '', '', '0', '0', '0', '0'])
            continue

        upstream = []
        downstream = []
        for i in range(best_idx - 1, max(best_idx - FLANK_N - 1, -1), -1):
            g = cgenes[i]
            upstream.append(classify_gene(contig, g[1], g[2]))
        upstream.reverse()

        for i in range(best_idx + 1, min(best_idx + FLANK_N + 1, len(cgenes))):
            g = cgenes[i]
            downstream.append(classify_gene(contig, g[1], g[2]))

        all_neighbors = upstream + downstream
        n_transposase = sum(1 for c in all_neighbors if c == 'transposase')
        n_amr = sum(1 for c in all_neighbors if c == 'AMR')
        n_stress = sum(1 for c in all_neighbors if c == 'stress')
        n_hypothetical = sum(1 for c in all_neighbors if c == 'hypothetical')

        results.append([
            gene, contig, str(start),
            ';'.join(upstream + ['[ARG]'] + downstream),
            str(len(all_neighbors)),
            str(n_transposase), str(n_amr), str(n_stress),
            str(n_hypothetical),
            str(sum(1 for c in all_neighbors if c == 'other')),
            str(sum(1 for c in all_neighbors if c == 'virulence')),
        ])

with open(out_tsv, 'w') as f:
    f.write("gene\tcontig\tstart\tsynteny_pattern\tn_neighbors\tn_transposase\tn_neighbor_amr\tn_stress\tn_hypothetical\tn_other\tn_virulence\n")
    for r in results:
        f.write('\t'.join(r) + '\n')
PYEOF
fi

# ════════════════════════════════════════════════════════════════════════════
# Cross-reference summary: join all 17 annotation results into one TSV
#   CRITICAL: Plasmid detection uses geNomad aggregated classification (PRIMARY),
#             falls back to MOB-recon only if geNomad unavailable
#   CRITICAL: IS elements use geNomad gene annotations (PRIMARY),
#             falls back to MEFinder only if geNomad unavailable
# ════════════════════════════════════════════════════════════════════════════
SUMMARY_TSV="${SAMPLE_DIR}/arg_context_summary.tsv"
if [[ -f "$AMR_TSV" ]] && [[ ! -f "$SUMMARY_TSV" ]]; then
    python3 - "$AMR_TSV" "$MEF_CSV" "$MOB_REPORT" "$VIRUS_SUMMARY" \
               "$PROM_TSV" "$RBS_TSV" "$SRNA_TBLOUT" "$OPERON_TSV" \
               "$DOSAGE_TSV" "$CAI_TSV" "$GC_TSV" "$MLST_TSV" \
               "$POINT_DIR" "$INTEGRON_DIR" "$SYNTENY_TSV" \
               "$SUMMARY_TSV" "$PROXIMITY_BP" <<'PYEOF' 2>/dev/null || true
import sys, csv, os, re
from collections import defaultdict

(amr_tsv, mef_csv, mob_report, virus_summary,
 prom_tsv, rbs_tsv, srna_tblout, operon_tsv,
 dosage_tsv, cai_tsv, gc_tsv, mlst_tsv,
 point_dir, integron_dir, synteny_tsv,
 out_tsv, prox_bp) = sys.argv[1:18]
prox_bp = int(prox_bp)

# ── Plasmid contigs (geNomad aggregated classification — PRIMARY) ──
plasmid_contigs = set()
genomad_class = os.path.join(os.path.dirname(out_tsv), "genomad",
    "assembly_aggregated_classification", "assembly_aggregated_classification.tsv")
if os.path.exists(genomad_class):
    with open(genomad_class) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            try:
                if float(row.get("plasmid_score", 0)) > 0.7:
                    plasmid_contigs.add(row.get("seq_name", "").split()[0])
            except (ValueError, TypeError):
                pass
# Fallback to MOB-recon if geNomad not available
if not plasmid_contigs and os.path.exists(mob_report):
    with open(mob_report) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("molecule_type") == "plasmid":
                cid = row.get("contig_id", "").split()[0]
                if cid:
                    plasmid_contigs.add(cid)

# ── Prophage regions (from geNomad virus summary) ──
prophage_regions = []
if os.path.exists(virus_summary):
    with open(virus_summary) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            seq_name = row.get("seq_name", "")
            contig = seq_name.split("|")[0] if "|" in seq_name else seq_name
            coords = row.get("coordinates", "")
            if "-" in coords:
                parts = coords.split("-")
                try:
                    prophage_regions.append((contig, int(parts[0]), int(parts[1])))
                except ValueError:
                    pass

# ── IS elements / transposases (geNomad gene annotations — PRIMARY) ──
# Stores (contig, start, end, name, strand) — strand as 1 or -1
is_elements = []
genomad_genes = os.path.join(os.path.dirname(out_tsv), "genomad",
    "assembly_annotate", "assembly_genes.tsv")
if os.path.exists(genomad_genes):
    with open(genomad_genes) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            desc = row.get("annotation_description", "").lower()
            if "transposase" in desc or "insertion element" in desc or "insertion sequence" in desc:
                gene_id = row.get("gene", "")
                contig = gene_id.rsplit("_", 1)[0] if "_" in gene_id else gene_id
                try:
                    is_strand = int(row.get("strand", 0))
                    is_elements.append((contig, int(row.get("start", 0)),
                                        int(row.get("end", 0)), desc[:50], is_strand))
                except (ValueError, TypeError):
                    pass
# Fallback to MEFinder if geNomad not available
if not is_elements and os.path.exists(mef_csv):
    with open(mef_csv) as f:
        try:
            for row in csv.DictReader(f):
                is_elements.append((row.get("contig",""), int(row.get("start",0)),
                                    int(row.get("end",0)), row.get("name", row.get("type", ""))))
        except Exception:
            pass

# ── Promoters ──
promoters = {}
if os.path.exists(prom_tsv):
    with open(prom_tsv) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            promoters[(row.get("gene",""), row.get("contig",""), row.get("start",""))] = row

# ── RBS ──
rbs_data = {}
if os.path.exists(rbs_tsv):
    with open(rbs_tsv) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rbs_data[(row.get("gene",""), row.get("contig",""), row.get("start",""))] = row

# ── sRNA hits from cmscan tblout ──
srna_hits = []
if os.path.exists(srna_tblout):
    with open(srna_tblout) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) < 18:
                continue
            target_name = fields[0]
            target_acc = fields[1]
            query_name = fields[2]
            seq_from = int(fields[7])
            seq_to = int(fields[8])
            strand = fields[9]
            score = float(fields[14])
            evalue = fields[15]
            clan = fields[17] if len(fields) > 17 else ""
            contig = query_name
            offset = 0
            if ':' in query_name:
                parts = query_name.rsplit(':', 1)
                contig = parts[0]
                range_parts = parts[1].split('-')
                if len(range_parts) == 2:
                    try:
                        offset = int(range_parts[0])
                    except ValueError:
                        pass
            s, e = min(seq_from, seq_to) + offset, max(seq_from, seq_to) + offset
            srna_hits.append((contig, s, e, target_name, target_acc, score, clan))

# ── Operon data ──
operon_data = {}
if os.path.exists(operon_tsv):
    with open(operon_tsv) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            operon_data[(row.get("gene",""), row.get("contig",""), row.get("start",""))] = row

# ── Gene dosage ──
dosage_data = {}
if os.path.exists(dosage_tsv):
    with open(dosage_tsv) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            dosage_data[(row.get("gene",""), row.get("contig",""), row.get("start",""))] = row

# ── CAI ──
cai_data = {}
if os.path.exists(cai_tsv):
    with open(cai_tsv) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            cai_data[(row.get("gene",""), row.get("contig",""), row.get("start",""))] = row

# ── GC deviation ──
gc_data = {}
if os.path.exists(gc_tsv):
    with open(gc_tsv) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            gc_data[(row.get("gene",""), row.get("contig",""), row.get("start",""))] = row

# ── MLST ──
mlst_scheme = ""
mlst_st = ""
if os.path.exists(mlst_tsv):
    with open(mlst_tsv) as f:
        line = f.readline().strip()
        parts = line.split('\t')
        if len(parts) >= 3:
            mlst_scheme = parts[1]
            mlst_st = parts[2]

# ── Point mutations ──
point_mutations = {}
point_results = os.path.join(point_dir, "PointFinder_results.txt")
if os.path.exists(point_results):
    with open(point_results) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            mutation = row.get("Mutation", row.get("mutation", ""))
            resistance = row.get("Resistance", row.get("resistance", ""))
            gene = row.get("Gene_name", row.get("gene_name", ""))
            if gene:
                if gene not in point_mutations:
                    point_mutations[gene] = []
                point_mutations[gene].append(f"{mutation}({resistance})")

# ── Integrons ──
integron_contigs = {}
for root, dirs, files in os.walk(integron_dir):
    for fname in files:
        if fname.endswith('.integrons'):
            fpath = os.path.join(root, fname)
            with open(fpath) as f:
                for line in f:
                    if line.startswith('#') or line.startswith('ID_integron'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 11:
                        contig = parts[1]
                        itype = parts[10]
                        try:
                            istart = int(parts[3])
                            iend = int(parts[4])
                        except (ValueError, IndexError):
                            continue
                        if contig not in integron_contigs:
                            integron_contigs[contig] = []
                        integron_contigs[contig].append((istart, iend, itype))

# ── Synteny ──
synteny_data = {}
if os.path.exists(synteny_tsv):
    with open(synteny_tsv) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            synteny_data[(row.get("gene",""), row.get("contig",""), row.get("start",""))] = row

# ── Build output ──
rows_out = []
with open(amr_tsv) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        etype = row.get("Element type", row.get("Type", ""))
        gene = row.get("Gene symbol", row.get("Element symbol", "unknown"))
        contig = row.get("Contig id", "")
        drug_class = row.get("Class", "")
        subclass = row.get("Subclass", "")
        if subclass and subclass != drug_class:
            drug_class = f"{drug_class}; {subclass}" if drug_class else subclass
        identity = row.get("% Identity to reference sequence", row.get("% Identity to reference", ""))
        coverage = row.get("% Coverage of reference sequence", row.get("% Coverage of reference", ""))
        method = row.get("Method", "")
        try:
            start = int(row.get("Start", "")); end = int(row.get("Stop", ""))
        except (ValueError, TypeError):
            start = end = None

        on_plasmid = contig in plasmid_contigs
        on_prophage = False
        if start is not None:
            for pc, ps, pe in prophage_regions:
                if contig == pc and ps <= start <= pe:
                    on_prophage = True; break

        nearest_is = ""; nearest_is_dist = ""; nearest_is_orient = ""
        # ARG strand: convert +/- to 1/-1
        arg_strand_str = row.get("Strand", "+")
        arg_strand = 1 if arg_strand_str == "+" else -1
        if start is not None:
            for ic, is_s, is_e, is_name, is_strand in is_elements:
                if ic != contig:
                    continue
                dist = 0 if (start <= is_e and end >= is_s) else min(abs(start - is_e), abs(is_s - end))
                if dist <= prox_bp and (not nearest_is_dist or dist < int(nearest_is_dist)):
                    nearest_is = is_name; nearest_is_dist = str(dist)
                    # Relative orientation: same=1, opposite=-1, unknown=0
                    if is_strand != 0:
                        nearest_is_orient = "same" if (arg_strand == is_strand) else "opposite"
                    else:
                        nearest_is_orient = ""

        contig_type = "prophage" if on_prophage else ("plasmid" if on_plasmid else "chromosome")

        # sRNA: nearest Rfam hit
        nearest_srna = ""; nearest_srna_dist = ""; nearest_srna_acc = ""
        if start is not None:
            for sc, ss, se, sname, sacc, sscore, sclan in srna_hits:
                if sc != contig:
                    continue
                dist = 0 if (start <= se and end >= ss) else min(abs(start - se), abs(ss - end))
                if dist <= prox_bp and (not nearest_srna_dist or dist < int(nearest_srna_dist)):
                    nearest_srna = sname
                    nearest_srna_acc = sacc
                    nearest_srna_dist = str(dist)

        # Integron proximity
        in_integron = "False"
        nearest_integron_dist = ""
        nearest_integron_type = ""
        if start is not None and contig in integron_contigs:
            for istart, iend, itype in integron_contigs[contig]:
                if start <= iend and end >= istart:
                    in_integron = "True"
                    nearest_integron_dist = "0"
                    nearest_integron_type = itype
                    break
                dist = min(abs(start - iend), abs(istart - end))
                if dist <= prox_bp and (not nearest_integron_dist or dist < int(nearest_integron_dist)):
                    nearest_integron_dist = str(dist)
                    nearest_integron_type = itype

        pkey = (gene, contig, str(start))
        p = promoters.get(pkey, {})
        r = rbs_data.get(pkey, {})
        op = operon_data.get(pkey, {})
        dos = dosage_data.get(pkey, {})
        ca = cai_data.get(pkey, {})
        gc = gc_data.get(pkey, {})
        sy = synteny_data.get(pkey, {})

        # Point mutations for this gene
        pmut = "; ".join(point_mutations.get(gene, []))

        rows_out.append([
            gene, etype, drug_class, contig, str(start or ""), str(end or ""),
            identity, coverage, method,
            contig_type, str(on_plasmid), str(on_prophage),
            nearest_is, nearest_is_dist, nearest_is_orient,
            p.get("ldf_score",""), p.get("tf_binding_sites",""),
            p.get("promoter_distance",""), p.get("up_element_at_ratio",""),
            r.get("expression",""), r.get("dg_total",""), r.get("dg_mrna",""),
            nearest_srna, nearest_srna_acc, nearest_srna_dist,
            op.get("operon_size",""), op.get("operon_position",""),
            dos.get("gene_copies",""),
            ca.get("cai",""), ca.get("rare_codon_pct",""), ca.get("rare_codon_clusters",""),
            gc.get("gene_gc",""), gc.get("genome_gc",""), gc.get("gc_deviation",""),
            mlst_scheme, mlst_st,
            pmut,
            in_integron, nearest_integron_dist, nearest_integron_type,
            sy.get("synteny_pattern",""),
            sy.get("n_transposase",""), sy.get("n_neighbor_amr",""),
            sy.get("n_stress",""), sy.get("n_hypothetical",""),
            sy.get("n_virulence",""),
        ])

with open(out_tsv, "w") as f:
    f.write("\t".join([
        "gene","element_type","drug_class","contig","start","end",
        "identity","coverage","method",
        "contig_type","on_plasmid","on_prophage",
        "nearest_IS","nearest_IS_distance_bp","nearest_IS_orientation",
        "promoter_ldf","promoter_tf_sites",
        "promoter_distance","promoter_up_at_ratio",
        "rbs_expression","rbs_dg_total","rbs_dg_mrna",
        "nearest_sRNA","nearest_sRNA_accession","nearest_sRNA_distance_bp",
        "operon_size","operon_position",
        "gene_copies",
        "cai","rare_codon_pct","rare_codon_clusters",
        "gene_gc","genome_gc","gc_deviation",
        "mlst_scheme","mlst_st",
        "point_mutations",
        "in_integron","nearest_integron_distance_bp","nearest_integron_type",
        "synteny_pattern",
        "synteny_n_transposase","synteny_n_amr",
        "synteny_n_stress","synteny_n_hypothetical","synteny_n_virulence",
    ]) + "\n")
    for r in rows_out:
        f.write("\t".join(r) + "\n")
PYEOF
fi

echo "[DONE] $SAMPLE_ID"
