#!/usr/bin/env python3
"""
Fast regeneration of synteny_results.tsv and arg_context_summary.tsv
using geNomad data for plasmid/IS/transposase detection.
Only reads existing output files — does NOT re-run any tools.
"""
import os, sys, csv, re, glob
from collections import defaultdict

RESULTS_DIR = "/home/unnot/WGS/phenotype/results"


def load_genomad_plasmid_contigs(sample_dir):
    """Get plasmid contig IDs from geNomad aggregated classification."""
    path = os.path.join(sample_dir, "genomad", "assembly_aggregated_classification",
                        "assembly_aggregated_classification.tsv")
    plasmid_contigs = set()
    if os.path.exists(path):
        with open(path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                try:
                    if float(row.get("plasmid_score", 0)) > 0.7:
                        plasmid_contigs.add(row.get("seq_name", "").split()[0])
                except (ValueError, TypeError):
                    pass
    return plasmid_contigs


def load_genomad_transposases(sample_dir):
    """Get IS/transposase locations from geNomad gene annotations."""
    # Try both possible paths
    for assembly_name in ["assembly"]:
        path = os.path.join(sample_dir, "genomad",
                            f"{assembly_name}_annotate", f"{assembly_name}_genes.tsv")
        if os.path.exists(path):
            break
    else:
        return []

    is_elements = []
    with open(path) as f:
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
    return is_elements


def load_prodigal_genes(gff_path):
    """Parse Prodigal GFF to get gene positions."""
    genes = []
    if not os.path.exists(gff_path):
        return genes
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or "\tCDS\t" not in line:
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            contig = parts[0]
            try:
                start, end = int(parts[3]), int(parts[4])
            except ValueError:
                continue
            strand = parts[6]
            gene_id = ""
            for attr in parts[8].split(";"):
                if attr.startswith("ID="):
                    gene_id = attr[3:]
                    break
            genes.append({"contig": contig, "start": start, "end": end,
                          "strand": strand, "gene_id": gene_id})
    return genes


def load_amr_results(amr_tsv):
    """Load AMRFinderPlus results."""
    results = []
    if not os.path.exists(amr_tsv):
        return results
    with open(amr_tsv) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            results.append(row)
    return results


def regen_synteny(sample_dir, is_elements):
    """Regenerate synteny_results.tsv using geNomad transposase data."""
    gff_path = os.path.join(sample_dir, "prodigal", "prodigal.gff")
    amr_tsv = os.path.join(sample_dir, "amr", "amrfinderplus.tsv")
    out_tsv = os.path.join(sample_dir, "synteny", "synteny_results.tsv")

    genes = load_prodigal_genes(gff_path)
    amr_results = load_amr_results(amr_tsv)

    if not genes or not amr_results:
        return

    # Build gene index by contig
    contig_genes = defaultdict(list)
    for g in genes:
        contig_genes[g["contig"]].append(g)
    for ctg in contig_genes:
        contig_genes[ctg].sort(key=lambda x: x["start"])

    # Build IS element index by contig
    is_by_contig = defaultdict(list)
    for (ctg, s, e, name) in is_elements:
        is_by_contig[ctg].append((s, e, name))

    # AMR gene positions
    amr_positions = set()
    for row in amr_results:
        contig = row.get("Contig id", "").split()[0]
        try:
            amr_positions.add((contig, int(row.get("Start", 0)), int(row.get("Stop", 0))))
        except (ValueError, TypeError):
            pass

    def is_amr(contig, start, end):
        for (c, s, e) in amr_positions:
            if c == contig and max(start, s) <= min(end, e):
                return True
        return False

    def is_transposase(contig, start, end):
        for (s, e, _) in is_by_contig.get(contig, []):
            if max(start, s) <= min(end, e):
                return True
        return False

    WINDOW = 5  # genes upstream/downstream

    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    with open(out_tsv, "w") as fout:
        fout.write("gene\tcontig\tstart\tend\tstrand\tn_transposase\tn_neighbor_amr\n")
        for row in amr_results:
            contig = row.get("Contig id", "").split()[0]
            try:
                amr_start = int(row.get("Start", 0))
                amr_end = int(row.get("Stop", 0))
            except (ValueError, TypeError):
                continue

            gene_name = row.get("Gene symbol", "unknown")
            strand = row.get("Strand", "+")

            # Find this gene's index in contig
            cgenes = contig_genes.get(contig, [])
            idx = -1
            for i, g in enumerate(cgenes):
                if max(g["start"], amr_start) <= min(g["end"], amr_end):
                    idx = i
                    break

            n_transposase = 0
            n_amr = 0
            if idx >= 0:
                lo = max(0, idx - WINDOW)
                hi = min(len(cgenes), idx + WINDOW + 1)
                for j in range(lo, hi):
                    if j == idx:
                        continue
                    ng = cgenes[j]
                    if is_transposase(contig, ng["start"], ng["end"]):
                        n_transposase += 1
                    if is_amr(contig, ng["start"], ng["end"]):
                        n_amr += 1

            fout.write(f"{gene_name}\t{contig}\t{amr_start}\t{amr_end}\t{strand}\t"
                       f"{n_transposase}\t{n_amr}\n")


def regen_summary(sample_dir, plasmid_contigs, is_elements, integron_contigs=None):
    """Regenerate arg_context_summary.tsv with fixed plasmid/IS data."""
    amr_tsv = os.path.join(sample_dir, "amr", "amrfinderplus.tsv")
    synteny_tsv = os.path.join(sample_dir, "synteny", "synteny_results.tsv")
    out_tsv = os.path.join(sample_dir, "arg_context_summary.tsv")

    # Read the existing summary to preserve all other columns
    if not os.path.exists(out_tsv):
        return False

    rows = []
    with open(out_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames
        for row in reader:
            rows.append(dict(row))

    if not rows or not fieldnames:
        return False

    # Build IS element index by contig for distance calculation
    is_by_contig = defaultdict(list)
    for elem in is_elements:
        ctg, s, e = elem[0], elem[1], elem[2]
        is_strand = elem[4] if len(elem) > 4 else 0
        is_by_contig[ctg].append((s, e, is_strand))

    # Load synteny data
    synteny_data = {}
    if os.path.exists(synteny_tsv):
        with open(synteny_tsv) as f:
            for srow in csv.DictReader(f, delimiter="\t"):
                key = (srow.get("contig", ""), srow.get("start", ""), srow.get("end", ""))
                synteny_data[key] = srow

    # Update rows
    for row in rows:
        cid = row.get("contig", "").strip()
        contig = cid.split()[0] if cid else ""

        # Fix plasmid detection
        row["on_plasmid"] = "True" if contig in plasmid_contigs else "False"
        row["contig_type"] = "plasmid" if contig in plasmid_contigs else "chromosome"

        # Fix IS distance
        try:
            gene_start = int(row.get("start", 0))
            gene_end = int(row.get("end", 0))
        except (ValueError, TypeError):
            gene_start = gene_end = 0

        min_is_dist = None
        nearest_is_orient = ""
        # ARG strand from summary (if available)
        arg_strand_raw = row.get("strand", "")
        arg_strand = 1 if arg_strand_raw == "+" else (-1 if arg_strand_raw == "-" else 0)
        for (is_s, is_e, is_strand) in is_by_contig.get(contig, []):
            if gene_start and is_s:
                dist = max(0, max(gene_start, is_s) - min(gene_end, is_e))
                if min_is_dist is None or dist < min_is_dist:
                    min_is_dist = dist
                    if is_strand != 0 and arg_strand != 0:
                        nearest_is_orient = "same" if arg_strand == is_strand else "opposite"
                    else:
                        nearest_is_orient = ""
        row["nearest_IS_distance_bp"] = str(min_is_dist) if min_is_dist is not None else ""
        if "nearest_IS_orientation" in (fieldnames or []):
            row["nearest_IS_orientation"] = nearest_is_orient

        # Fix integron proximity
        if integron_contigs and contig in integron_contigs and gene_start:
            in_integron = False
            nearest_int_dist = None
            nearest_int_type = ""
            for (istart, iend, itype) in integron_contigs[contig]:
                if gene_start <= iend and gene_end >= istart:
                    in_integron = True
                    nearest_int_dist = 0
                    nearest_int_type = itype
                    break
                dist = min(abs(gene_start - iend), abs(istart - gene_end))
                if dist <= 5000 and (nearest_int_dist is None or dist < nearest_int_dist):
                    nearest_int_dist = dist
                    nearest_int_type = itype
            row["in_integron"] = "True" if in_integron else row.get("in_integron", "False")
            if nearest_int_dist is not None:
                row["nearest_integron_distance_bp"] = str(nearest_int_dist)
                row["nearest_integron_type"] = nearest_int_type

        # Fix synteny transposase
        key = (contig, str(gene_start), str(gene_end))
        if key in synteny_data:
            row["synteny_n_transposase"] = synteny_data[key].get("n_transposase", "0")
            row["synteny_n_amr_neighbor"] = synteny_data[key].get("n_neighbor_amr", "0")

    # Write back
    with open(out_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    return True


def load_integrons(sample_dir):
    """Parse IntegronFinder results into per-contig integron regions."""
    integron_contigs = {}
    integron_dir = os.path.join(sample_dir, "integron")
    if not os.path.exists(integron_dir):
        return integron_contigs
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
    return integron_contigs


def process_sample(sample_dir):
    """Process a single sample."""
    sample_id = os.path.basename(sample_dir)

    # Skip if no AMR results
    if not os.path.exists(os.path.join(sample_dir, "amr", "amrfinderplus.tsv")):
        return

    # Skip if no summary exists
    if not os.path.exists(os.path.join(sample_dir, "arg_context_summary.tsv")):
        return

    # Load geNomad data
    plasmid_contigs = load_genomad_plasmid_contigs(sample_dir)
    is_elements = load_genomad_transposases(sample_dir)
    integron_contigs = load_integrons(sample_dir)

    # Regenerate synteny with geNomad transposase data
    regen_synteny(sample_dir, is_elements)

    # Update summary with fixed plasmid/IS/transposase/integron data
    regen_summary(sample_dir, plasmid_contigs, is_elements, integron_contigs)


if __name__ == "__main__":
    import multiprocessing as mp

    sample_dirs = sorted(glob.glob(os.path.join(RESULTS_DIR, "SRR*")))
    # Filter to those with summaries
    sample_dirs = [d for d in sample_dirs
                   if os.path.exists(os.path.join(d, "arg_context_summary.tsv"))]

    print(f"Processing {len(sample_dirs)} samples...")

    with mp.Pool(20) as pool:
        for i, _ in enumerate(pool.imap_unordered(process_sample, sample_dirs), 1):
            if i % 500 == 0:
                print(f"  {i}/{len(sample_dirs)} done")

    print("Done.")
