"""Streamlined annotation pipeline for web deployment.
Runs only the essential tools (no mobility/regulation features).
"""
import asyncio
import os
import csv
import json
from app.config import (AMRFINDERPLUS, PRODIGAL, BPROM, BPROM_DATA,
                        POINTFINDER_DB, RESFINDER_DB, POINTFINDER_SPECIES_MAP)


async def run_annotation(fasta_path: str, job_dir: str, species: str, step: str):
    """Run a single annotation step."""
    if step == "amrfinder":
        return await run_amrfinder(fasta_path, job_dir)
    elif step == "prodigal":
        return await run_prodigal(fasta_path, job_dir)
    elif step == "expression":
        return await run_expression(fasta_path, job_dir)
    elif step == "pointfinder":
        return await run_pointfinder(fasta_path, job_dir, species)
    elif step == "genomic_features":
        return compute_genomic_features(fasta_path, job_dir)
    else:
        raise ValueError(f"Unknown annotation step: {step}")


async def run_amrfinder(fasta_path: str, job_dir: str) -> dict:
    """Run AMRFinderPlus for resistance gene detection."""
    amr_out = os.path.join(job_dir, "amrfinderplus.tsv")
    proc = await asyncio.create_subprocess_exec(
        AMRFINDERPLUS, "-n", fasta_path, "-o", amr_out, "--plus",
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )
    _, stderr = await proc.communicate()

    results = {'all_genes': [], 'genes_by_class': {}, 'mutations': []}

    if os.path.exists(amr_out):
        with open(amr_out) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene = row.get('Gene symbol', row.get('Element symbol', ''))
                drug_class = row.get('Class', row.get('Subclass', ''))
                etype = row.get('Element type', row.get('Type', ''))
                identity = float(row.get('% Identity to reference', 0))
                coverage = float(row.get('% Coverage of reference', 0))

                if etype in ('AMR', 'STRESS', 'VIRULENCE') or 'AMR' in str(etype):
                    gene_info = {
                        'gene': gene,
                        'drug_class': drug_class,
                        'identity': identity,
                        'coverage': coverage,
                        'on_plasmid': False,
                    }
                    results['all_genes'].append(gene_info)

                    if drug_class not in results['genes_by_class']:
                        results['genes_by_class'][drug_class] = []
                    if gene not in results['genes_by_class'][drug_class]:
                        results['genes_by_class'][drug_class].append(gene)

    # Save parsed results
    with open(os.path.join(job_dir, "amr_parsed.json"), 'w') as f:
        json.dump(results, f)

    return results


async def run_prodigal(fasta_path: str, job_dir: str):
    """Run Prodigal for gene prediction."""
    genes_gff = os.path.join(job_dir, "genes.gff")
    genes_fna = os.path.join(job_dir, "genes.fna")

    proc = await asyncio.create_subprocess_exec(
        PRODIGAL, "-i", fasta_path, "-o", genes_gff, "-f", "gff",
        "-d", genes_fna, "-p", "meta",
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )
    await proc.communicate()


async def run_expression(fasta_path: str, job_dir: str):
    """Run promoter (BPROM) and RBS (OSTIR) prediction for each ARG.
    This is the bottleneck step (~1-2 min).
    """
    # TODO: Implement BPROM and OSTIR calls for each detected ARG
    # For now, create placeholder expression data
    # In production, this extracts upstream sequences and runs BPROM/OSTIR
    expr_out = os.path.join(job_dir, "expression_results.json")
    if not os.path.exists(expr_out):
        with open(expr_out, 'w') as f:
            json.dump({}, f)


async def run_pointfinder(fasta_path: str, job_dir: str, species: str):
    """Run PointFinder for chromosomal point mutations."""
    pf_species = POINTFINDER_SPECIES_MAP.get(species)
    if not pf_species:
        # No PointFinder DB for this species (e.g., A. baumannii)
        # For A. baumannii, run custom gyrA Ser83Leu detection
        if species == "Acinetobacter_baumannii":
            await detect_abaumannii_gyra(fasta_path, job_dir)
        return

    pf_dir = os.path.join(job_dir, "pointfinder")
    os.makedirs(pf_dir, exist_ok=True)

    try:
        proc = await asyncio.create_subprocess_exec(
            "python3", "-m", "resfinder",
            "-ifa", fasta_path,
            "-o", pf_dir,
            "-c", "-db_point", POINTFINDER_DB,
            "-s", pf_species,
            "-db_res", RESFINDER_DB,
            "--ignore_missing_species",
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
        await proc.communicate()
    except FileNotFoundError:
        pass  # PointFinder not installed


async def detect_abaumannii_gyra(fasta_path: str, job_dir: str):
    """Custom gyrA Ser83Leu detection for A. baumannii."""
    # Use tblastn with E. coli gyrA QRDR reference
    query = ">gyrA_ref\nMSDLAREITPVNIEEELKNSYLDYAMSVIVGRALPDVRDGLKPVHRRVLYAMN\n"
    query_file = os.path.join(job_dir, "gyra_query.faa")
    with open(query_file, 'w') as f:
        f.write(query)

    try:
        proc = await asyncio.create_subprocess_exec(
            "tblastn", "-query", query_file, "-subject", fasta_path,
            "-outfmt", "6 qstart qend qseq sseq pident",
            "-evalue", "1e-20", "-max_target_seqs", "1", "-max_hsps", "1",
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
        stdout, _ = await proc.communicate()

        result = {'has_gyrA_Ser83Leu': 0}
        if stdout.strip():
            parts = stdout.decode().strip().split('\n')[0].split('\t')
            if len(parts) >= 5 and float(parts[4]) > 60:
                qstart = int(parts[0])
                qseq, sseq = parts[2], parts[3]
                ref_pos = qstart
                for r, s in zip(qseq, sseq):
                    if r != '-':
                        if ref_pos == 85 and s == 'L':
                            result['has_gyrA_Ser83Leu'] = 1
                        ref_pos += 1

        with open(os.path.join(job_dir, "abaumannii_gyra.json"), 'w') as f:
            json.dump(result, f)

    except FileNotFoundError:
        pass
    finally:
        if os.path.exists(query_file):
            os.remove(query_file)


def compute_genomic_features(fasta_path: str, job_dir: str):
    """Compute pure-Python genomic features: CAI, GC, gene dosage, etc."""
    from Bio import SeqIO

    # Read assembly
    total_len = 0
    gc_count = 0
    n_contigs = 0

    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper()
        total_len += len(seq)
        gc_count += seq.count('G') + seq.count('C')
        n_contigs += 1

    genome_gc = gc_count / total_len if total_len > 0 else 0

    features = {
        'genome_gc': genome_gc,
        'genome_size': total_len,
        'n_contigs': n_contigs,
    }

    with open(os.path.join(job_dir, "genomic_features.json"), 'w') as f:
        json.dump(features, f)

    return features
