"""Pipeline orchestrator - runs annotate_minimal.sh then predicts."""
import os
import csv
import json
import asyncio
import traceback

from app.models.job_store import get_job
from app.models.schemas import JobStage
from app.pipeline.features import extract_features
from app.pipeline.predict import predict_phenotypes
from app.config import MODELS_DIR


SCRIPT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ANNOTATE_SCRIPT = os.path.join(SCRIPT_DIR, "annotate_minimal.sh")


async def run_pipeline(job_id: str, fasta_path: str):
    """Run the minimal annotation pipeline and predict phenotypes."""
    job = get_job(job_id)
    if not job:
        return

    job_dir = os.path.dirname(fasta_path)

    try:
        # Step 1: Validate input
        job.update(JobStage.VALIDATING)
        validate_fasta(fasta_path)

        # Step 2-6: Run annotate_minimal.sh (MLST + AMRFinderPlus + Prodigal + PointFinder + BPROM + OSTIR + genomic features)
        job.update(JobStage.IDENTIFYING_SPECIES)

        proc = await asyncio.create_subprocess_exec(
            "bash", ANNOTATE_SCRIPT, fasta_path, job_dir, "4",
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )

        # Monitor output for progress updates
        while True:
            line = await proc.stdout.readline()
            if not line:
                break
            text = line.decode().strip()
            if "Step 1" in text or "MLST" in text:
                job.update(JobStage.IDENTIFYING_SPECIES)
            elif "Step 2" in text or "AMRFinder" in text:
                job.update(JobStage.DETECTING_GENES)
                # Parse species from MLST
                if "Species scheme:" in text:
                    parts = text.split("Species scheme:")[1].strip()
                    scheme = parts.split(",")[0].strip()
                    st = parts.split("ST:")[1].strip() if "ST:" in parts else "unknown"
                    from app.config import MLST_SCHEME_MAP, SPECIES_DISPLAY
                    species = MLST_SCHEME_MAP.get(scheme.lower())
                    if species:
                        job.update(JobStage.DETECTING_GENES,
                                   species=SPECIES_DISPLAY.get(species, species),
                                   mlst_st=st)
            elif "Step 3" in text or "Prodigal" in text:
                job.update(JobStage.PREDICTING_STRUCTURE)
            elif "Step 4" in text or "PointFinder" in text:
                job.update(JobStage.DETECTING_MUTATIONS)
            elif "Step 5" in text or "Promoter" in text or "RBS" in text:
                job.update(JobStage.ANALYZING_EXPRESSION)
            elif "Step" in text and ("6" in text or "7" in text):
                job.update(JobStage.COMPUTING_FEATURES)

        await proc.wait()

        # Check if annotation succeeded
        summary_path = os.path.join(job_dir, "arg_context_summary.tsv")
        if not os.path.exists(summary_path):
            stderr_out = (await proc.stderr.read()).decode()
            job.update(JobStage.FAILED,
                       error=f"Annotation pipeline failed. {stderr_out[:500]}")
            return

        # Determine species from MLST results
        mlst_path = os.path.join(job_dir, "mlst_results.tsv")
        species, mlst_st = parse_mlst(mlst_path)

        if not species:
            job.update(JobStage.FAILED,
                       error="Unsupported or unrecognized species. "
                             "Supported: Salmonella enterica, Escherichia coli, "
                             "Klebsiella pneumoniae, Staphylococcus aureus, "
                             "Acinetobacter baumannii.")
            return

        job.update(JobStage.RUNNING_PREDICTIONS, species=species, mlst_st=mlst_st)

        # Load AMR gene detection results
        amr_results = parse_amrfinderplus(os.path.join(job_dir, "amrfinderplus.tsv"))

        # Extract features and predict
        feature_data = extract_features(job_dir, species, mlst_st)
        results = predict_phenotypes(species, feature_data, mlst_st, amr_results=amr_results)

        job.update(JobStage.COMPLETE, results=results)

    except Exception as e:
        traceback.print_exc()
        job.update(JobStage.FAILED, error=str(e))


def parse_amrfinderplus(tsv_path: str):
    """Parse AMRFinderPlus TSV output into structured results."""
    results = {'all_genes': [], 'genes_by_class': {}, 'mutations': []}
    if not os.path.exists(tsv_path):
        return results

    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene = row.get('Element symbol', '')
            drug_class = row.get('Class', row.get('Subclass', ''))
            etype = row.get('Type', '')
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

    return results


def parse_mlst(mlst_path: str):
    """Parse MLST results to determine species."""
    from app.config import MLST_SCHEME_MAP
    if not os.path.exists(mlst_path):
        return None, None
    with open(mlst_path) as f:
        line = f.readline().strip()
    parts = line.split('\t')
    if len(parts) < 3:
        return None, None
    scheme = parts[1].strip().lower()
    st = parts[2].strip()
    if st == '-':
        st = 'unknown'
    species = MLST_SCHEME_MAP.get(scheme)
    return species, st


def validate_fasta(fasta_path: str):
    """Basic FASTA validation."""
    with open(fasta_path) as f:
        first_line = f.readline().strip()
    if not first_line.startswith('>'):
        raise ValueError("Invalid FASTA format: file must start with '>'")
    n_contigs = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                n_contigs += 1
    if n_contigs == 0:
        raise ValueError("No sequences found in FASTA file")
    if n_contigs > 5000:
        raise ValueError(f"Too many contigs ({n_contigs}). Please provide an assembled genome, not raw reads.")
