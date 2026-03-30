"""Species identification using MLST."""
import asyncio
import os
from app.config import MLST, MLST_SCHEME_MAP, SPECIES_DISPLAY


async def identify_species(fasta_path: str, job_dir: str):
    """Identify species using mlst tool. Returns (species_code, ST, scheme)."""
    proc = await asyncio.create_subprocess_exec(
        MLST, fasta_path,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )
    stdout, stderr = await proc.communicate()

    if proc.returncode != 0:
        raise RuntimeError(f"MLST failed: {stderr.decode()}")

    # Parse mlst output: filename\tscheme\tST\talleles...
    line = stdout.decode().strip()
    parts = line.split('\t')
    if len(parts) < 3:
        return None, None, "unknown"

    scheme = parts[1].strip().lower().replace(' ', '_')
    st = parts[2].strip()
    if st == '-':
        st = 'unknown'

    # Map scheme to species
    species = MLST_SCHEME_MAP.get(scheme)

    # Save MLST result
    mlst_out = os.path.join(job_dir, "mlst_result.tsv")
    with open(mlst_out, 'w') as f:
        f.write(line + '\n')

    return species, st, scheme
