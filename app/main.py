"""AMR Predictor - FastAPI Application."""
import os
import uuid
import asyncio
import shutil
from contextlib import asynccontextmanager

from fastapi import FastAPI, UploadFile, File, HTTPException
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse, JSONResponse

from app.config import UPLOAD_DIR, MAX_FILE_SIZE_MB, MAX_CONCURRENT_JOBS, CLEANUP_HOURS
from app.models.job_store import create_job, get_job, cleanup_old_jobs
from app.models.schemas import JobStage
from app.pipeline.runner import run_pipeline

# Semaphore for concurrent job limit
job_semaphore = asyncio.Semaphore(MAX_CONCURRENT_JOBS)


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Startup/shutdown."""
    os.makedirs(UPLOAD_DIR, exist_ok=True)
    # Start periodic cleanup
    cleanup_task = asyncio.create_task(periodic_cleanup())
    # Load models at startup
    from app.pipeline.predict import load_all_models
    load_all_models()
    yield
    cleanup_task.cancel()


app = FastAPI(
    title="AMR Predictor",
    description="Antimicrobial resistance phenotype prediction from WGS",
    version="1.0.0",
    lifespan=lifespan,
)

# Static files
static_dir = os.path.join(os.path.dirname(__file__), "static")
app.mount("/static", StaticFiles(directory=static_dir), name="static")


@app.get("/")
async def index():
    return FileResponse(os.path.join(static_dir, "index.html"))


@app.post("/api/upload")
async def upload_genome(file: UploadFile = File(...)):
    """Upload a FASTA file and start analysis."""
    # Validate filename
    valid_exts = ('.fasta', '.fa', '.fna', '.fasta.gz')
    if not file.filename.lower().endswith(valid_exts):
        raise HTTPException(400, "Invalid file type. Please upload a FASTA file.")

    # Read and validate size
    contents = await file.read()
    if len(contents) > MAX_FILE_SIZE_MB * 1024 * 1024:
        raise HTTPException(400, f"File too large. Maximum size is {MAX_FILE_SIZE_MB} MB.")
    if len(contents) < 100:
        raise HTTPException(400, "File too small. Please upload a valid assembled genome.")

    # Create job
    job_id = str(uuid.uuid4())[:12]
    job_dir = os.path.join(UPLOAD_DIR, job_id)
    os.makedirs(job_dir, exist_ok=True)

    # Save file
    fasta_path = os.path.join(job_dir, "assembly.fasta")
    with open(fasta_path, "wb") as f:
        f.write(contents)

    # Create job state
    job = create_job(job_id, file.filename)

    # Launch pipeline in background
    asyncio.create_task(run_with_semaphore(job_id, fasta_path))

    return {"job_id": job_id, "status": "queued"}


@app.get("/api/job/{job_id}")
async def get_job_status(job_id: str):
    """Poll job progress."""
    job = get_job(job_id)
    if not job:
        raise HTTPException(404, "Job not found")
    return job.to_status().model_dump()


@app.get("/api/results/{job_id}")
async def get_results(job_id: str):
    """Get final prediction results."""
    job = get_job(job_id)
    if not job:
        raise HTTPException(404, "Job not found")
    if job.stage != JobStage.COMPLETE:
        raise HTTPException(400, "Analysis not yet complete")
    if not job.results:
        raise HTTPException(500, "Results not available")
    return job.results


async def run_with_semaphore(job_id: str, fasta_path: str):
    """Run pipeline with concurrency limit."""
    async with job_semaphore:
        await run_pipeline(job_id, fasta_path)


async def periodic_cleanup():
    """Clean up old jobs periodically."""
    while True:
        await asyncio.sleep(3600)
        cleanup_old_jobs(CLEANUP_HOURS)
        # Clean up old upload directories
        if os.path.exists(UPLOAD_DIR):
            import time
            cutoff = time.time() - CLEANUP_HOURS * 3600
            for d in os.listdir(UPLOAD_DIR):
                dp = os.path.join(UPLOAD_DIR, d)
                if os.path.isdir(dp) and os.path.getmtime(dp) < cutoff:
                    shutil.rmtree(dp, ignore_errors=True)
