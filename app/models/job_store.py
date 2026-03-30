"""In-memory job state store."""
import time
from typing import Optional
from app.models.schemas import JobStage, STAGE_LABELS, STAGE_PROGRESS, JobStatus


class JobState:
    def __init__(self, job_id: str, filename: str):
        self.job_id = job_id
        self.filename = filename
        self.stage = JobStage.QUEUED
        self.species: Optional[str] = None
        self.mlst_st: Optional[str] = None
        self.error: Optional[str] = None
        self.results: Optional[dict] = None
        self.created_at = time.time()

    def update(self, stage: JobStage, **kwargs):
        self.stage = stage
        for k, v in kwargs.items():
            setattr(self, k, v)

    def to_status(self) -> JobStatus:
        return JobStatus(
            job_id=self.job_id,
            stage=self.stage,
            stage_label=STAGE_LABELS[self.stage],
            progress=STAGE_PROGRESS[self.stage],
            species=self.species,
            mlst_st=self.mlst_st,
            error=self.error,
        )


# Global store
_jobs: dict[str, JobState] = {}


def create_job(job_id: str, filename: str) -> JobState:
    job = JobState(job_id, filename)
    _jobs[job_id] = job
    return job


def get_job(job_id: str) -> Optional[JobState]:
    return _jobs.get(job_id)


def cleanup_old_jobs(max_age_hours: int = 24):
    cutoff = time.time() - max_age_hours * 3600
    expired = [jid for jid, j in _jobs.items() if j.created_at < cutoff]
    for jid in expired:
        del _jobs[jid]
