"""Pydantic models for API request/response."""
from pydantic import BaseModel
from typing import Optional
from enum import Enum


class JobStage(str, Enum):
    QUEUED = "queued"
    VALIDATING = "validating"
    IDENTIFYING_SPECIES = "identifying_species"
    DETECTING_GENES = "detecting_genes"
    PREDICTING_STRUCTURE = "predicting_structure"
    ANALYZING_EXPRESSION = "analyzing_expression"
    DETECTING_MUTATIONS = "detecting_mutations"
    COMPUTING_FEATURES = "computing_features"
    RUNNING_PREDICTIONS = "running_predictions"
    COMPLETE = "complete"
    FAILED = "failed"


STAGE_LABELS = {
    JobStage.QUEUED: "Queued",
    JobStage.VALIDATING: "Validating input",
    JobStage.IDENTIFYING_SPECIES: "Identifying species",
    JobStage.DETECTING_GENES: "Detecting resistance genes",
    JobStage.PREDICTING_STRUCTURE: "Predicting gene structure",
    JobStage.ANALYZING_EXPRESSION: "Analyzing gene expression context",
    JobStage.DETECTING_MUTATIONS: "Detecting point mutations",
    JobStage.COMPUTING_FEATURES: "Computing genomic features",
    JobStage.RUNNING_PREDICTIONS: "Running ML predictions",
    JobStage.COMPLETE: "Analysis complete",
    JobStage.FAILED: "Analysis failed",
}

STAGE_PROGRESS = {
    JobStage.QUEUED: 0,
    JobStage.VALIDATING: 5,
    JobStage.IDENTIFYING_SPECIES: 10,
    JobStage.DETECTING_GENES: 25,
    JobStage.PREDICTING_STRUCTURE: 40,
    JobStage.ANALYZING_EXPRESSION: 55,
    JobStage.DETECTING_MUTATIONS: 70,
    JobStage.COMPUTING_FEATURES: 80,
    JobStage.RUNNING_PREDICTIONS: 90,
    JobStage.COMPLETE: 100,
    JobStage.FAILED: 0,
}


class JobStatus(BaseModel):
    job_id: str
    stage: JobStage
    stage_label: str
    progress: int
    species: Optional[str] = None
    mlst_st: Optional[str] = None
    error: Optional[str] = None


class GeneResult(BaseModel):
    gene: str
    drug_class: str
    identity: float
    coverage: float
    on_plasmid: bool = False


class AntibioticPrediction(BaseModel):
    antibiotic: str
    drug_class: str
    prediction: str  # "Resistant" or "Susceptible"
    probability: float
    confidence: str  # "High", "Moderate", "Low"
    key_genes: list[str]
    key_mutations: list[str]


class PredictionResults(BaseModel):
    job_id: str
    species: str
    species_display: str
    mlst_st: str
    n_antibiotics: int
    n_resistant: int
    n_susceptible: int
    predictions: list[AntibioticPrediction]
    detected_genes: list[GeneResult]
    pipeline_version: str = "1.0.0"
