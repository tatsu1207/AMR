"""Application configuration."""
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODELS_DIR = os.environ.get("AMR_MODELS_DIR", os.path.join(BASE_DIR, "models"))
UPLOAD_DIR = os.environ.get("AMR_UPLOAD_DIR", os.path.join(BASE_DIR, "uploads"))
MAX_CONCURRENT_JOBS = int(os.environ.get("MAX_CONCURRENT_JOBS", "5"))
CLEANUP_HOURS = int(os.environ.get("CLEANUP_HOURS", "24"))
MAX_FILE_SIZE_MB = int(os.environ.get("MAX_FILE_SIZE_MB", "100"))

# Tool paths (overridable via env for Docker)
AMRFINDERPLUS = os.environ.get("AMRFINDERPLUS", "amrfinder")
PRODIGAL = os.environ.get("PRODIGAL", "prodigal")
MLST = os.environ.get("MLST", "mlst")
BPROM = os.environ.get("BPROM", os.path.join(BASE_DIR, "bin", "bprom"))
BPROM_DATA = os.environ.get("BPROM_DATA", os.path.join(BASE_DIR, "databases", "bprom_data"))
POINTFINDER_DB = os.environ.get("POINTFINDER_DB", os.path.join(BASE_DIR, "databases", "pointfinder_db"))
RESFINDER_DB = os.environ.get("RESFINDER_DB", os.path.join(BASE_DIR, "databases", "resfinder_db"))

SUPPORTED_SPECIES = {
    "senterica": "Salmonella_enterica",
    "ecoli": "Escherichia_coli",
    "klebsiella": "Klebsiella_pneumoniae",
    "saureus": "Staphylococcus_aureus",
    "abaumannii": "Acinetobacter_baumannii",
}

SPECIES_DISPLAY = {
    "Salmonella_enterica": "Salmonella enterica",
    "Escherichia_coli": "Escherichia coli",
    "Klebsiella_pneumoniae": "Klebsiella pneumoniae",
    "Staphylococcus_aureus": "Staphylococcus aureus",
    "Acinetobacter_baumannii": "Acinetobacter baumannii",
}

# MLST scheme to species mapping
MLST_SCHEME_MAP = {
    "senterica": "Salmonella_enterica",
    "ecoli": "Escherichia_coli",
    "ecoli_achtman_4": "Escherichia_coli",
    "klebsiella": "Klebsiella_pneumoniae",
    "kpneumoniae": "Klebsiella_pneumoniae",
    "saureus": "Staphylococcus_aureus",
    "abaumannii": "Acinetobacter_baumannii",
    "abaumannii_2": "Acinetobacter_baumannii",
}

# PointFinder species mapping
POINTFINDER_SPECIES_MAP = {
    "Salmonella_enterica": "salmonella",
    "Escherichia_coli": "escherichia_coli",
    "Klebsiella_pneumoniae": "klebsiella",
    "Staphylococcus_aureus": "staphylococcus_aureus",
}
