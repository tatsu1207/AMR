FROM mambaorg/micromamba:1.5-jammy

USER root

# Install utilities; on x86, also install 32-bit libs for BPROM
ARG TARGETARCH
RUN apt-get update && \
    if [ "$TARGETARCH" = "amd64" ]; then \
        dpkg --add-architecture i386 && apt-get update && \
        apt-get install -y --no-install-recommends libc6:i386 libgcc-s1:i386; \
    fi && \
    apt-get install -y --no-install-recommends wget curl git procps && \
    rm -rf /var/lib/apt/lists/*

# Copy environment definition and create conda env
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba create -n amr -f /tmp/environment.yml --channel-priority flexible && \
    micromamba clean --all --yes

# Set PATH so all conda tools are available
ENV PATH="/opt/conda/envs/amr/bin:$PATH"
ENV MAMBA_DOCKERFILE_ACTIVATE_ENVIRONMENT=amr

# Install ViennaRNA (for OSTIR) - separate step to avoid solver conflicts
RUN micromamba install -n amr -c bioconda -c conda-forge viennarna --channel-priority flexible -y && \
    micromamba clean --all --yes || true

# Update AMRFinderPlus database
RUN amrfinder --update || true

# Update MLST schemes
RUN mlst --check || true

# Set working directory
WORKDIR /opt/amr

# Copy application code and resources
COPY app/ app/
COPY models/ models/
COPY databases/ databases/
COPY bin/ bin/
COPY annotate_minimal.sh .
COPY scripts/ scripts/
COPY run.sh .

# Make scripts executable; BPROM binary only works on x86
RUN chmod +x annotate_minimal.sh run.sh && \
    if [ "$TARGETARCH" = "amd64" ]; then chmod +x bin/bprom; fi

# Create upload directory
RUN mkdir -p uploads

# Environment variables for tool/database paths
# CONDA_PREFIX must point to the amr env so AMRFinderPlus finds its database
ENV CONDA_PREFIX=/opt/conda/envs/amr \
    AMR_MODELS_DIR=/opt/amr/models \
    AMR_UPLOAD_DIR=/opt/amr/uploads \
    POINTFINDER_DB=/opt/amr/databases/pointfinder_db \
    RESFINDER_DB=/opt/amr/databases/resfinder_db \
    BPROM=/opt/amr/bin/bprom \
    BPROM_DATA=/opt/amr/databases/bprom_data \
    TSS_DATA=/opt/amr/databases/bprom_data

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=5s --retries=3 \
    CMD curl -f http://localhost:8000/ || exit 1

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
