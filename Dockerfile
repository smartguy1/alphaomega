# AlphaOmega — Multi-stage Dockerfile
# Phase 1: Core pipeline (BAM → VCF)

FROM python:3.11-slim AS base

# Install bioinformatics tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    bcftools \
    samtools \
    tabix \
    wget \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install Python dependencies
COPY requirements.txt .
RUN pip install --default-timeout=1000 --no-cache-dir -r requirements.txt

# Copy source code
COPY src/ ./src/
COPY config/ ./config/
COPY scripts/ ./scripts/

# Make scripts executable
RUN chmod +x scripts/*.sh

# Create data directories
RUN mkdir -p profiles db reference data output

ENTRYPOINT ["python", "-m", "src.cli"]
