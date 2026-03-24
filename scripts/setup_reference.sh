#!/usr/bin/env bash
# AlphaOmega — Reference Genome Setup
#
# Downloads and indexes reference genomes for variant calling.
# Run once before first analysis.
#
# Usage:
#   ./scripts/setup_reference.sh [--grch37] [--grch38] [--both]
#   Default: --grch38 only

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
REF_DIR="${PROJECT_DIR}/reference"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

log_info()  { echo -e "${GREEN}[INFO]${NC} $*"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC} $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*"; }

# GRCh38 reference (recommended)
GRCH38_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
GRCH38_FILE="GRCh38.fa"

# GRCh37 reference (fallback)
GRCH37_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"
GRCH37_FILE="GRCh37.fa"

check_dependencies() {
    local missing=()
    for cmd in samtools wget; do
        if ! command -v "$cmd" &> /dev/null; then
            missing+=("$cmd")
        fi
    done
    if [ ${#missing[@]} -gt 0 ]; then
        log_error "Missing dependencies: ${missing[*]}"
        log_error "Install with: apt-get install ${missing[*]} (or conda)"
        exit 1
    fi
}

download_grch38() {
    local target="${REF_DIR}/${GRCH38_FILE}"

    if [ -f "$target" ]; then
        log_info "GRCh38 already exists: $target"
        return
    fi

    log_info "Downloading GRCh38 reference (~900MB compressed)..."
    log_info "This will take several minutes depending on your connection."

    wget -c -O "$target" "$GRCH38_URL"

    log_info "Indexing GRCh38..."
    samtools faidx "$target"

    log_info "GRCh38 ready: $target"
}

download_grch37() {
    local target="${REF_DIR}/${GRCH37_FILE}"
    local gz_target="${target}.gz"

    if [ -f "$target" ]; then
        log_info "GRCh37 already exists: $target"
        return
    fi

    log_info "Downloading GRCh37 reference (~900MB compressed)..."
    wget -c -O "$gz_target" "$GRCH37_URL"

    log_info "Decompressing GRCh37..."
    gunzip -k "$gz_target"
    mv "${REF_DIR}/human_g1k_v37.fasta" "$target" 2>/dev/null || \
    mv "$gz_target" "$target" 2>/dev/null || true

    log_info "Indexing GRCh37..."
    samtools faidx "$target"

    log_info "GRCh37 ready: $target"
}

show_usage() {
    echo "Usage: $0 [--grch37] [--grch38] [--both]"
    echo ""
    echo "Options:"
    echo "  --grch38    Download GRCh38 only (default)"
    echo "  --grch37    Download GRCh37 only"
    echo "  --both      Download both references"
    echo ""
    echo "Reference genomes are stored in: ${REF_DIR}/"
    echo "Total disk space needed: ~3GB per build"
}

main() {
    local do_grch37=false
    local do_grch38=false

    # Parse args
    if [ $# -eq 0 ]; then
        do_grch38=true
    fi

    for arg in "$@"; do
        case $arg in
            --grch37) do_grch37=true ;;
            --grch38) do_grch38=true ;;
            --both)   do_grch37=true; do_grch38=true ;;
            --help|-h) show_usage; exit 0 ;;
            *) log_error "Unknown option: $arg"; show_usage; exit 1 ;;
        esac
    done

    check_dependencies

    # Create reference directory
    mkdir -p "$REF_DIR"

    if $do_grch38; then
        download_grch38
    fi

    if $do_grch37; then
        download_grch37
    fi

    # Show disk usage
    echo ""
    log_info "Reference directory contents:"
    ls -lh "$REF_DIR/" 2>/dev/null || true
    echo ""
    log_info "Setup complete!"
}

main "$@"
