#!/usr/bin/env bash

# Search protein sequences for DNA repair proteins using HMM profiles
#
# Usage:
#   ./search_dna_repair_hmms.sh \
#       --subjects /path/to/proteins \
#       --hmms /path/to/hmms \
#       --output /path/to/results

set -euo pipefail


# ============================================================
# Default values
# ============================================================

SUBJECT_DIR=""
HMM_PROFILE_DIR=""
OUTPUT_DIR=""


# ============================================================
# Help message
# ============================================================

usage() {
    echo ""
    echo "DNA Repair HMM Search"
    echo ""
    echo "Usage:"
    echo "  $0 --subjects DIR --hmms DIR --output DIR"
    echo ""
    echo "Options:"
    echo "  --subjects DIR    Directory containing protein FASTA files (.faa)"
    echo "  --hmms DIR        Directory containing HMM profiles (.hmm)"
    echo "  --output DIR      Directory where results will be written"
    echo "  -h, --help        Show this help message"
    echo ""
    exit 0
}


# ============================================================
# Parse command-line arguments
# ============================================================

while [[ $# -gt 0 ]]; do

    case "$1" in

        --subjects)
            SUBJECT_DIR="$2"
            shift 2
            ;;

        --hmms)
            HMM_PROFILE_DIR="$2"
            shift 2
            ;;

        --output)
            OUTPUT_DIR="$2"
            shift 2
            ;;

        -h|--help)
            usage
            ;;

        *)
            echo "Error: Unknown option: $1"
            usage
            ;;

    esac

done


# ============================================================
# Check required arguments
# ============================================================

if [[ -z "$SUBJECT_DIR" || -z "$HMM_PROFILE_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: --subjects, --hmms and --output are required."
    echo ""
    usage
fi


# ============================================================
# Check input directories
# ============================================================

if [[ ! -d "$SUBJECT_DIR" ]]; then
    echo "Error: Subject directory does not exist:"
    echo "$SUBJECT_DIR"
    exit 1
fi

if [[ ! -d "$HMM_PROFILE_DIR" ]]; then
    echo "Error: HMM directory does not exist:"
    echo "$HMM_PROFILE_DIR"
    exit 1
fi


# ============================================================
# Check HMMER
# ============================================================

if ! command -v hmmsearch &> /dev/null; then
    echo "Error: hmmsearch was not found."
    echo "Please install HMMER or activate the appropriate environment."
    exit 1
fi


# ============================================================
# Create output directory
# ============================================================

mkdir -p "$OUTPUT_DIR"


# ============================================================
# Count input files
# ============================================================

FAA_COUNT=$(find "$SUBJECT_DIR" -maxdepth 1 -type f -name "*.faa" | wc -l)
HMM_COUNT=$(find "$HMM_PROFILE_DIR" -maxdepth 1 -type f -name "*.hmm" | wc -l)

echo ""
echo "=========================================="
echo " DNA Repair HMM Search"
echo "=========================================="
echo ""
echo "Protein directory : $SUBJECT_DIR"
echo "HMM directory     : $HMM_PROFILE_DIR"
echo "Output directory  : $OUTPUT_DIR"
echo ""
echo "Protein FASTA files: $FAA_COUNT"
echo "HMM profiles       : $HMM_COUNT"
echo ""


# ============================================================
# Main analysis
# ============================================================

for faaFile in "$SUBJECT_DIR"/*.faa; do

    # Skip if no FASTA files were found
    [[ -e "$faaFile" ]] || continue

    FILENAME=$(basename "$faaFile")
    BASENAME="${FILENAME%.faa}"

    echo "------------------------------------------"
    echo "Analyzing: $FILENAME"
    echo "------------------------------------------"

    # Create directory for this genome/sample
    SAMPLE_OUTPUT="$OUTPUT_DIR/$BASENAME"
    mkdir -p "$SAMPLE_OUTPUT"


    # ------------------------------------------
    # Search every HMM against this protein set
    # ------------------------------------------

    for hmmFile in "$HMM_PROFILE_DIR"/*.hmm; do

        [[ -e "$hmmFile" ]] || continue

        HMM_FILENAME=$(basename "$hmmFile")
        HMM_PREFIX="${HMM_FILENAME%.hmm}"

        echo "  Searching HMM: $HMM_PREFIX"

        hmmsearch \
            --tblout "$SAMPLE_OUTPUT/${HMM_PREFIX}.tblout" \
            --domtblout "$SAMPLE_OUTPUT/${HMM_PREFIX}.domtblout" \
            "$hmmFile" \
            "$faaFile" \
            > "$SAMPLE_OUTPUT/${HMM_PREFIX}.txt"

    done

    echo "Finished: $BASENAME"
    echo ""

done


echo "=========================================="
echo "Analysis completed!"
echo "Results: $OUTPUT_DIR"
echo "=========================================="
