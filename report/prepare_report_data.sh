#!/usr/bin/env bash

set -e  # exit on error

ROOT_DIR="$PWD"
BIOPROJECT=""
SAMPLES=()
REPORT_NAME=""

if [ "$(basename $PWD)" == "mgs-pipeline" ]; then
    S3_DIR="s3://nao-mgs/"
    MGS_PIPELINE_DIR="$ROOT_DIR"
elif [ "$(basename $PWD)" == "mgs-restricted" ]; then
    S3_DIR="s3://nao-restricted/"
    MGS_PIPELINE_DIR=$ROOT_DIR/../mgs-pipeline
else
    echo "Run this from either mgs-pipeline or mgs-restricted"
    exit 1
fi

# parse arguments. Use -b to specify the bioproject, -s for the samples (comma-separated), 
# and -r for the report name.
while getopts "b:s:r:" opt; do
    case "$opt" in
    b)  BIOPROJECT="$OPTARG"
        ;;
    s)  IFS=',' read -ra ADDR <<< "$OPTARG"
        for sample in "${ADDR[@]}"; do
            SAMPLES+=("$sample")
        done
        ;;
    r)  REPORT_NAME="$OPTARG"
        ;;
*)  echo "Usage: $0 -b bioproject [-s samples (comma separated, optional)] [-r reportname (optional)]"
        exit 1
    esac
done

# Check if the BIOPROJECT variable is empty
if [ -z "$BIOPROJECT" ]; then
    echo "Error: A bioproject must be provided using the -b option."
    echo "Usage: $0 -b bioproject [-s samples (comma separated, optional)] [-r reportname (optional)]"
    exit 1
fi 

# Set default report name if it's not provided by the user.
if [ -z "$REPORT_NAME" ]; then
    DATE=$(date +"%Y%m%d")
    REPORT_NAME="${DATE}"
fi

# Determine the list of samples.
# If SAMPLES array is empty, grab all sample files from the bioproject.
# We are just using the cladecounts file to get sample IDs
SAMPLE_LIST=()
if [ ${#SAMPLES[@]} -eq 0 ]; then
    S3_BIOPROJECT_PATH="${S3_DIR}${BIOPROJECT}/cladecounts/"
    for cc in $(aws s3 ls $S3_BIOPROJECT_PATH | awk '{print $NF}' | sed 's/.tsv.gz//'); do
        SAMPLE_LIST+=("$cc")
    done
else
    for sample in "${SAMPLES[@]}"; do
        SAMPLE_LIST+=("$sample")
    done
fi

# Create directory for the report
REPORT_DIR="${MGS_PIPELINE_DIR}/bioprojects/${BIOPROJECT}/reports/${REPORT_NAME}"
mkdir -p "$REPORT_DIR"

# Fetch cladecounts for the samples in SAMPLE_LIST
for sample in "${SAMPLE_LIST[@]}"; do
    S3_SAMPLE_PATH="${S3_DIR}${BIOPROJECT}/cladecounts/${sample}.tsv.gz"
    if [ ! -s "$REPORT_DIR/cladecounts/${sample}.tsv.gz" ]; then
        echo $S3_SAMPLE_PATH
    fi
done | xargs -I {} -P 32 aws s3 cp {} "${REPORT_DIR}/cladecounts/"

$MGS_PIPELINE_DIR/report/prepare_report_data.py $ROOT_DIR $MGS_PIPELINE_DIR $BIOPROJECT "${SAMPLE_LIST[*]}" $REPORT_NAME

