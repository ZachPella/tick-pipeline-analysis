#!/bin/bash
#SBATCH --job-name=remove_dups
#SBATCH --time=2-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --mem=35G
#SBATCH --array=1-62
#SBATCH --partition=batch

## Record relevant job info
START_DIR=$(pwd)
HOST_NAME=$(hostname)
RUN_DATE=$(date)
echo "Starting working directory: ${START_DIR}"
echo "Host name: ${HOST_NAME}"
echo "Run date: ${RUN_DATE}"
printf "\n"

## Set working directory and variables
BASEDIR=/work/fauverlab/zachpella/scripts_ticksJune2025_10_scatter
INPUTDIR=${BASEDIR}/readgroups  # Input directory for BAM files with read groups
WORKDIR=${BASEDIR}/dedup        # Output directory for deduplicated BAM files
SAMPLE_LIST=${BASEDIR}/sample_list.txt

mkdir -p ${WORKDIR}

## Check if sample list exists
if [ ! -f "$SAMPLE_LIST" ]; then
    echo "Error: Sample list file not found: $SAMPLE_LIST"
    exit 1
fi

## Get total number of samples
TOTAL_SAMPLES=$(wc -l < "$SAMPLE_LIST")

## Check if array task ID is valid
if [ ${SLURM_ARRAY_TASK_ID} -gt ${TOTAL_SAMPLES} ]; then
    echo "Error: Array task ID ${SLURM_ARRAY_TASK_ID} exceeds number of samples (${TOTAL_SAMPLES})"
    exit 1
fi

## Get sample name from list - bulletproof method!
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

## Verify sample name is not empty
if [[ -z "$SAMPLE" ]]; then
    echo "Error: Empty sample name for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

## Set file paths
INPUT_BAM="${INPUTDIR}/${SAMPLE}.rg.sorted.bam"
OUTPUT_BAM=${WORKDIR}/${SAMPLE}.dedup.rg.sorted.bam
METRICS_FILE=${WORKDIR}/${SAMPLE}.dedup_metrics.txt

## Check if input BAM file exists
if [ ! -f "${INPUT_BAM}" ]; then
    echo "Error: Input BAM file not found: ${INPUT_BAM}"
    exit 1
fi

## Confirm that variables are properly assigned
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Sample: ${SAMPLE}"
echo "Input BAM: ${INPUT_BAM}"
echo "Output BAM: ${OUTPUT_BAM}"
echo "Metrics file: ${METRICS_FILE}"
echo "Marking duplicates for ${SAMPLE}..."

## Load modules
module purge
module load picard
module load samtools/1.19

## Move into working directory
cd ${WORKDIR}

## Mark duplicates using Picard
echo "Marking duplicates in ${INPUT_BAM}..."
picard MarkDuplicates \
    I=${INPUT_BAM} \
    O=${OUTPUT_BAM} \
    M=${METRICS_FILE} \
    REMOVE_DUPLICATES=true \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
    CREATE_INDEX=true

