#!/bin/bash
#SBATCH --job-name=g1.5_fastqc
#SBATCH --mail-user=zpella@unmc.edu
#SBATCH --mail-type=ALL
#SBATCH --time=0-06:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G
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
WORKDIR=${BASEDIR}/trimmed_reads
QCDIR=${BASEDIR}/fastqc_results
SAMPLE_LIST=${BASEDIR}/sample_list.txt

mkdir -p ${QCDIR}

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

## Set array-specific variables
READS1_TRIMMED=${SAMPLE}_R1_trimmed.fastq.gz
READS2_TRIMMED=${SAMPLE}_R2_trimmed.fastq.gz

## Confirm that variables are properly assigned
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "SAMPLE NAME: ${SAMPLE}"
echo "TRIMMED FORWARD READS: ${READS1_TRIMMED}"
echo "TRIMMED REVERSE READS: ${READS2_TRIMMED}"
echo "Starting fastqc for ${SAMPLE}..."

## Check if input files exist
if [[ ! -f "${WORKDIR}/${READS1_TRIMMED}" || ! -f "${WORKDIR}/${READS2_TRIMMED}" ]]; then
    echo "ERROR: Trimmed files not found for ${SAMPLE}"
    echo "Looking for: ${WORKDIR}/${READS1_TRIMMED}"
    echo "Looking for: ${WORKDIR}/${READS2_TRIMMED}"
    exit 1
fi

## Load modules
module purge
module load fastqc/0.12

## Move into working directory
cd ${WORKDIR}

## Run fastqc on individual files
echo "Running FastQC on R1..."
fastqc --threads 4 --outdir=${QCDIR} ${READS1_TRIMMED}

echo "Running FastQC on R2..."
fastqc --threads 4 --outdir=${QCDIR} ${READS2_TRIMMED}

## Verify output files were created
R1_HTML="${QCDIR}/${SAMPLE}_R1_trimmed_fastqc.html"
R2_HTML="${QCDIR}/${SAMPLE}_R2_trimmed_fastqc.html"

if [[ -f "${R1_HTML}" && -f "${R2_HTML}" ]]; then
    echo "✓ FastQC completed successfully for ${SAMPLE}"
    echo "  R1 report: ${R1_HTML}"
    echo "  R2 report: ${R2_HTML}"
else
    echo "✗ Error: FastQC output files not created for ${SAMPLE}"
    exit 1
fi

echo "Completed at: $(date)"
