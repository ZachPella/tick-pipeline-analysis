#!/bin/bash
#SBATCH --job-name=fastp_ticks_processing
#SBATCH --mail-user=zpella@unmc.edu
#SBATCH --mail-type=ALL
#SBATCH --time=0-06:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --ntasks=1
#SBATCH --mem=15G
#SBATCH --array=1-62
#SBATCH --partition=batch

## Set up directories and variables
BASEDIR=/work/fauverlab/zachpella/scripts_ticksJune2025_10_scatter
WORKDIR=${BASEDIR}/trimmed_reads
READSDIR=${BASEDIR}/concatenated_fastq
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

## Set array-specific variables for concatenated files
READS1=${SAMPLE}_R1_merged.fastq.gz
READS2=${SAMPLE}_R2_merged.fastq.gz
READS1_TRIMMED=${SAMPLE}_R1_trimmed.fastq.gz
READS2_TRIMMED=${SAMPLE}_R2_trimmed.fastq.gz

## Confirm that variables are properly assigned
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "SAMPLE NAME: ${SAMPLE}"
echo "FORWARD READS: ${READS1}"
echo "REVERSE READS: ${READS2}"
echo "TRIMMED FORWARD READS TO BE NAMED: ${READS1_TRIMMED}"
echo "TRIMMED REVERSE READS TO BE NAMED: ${READS2_TRIMMED}"
echo "Starting fastp for ${SAMPLE}..."

## Check if input files exist
if [ ! -f "${READSDIR}/${READS1}" ] || [ ! -f "${READSDIR}/${READS2}" ]; then
    echo "Error: Input read files not found:"
    echo "  ${READSDIR}/${READS1}"
    echo "  ${READSDIR}/${READS2}"
    exit 1
fi

## Load modules
module purge
module load fastp

## Move into working directory
cd ${WORKDIR}

## Run fastp
fastp \
    --in1 ${READSDIR}/${READS1} \
    --in2 ${READSDIR}/${READS2} \
    --out1 ${WORKDIR}/${READS1_TRIMMED} \
    --out2 ${WORKDIR}/${READS2_TRIMMED} \
    -l 50 \
    -h ${SAMPLE}.fastp.html

## Verify output files were created
if [[ -f "${WORKDIR}/${READS1_TRIMMED}" && -f "${WORKDIR}/${READS2_TRIMMED}" ]]; then
    echo "✓ fastp completed successfully for ${SAMPLE}"
    echo "  Output R1: ${WORKDIR}/${READS1_TRIMMED}"
    echo "  Output R2: ${WORKDIR}/${READS2_TRIMMED}"
else
    echo "✗ Error: fastp output files not created for ${SAMPLE}"
    exit 1
fi
