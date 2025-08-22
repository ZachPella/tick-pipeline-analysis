#!/bin/bash
#SBATCH --job-name=g2_make_sam
#SBATCH --mail-user=zpella@unmc.edu
#SBATCH --mail-type=ALL
#SBATCH --time=6-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=45G
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
WORKDIR=${BASEDIR}/sam_files
READSDIR=${BASEDIR}/trimmed_reads
REFERENCEDIR=/work/fauverlab/zachpella/practice_pop_gen/reference
REFERENCE=masked_ixodes_ref_genome.fasta
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

## Set file variables
READS1_TRIMMED=${SAMPLE}_R1_trimmed.fastq.gz
READS2_TRIMMED=${SAMPLE}_R2_trimmed.fastq.gz

## Confirm variables are assigned correctly
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "SAMPLE NAME: ${SAMPLE}"
echo "TRIMMED FORWARD READS: ${READS1_TRIMMED}"
echo "TRIMMED REVERSE READS: ${READS2_TRIMMED}"
echo "REFERENCE: ${REFERENCE}"
echo "Starting bwa mem for ${SAMPLE}..."
printf "\n"

## Check if input files exist
if [ ! -f "${READSDIR}/${READS1_TRIMMED}" ] || [ ! -f "${READSDIR}/${READS2_TRIMMED}" ]; then
    echo "Error: Input read files not found:"
    echo "  ${READSDIR}/${READS1_TRIMMED}"
    echo "  ${READSDIR}/${READS2_TRIMMED}"
    exit 1
fi

## Check if reference exists
if [ ! -f "${REFERENCEDIR}/${REFERENCE}" ]; then
    echo "Error: Reference file not found: ${REFERENCEDIR}/${REFERENCE}"
    exit 1
fi

## Load modules
module purge
module load bwa

## Move into working directory
cd ${WORKDIR}

## Align short reads to the reference assembly with BWA MEM
echo "Running BWA MEM alignment..."
bwa mem \
    -t 4 \
    -M \
    ${REFERENCEDIR}/${REFERENCE} \
    ${READSDIR}/${READS1_TRIMMED} ${READSDIR}/${READS2_TRIMMED} \
    -o ${WORKDIR}/${SAMPLE}.sam

## Verify output file was created
if [[ -f "${WORKDIR}/${SAMPLE}.sam" && -s "${WORKDIR}/${SAMPLE}.sam" ]]; then
    echo "✓ BWA alignment completed successfully for ${SAMPLE}"
    echo "  Output SAM: ${WORKDIR}/${SAMPLE}.sam"
    echo "  File size: $(ls -lh ${WORKDIR}/${SAMPLE}.sam | awk '{print $5}')"
else
    echo "✗ Error: SAM file not created or is empty for ${SAMPLE}"
    exit 1
fi

echo "Completed at: $(date)"
