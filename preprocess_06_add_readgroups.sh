#!/bin/bash
#SBATCH --job-name=add_read_groups
#SBATCH --mail-user=zpella@unmc.edu
#SBATCH --mail-type=ALL
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
INPUTDIR=${BASEDIR}/bam_files
WORKDIR=${BASEDIR}/readgroups
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

## Set file paths - Updated to match actual file structure
INPUT_BAM=${INPUTDIR}/${SAMPLE}.sorted.bam
OUTPUT_BAM=${WORKDIR}/${SAMPLE}.rg.sorted.bam


## Load modules
module purge
module load picard
module load samtools/1.19

## Move to working directory
cd ${WORKDIR}

## Add read groups using Picard
echo "Adding read groups to ${INPUT_BAM}..."
picard AddOrReplaceReadGroups \
    I=${INPUT_BAM} \
    O=${OUTPUT_BAM} \
    RGID=${SAMPLE} \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=H5VFNDMX2.${SAMPLE} \
    RGSM=${SAMPLE} \
    SORT_ORDER=coordinate


