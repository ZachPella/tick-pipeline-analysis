#!/bin/bash
#SBATCH --job-name=g3_sam_to_bam
#SBATCH --mail-user=zpella@unmc.edu
#SBATCH --mail-type=ALL
#SBATCH --time=4-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
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
SAMDIR=${BASEDIR}/sam_files  # SAM files from BWA (updated directory)
WORKDIR=${BASEDIR}/bam_files    # Output directory for BAM and stats
REFERENCEDIR=/work/fauverlab/zachpella/practice_pop_gen/reference
REFERENCE=masked_ixodes_ref_genome.fasta
TARGETS=${REFERENCEDIR}/${REFERENCE}.bed
SAMPLE_LIST=${BASEDIR}/sample_list.txt

# Create working directory if it doesn't exist
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

## Set SAM file path
SAMPLE_FILE="${SAMDIR}/${SAMPLE}.sam"

## Confirm that variables are properly assigned
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Sample: ${SAMPLE}"
echo "SAM file: ${SAMPLE_FILE}"
echo "Starting SAM-to-BAM conversion and stats summaries for ${SAMPLE}..."

## Check if SAM file exists
if [ ! -f "${SAMPLE_FILE}" ]; then
    echo "Error: SAM file not found: ${SAMPLE_FILE}"
    exit 1
fi

## Load modules
module purge
module load samtools/1.19

## Move into working directory
cd ${WORKDIR}

## Compress SAM to BAM and sort
echo "Converting SAM to BAM..."
samtools view -b "${SAMPLE_FILE}" > "${WORKDIR}/${SAMPLE}.bam"

echo "Sorting BAM file..."
samtools sort -m 8G "${WORKDIR}/${SAMPLE}.bam" -o "${WORKDIR}/${SAMPLE}.sorted.bam" -T "${WORKDIR}/${SAMPLE}.reads.tmp"

# Only remove SAM file if BAM conversion was successful
if [ -s "${WORKDIR}/${SAMPLE}.sorted.bam" ]; then
    rm "${SAMPLE_FILE}"
    rm "${WORKDIR}/${SAMPLE}.bam"
    echo "✓ SAM-to-BAM conversion successful for ${SAMPLE}"
else
    echo "Error: BAM file is empty or conversion failed"
    exit 1
fi

## Generate basic samtools summary stats
echo "Generating statistics..."
samtools flagstat "${WORKDIR}/${SAMPLE}.sorted.bam" > "${WORKDIR}/flagstats.${SAMPLE}.out"

# Run samtools stats for coverage thresholds 5 and 10
for cov in 5 10; do
    samtools stats -t "${TARGETS}" --cov-threshold ${cov} "${WORKDIR}/${SAMPLE}.sorted.bam" > "${WORKDIR}/stats.${cov}x.${SAMPLE}.out"
done

# Calculate average depth of coverage to reference
echo "Calculating depth of coverage..."
samtools depth -a "${WORKDIR}/${SAMPLE}.sorted.bam" > "${WORKDIR}/${SAMPLE}.depth"
AVGDOC=$(awk '{ total += $3; count++ } END { print total/count }' "${WORKDIR}/${SAMPLE}.depth")
echo "Average depth of coverage: ${AVGDOC}" > "${WORKDIR}/averageDOC.${SAMPLE}.out"

# Calculate by-contig coverage stats
samtools coverage -o "${WORKDIR}/coverage.${SAMPLE}.out" "${WORKDIR}/${SAMPLE}.sorted.bam"
samtools coverage --plot-depth -o "${WORKDIR}/hist.coverage.${SAMPLE}.out" "${WORKDIR}/${SAMPLE}.sorted.bam"

# Count number of reads that align to reference
samtools view -c -F 4 -F 2048 "${WORKDIR}/${SAMPLE}.sorted.bam" > "${WORKDIR}/countmappedreads.${SAMPLE}.out"

# Generate a sorted BAM file with only mapped reads
samtools view -b -F 4 "${WORKDIR}/${SAMPLE}.sorted.bam" > "${WORKDIR}/${SAMPLE}.sorted.mapped.bam"

# Index the final BAM files
samtools index "${WORKDIR}/${SAMPLE}.sorted.bam"
samtools index "${WORKDIR}/${SAMPLE}.sorted.mapped.bam"

## Print confirmation statement
echo "✓ Processing for ${SAMPLE} completed successfully."
echo "Output files are in: ${WORKDIR}"
printf "\n"
