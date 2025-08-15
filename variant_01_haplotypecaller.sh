#!/bin/bash
#SBATCH --job-name=haplotype_scatter_gather
#SBATCH --time=7-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20 # 10 chunks * 2 cpus/chunk = 20 cpus total
#SBATCH --mem=110G # 10 chunks * 10 GB/chunk + 10 GB buffer = 110 GB total
#SBATCH --array=1-61
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
WORKDIR=${BASEDIR}/genotyping
BAMDIR=${BASEDIR}/dedup
REFERENCEDIR=/work/fauverlab/zachpella/practice_pop_gen/reference
REFERENCE=masked_ixodes_ref_genome.fasta
INTERVAL_LIST=${REFERENCEDIR}/${REFERENCE}.interval_list
SAMPLE_LIST=${BASEDIR}/sample_list.txt
SCATTER_COUNT=10
CPUS_PER_CHUNK=2

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
if [[-z "$SAMPLE"]]; then
    echo "Error: Empty sample name for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

## Create sample-specific directories
mkdir -p ${WORKDIR}/scattered_intervals/${SAMPLE}
mkdir -p ${WORKDIR}/scattered_gvcfs/${SAMPLE}

## Set file paths
BAM_FILE="${BAMDIR}/${SAMPLE}.dedup.rg.sorted.bam"

## Check if input BAM file exists
if [ ! -f "${BAM_FILE}" ]; then
    echo "Error: Input BAM file not found: ${BAM_FILE}"
    exit 1
fi

echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Starting Scatter-Gather HaplotypeCaller for ${SAMPLE}"
echo "Scatter count: ${SCATTER_COUNT}"
echo "BAM file: ${BAM_FILE}"

## Load modules
module purge
module load gatk4/4.6

cd ${WORKDIR}

## Step 1: Split intervals
echo "Splitting intervals into ${SCATTER_COUNT} chunks..."
gatk SplitIntervals \
    -R ${REFERENCEDIR}/${REFERENCE} \
    -L ${INTERVAL_LIST} \
    --scatter-count ${SCATTER_COUNT} \
    -O scattered_intervals/${SAMPLE}/

## Step 2: Run HaplotypeCaller chunks in parallel
echo "Running ${SCATTER_COUNT} parallel HaplotypeCaller jobs..."
for i in $(seq 0 $((SCATTER_COUNT-1))); do
    CHUNK=$(printf "%04d" $i)
    SCATTERED_INTERVAL="scattered_intervals/${SAMPLE}/${CHUNK}-scattered.interval_list"
    CHUNK_GVCF="scattered_gvcfs/${SAMPLE}/${SAMPLE}.${CHUNK}.g.vcf"
    gatk --java-options "-Xmx10G" HaplotypeCaller \
        -R ${REFERENCEDIR}/${REFERENCE} \
        -I ${BAM_FILE} \
        -L ${SCATTERED_INTERVAL} \
        -ploidy 2 \
        -O ${CHUNK_GVCF} \
        --ERC GVCF &
done

## Wait for all jobs to finish
wait

echo "✓ Scatter-Gather HaplotypeCaller completed successfully for ${SAMPLE}"
echo "  Individual GVCFs are located in: scattered_gvcfs/${SAMPLE}/"
echo "Completed at: $(date)"
printf "\n"
