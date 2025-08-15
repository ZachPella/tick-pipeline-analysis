#!/bin/bash
#SBATCH --job-name=genotype_gvcfs_ticks_na
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=115G
#SBATCH --time=4-00:00:00
#SBATCH --array=0-9
#SBATCH --error=%x_%A_%a.err
#SBATCH --output=%x_%A_%a.out

module purge
module load gatk4/4.6

# Set up working directories and variables
START_DIR=$(pwd)
BASEDIR=/work/fauverlab/zachpella/scripts_ticksJune2025_10_scatter
WORKDIR=${BASEDIR}/genotyping
REFERENCEDIR=/work/fauverlab/zachpella/practice_pop_gen/reference
REFERENCE=masked_ixodes_ref_genome.fasta
GENOMICSDB_DIR=${WORKDIR}/genomicsdb_chunks

# The interval list and GenomicsDB workspace for this specific job will be selected based on the array task ID
CHUNK=$(printf "%04d" ${SLURM_ARRAY_TASK_ID})
INTERVAL_LIST_CHUNK="${WORKDIR}/scattered_intervals/$(sed -n "1p" ${BASEDIR}/sample_list.txt)/${CHUNK}-scattered.interval_list"
GENOMICSDB_CHUNK_PATH="${GENOMICSDB_DIR}/chunk_${CHUNK}"

# Set output paths
FINAL_OUTPUT_DIR=${WORKDIR}/genotyped_vcfs
OUTPUT_VCF_NAME="chunk_${CHUNK}_genotyped.vcf.gz"

# Create job-specific directories on scratch
mkdir -p /scratch/$SLURM_JOBID/tmp
mkdir -p /scratch/$SLURM_JOBID/output

# Copy reference and interval list to scratch
echo "Copying reference and interval list to scratch for chunk ${SLURM_ARRAY_TASK_ID}..."
cp "${REFERENCEDIR}/${REFERENCE}" /scratch/$SLURM_JOBID/
cp "${REFERENCEDIR}/${REFERENCE}.fai" /scratch/$SLURM_JOBID/
cp "${REFERENCEDIR}/${REFERENCE%.*}.dict" /scratch/$SLURM_JOBID/
cp "${INTERVAL_LIST_CHUNK}" /scratch/$SLURM_JOBID/
cp -r "${GENOMICSDB_CHUNK_PATH}" /scratch/$SLURM_JOBID/

# Run GenotypeGVCFs for the current chunk
echo "Starting GenotypeGVCFs for chunk ${SLURM_ARRAY_TASK_ID}..."
gatk --java-options "-Djava.io.tmpdir=/scratch/$SLURM_JOBID/tmp -Xms2G -Xmx90G -XX:ParallelGCThreads=20" \
    GenotypeGVCFs \
    -R /scratch/$SLURM_JOBID/${REFERENCE} \
    -V "gendb:///scratch/$SLURM_JOBID/chunk_${CHUNK}" \
    -O /scratch/$SLURM_JOBID/output/"${OUTPUT_VCF_NAME}" \
    -L /scratch/$SLURM_JOBID/$(basename "${INTERVAL_LIST_CHUNK}")

echo "Copying final VCF and index to ${FINAL_OUTPUT_DIR} for chunk ${SLURM_ARRAY_TASK_ID}..."
mkdir -p "${FINAL_OUTPUT_DIR}"
cp /scratch/$SLURM_JOBID/output/"${OUTPUT_VCF_NAME}" "${FINAL_OUTPUT_DIR}/"

echo "Cleaning up scratch directory: /scratch/${SLURM_JOBID}"
rm -rf /scratch/${SLURM_JOBID}

echo "Job finished for chunk ${SLURM_ARRAY_TASK_ID}."
