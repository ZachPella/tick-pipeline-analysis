#!/bin/bash
#SBATCH --job-name=gather_vcfs_ticks
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --time=1-00:00:00
#SBATCH --error=%x_%j.err
#SBATCH --output=%x_%j.out

module purge
module load gatk4/4.6

# Set up directories
INPUT_DIR=/work/fauverlab/zachpella/scripts_ticksJune2025_10_scatter/genotyping/genotyped_vcfs
FINAL_OUTPUT_DIR=/work/fauverlab/zachpella/scripts_ticksJune2025_10_scatter/genotyping
FINAL_VCF="cohort_ticks_june2025_final.vcf.gz"

# Create a list of input VCFs
INPUT_VCFS=()
for file in "${INPUT_DIR}"/chunk_*.vcf.gz; do
    INPUT_VCFS+=("-I" "${file}")
done

# Run GatherVcfs
echo "Starting GatherVcfs to combine all genotyped VCFs..."
gatk --java-options "-Xms10G -Xmx50G" GatherVcfs \
    "${INPUT_VCFS[@]}" \
    -O "${FINAL_OUTPUT_DIR}/${FINAL_VCF}"

echo "GatherVcfs completed. Final VCF is located at ${FINAL_OUTPUT_DIR}/${FINAL_VCF}"
