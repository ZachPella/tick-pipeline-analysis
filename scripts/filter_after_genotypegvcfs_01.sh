#!/bin/bash
#SBATCH --job-name=select_snps
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1  # SelectVariants is not highly multi-threaded for simple selections
#SBATCH --mem=60G          # Adjust as needed, depends on size of input VCF
#SBATCH --time=0-10:00:00  # Sufficient time, adjust if very large VCF
#SBATCH --error=%x_%j.err
#SBATCH --output=%x_%j.out

module purge
module load gatk4/4.6

# --- Inherit from the previous script's output ---
# Define your input and output paths
BASEDIR=/work/fauverlab/zachpella/scripts_ticksJune2025_10_scatter
INPUT_VCF="${BASEDIR}/genotyping/cohort_ticks_june2025_final.vcf.gz"
OUTPUT_SNP_VCF="${BASEDIR}/genotyping/cohort_ticks_june2025_snps_only.vcf.gz"

echo "Starting to index the input VCF file..."

# Index the input VCF file
gatk --java-options "-Xms2G -Xmx15G" IndexFeatureFile \
    --input "${INPUT_VCF}"

echo "Indexing completed."

echo "Starting SelectVariants to extract SNPs..."

gatk --java-options "-Xms2G -Xmx35G" SelectVariants \
    --variant "${INPUT_VCF}" \
    --select-type-to-include SNP \
    --output "${OUTPUT_SNP_VCF}"

echo "SelectVariants completed. SNPs-only VCF: ${OUTPUT_SNP_VCF}"
