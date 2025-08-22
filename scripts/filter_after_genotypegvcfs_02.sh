#!/bin/bash
#SBATCH --job-name=variant_filtration
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1  # VariantFiltration can benefit from more cores for larger VCFs, but 1 is a safe start
#SBATCH --mem=40G          # Adjust as needed, depends on size of input VCF
#SBATCH --time=0-4:00:00   # Sufficient time, adjust if very large VCF
#SBATCH --error=%x_%j.err
#SBATCH --output=%x_%j.out

module purge
module load gatk4/4.6

# --- Inherit from the previous script's output ---
# Define your input and output paths based on the previous script's output
BASEDIR=/work/fauverlab/zachpella/scripts_ticksJune2025_10_scatter
INPUT_SNP_VCF="${BASEDIR}/genotyping/cohort_ticks_june2025_snps_only.vcf.gz"
OUTPUT_FILTERED_VCF="${BASEDIR}/genotyping/cohort_ticks_june2025_snps_filtered.vcf.gz"

echo "Starting VariantFiltration for SNPs..."

gatk --java-options "-Xms4G -Xmx35G" VariantFiltration \
    --variant "${INPUT_SNP_VCF}" \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
    --filter-expression "SOR > 3.0" --filter-name "SOR3" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    --output "${OUTPUT_FILTERED_VCF}"

echo "VariantFiltration completed. Filtered VCF: ${OUTPUT_FILTERED_VCF}"
