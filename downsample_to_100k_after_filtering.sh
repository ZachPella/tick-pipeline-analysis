#!/bin/bash
#SBATCH --job-name=b4_joint_vcf_generation_prelim_iscap_popgen
#SBATCH --time=3-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --partition=batch

START_DIR=$(pwd)
HOST_NAME=$(hostname)
RUN_DATE=$(date)
echo "Starting working directory: ${START_DIR}"
echo "Host name: ${HOST_NAME}"
echo "Run date: ${RUN_DATE}"
printf "\n"

WORKDIR=/work/fauverlab/zachpella/scripts_ticksJune2025_10_scatter/genotyping
MERGEDVCFNAME=cohort_ticks_june2025_snps_filtered_only

module purge
module load bcftools
module load vcflib

cd ${WORKDIR}

echo "Downsampling joint VCF to ~100K variants..."
bcftools view ${MERGEDVCFNAME}.vcf.gz | vcfrandomsample -r 0.00037 > 100Ksubset_after.${MERGEDVCFNAME}.vcf

echo "Gzipping & indexing downsampled VCF..."
bgzip 100Ksubset_after.${MERGEDVCFNAME}.vcf
bcftools index 100Ksubset_after.${MERGEDVCFNAME}.vcf.gz

echo "No. sites in downsampled VCF:"
bcftools view -H 100Ksubset_after.${MERGEDVCFNAME}.vcf.gz | wc -l
