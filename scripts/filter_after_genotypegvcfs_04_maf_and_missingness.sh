#!/bin/bash
#SBATCH --job-name=f1_VCFtools_filter
#SBATCH --time=2-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --mem=64G
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
# Assuming a consistent base directory for your entire workflow
BASEDIR=/work/fauverlab/zachpella/scripts_ticksJune2025_10_scatter
WORKDIR=${BASEDIR}/genotyping
INPUT_VCF=${WORKDIR}/cohort_ticks_june2025_snps_passing_only.vcf.gz
OUTPUT_VCF_PREFIX=cohort_ticks_june2025_snps_passing_only

## Move into working directory
cd ${WORKDIR}

## Load modules
module purge
module load vcftools

## Check if the input VCF exists
if [ ! -f "${INPUT_VCF}" ]; then
    echo "Error: Input VCF file not found: ${INPUT_VCF}"
    exit 1
fi

echo "Filtering input VCF: ${INPUT_VCF}"

## Filter joint VCF by minor allele frequency & site-level missingness
vcftools --gzvcf "${INPUT_VCF}" \
    --maf 0.05 \
    --max-missing 0.7 \
    --recode \
    --recode-INFO-all \
    --stdout | gzip -c > ${OUTPUT_VCF_PREFIX}.maf005.miss03.vcf.gz

echo "VCFtools filtering completed."
echo "Filtered VCF saved to: ${WORKDIR}/${OUTPUT_VCF_PREFIX}.maf005.miss03.vcf.gz"
echo "Completed at: $(date)"
printf "\n"
