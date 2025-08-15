#!/bin/bash
#SBATCH --job-name=vcf_stats_zachpella
#SBATCH --output=vcf_stats_%j.out
#SBATCH --error=vcf_stats_%j.err
#SBATCH --time=10:00:00
#SBATCH --mem=35G
#SBATCH --nodes=1
#SBATCH --ntasks=1

WORKDIR=/work/fauverlab/zachpella/scripts_ticksJune2025/genotyping
SUBSETVCF=100Ksubset_before.cohort_ticks_june2025_combined

cd ${WORKDIR}

module load vcftools

if [ ! -f "${SUBSETVCF}.vcf.gz" ]; then
    echo "Error: VCF file not found at ${WORKDIR}/${SUBSETVCF}.vcf.gz."
    exit 1
fi

vcftools --gzvcf ${SUBSETVCF}.vcf.gz --freq2 --out out.${SUBSETVCF} --max-alleles 2
vcftools --gzvcf ${SUBSETVCF}.vcf.gz --depth --out out.${SUBSETVCF}
vcftools --gzvcf ${SUBSETVCF}.vcf.gz --site-mean-depth --out out.${SUBSETVCF}
vcftools --gzvcf ${SUBSETVCF}.vcf.gz --site-quality --out out.${SUBSETVCF}
vcftools --gzvcf ${SUBSETVCF}.vcf.gz --missing-indv --out out.${SUBSETVCF}
vcftools --gzvcf ${SUBSETVCF}.vcf.gz --missing-site --out out.${SUBSETVCF}
vcftools --gzvcf ${SUBSETVCF}.vcf.gz --het --out out.${SUBSETVCF}
