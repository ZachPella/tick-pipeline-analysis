#!/bin/bash
#SBATCH --job-name=plink_10scatter_matching_5scatter_approach
#SBATCH --time=0-10:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --partition=batch

## record relevant job info
START_DIR=$(pwd)
HOST_NAME=$(hostname)
RUN_DATE=$(date)
echo "Starting working directory: ${START_DIR}"
echo "Host name: ${HOST_NAME}"
echo "Run date: ${RUN_DATE}"
printf "\n"

## set working directory and variables - UPDATED FOR 10-SCATTER DATA
BASEDIR=/work/fauverlab/zachpella/scripts_ticksJune2025_10_scatter
WORKDIR=${BASEDIR}/genotyping
JOINTVCF=cohort_ticks_june2025_snps_passing_only.maf005.miss03

## load modules
module purge
module load plink2

## move into VCF directory
cd ${WORKDIR}

## filter for linkage disequilibrium
plink2 \
        --vcf ${JOINTVCF}.vcf.gz \
        --double-id \
        --allow-extra-chr \
        --set-missing-var-ids @:# \
        --indep-pairwise 50 10 0.1 \
        --out ${JOINTVCF}

## prune and create pca
plink2 \
        --vcf ${JOINTVCF}.vcf.gz \
        --double-id \
        --allow-extra-chr \
        --set-missing-var-ids @:# \
        --extract ${JOINTVCF}.prune.in \
        --make-bed \
        --pca \
        --out ${JOINTVCF}

echo "✓ PLINK analysis completed using identical approach to 5-scatter"
echo "  LD pruning results: ${JOINTVCF}.prune.in/out"
echo "  PCA results: ${JOINTVCF}.eigenvec/eigenval"
echo "  PLINK files: ${JOINTVCF}.bed/bim/fam"
echo "Completed at: $(date)"
printf "\n"

## Expected input: cohort_ticks_june2025_snps_passing_only.maf005.miss03.vcf.gz
