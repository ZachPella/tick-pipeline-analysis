#!/bin/bash
#SBATCH --job-name=pop_structure
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=90G
#SBATCH --partition=batch

## Record relevant job info
START_DIR=$(pwd)
HOST_NAME=$(hostname)
RUN_DATE=$(date)
echo "Starting working directory: ${START_DIR}"
echo "Host name: ${HOST_NAME}"
echo "Run date: ${RUN_DATE}"
printf "\n"

## Set working directory and variables - UPDATED FOR NEW PIPELINE
BASEDIR=/work/fauverlab/zachpella/scripts_ticksJune2025_10_scatter
WORKDIR=${BASEDIR}
INPUT_VCF=${BASEDIR}/genotyping/cohort_ticks_june2025_snps_passing_only.maf005.miss03.vcf.gz
cd ${WORKDIR}

## Create output directory
mkdir -p population_structure

## Check if input VCF exists
if [ ! -f "${INPUT_VCF}" ]; then
    echo "Error: Input VCF file not found: ${INPUT_VCF}"
    exit 1
fi

echo "Using input VCF: ${INPUT_VCF}"
echo "Output directory: ${WORKDIR}/population_structure"

## Load modules
module purge
module load plink/1.90  # Load PLINK1.9
module load plink2      # Load PLINK2

## Step 1: Convert VCF to PLINK format using PLINK1.9 and filter for biallelic variants
echo "Converting VCF to PLINK format and filtering for biallelic variants..."
plink --vcf ${INPUT_VCF} \
      --const-fid 0\
      --make-bed \
      --out population_structure/tick_population \
      --allow-extra-chr \
      --biallelic-only strict

## Verify PLINK files were created
if [ ! -f "population_structure/tick_population.bed" ]; then
    echo "Error: PLINK bed file not created"
    exit 1
fi

## Step 2: Perform linkage disequilibrium (LD) pruning
echo "Performing linkage disequilibrium pruning..."
plink2 --bfile population_structure/tick_population \
       --indep-pairwise 50 10 0.1 \
       --out population_structure/tick_population_ld \
       --allow-extra-chr

## Verify LD pruning output was created
if [ ! -f "population_structure/tick_population_ld.prune.in" ]; then
    echo "Error: LD pruning output files not created"
    exit 1
fi

## Step 3: Extract LD-pruned variants
echo "Extracting LD-pruned variants..."
plink2 --bfile population_structure/tick_population \
       --extract population_structure/tick_population_ld.prune.in \
       --make-bed \
       --out population_structure/tick_population_pruned \
       --allow-extra-chr

## Verify pruned PLINK files were created
if [ ! -f "population_structure/tick_population_pruned.bed" ]; then
    echo "Error: Pruned PLINK bed file not created"
    exit 1
fi

## Step 4: Perform PCA using PLINK2 on the LD-pruned dataset
echo "Performing PCA analysis on LD-pruned data..."
plink2 --bfile population_structure/tick_population_pruned \
       --pca 20 \
       --out population_structure/tick_pca_pruned \
       --allow-extra-chr

## Verify PCA output was created
if [ ! -f "population_structure/tick_pca_pruned.eigenvec" ]; then
    echo "Error: PCA output files not created"
    exit 1
fi

echo "✓ Population structure analysis with LD pruning completed successfully"
echo "  Original PLINK files: population_structure/tick_population.*"
echo "  LD pruning results: population_structure/tick_population_ld.*"
echo "  Pruned PLINK files: population_structure/tick_population_pruned.*"
echo "  PCA results: population_structure/tick_pca_pruned.*"
echo "Completed at: $(date)"
printf "\n"
