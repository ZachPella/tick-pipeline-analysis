#!/bin/bash
#SBATCH --job-name=concatenate_lanes
#SBATCH --time=0-06:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1-62
#SBATCH --partition=batch

# Set directories
BASE_DIR="/work/fauverlab/shared/iscapularis_NovaSeq_S2/250610"
OUTPUT_DIR="/work/fauverlab/zachpella/scripts_ticksJune2025_10_scatter/concatenated_fastq"
SAMPLE_LIST="/work/fauverlab/zachpella/scripts_ticksJune2025_10_scatter/sample_list.txt"

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Check if sample list exists
if [ ! -f "$SAMPLE_LIST" ]; then
    echo "Error: Sample list file not found: $SAMPLE_LIST"
    exit 1
fi

# Get total number of samples
TOTAL_SAMPLES=$(wc -l < "$SAMPLE_LIST")

# Check if array task ID is valid
if [ ${SLURM_ARRAY_TASK_ID} -gt ${TOTAL_SAMPLES} ]; then
    echo "Error: Array task ID ${SLURM_ARRAY_TASK_ID} exceeds number of samples (${TOTAL_SAMPLES})"
    exit 1
fi

# Get sample name from list - this is bulletproof!
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
SAMPLE_DIR="${BASE_DIR}/${SAMPLE_NAME}"

# Verify sample name is not empty
if [[ -z "$SAMPLE_NAME" ]]; then
    echo "Error: Empty sample name for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Processing sample: ${SAMPLE_NAME}"
echo "Sample directory: ${SAMPLE_DIR}"
echo "Started at: $(date)"

# Check for required files
L001_R1=(${SAMPLE_DIR}/*_L001_R1_*.fastq.gz)
L001_R2=(${SAMPLE_DIR}/*_L001_R2_*.fastq.gz)
L002_R1=(${SAMPLE_DIR}/*_L002_R1_*.fastq.gz)
L002_R2=(${SAMPLE_DIR}/*_L002_R2_*.fastq.gz)

# Verify files exist and are readable
if [[ -f "${L001_R1[0]}" && -f "${L001_R2[0]}" && -f "${L002_R1[0]}" && -f "${L002_R2[0]}" ]]; then
    echo "Found all required files:"
    echo "  L001_R1: ${L001_R1[0]}"
    echo "  L001_R2: ${L001_R2[0]}"
    echo "  L002_R1: ${L002_R1[0]}"
    echo "  L002_R2: ${L002_R2[0]}"
    
    echo "Concatenating R1 files..."
    cat "${L001_R1[0]}" "${L002_R1[0]}" > "${OUTPUT_DIR}/${SAMPLE_NAME}_R1_merged.fastq.gz"
    
    echo "Concatenating R2 files..."
    cat "${L001_R2[0]}" "${L002_R2[0]}" > "${OUTPUT_DIR}/${SAMPLE_NAME}_R2_merged.fastq.gz"
    
    # Verify output files were created successfully
    if [[ -f "${OUTPUT_DIR}/${SAMPLE_NAME}_R1_merged.fastq.gz" && -f "${OUTPUT_DIR}/${SAMPLE_NAME}_R2_merged.fastq.gz" ]]; then
        echo "✓ Successfully concatenated ${SAMPLE_NAME}"
        echo "  Output R1: ${OUTPUT_DIR}/${SAMPLE_NAME}_R1_merged.fastq.gz"
        echo "  Output R2: ${OUTPUT_DIR}/${SAMPLE_NAME}_R2_merged.fastq.gz"
    else
        echo "✗ Error: Output files not created for ${SAMPLE_NAME}"
        exit 1
    fi
else
    echo "✗ Missing required files for ${SAMPLE_NAME}"
    echo "Looking for files in: ${SAMPLE_DIR}"
    ls -la "${SAMPLE_DIR}" 2>/dev/null || echo "Directory does not exist"
    exit 1
fi

echo "Completed at: $(date)"
