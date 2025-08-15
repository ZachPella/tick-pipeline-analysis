#!/bin/bash
#SBATCH --job-name=comprehensive_stats_summary
#SBATCH --mail-user=zpella@unmc.edu
#SBATCH --mail-type=ALL
#SBATCH --time=0-01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --partition=batch

## Set working directories to match g3 script paths
BASEDIR=/work/fauverlab/zachpella/scripts_ticksJune2025_10_scatter
ALIGNDIR=${BASEDIR}/bam_files
FASTPDIR=${BASEDIR}/err_and_out
OUTDIR=${BASEDIR}

## Set reference and targets files
REFERENCEDIR=/work/fauverlab/zachpella/practice_pop_gen/reference
REFERENCE=masked_ixodes_ref_genome.fasta
TARGETS=${REFERENCEDIR}/${REFERENCE}.bed

## Load samtools module
module load samtools/1.19

## Create a header for the comprehensive summary file
echo "Sample,Total_Reads_Before_Fastp,Total_Reads_After_Fastp,Percent_Kept,Duplication_Rate_Fastp,Primary_Reads,Primary_Mapped,Perc_Primary_Mapped,Properly_Paired,Perc_Properly_Paired,Mean_Coverage,Coverage_5x_Percent,Coverage_10x_Percent" > ${OUTDIR}/comprehensive_stats_summary.csv

## First, create a dictionary of fastp metrics
declare -A fastp_before
declare -A fastp_after
declare -A fastp_dup_rate

echo "Processing fastp metrics..."
for i in ${FASTPDIR}/fastp_ticks_processing_*.err; do
    SAMPLE=$(grep "HTML report:" $i | sed 's/.*HTML report: \(.*\)\.fastp\.html.*/\1/')
    READ_BEFORE=$(grep -A 1 "Read1 before filtering:" $i | grep "total reads:" | awk '{print $3}')
    READ_AFTER=$(grep -A 1 "Read1 after filtering:" $i | grep "total reads:" | awk '{print $3}')
    DUP_RATE=$(grep "Duplication rate:" $i | awk '{print $3}' | sed 's/%//')

    # Double the read counts for paired-end data
    if [ -n "$READ_BEFORE" ]; then
        READ_BEFORE=$((READ_BEFORE * 2))
    fi
    if [ -n "$READ_AFTER" ]; then
        READ_AFTER=$((READ_AFTER * 2))
    fi

    # Store in associative arrays
    fastp_before["$SAMPLE"]=$READ_BEFORE
    fastp_after["$SAMPLE"]=$READ_AFTER
    fastp_dup_rate["$SAMPLE"]=$DUP_RATE
done


## Now process flagstat files and combine with fastp data and coverage stats
echo "Processing alignment statistics and creating comprehensive summary..."

cd ${ALIGNDIR} 

for i in flagstats.*.out; do
    BASENAME=$(basename $i .out | sed 's/flagstats\.//')
    echo "Processing $BASENAME..."

    # Extract alignment metrics
    TOTAL=$(grep "in total" $i | head -1 | awk '{print $1}')
    PRIMARY=$(grep "^[0-9]* + [0-9]* primary$" $i | awk '{print $1}')

    # Get primary mapped stats specifically
    PRIMARY_MAPPED_LINE=$(grep "primary mapped" $i)
    PRIMARY_MAPPED=$(echo "$PRIMARY_MAPPED_LINE" | awk '{print $1}')
    PERC_PRIMARY_MAPPED=$(echo "$PRIMARY_MAPPED_LINE" | awk -F'[()]' '{print $2}' | awk '{print $1}')

    # Get properly paired stats
    PAIRED_LINE=$(grep "properly paired" $i)
    PROPERLY_PAIRED=$(echo "$PAIRED_LINE" | awk '{print $1}')
    PERC_PROPERLY_PAIRED=$(echo "$PAIRED_LINE" | awk -F'[()]' '{print $2}' | awk '{print $1}')

    # Handle special case for NTC or any file with 0 reads
    if [ -z "$PERC_PRIMARY_MAPPED" ]; then PERC_PRIMARY_MAPPED="0"; fi
    if [ -z "$PERC_PROPERLY_PAIRED" ]; then PERC_PROPERLY_PAIRED="0"; fi
    if [ -z "$PROPERLY_PAIRED" ]; then PROPERLY_PAIRED="0"; fi
    if [ -z "$PRIMARY_MAPPED" ]; then PRIMARY_MAPPED="0"; fi

    # Get fastp data for this sample
    READ_BEFORE=${fastp_before["$BASENAME"]}
    READ_AFTER=${fastp_after["$BASENAME"]}
    DUP_RATE=${fastp_dup_rate["$BASENAME"]}

    # Get mean coverage from averageDOC file
    AVGDOC_FILE="${ALIGNDIR}/averageDOC.${BASENAME}.out"
    if [ -f "${AVGDOC_FILE}" ]; then
        # Handle both formats:
        # Format 1: "Average depth of coverage: 8.03604" (field $5)
        # Format 2: "Average depth of coverage for SampleName: 6.73386" (field $7)

        # Try format 1 first (original format)
        MEAN_COVERAGE=$(grep "Average depth of coverage:" "${AVGDOC_FILE}" | awk '{print $5}')

        # If that didn't work, try format 2 (reprocessed format)
        if [ -z "$MEAN_COVERAGE" ]; then
            MEAN_COVERAGE=$(grep "Average depth of coverage for" "${AVGDOC_FILE}" | awk '{print $7}')
        fi

        # Check if we got a valid number and format it
        if [ -n "$MEAN_COVERAGE" ] && [[ "$MEAN_COVERAGE" =~ ^[0-9]+\.?[0-9]*$ ]]; then
            MEAN_COVERAGE=$(awk "BEGIN {printf \"%.2f\", ${MEAN_COVERAGE}}")
        else
            MEAN_COVERAGE="NA"
        fi
    else
        MEAN_COVERAGE="NA"
    fi

    # Get coverage depth percentages for 5x and 10x from NEW samtools stats files
    COVERAGE_5X="NA"
    COVERAGE_10X="NA"

    # Extract from stats.5x.${SAMPLE}.out file (new method)
    STATS_5X_FILE="${ALIGNDIR}/stats.5x.${BASENAME}.out"
    if [ -f "${STATS_5X_FILE}" ]; then
        COVERAGE_5X=$(grep "percentage of target genome with coverage >" "${STATS_5X_FILE}" | awk '{print $NF}' | sed 's/%//')
    fi

    # Extract from stats.10x.${SAMPLE}.out file (new method)
    STATS_10X_FILE="${ALIGNDIR}/stats.10x.${BASENAME}.out"
    if [ -f "${STATS_10X_FILE}" ]; then
        COVERAGE_10X=$(grep "percentage of target genome with coverage >" "${STATS_10X_FILE}" | awk '{print $NF}' | sed 's/%//')
    fi

    # Ensure we have valid values
    if [ -z "$COVERAGE_5X" ] || [ "$COVERAGE_5X" = "" ]; then COVERAGE_5X="NA"; fi
    if [ -z "$COVERAGE_10X" ] || [ "$COVERAGE_10X" = "" ]; then COVERAGE_10X="NA"; fi

    # Calculate percent kept
    if [ -n "$READ_BEFORE" ] && [ -n "$READ_AFTER" ] && [ "$READ_BEFORE" -ne 0 ]; then
        PERCENT_KEPT=$(awk "BEGIN {printf \"%.2f\", ($READ_AFTER/$READ_BEFORE)*100}")
    else
        # If we don't have fastp data, use 0 as placeholders
        READ_BEFORE="NA"
        READ_AFTER="NA"
        PERCENT_KEPT="NA"
        DUP_RATE="NA"
    fi

    # Output to comprehensive summary
    echo "$BASENAME,$READ_BEFORE,$READ_AFTER,$PERCENT_KEPT,$DUP_RATE,$PRIMARY,$PRIMARY_MAPPED,$PERC_PRIMARY_MAPPED,$PROPERLY_PAIRED,$PERC_PROPERLY_PAIRED,$MEAN_COVERAGE,$COVERAGE_5X,$COVERAGE_10X" >> ${OUTDIR}/comprehensive_stats_summary.csv
done

echo "Comprehensive summary created: ${OUTDIR}/comprehensive_stats_summary.csv"

## Create a tab-separated version for command line viewing
sed 's/,/\t/g' ${OUTDIR}/comprehensive_stats_summary.csv > ${OUTDIR}/comprehensive_stats_summary.tsv
echo "Tab-separated version created: ${OUTDIR}/comprehensive_stats_summary.tsv"
