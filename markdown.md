# Tick Population Genomics Pipeline - Comprehensive Analysis

## Table of Contents
- [Pipeline Overview](#pipeline-overview)
- [Preprocessing Steps (1-8)](#preprocessing-steps-1-8)
  - [Step 1: Lane Concatenation](#step-1-lane-concatenation)
  - [Step 2: Quality Control and Read Trimming](#step-2-quality-control-and-read-trimming)
  - [Step 3: Post-Cleaning Quality Assessment](#step-3-post-cleaning-quality-assessment)
  - [Step 4: Genome Alignment](#step-4-genome-alignment)
  - [Step 5: SAM to BAM Conversion and Sorting](#step-5-sam-to-bam-conversion-and-sorting)
  - [Step 6: Read Group Addition](#step-6-read-group-addition)
  - [Step 7: Duplicate Removal](#step-7-duplicate-removal)
  - [Step 8: Comprehensive Alignment Statistics Summary](#step-8-comprehensive-alignment-statistics-summary)
- [GATK Variant Calling Steps (9-11)](#gatk-variant-calling-steps-9-11)
  - [Step 9: Individual Sample Variant Calling](#step-9-individual-sample-variant-calling)
  - [Step 10: GenomicsDB Consolidation](#step-10-genomicsdb-consolidation)
  - [Step 11: Joint Genotyping and VCF Assembly](#step-11-joint-genotyping-and-vcf-assembly)
- [Post-VCF Filtering and Analysis Steps (12-17)](#post-vcf-filtering-and-analysis-steps-12-17)
  - [Step 12: SNP Extraction and Indexing](#step-12-snp-extraction-and-indexing)
  - [Step 13: Hard Quality Filtering](#step-13-hard-quality-filtering)
  - [Step 14: High-Quality Variant Selection](#step-14-high-quality-variant-selection)
  - [Step 15: Population-Level Filtering](#step-15-population-level-filtering)
  - [Step 16: Population Structure Analysis](#step-16-population-structure-analysis)
  - [Step 17: Advanced PCA Visualization and Population Analysis](#step-17-advanced-pca-visualization-and-population-analysis)
- [Pipeline Completion Summary](#pipeline-completion-summary)

---

## Pipeline Overview

This pipeline processes 61 tick samples through a comprehensive genomics workflow covering:
1. **Preprocessing (Steps 1-8)**: Raw read processing, quality control, alignment, and preparation
2. **Variant Calling (Steps 9-11)**: GATK-based variant discovery and joint genotyping
3. **Filtering & Analysis (Steps 12-17)**: Post-processing, filtering, and population genetic analysis

**Dataset**: NovaSeq sequencing data from 61 *Ixodes scapularis* samples
**Output Goal**: Population structure analysis via PCA and genetic diversity metrics



<img width="801" height="661" alt="pipeline_diagram" src="https://github.com/user-attachments/assets/f11c1b42-fd5a-41bf-b45a-118bb8a1a209" />
---

## Preprocessing Steps (1-8)

### Step 1: Lane Concatenation
**Script**: `preprocess_01_concatenate.sh`  
**Tools**: bash commands  
**Runtime**: 30 seconds, 61 parallel jobs  

#### Purpose and Context
This script addresses a common sequencing scenario where individual samples were sequenced across multiple lanes of the same sequencing run to increase coverage depth or optimize sequencing performance. Each tick sample was sequenced on two lanes (L001 and L002) of a NovaSeq run, requiring consolidation before downstream processing.

#### Biological Rationale
When samples are sequenced across multiple lanes, you get separate FASTQ files for each lane containing reads from the same DNA library. These represent identical genetic material processed through different physical lanes of the sequencer. The reads in L001 and L002 files sample the same underlying tick genome - they're technical replicates from the sequencing process rather than biological replicates.

#### Input Data Structure
For each of the 61 tick samples, the script expects four separate files:
- `SAMPLE_L001_R1_*.fastq.gz` (forward reads from lane 1)  
- `SAMPLE_L001_R2_*.fastq.gz` (reverse reads from lane 1)
- `SAMPLE_L002_R1_*.fastq.gz` (forward reads from lane 2)
- `SAMPLE_L002_R2_*.fastq.gz` (reverse reads from lane 2)

#### Processing Details
The script performs simple concatenation (not interleaving) where:
- All R1 reads from both lanes combine into one forward read file
- All R2 reads from both lanes combine into one reverse read file
- Maintains proper paired-end read organization
- Preserves gzip compression for storage efficiency

#### Why This Step is Essential
1. **Downstream Compatibility**: Most alignment and analysis tools expect one R1 and one R2 file per sample, not multiple lane-split files
2. **Coverage Consolidation**: Combines sequencing depth from both lanes into unified files for more robust downstream analysis  
3. **Simplified Processing**: Reduces file management complexity for the remaining pipeline steps
4. **Maintains Pair Information**: Keeps forward and reverse reads properly organized for paired-end analysis

#### Output Structure
After processing, each sample has exactly two files:
- `SAMPLE_NAME_R1_merged.fastq.gz` (all forward reads)
- `SAMPLE_NAME_R2_merged.fastq.gz` (all reverse reads)

This creates a clean, standardized input format for the subsequent quality control step where read quality will be assessed and improved before genome alignment.

#### Quality Control Notes
- Script includes robust error checking for missing files
- Verifies successful output file creation
- Provides detailed logging for troubleshooting
- Uses SLURM array jobs for efficient parallel processing of all 61 samples

---

### Step 2: Quality Control and Read Trimming
**Script**: `preprocess_02_qc_fastp.sh`  
**Tools**: fastp  
**Runtime**: 1 hour, 61 parallel jobs  

#### Purpose and Context
This step performs comprehensive quality control on the concatenated raw sequencing reads. After combining lane data in Step 1, we now need to assess and improve read quality by removing low-quality bases, adapter sequences, and reads that don't meet minimum quality standards. This is critical for accurate downstream alignment and variant calling.

#### Biological Rationale
Raw sequencing reads contain various artifacts that can compromise analysis accuracy:
- **Quality degradation**: Base calling accuracy typically decreases toward read ends due to sequencing chemistry limitations
- **Adapter contamination**: Residual adapter sequences from library preparation that don't represent tick genomic DNA
- **Short reads**: Reads shorter than useful alignment length that provide little mapping information
- **Low-quality reads**: Reads with excessive sequencing errors that could introduce false variants

#### Processing Parameters
The fastp tool is configured with specific parameters optimized for this dataset:
- **Length filtering (`-l 50`)**: Removes reads shorter than 50 base pairs, ensuring sufficient length for reliable genome alignment
- **Automatic adapter detection**: fastp automatically identifies and removes common Illumina adapter sequences
- **Quality trimming**: Trims low-quality bases from read ends (default quality threshold)
- **HTML report generation**: Creates detailed quality metrics for each sample

#### Input Data Structure
Takes the concatenated files from Step 1:
- `SAMPLE_NAME_R1_merged.fastq.gz` (forward reads)
- `SAMPLE_NAME_R2_merged.fastq.gz` (reverse reads)

#### What fastp Accomplishes
1. **Adapter Removal**: Automatically detects and removes adapter sequences that would interfere with alignment
2. **Quality Trimming**: Removes low-quality bases from read ends where sequencing errors are most common
3. **Length Filtering**: Discards reads too short for reliable mapping (< 50bp)
4. **Pair Synchronization**: Maintains proper pairing between R1 and R2 reads after filtering
5. **Quality Assessment**: Generates comprehensive quality metrics and visualizations

#### Quality Control Metrics Generated
The HTML report provides essential metrics including:
- Read count before/after filtering
- Base quality score distributions
- GC content analysis  
- Adapter contamination levels
- Length distribution of reads
- Overall data quality assessment

#### Output Structure
After processing, each sample has:
- `SAMPLE_NAME_R1_trimmed.fastq.gz` (cleaned forward reads)
- `SAMPLE_NAME_R2_trimmed.fastq.gz` (cleaned reverse reads)  
- `SAMPLE_NAME.fastp.html` (detailed quality report)

#### Why This Step is Essential
1. **Improves Alignment Accuracy**: Clean reads map more accurately to the reference genome
2. **Reduces False Variants**: Removes sequencing artifacts that could be miscalled as genetic variants
3. **Optimizes Coverage**: Retains high-quality reads while removing uninformative data
4. **Standardizes Input**: Ensures consistent read quality across all tick samples
5. **Enables Quality Assessment**: Provides metrics to identify problematic samples early

#### Expected Outcomes
Well-performing samples typically show:
- 85-95% of reads retained after filtering
- Improved average quality scores
- Minimal adapter contamination (< 1%)
- Consistent read length distributions

This cleaned, high-quality read data serves as optimal input for the subsequent FastQC quality assessment and BWA genome alignment steps.

---

### Step 3: Post-Cleaning Quality Assessment
**Script**: `preprocess_03_qc_fastqc.sh`  
**Tools**: FastQC v0.12  
**Runtime**: 1 minute, 61 parallel jobs, 4 threads each  

#### Purpose and Context
This step provides comprehensive quality assessment of the fastp-cleaned reads to verify that the trimming and filtering process was effective. While fastp performs quality control, FastQC provides detailed diagnostic reports that allow visual inspection of read quality metrics and identification of any remaining issues before proceeding to genome alignment.

#### Biological Rationale
After fastp processing, it's essential to verify that:
- Quality trimming was effective and reads now meet alignment standards
- No systematic biases were introduced during processing
- Read quality is consistent across all 61 tick samples
- Any remaining technical artifacts are identified before expensive alignment steps

This quality assessment serves as a checkpoint to catch problems early rather than discovering issues after time-intensive downstream processing.

#### What FastQC Analyzes
FastQC generates comprehensive quality reports examining multiple aspects of the cleaned reads:

**Per-base sequence quality**: Confirms that fastp successfully removed low-quality bases and that remaining bases meet high-quality thresholds across read positions.

**Per-sequence quality scores**: Verifies overall read quality distribution and ensures the majority of reads have high average quality scores.

**Per-base sequence content**: Detects compositional biases that might indicate adapter contamination or other technical artifacts not caught by fastp.

**GC content distribution**: Assesses whether GC content matches expected patterns for tick genomic DNA and identifies potential contamination.

**Sequence length distribution**: Confirms that length filtering retained reads of appropriate sizes for reliable genome mapping.

**Sequence duplication levels**: Identifies excessive PCR duplication that might affect downstream variant calling (though formal deduplication occurs in Step 7).

**Adapter content**: Verifies that fastp successfully removed adapter sequences and no residual contamination remains.

#### Processing Details
- **Multi-threading**: Uses 4 CPU threads per sample for faster processing
- **Separate analysis**: Processes R1 and R2 reads independently to identify any read-direction-specific issues
- **Standardized output**: All reports go to a centralized directory for easy comparison across samples

#### Input Data Structure
Analyzes the fastp-cleaned files from Step 2:
- `SAMPLE_NAME_R1_trimmed.fastq.gz` (cleaned forward reads)
- `SAMPLE_NAME_R2_trimmed.fastq.gz` (cleaned reverse reads)

#### Output Structure
For each sample, generates:
- `SAMPLE_NAME_R1_trimmed_fastqc.html` (detailed R1 quality report)
- `SAMPLE_NAME_R2_trimmed_fastqc.html` (detailed R2 quality report)
- Associated `.zip` files containing raw data for programmatic analysis

#### Why This Step is Essential
1. **Quality Verification**: Confirms fastp cleaning was effective and reads are suitable for alignment
2. **Problem Detection**: Identifies samples with persistent quality issues that might need special handling
3. **Consistency Assessment**: Ensures all tick samples have comparable data quality
4. **Downstream Optimization**: Informs alignment parameter choices based on actual read characteristics
5. **Documentation**: Provides quality metrics for methods sections and supplementary materials

#### Expected Quality Metrics for Good Samples
- **Per-base quality**: Phred scores > 30 across most of read length
- **Sequence quality**: >90% of reads with mean quality > 30
- **GC content**: Consistent with expected tick genome composition (~30-35%)
- **Adapter content**: <0.1% residual adapter sequences
- **Length distribution**: Tight distribution around expected read lengths

#### Quality Control Decision Points
This assessment helps determine:
- Whether any samples need re-processing with different fastp parameters
- If additional filtering steps are needed before alignment
- Whether any samples should be excluded from downstream analysis due to persistent quality issues
- Optimal alignment parameters based on actual read characteristics

The FastQC reports provide the final quality checkpoint before proceeding to the computationally expensive BWA alignment step, ensuring that only high-quality, clean reads are used for genome mapping.

---

### Step 4: Genome Alignment
**Script**: `preprocess_04_align_bwa.sh`  
**Tools**: BWA-MEM  
**Runtime**: 13 hours, 61 parallel jobs, 4 threads each, 45GB memory  

#### Purpose and Context
This step performs the core genomic alignment, mapping the cleaned sequencing reads to the *Ixodes scapularis* reference genome. This is the most computationally intensive preprocessing step, transforming millions of short DNA sequences into genomic coordinates that reveal where each read originated in the tick genome. The alignment results form the foundation for all subsequent variant calling and population genetic analyses.

#### Biological Rationale
Genome alignment addresses the fundamental question: "Where in the tick genome did each sequencing read originate?" By mapping reads to their genomic locations, we can:
- Identify positions where individuals differ from the reference genome (variants)
- Quantify coverage depth across genomic regions
- Detect structural variations and copy number changes
- Enable accurate variant calling in subsequent steps

The alignment process essentially reconstructs the relationship between the fragmented sequencing reads and the complete tick genome.

#### Reference Genome Details
Uses a **masked reference genome** (`masked_ixodes_ref_genome.fasta`):
- **Masking significance**: Repetitive regions and low-complexity sequences are masked (replaced with N's) to prevent spurious alignments
- **Improved specificity**: Reduces false positive alignments in repetitive genomic regions
- **Population genomics optimization**: Focuses analysis on unique, informative genomic regions suitable for variant detection

#### BWA-MEM Algorithm Choice
BWA-MEM is specifically chosen for this analysis because:
- **Paired-end optimization**: Handles paired-end reads effectively, using mate-pair information to improve alignment accuracy
- **Long read support**: Performs well with reads >70bp, ideal for NovaSeq data
- **Speed and accuracy**: Balances computational efficiency with alignment precision
- **Population genomics standard**: Widely used and validated for variant calling workflows

#### Processing Parameters
The BWA-MEM command uses optimized parameters:
- **Threading (`-t 4`)**: Uses 4 CPU cores per sample for faster processing
- **Mark duplicates flag (`-M`)**: Marks shorter split reads as secondary alignments (required for GATK compatibility)
- **Paired-end mode**: Automatically detects and processes R1/R2 read pairs together
- **Memory allocation**: 45GB per job to handle large genome alignment efficiently

#### Input Data Structure
Takes the high-quality cleaned reads from Steps 2-3:
- `SAMPLE_NAME_R1_trimmed.fastq.gz` (forward reads)
- `SAMPLE_NAME_R2_trimmed.fastq.gz` (reverse reads)

#### Alignment Process Details
1. **Read pair processing**: BWA-MEM analyzes both R1 and R2 reads simultaneously
2. **Seed finding**: Identifies exact matches between reads and reference genome
3. **Extension**: Extends seeds using dynamic programming to find optimal alignments
4. **Scoring**: Calculates alignment scores considering mismatches, gaps, and mate-pair information
5. **Mapping quality**: Assigns confidence scores to each alignment based on uniqueness

#### Output Structure
Generates SAM (Sequence Alignment/Map) files:
- `SAMPLE_NAME.sam` (text-based alignment file for each sample)
- Contains alignment coordinates, mapping quality scores, and read pair information
- Human-readable format enabling quality assessment and troubleshooting

#### Expected Alignment Metrics
Well-performing tick samples typically show:
- **Mapping rate**: 85-95% of reads successfully aligned
- **Proper pairs**: >90% of read pairs align in expected orientation and distance
- **Mapping quality**: Most reads have MAPQ scores >20 (high confidence alignments)
- **Coverage distribution**: Relatively uniform coverage across accessible genomic regions

#### Why This Step is Essential
1. **Genomic Context**: Transforms sequence data into genomic coordinates
2. **Variant Detection Foundation**: Enables identification of differences from reference genome
3. **Coverage Assessment**: Reveals sequencing depth across genomic regions
4. **Quality Control**: Alignment statistics indicate sample and sequencing quality
5. **Population Analysis Prep**: Creates standardized coordinate system for cross-sample comparisons

#### Computational Considerations
- **Runtime**: 6-day allocation reflects the computational intensity of aligning millions of reads
- **Memory requirements**: 45GB accommodates large genome and alignment index in memory
- **Parallel processing**: 61 simultaneous jobs maximize cluster efficiency
- **I/O optimization**: Direct output to SAM format minimizes intermediate file handling

The SAM files generated in this step contain the complete alignment information needed for downstream BAM conversion, read group addition, and duplicate removal before variant calling.

---

### Step 5: SAM to BAM Conversion and Sorting
**Script**: `preprocess_05_sam_to_bam.sh`  
**Tools**: samtools v1.19  
**Runtime**: 1 hour, 61 parallel jobs, 2 threads each, 32GB memory  

#### Purpose and Context
This step converts the text-based SAM alignment files to the binary BAM format and performs coordinate-based sorting. While SAM files are human-readable, they are inefficient for computational processing. BAM files provide the same information in a compressed, indexed format that enables rapid access and analysis. Additionally, this step generates comprehensive alignment statistics that are crucial for quality assessment.

#### Biological Rationale
The conversion from SAM to BAM format is essential for downstream processing because:
- **Storage efficiency**: BAM files are typically 3-5x smaller than equivalent SAM files
- **Processing speed**: Binary format enables much faster random access to genomic regions
- **Tool compatibility**: GATK and other variant calling tools require sorted, indexed BAM files
- **Memory optimization**: Compressed format reduces memory requirements for large-scale analysis

#### Technical Processing Details
The script performs several critical operations:

**SAM to BAM conversion**: Transforms text-based alignment data to compressed binary format while preserving all alignment information.

**Coordinate sorting**: Reorganizes reads by genomic position rather than the original read order. This sorting is essential because:
- Enables efficient random access to any genomic region
- Required for duplicate detection in subsequent steps
- Optimizes performance for variant calling algorithms
- Allows for indexed access patterns

**Memory management**: Uses 8GB sorting buffer to handle large datasets efficiently while maintaining system stability.

#### Comprehensive Statistics Generation
This step produces extensive quality metrics for each sample:

**Basic alignment statistics (`flagstat`)**:
- Total reads processed
- Properly paired reads percentage
- Mapping rate and mapping quality distribution
- Duplicate rates (preliminary assessment)

**Coverage analysis at multiple thresholds**:
- Coverage statistics at 5x and 10x depth thresholds
- Identifies regions with sufficient depth for reliable variant calling
- Assesses uniformity of coverage across the genome

**Depth of coverage analysis**:
- Per-base depth calculation across the entire reference genome
- Average depth computation for overall coverage assessment
- Depth distribution analysis to identify coverage patterns

**Contig-level coverage statistics**:
- Coverage metrics for each chromosome/scaffold
- Identifies potential coverage biases across genomic regions
- Generates coverage histograms for visual assessment

#### File Management and Optimization
The script implements several efficiency measures:
- **Temporary file cleanup**: Removes intermediate unsorted BAM files after successful sorting
- **SAM file removal**: Deletes large SAM files after successful BAM conversion to save storage
- **Mapped-only BAM creation**: Generates separate BAM files containing only successfully mapped reads
- **Indexing**: Creates BAM index files (.bai) for rapid random access

#### Input Data Structure
Takes SAM files from Step 4:
- `SAMPLE_NAME.sam` (text-based alignment files)

#### Output Structure
Generates multiple files per sample:
- `SAMPLE_NAME.sorted.bam` (primary sorted alignment file)
- `SAMPLE_NAME.sorted.mapped.bam` (mapped reads only)
- `SAMPLE_NAME.sorted.bam.bai` (BAM index file)
- Multiple statistics files with various metrics

#### Quality Control Metrics Generated
The comprehensive statistics enable assessment of:
- **Alignment success**: Overall mapping rates and read pair concordance
- **Coverage adequacy**: Depth distribution and coverage uniformity
- **Sample quality**: Comparison of metrics across all tick samples
- **Reference coverage**: How well the sequencing covers the tick genome

#### Expected Performance Metrics
Well-performing tick samples typically show:
- **Mapping rate**: 85-95% of reads successfully mapped
- **Proper pairs**: >90% of read pairs in expected orientation
- **Average depth**: 10-30x coverage across accessible genome regions
- **Coverage breadth**: >85% of reference genome covered at ≥5x depth

#### Why This Step is Essential
1. **Format optimization**: Converts to efficient binary format required for downstream tools
2. **Access efficiency**: Coordinate sorting enables rapid genomic region access
3. **Quality assessment**: Generates comprehensive metrics for sample evaluation
4. **Storage optimization**: Compressed format reduces storage requirements significantly
5. **Processing preparation**: Creates properly formatted input for read group addition and deduplication

#### Computational Considerations
- **Memory allocation**: 32GB accommodates large tick genome datasets and sorting operations
- **Temporary storage**: Uses dedicated temporary directories to avoid I/O conflicts
- **Parallel processing**: 61 simultaneous conversions maximize cluster efficiency
- **Quality control**: Multiple verification steps ensure successful file conversion

The sorted BAM files and comprehensive statistics generated in this step provide the foundation for subsequent read group addition and duplicate removal, while the quality metrics enable early identification of problematic samples before proceeding to computationally expensive variant calling steps.

---

### Step 6: Read Group Addition
**Script**: `preprocess_06_add_readgroups.sh`  
**Tools**: Picard AddOrReplaceReadGroups, samtools v1.19  
**Runtime**: 25 minutes, 61 parallel jobs, 35GB memory  

#### Purpose and Context
This step adds essential metadata (read groups) to each BAM file that identifies the source, library preparation, and sequencing details for every read. Read groups are mandatory for GATK variant calling tools and enable proper handling of technical artifacts, batch effects, and multi-sample analysis. Without proper read group information, GATK will reject the BAM files in subsequent variant calling steps.

#### Biological and Technical Rationale
Read groups serve multiple critical functions in genomics analysis:

**Sample identification**: Each read is tagged with its sample of origin, enabling GATK to properly handle multi-sample variant calling while maintaining sample-specific information.

**Technical artifact tracking**: Read group metadata allows GATK to model and correct for batch effects, lane-specific biases, and sequencing platform differences that could otherwise introduce false positive variants.

**Library preparation tracking**: Links reads to their original library preparation, which is important for identifying PCR artifacts and other library-specific technical issues.

**Sequencing platform information**: Records the sequencing technology used, allowing GATK to apply platform-specific error models and quality score interpretations.

#### Read Group Metadata Components
The script adds five essential read group fields to each read:

**RGID (Read Group Identifier)**: Set to the sample name, providing a unique identifier for this particular sequencing run of this sample.

**RGLB (Read Group Library)**: Set to "lib1", indicating the library preparation batch. While simplified here, this becomes important when samples have multiple libraries.

**RGPL (Read Group Platform)**: Set to "ILLUMINA", informing GATK about the sequencing technology used and enabling platform-specific quality score recalibration.

**RGPU (Read Group Platform Unit)**: Set to "H5VFNDMX2.${SAMPLE}", identifying the specific flow cell and sample combination. This tracks technical replicates and enables flow cell-specific bias correction.

**RGSM (Read Group Sample)**: Set to the sample name, linking reads back to their biological sample of origin for multi-sample variant calling.

#### Why This Metadata is Essential for GATK
GATK's sophisticated algorithms rely on read group information for:
- **Base quality score recalibration**: Platform and batch-specific quality score corrections
- **Variant calling confidence**: Sample-aware statistical models for variant detection
- **Multi-sample processing**: Proper attribution of reads to samples in joint calling
- **Technical artifact filtering**: Identification and correction of sequencing-specific biases

#### Processing Details
The Picard AddOrReplaceReadGroups tool:
- **Preserves alignment data**: Maintains all mapping information while adding metadata
- **Maintains sorting**: Keeps BAM files in coordinate-sorted order required for downstream tools
- **Validates format**: Ensures output BAM files meet GATK requirements
- **Handles large files**: Efficiently processes large tick genome BAM files

#### Input Data Structure
Takes sorted BAM files from Step 5:
- `SAMPLE_NAME.sorted.bam` (coordinate-sorted alignment files)

#### Output Structure
Generates read group-annotated BAM files:
- `SAMPLE_NAME.rg.sorted.bam` (BAM files with read group metadata)

#### Quality Control Considerations
Proper read group addition is verified by:
- Successful Picard execution without errors
- Output file size comparable to input (minimal overhead)
- Read group header presence in BAM file
- Maintained coordinate sorting order

#### Expected Performance
Well-executed read group addition should show:
- **Complete metadata**: All five read group fields properly populated
- **Preserved data**: No loss of reads or alignment information
- **Format compliance**: BAM files pass GATK validation checks
- **Efficient processing**: Reasonable runtime given file sizes

#### Why This Step is Essential
1. **GATK Compatibility**: Mandatory requirement for all downstream GATK tools
2. **Sample Tracking**: Enables proper multi-sample variant calling workflows
3. **Quality Control**: Allows platform and batch-specific bias correction
4. **Technical Artifact Correction**: Enables sophisticated error modeling in variant calling
5. **Metadata Preservation**: Maintains traceability from reads back to original samples

#### Computational Considerations
- **Memory allocation**: 35GB handles large BAM files and Picard's memory requirements
- **Processing time**: 2-day allocation accommodates potential I/O bottlenecks
- **File handling**: Input and output in separate directories for organization
- **Format validation**: Picard ensures GATK-compatible output format

#### Integration with Downstream Steps
The read group information added here will be utilized by:
- **Step 7**: Duplicate detection algorithms can distinguish technical vs. biological duplicates
- **Steps 9-11**: GATK variant calling tools require this metadata for proper function
- **Quality control**: Base quality score recalibration uses platform and batch information

The read group-annotated BAM files generated in this step are now properly formatted for GATK processing and ready for the crucial duplicate removal step that precedes variant calling.

---

### Step 7: Duplicate Removal
**Script**: `preprocess_07_dedup.sh`  
**Tools**: Picard MarkDuplicates, samtools v1.19  
**Runtime**: 30 minutes, 61 parallel jobs, 35GB memory  

#### Purpose and Context
This step identifies and removes PCR duplicate reads that artificially inflate coverage at specific genomic positions. PCR duplicates arise during library preparation when the same DNA template molecule is amplified multiple times, creating multiple sequencing reads that represent the same original genomic location. Removing these duplicates is essential for accurate variant calling, as they can create false signals that bias variant frequency estimates and quality scores.

#### Biological Rationale for Duplicate Removal
PCR duplicates pose several problems for genomics analysis:

**Coverage inflation**: Multiple reads from the same original DNA molecule artificially increase apparent sequencing depth at specific positions, potentially leading to overconfident variant calls.

**Bias in allele frequencies**: If one allele is preferentially amplified during PCR, duplicates can skew the apparent ratio of alternative to reference alleles, affecting variant calling accuracy.

**Quality score inflation**: Duplicate reads contribute inflated quality scores since they represent the same original sequencing event, not independent observations of the genomic position.

**Statistical model violation**: Variant calling algorithms assume independent sampling of genomic positions, which duplicates violate.

#### Technical Background on Duplicate Types
Picard MarkDuplicates identifies two types of duplicates:

**PCR duplicates**: Reads with identical 5' and 3' coordinates from both read pairs, indicating they originated from the same DNA template molecule during library preparation.

**Optical duplicates**: Reads that appear duplicated due to optical artifacts during sequencing (e.g., signal bleeding between adjacent clusters on the flow cell). The `OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500` parameter defines the pixel distance threshold for identifying these artifacts.

#### Processing Parameters and Their Significance

**REMOVE_DUPLICATES=true**: Actually removes duplicate reads from the output rather than just marking them. This creates cleaner data for variant calling while the metrics file preserves information about duplication rates.

**OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500**: Sets the pixel distance threshold for optical duplicate detection. This value is optimized for NovaSeq flow cell geometry and cluster density.

**CREATE_INDEX=true**: Automatically generates BAM index files (.bai) required for efficient random access by downstream GATK tools.

#### Quality Metrics Generated
The deduplication process produces detailed metrics including:
- **Total reads processed**: Input read count for the sample
- **Duplicate reads identified**: Number and percentage of reads marked as duplicates
- **Optical duplicates**: Subset of duplicates due to sequencing artifacts
- **Library complexity**: Estimation of original library diversity before PCR amplification
- **Duplication rate**: Overall percentage of reads removed

#### Expected Duplication Rates
Well-prepared tick genomic libraries typically show:
- **Total duplication rate**: 5-20% depending on library complexity and PCR cycles
- **Optical duplication rate**: <1% for properly prepared libraries
- **Library complexity**: High complexity indicates good library preparation

Excessively high duplication rates (>30%) may indicate:
- Over-amplification during library preparation
- Low input DNA amounts requiring excessive PCR cycles
- Technical problems during library construction

#### Input Data Structure
Takes read group-annotated BAM files from Step 6:
- `SAMPLE_NAME.rg.sorted.bam` (BAM files with read group metadata)

#### Output Structure
Generates deduplicated files:
- `SAMPLE_NAME.dedup.rg.sorted.bam` (duplicate-free BAM file)
- `SAMPLE_NAME.dedup.rg.sorted.bai` (BAM index file)
- `SAMPLE_NAME.dedup_metrics.txt` (detailed duplication statistics)

#### Why This Step is Essential
1. **Accurate variant calling**: Removes artificial coverage inflation that could bias variant detection
2. **Proper allele frequencies**: Ensures variant allele frequencies reflect true biological ratios
3. **Statistical validity**: Maintains independence assumptions required by variant calling algorithms
4. **Quality assessment**: Duplication metrics reveal library preparation quality
5. **Storage efficiency**: Reduces file sizes by removing redundant data

#### Quality Control Assessment
The duplication metrics enable evaluation of:
- **Library quality**: Low duplication rates indicate good library preparation
- **Sample comparability**: Consistent duplication rates across samples suggest uniform processing
- **Technical issues**: Extreme duplication rates may indicate problematic samples
- **Coverage estimation**: Effective coverage after duplicate removal

#### Impact on Downstream Analysis
Duplicate removal affects:
- **Coverage depth**: Reduces apparent coverage to reflect true independent observations
- **Variant calling sensitivity**: Improves accuracy by removing biased coverage
- **Population genetics**: Ensures variant frequencies accurately represent population patterns
- **Statistical power**: Maintains proper statistical assumptions for association tests

#### Computational Considerations
- **Memory requirements**: 35GB accommodates large tick genome datasets and Picard's sorting operations
- **I/O efficiency**: Creates indexed output files for rapid downstream access
- **Quality control**: Comprehensive metrics enable sample quality assessment
- **File management**: Organized output structure facilitates downstream processing

#### Expected Performance for Tick Samples
Well-performing samples should demonstrate:
- **Reasonable duplication rates**: 5-20% total duplicates
- **Low optical duplicates**: <1% of total reads
- **Successful processing**: Complete metrics files and indexed BAM output
- **Maintained data quality**: Preserved alignment quality after deduplication

The deduplicated BAM files generated in this step represent the final, analysis-ready preprocessing output. These files contain high-quality, unique alignments with proper metadata formatting required for GATK variant calling. The comprehensive metrics generated also provide important quality control information for assessing library preparation success and sample comparability across the 61 tick samples.

---

### Step 8: Comprehensive Alignment Statistics Summary
**Script**: `preprocess_08_flagstat_summary.sh`  
**Tools**: bash commands, samtools v1.19  
**Runtime**: 1 second, single job, 50GB memory  

#### Purpose and Context
This step consolidates quality metrics from all previous preprocessing steps into a comprehensive summary table that enables systematic evaluation of the entire dataset. Rather than examining individual sample files scattered across multiple directories, this creates a single unified view of preprocessing success across all 61 tick samples plus the negative control. This summary is essential for identifying problematic samples and assessing overall dataset quality before proceeding to computationally expensive variant calling.

#### Biological and Technical Rationale
After processing 61 samples through 7 preprocessing steps, it's crucial to:
- **Identify outlier samples**: Detect samples with unusual quality metrics that might need special handling
- **Assess processing uniformity**: Ensure consistent processing success across all samples
- **Quality gate checkpoint**: Make informed decisions about which samples to include in variant calling
- **Document pipeline performance**: Create records for methods sections and quality reporting

#### Comprehensive Metrics Integrated
The script consolidates metrics from multiple preprocessing stages:

**Fastp Quality Control Metrics**:
- Total reads before and after quality filtering
- Percentage of reads retained after cleaning
- Duplication rate detected during quality control

**Alignment Performance Metrics**:
- Primary reads processed through alignment
- Successfully mapped reads and mapping percentage
- Properly paired reads indicating good insert size distribution

**Coverage Analysis Metrics**:
- Mean coverage depth across the tick genome
- Percentage of genome covered at 5x depth (sufficient for variant detection)
- Percentage of genome covered at 10x depth (high-confidence variant calling)

#### Data Integration Challenges Addressed
The script handles several technical complexities:

**Multiple file formats**: Integrates data from fastp error logs, samtools flagstat outputs, and coverage analysis files across different directories.

**Paired-end read counting**: Correctly doubles read counts from fastp (which reports per-direction) to match total read counts from alignment statistics.

**Format variations**: Handles different output formats that may have been generated during reprocessing or software updates.

**Missing data handling**: Gracefully manages samples with missing files or failed processing steps.

#### Output Structure and Format
Generates two complementary summary files:

**CSV format** (`comprehensive_stats_summary.csv`): Machine-readable format suitable for statistical analysis, plotting, and import into spreadsheet applications.

**Tab-separated format** (`comprehensive_stats_summary.tsv`): Human-readable format optimized for command-line viewing and Unix text processing tools.

Both files contain identical data with 13 key metrics per sample:
1. Sample identifier
2. Total reads before fastp filtering
3. Total reads after fastp filtering  
4. Percentage of reads retained
5. Duplication rate from fastp
6. Primary reads from alignment
7. Primary mapped reads
8. Percentage of primary reads mapped
9. Properly paired reads
10. Percentage of reads properly paired
11. Mean coverage depth
12. Genome coverage at 5x depth
13. Genome coverage at 10x depth

#### Quality Assessment Benchmarks
The consolidated metrics enable assessment against expected performance standards:

**Read retention**: Well-processed samples typically retain 85-95% of reads after quality filtering.

**Mapping success**: Good tick samples show 85-95% mapping rates to the reference genome.

**Proper pairing**: >90% of reads should align in proper paired-end orientation.

**Coverage metrics**: Adequate samples show >10x mean coverage with >85% of genome covered at 5x depth.

**Sample consistency**: Metrics should be relatively consistent across samples, with outliers flagged for investigation.

#### Why This Step is Essential
1. **Quality control gateway**: Final checkpoint before expensive variant calling steps
2. **Sample selection**: Enables informed decisions about sample inclusion/exclusion
3. **Troubleshooting**: Identifies specific processing failures requiring attention
4. **Documentation**: Provides comprehensive quality metrics for publication
5. **Comparative analysis**: Enables identification of technical vs. biological sample differences

#### Expected Outcomes and Decision Points
The summary enables several key decisions:

**Sample inclusion**: Samples with poor mapping rates or low coverage may need exclusion from variant calling.

**Processing success**: Identifies samples requiring reprocessing due to technical failures.

**Dataset quality**: Overall assessment of whether the dataset meets standards for population genomics analysis.

**Parameter optimization**: Informs whether variant calling parameters need adjustment based on actual coverage distributions.

#### Computational Efficiency
- **Single job processing**: Consolidates metrics without requiring array jobs
- **Memory allocation**: 50GB accommodates processing of all sample files simultaneously
- **I/O optimization**: Minimizes file system access by processing all samples in one pass
- **Output formats**: Provides both machine and human-readable formats for different use cases

This comprehensive summary provides the critical quality assessment needed before transitioning from preprocessing to the computationally intensive GATK variant calling phase. The consolidated metrics enable informed decisions about sample quality and processing success across the entire tick population dataset.

---

## GATK Variant Calling Steps (9-11)

### Step 9: Individual Sample Variant Calling
**Script**: `variant_01_haplotypecaller.sh`  
**Tools**: GATK4 v4.6 HaplotypeCaller  
**Runtime**: 10 hours, 61 parallel jobs, 20 CPUs each, 110GB memory  

#### Purpose and Context
This step performs the core variant discovery process, identifying SNPs and small indels in each tick sample individually using GATK's HaplotypeCaller algorithm. Rather than calling variants across all samples simultaneously, this approach first generates comprehensive variant data for each sample in GVCF (Genomic Variant Call Format) files. The GVCF format captures both variant and non-variant sites with confidence information, enabling sophisticated joint genotyping in subsequent steps.

#### Biological Rationale for Individual Variant Calling
Individual sample variant calling offers several advantages for population genomics:

**Comprehensive site coverage**: GVCF format records confidence information for all genomic positions, not just variant sites, enabling accurate joint genotyping across the population.

**Sample-specific optimization**: HaplotypeCaller can optimize variant detection parameters for each sample's unique coverage and quality characteristics.

**Scalability**: Processing samples individually allows for easy addition of new samples to the dataset without reprocessing existing data.

**Quality control**: Individual GVCFs enable sample-specific quality assessment before joint analysis.

#### GATK HaplotypeCaller Algorithm
HaplotypeCaller uses a sophisticated approach fundamentally different from simpler variant callers:

**Local reassembly**: Rather than examining individual read positions, HaplotypeCaller performs local de novo assembly in regions showing variation, reconstructing possible haplotypes that could explain the observed read data.

**Haplotype-aware calling**: Simultaneously considers all variants within a local region, properly handling linked variants and complex polymorphisms that simpler callers might miss.

**Statistical modeling**: Uses advanced statistical models to calculate variant confidence scores based on read support, mapping quality, and base quality information.

**Diploid genome modeling**: Explicitly models diploid genomes, properly handling heterozygous sites and complex genotype combinations.

#### Scatter-Gather Processing Strategy
The script implements an efficient parallel processing approach:

**Genome partitioning**: Splits the tick reference genome into 10 independent chunks using GATK SplitIntervals, enabling parallel processing within each sample.

**Resource optimization**: Each chunk uses 2 CPU threads and 10GB memory, with 10 chunks running simultaneously per sample (20 CPUs, 110GB total).

**Interval-based processing**: Uses pre-defined interval lists to focus variant calling on accessible genomic regions while avoiding repetitive or problematic sequences.

**Parallel execution**: All 10 chunks for a sample run simultaneously, dramatically reducing processing time from days to hours per sample.

#### Processing Parameters and Optimization

**Ploidy specification (`-ploidy 2`)**: Explicitly sets diploid genome model appropriate for tick genetics.

**GVCF output (`--ERC GVCF`)**: Generates comprehensive variant files including non-variant sites with confidence information.

**Native threading (`--native-pair-hmm-threads 2`)**: Optimizes computational performance using multiple CPU threads for intensive calculations.

**Memory allocation (`-Xmx10G`)**: Allocates sufficient memory per chunk for complex local assembly operations.

#### GVCF Format Advantages
The Genomic VCF format provides several benefits over standard VCF:

**Non-variant site information**: Records confidence levels for sites that appear homozygous reference, enabling accurate joint genotyping.

**Scalable joint calling**: Allows efficient combination of samples without losing individual sample information.

**Quality preservation**: Maintains detailed quality metrics for downstream filtering and analysis.

**Flexible analysis**: Enables both population-level and individual-level variant analysis from the same data.

#### Input Data Structure
Takes final preprocessed BAM files from Step 7:
- `SAMPLE_NAME.dedup.rg.sorted.bam` (deduplicated, analysis-ready alignments)

#### Output Structure
Generates scattered GVCF files for each sample:
- `scattered_intervals/SAMPLE_NAME/` (interval chunks for the sample)
- `scattered_gvcfs/SAMPLE_NAME/SAMPLE_NAME.XXXX.g.vcf` (10 GVCF chunks per sample)

#### Expected Variant Discovery Metrics
Well-performing tick samples typically show:
- **Variant density**: 1-3 variants per kilobase in accessible genomic regions
- **Ti/Tv ratio**: Transition/transversion ratio of approximately 2.0-2.5 for high-quality variants
- **Het/Hom ratio**: Heterozygous to homozygous variant ratio reflecting population genetics
- **Quality distribution**: Most variants with QUAL scores >30 indicating high confidence

#### Why This Step is Essential
1. **Comprehensive variant discovery**: Identifies all genetic differences between samples and reference
2. **Population genomics foundation**: Creates the variant dataset for subsequent population analysis
3. **Quality-aware calling**: Provides confidence metrics for downstream filtering
4. **Scalable approach**: Enables efficient joint analysis across the tick population
5. **Statistical rigor**: Uses sophisticated algorithms for accurate variant detection

#### Computational Considerations
- **Intensive processing**: 7-day allocation reflects the computational complexity of local reassembly
- **Memory requirements**: 110GB accommodates simultaneous processing of multiple genomic chunks
- **Parallel efficiency**: Scatter-gather approach maximizes computational resource utilization
- **I/O optimization**: Organized output structure facilitates downstream processing

#### Integration with Population Analysis
The individual GVCFs generated here serve as input for:
- **Step 10**: Consolidation into GenomicsDB for efficient joint processing
- **Step 11**: Joint genotyping across all samples to create the final population variant dataset
- **Quality control**: Individual sample variant metrics for sample assessment

This step transforms the analysis-ready BAM files into comprehensive variant data, marking the transition from alignment-based processing to variant-based population genomics analysis. The GVCF format ensures that subsequent joint calling can leverage information from all samples while maintaining the detailed variant information discovered in each individual tick.

---

### Step 10: GenomicsDB Consolidation
**Script**: `variant_02_genomicsdb.sh`  
**Tools**: GATK4 v4.3 GenomicsDBImport  
**Runtime**: 8 hours, 10 parallel jobs (one per genomic chunk), 2 CPUs each, 115GB memory  

#### Purpose and Context
This step consolidates the scattered GVCF files from all 61 tick samples into an efficient GenomicsDB workspace that enables scalable joint genotyping. While Step 9 generated 610 individual GVCF chunks (61 samples × 10 chunks each), this step organizes them into 10 GenomicsDB workspaces corresponding to the genomic intervals. This data structure transformation is essential for the efficient joint variant calling that follows.

#### Biological Rationale for GenomicsDB
GenomicsDB addresses critical challenges in population genomics:

**Scalable joint analysis**: Traditional VCF-based approaches become computationally prohibitive with large sample sets. GenomicsDB enables efficient joint calling across dozens or hundreds of samples without memory limitations.

**Preserves individual sample information**: Unlike simple VCF merging, GenomicsDB maintains the detailed variant and non-variant site information from each sample's GVCF, enabling sophisticated population-level statistical modeling.

**Optimized data access**: The columnar storage format allows rapid access to all samples' data at any genomic position, which is exactly what joint genotyping algorithms require.

**Memory efficiency**: Compressed storage format reduces memory requirements compared to loading multiple large GVCF files simultaneously.

#### Technical Architecture and Data Organization
GenomicsDB implements a sophisticated storage strategy:

**Columnar storage**: Organizes variant data by genomic position rather than by sample, enabling efficient access to all samples' information at each genomic site.

**Compressed representation**: Uses advanced compression algorithms optimized for genomic data patterns, significantly reducing storage requirements.

**Indexed access**: Provides rapid random access to any genomic interval without scanning entire datasets.

**Parallel processing compatibility**: Designed to support efficient parallel access patterns required by joint genotyping algorithms.

#### Processing Strategy and Workflow
The script implements a chunk-based consolidation approach:

**Interval-based organization**: Creates separate GenomicsDB workspaces for each of the 10 genomic chunks, maintaining the same partitioning strategy used in individual variant calling.

**Cross-sample integration**: For each genomic chunk, consolidates the corresponding GVCF files from all 61 tick samples into a single GenomicsDB workspace.

**Scratch space optimization**: Uses high-performance scratch storage during processing to minimize I/O bottlenecks, then copies completed workspaces to permanent storage.

**File discovery automation**: Automatically identifies and imports all relevant GVCF chunks for each genomic interval, handling the complex file organization from scattered variant calling.

#### Input Data Structure
Takes scattered GVCF files from Step 9:
- `scattered_gvcfs/SAMPLE_NAME/SAMPLE_NAME.XXXX.g.vcf` (610 total GVCF chunks across all samples)

#### Processing Parameters and Optimization

**Memory allocation (`-Xmx90G`)**: Provides substantial memory for processing large numbers of samples and complex genomic regions.

**Temporary directory management**: Uses dedicated scratch space to optimize I/O performance during intensive database operations.

**POSIX optimization (`--genomicsdb-shared-posixfs-optimizations true`)**: Enables file system optimizations for improved performance on shared storage systems.

**Parallel garbage collection**: Configures Java memory management for optimal performance with large datasets.

#### Output Structure
Generates genomic interval-specific databases:
- `genomicsdb_chunks/chunk_XXXX/` (10 GenomicsDB workspaces, one per genomic interval)

Each workspace contains:
- Compressed variant data from all 61 samples for that genomic region
- Metadata and indexing information for rapid access
- Sample mapping and reference information

#### Quality Control and Validation
The consolidation process includes several validation steps:
- **File discovery verification**: Ensures all expected GVCF chunks are found and imported
- **Sample count validation**: Confirms the correct number of samples are included in each workspace
- **Workspace integrity**: Verifies successful database creation and indexing
- **Storage efficiency**: Monitors compression ratios and storage optimization

#### Why This Step is Essential
1. **Scalable joint calling**: Enables efficient joint genotyping across the entire tick population
2. **Memory optimization**: Reduces computational requirements for population-level analysis
3. **Data integration**: Combines individual sample information while preserving detailed variant data
4. **Performance enhancement**: Optimizes data access patterns for downstream joint calling
5. **Storage efficiency**: Compressed format reduces storage requirements significantly

#### Expected Performance Metrics
Well-executed GenomicsDB consolidation should show:
- **Complete sample integration**: All 61 samples successfully imported into each workspace
- **Compression efficiency**: Significant storage reduction compared to raw GVCF files
- **Access performance**: Rapid query response times for genomic intervals
- **Data integrity**: Perfect preservation of variant and non-variant site information

#### Computational Considerations
- **Memory requirements**: 115GB accommodates large sample sets and complex genomic regions
- **I/O optimization**: Scratch space utilization minimizes file system bottlenecks
- **Parallel efficiency**: 10 simultaneous jobs process different genomic regions independently
- **Storage management**: Automated copying from scratch to permanent storage

#### Integration with Joint Calling
The GenomicsDB workspaces created here serve as optimized input for:
- **Step 11**: Joint genotyping across all samples using GATK GenotypeGVCFs
- **Population analysis**: Efficient access to population-level variant data
- **Quality control**: Sample-level and site-level filtering based on joint statistics

This step represents a critical data transformation that enables population-scale analysis. By converting scattered individual variant files into an integrated, optimized database format, it creates the foundation for sophisticated joint variant calling that leverages information from all 61 tick samples simultaneously.

---

### Step 11: Joint Genotyping and VCF Assembly
**Scripts**: `variant_03_genotype_gvcfs.sh` + `variant_03b_gather_vcfs.sh`  
**Tools**: GATK4 v4.6 GenotypeGVCFs, GatherVcfs  
**Runtime**: 1 day (Step 11a) + 1 second (Step 11b), distributed processing  

#### Purpose and Context
This final variant calling step performs joint genotyping across all 61 tick samples simultaneously, leveraging information from the entire population to make more accurate variant calls. The process occurs in two phases: first, joint genotyping is performed on each genomic chunk using the GenomicsDB workspaces, then all chunks are assembled into a single comprehensive population VCF file. This approach dramatically improves variant calling accuracy compared to individual sample calling by using population-level information to distinguish true variants from sequencing artifacts.

#### Biological Rationale for Joint Genotyping
Joint genotyping provides several critical advantages for population genomics:

**Enhanced variant detection**: Low-frequency variants that might be missed in individual samples become detectable when supported by evidence across multiple samples.

**Improved genotype accuracy**: Sites that appear homozygous reference in individual calling may be revealed as heterozygous when population context indicates the presence of alternative alleles.

**Artifact discrimination**: Technical artifacts tend to occur randomly across samples, while true variants show consistent patterns across genetically related individuals.

**Allele frequency estimation**: Joint calling enables accurate estimation of population allele frequencies, essential for population genetic analyses.

**Statistical power**: Larger sample sizes provide increased statistical power for distinguishing true variants from noise.

#### GATK GenotypeGVCFs Algorithm
The joint genotyping process uses sophisticated statistical modeling:

**Population-aware calling**: Considers variant evidence from all samples simultaneously when making genotype calls at each position.

**Bayesian inference**: Uses prior probability distributions informed by population genetics principles to improve genotype accuracy.

**Quality score recalibration**: Adjusts confidence scores based on population-level evidence and known variant properties.

**Allele frequency estimation**: Calculates population-level allele frequencies while accounting for sample size and calling uncertainty.

#### Step 11a: Parallel Joint Genotyping
The first script performs joint genotyping on each genomic chunk:

**GenomicsDB input**: Uses the optimized database format from Step 10 for efficient access to all samples' data.

**Scratch optimization**: Copies reference genomes and databases to high-performance scratch storage to minimize I/O bottlenecks during intensive computations.

**Memory management**: Allocates 90GB per job to handle complex genotyping calculations across all 61 samples simultaneously.

**Parallel processing**: Runs 10 simultaneous jobs, each processing one genomic chunk independently.

#### Step 11b: VCF Assembly and Consolidation
The second script combines all genomic chunks into the final population dataset:

**Coordinate ordering**: Ensures proper genomic coordinate ordering when combining chunks from different chromosomal regions.

**Index generation**: Creates comprehensive index files for the final VCF to enable efficient downstream analysis.

**Quality preservation**: Maintains all variant quality information and population statistics from the joint calling process.

**Format standardization**: Produces a standard VCF format compatible with downstream population genetics tools.

#### Input Data Structure
Step 11a takes GenomicsDB workspaces from Step 10:
- `genomicsdb_chunks/chunk_XXXX/` (10 integrated databases)

Step 11b takes intermediate VCF files:
- `chunk_XXXX_genotyped.vcf.gz` (10 genotyped VCF chunks)

#### Output Structure
Step 11a generates:
- `genotyped_vcfs/chunk_XXXX_genotyped.vcf.gz` (10 population-called VCF chunks)

Step 11b produces the final output:
- `cohort_ticks_june2025_final.vcf.gz` (complete population variant dataset)

#### Expected Joint Calling Improvements
Joint genotyping typically shows substantial improvements over individual calling:

**Variant discovery**: 10-30% increase in total variants detected, particularly rare and low-frequency variants.

**Genotype accuracy**: Improved accuracy for heterozygous calls, especially in low-coverage regions.

**Quality scores**: More accurate variant quality scores reflecting population evidence.

**Allele frequencies**: Population-level allele frequency estimates enabling downstream analyses.

#### Why This Step is Essential
1. **Population-scale accuracy**: Leverages population information for more accurate variant calling
2. **Comprehensive variant discovery**: Detects variants missed by individual sample analysis
3. **Statistical rigor**: Provides population-level statistics required for genetic analyses
4. **Downstream compatibility**: Creates standardized VCF format for population genetics tools
5. **Quality optimization**: Improves variant quality scores using population evidence

#### Computational Considerations
**Step 11a considerations**:
- High memory requirements (115GB) for joint statistical calculations
- Scratch space optimization for I/O intensive operations
- Parallel processing across genomic chunks for efficiency

**Step 11b considerations**:
- Lower computational requirements for file concatenation
- Memory allocation sufficient for handling large VCF files
- Single job processing for maintaining coordinate order

#### Expected Population Variant Dataset
The final VCF should contain:
- **Total variants**: Hundreds of thousands to millions of SNPs and indels across the tick genome
- **Sample coverage**: All 61 tick samples with complete genotype information
- **Quality metrics**: Population-informed quality scores and statistics
- **Allele frequencies**: Population-level allele frequency estimates
- **Format compliance**: Standard VCF format compatible with downstream tools

#### Quality Control Metrics
Key metrics for assessing joint calling success:
- **Ti/Tv ratio**: Should be ~2.0-2.5 for high-quality variant sets
- **Het/Hom ratio**: Should reflect population genetic expectations
- **Singleton rate**: Proportion of variants found in only one sample
- **Allele frequency spectrum**: Distribution of variant frequencies across the population

#### Integration with Downstream Analysis
The final population VCF serves as input for:
- **Steps 12-17**: Variant filtering, quality control, and population genetic analysis
- **PLINK conversion**: Format conversion for population genetics software
- **PCA analysis**: Principal component analysis for population structure
- **Diversity calculations**: Population genetic diversity metrics

This two-phase joint calling process transforms the individual sample variant data into a comprehensive population-level dataset that captures the genetic diversity present in the 61 tick samples. The final VCF represents the complete variant catalog for subsequent population genetic analyses, filtering, and biological interpretation.

---

## Post-VCF Filtering and Analysis Steps (12-17)

### Step 12: SNP Extraction and Indexing
**Script**: `filter_after_genotypegvcfs_01.sh`  
**Tools**: GATK4 v4.6 SelectVariants, IndexFeatureFile  
**Runtime**: 2 hours, single job, 60GB memory  

#### Purpose and Context
This step initiates the post-variant calling filtering pipeline by extracting only single nucleotide polymorphisms (SNPs) from the comprehensive population VCF while excluding insertions, deletions, and complex variants. This separation is essential because SNPs and indels require different filtering strategies and quality thresholds. Most population genetic analyses focus primarily on SNPs due to their higher calling accuracy, simpler inheritance patterns, and extensive tool support for downstream analysis.

#### Biological Rationale for SNP-Focused Analysis
SNPs offer several advantages for population genomics studies:

**Calling accuracy**: SNPs generally have higher accuracy rates than indels because they involve simpler sequence changes that are easier to detect reliably.

**Inheritance simplicity**: Single base changes follow straightforward Mendelian inheritance patterns without the alignment complexities associated with length variants.

**Population genetics compatibility**: Most population genetic theory and software tools are optimized for bi-allelic SNP markers.

**Evolutionary significance**: SNPs represent the most common type of genetic variation and provide excellent resolution for population structure and demographic analyses.

**Technical uniformity**: SNPs have more consistent quality score distributions and filtering characteristics compared to the heterogeneous nature of indel variants.

#### Variant Type Separation Strategy
The separation of variant types enables optimized processing:

**Different quality distributions**: SNPs and indels have different quality score patterns that require distinct filtering thresholds.

**Algorithmic considerations**: Some population genetic analyses work exclusively with SNPs and cannot handle mixed variant types.

**Quality control focus**: Separating variant types allows for type-specific quality assessment and filtering optimization.

**Downstream tool compatibility**: Many population genetics tools expect SNP-only datasets in specific formats.

#### Processing Steps and Technical Details

**VCF indexing**: Creates an index file (.tbi) that enables rapid random access to the large population VCF file, essential for efficient variant extraction.

**SNP selection**: Uses GATK SelectVariants with `--select-type-to-include SNP` to extract only single nucleotide changes while excluding:
- Insertions and deletions (indels)
- Mixed variant types (e.g., complex substitutions)
- Symbolic variants (e.g., structural variants)
- Multi-allelic complex sites

**Memory optimization**: Allocates sufficient memory (35GB for selection) to handle the large population VCF efficiently without memory bottlenecks.

#### Input Data Structure
Takes the final population VCF from Step 11:
- `cohort_ticks_june2025_combined.vcf.gz` (complete population variant dataset)

#### Output Structure
Generates SNP-focused dataset:
- `cohort_ticks_june2025_snps_only.vcf.gz` (SNPs extracted from population data)
- Associated index files for efficient access

#### Expected SNP Dataset Characteristics
The extracted SNP dataset should show:
- **Variant count**: Typically 70-90% of total variants are SNPs in well-called datasets
- **Quality distribution**: More uniform quality scores compared to mixed variant types
- **Ti/Tv ratio**: Transition/transversion ratio should be ~2.0-2.5 for high-quality SNP sets
- **Allele frequency spectrum**: Full range of population allele frequencies preserved

#### Why This Step is Essential
1. **Optimized filtering**: Enables SNP-specific filtering thresholds in subsequent steps
2. **Tool compatibility**: Creates datasets compatible with population genetics software
3. **Quality focus**: Concentrates analysis on the most reliable variant type
4. **Processing efficiency**: Reduces computational burden for downstream analyses
5. **Methodological consistency**: Aligns with standard practices in population genomics

#### Quality Control Considerations
This step enables assessment of:
- **Variant type composition**: Proportion of SNPs vs. other variant types in the original dataset
- **Quality preservation**: Maintenance of variant quality information during extraction
- **Sample completeness**: Verification that all 61 samples retain data in the SNP subset
- **Coordinate integrity**: Proper genomic ordering and indexing of extracted variants

#### Expected Performance Metrics
Well-executed SNP extraction should demonstrate:
- **High SNP retention**: >90% of high-quality SNPs preserved from original dataset
- **Proper indexing**: Successful creation of index files for rapid access
- **Quality consistency**: Maintained variant quality scores and population statistics
- **Format compliance**: Standard VCF format suitable for downstream filtering

#### Integration with Downstream Filtering
The SNP-only dataset serves as input for:
- **Step 13**: Hard filtering based on GATK best practices for SNP quality thresholds
- **Step 14**: Selection of only variants passing all quality filters
- **Step 15**: Population-level filtering for missingness and minor allele frequency
- **Steps 16-17**: Format conversion and population genetic analysis

#### Computational Considerations
- **Memory allocation**: 60GB accommodates large population VCF files and indexing operations
- **Processing time**: 10-hour allocation provides buffer for potential I/O limitations
- **Storage efficiency**: Compressed output maintains manageable file sizes
- **Index generation**: Automatic indexing enables efficient downstream processing

This step marks the transition from comprehensive variant calling to focused population genetic analysis. By extracting and indexing the SNP subset, it creates an optimized dataset for the rigorous quality filtering and population structure analyses that follow in the remaining pipeline steps.

---

### Step 13: Hard Quality Filtering
**Script**: `filter_after_genotypegvcfs_02.sh`  
**Tools**: GATK4 v4.6 VariantFiltration  
**Runtime**: 1 second, single job, 40GB memory  

#### Purpose and Context
This step applies rigorous quality filters to the SNP dataset using GATK best practices to identify and flag variants that fail to meet statistical quality thresholds. Rather than removing variants immediately, this approach tags failing variants with filter flags, preserving the complete dataset while clearly marking questionable sites. This filtering strategy is based on extensive validation studies and represents the gold standard for variant quality control in population genomics.

#### Biological Rationale for Hard Filtering
Quality filtering addresses fundamental challenges in variant calling:

**Sequencing artifacts**: Technical errors during sequencing, library preparation, or base calling can create false positive variant calls that mimic real genetic variation.

**Alignment artifacts**: Mapping errors, particularly in repetitive regions, can generate spurious variant calls where reads are incorrectly aligned.

**Coverage biases**: Regions with very low or extremely high coverage produce unreliable variant calls due to insufficient evidence or technical artifacts.

**Statistical confidence**: Variants with low statistical support may represent noise rather than true genetic differences.

#### GATK Best Practices Filter Annotations
Each filter targets specific types of technical artifacts:

**QD (Quality by Depth) < 2.0**: Identifies variants where the variant confidence is low relative to the depth of coverage. Low QD suggests the variant call is not well-supported by the available reads.

**QUAL < 30.0**: Filters variants with overall low confidence scores. QUAL represents the Phred-scaled probability that the variant call is incorrect, so QUAL < 30 indicates >10% chance of error.

**SOR (Strand Odds Ratio) > 3.0**: Detects variants showing significant strand bias, where variant-supporting reads come predominantly from one DNA strand, indicating potential technical artifacts.

**FS (Fisher Strand) > 60.0**: Alternative measure of strand bias using Fisher's exact test. High FS values indicate significant imbalance in strand representation.

**MQ (Mapping Quality) < 40.0**: Flags variants in regions where reads have poor mapping quality, indicating potential alignment errors or repetitive sequences.

**MQRankSum < -12.5**: Identifies sites where variant-supporting reads have significantly lower mapping quality than reference-supporting reads, suggesting mapping artifacts.

**ReadPosRankSum < -8.0**: Detects variants where variant-supporting bases occur predominantly at read ends, where sequencing errors are more common.

#### Filter Strategy and Implementation
The VariantFiltration tool implements a tagging approach:

**Preservation of data**: Variants failing filters are marked with descriptive tags (e.g., "QD2", "QUAL30") rather than being removed from the dataset.

**Multiple filter support**: Variants can fail multiple filters simultaneously, with all applicable filter tags recorded.

**Downstream flexibility**: Researchers can choose different stringency levels by selecting which filter flags to exclude in subsequent analyses.

**Quality assessment**: Filter statistics provide insights into dataset quality and potential technical issues.

#### Input Data Structure
Takes SNP-only dataset from Step 12:
- `cohort_ticks_june2025_snps_only.vcf.gz` (extracted SNP variants)

#### Output Structure
Generates quality-annotated dataset:
- `cohort_ticks_june2025_snps_filtered.vcf.gz` (SNPs with quality filter annotations)

#### Expected Filtering Outcomes
Well-performing datasets typically show:
- **Filter rates**: 5-15% of variants failing one or more quality filters
- **Quality distribution**: Most high-quality variants passing all filters
- **Filter correlation**: Some correlation between different filter failures
- **Biological relevance**: Filtered variants concentrated in problematic genomic regions

#### Why This Step is Essential
1. **Technical artifact removal**: Identifies and flags likely false positive variant calls
2. **Statistical rigor**: Ensures downstream analyses use only high-confidence variants
3. **Reproducibility**: Applies standardized, validated filtering criteria
4. **Quality assessment**: Provides metrics for evaluating dataset quality
5. **Flexibility preservation**: Maintains complete dataset while marking questionable variants

#### Filter Interpretation and Biological Context
Each filter addresses specific biological scenarios:

**Coverage-related filters (QD, QUAL)**: Target regions with insufficient sequencing evidence for reliable variant calling, common in repetitive or difficult-to-sequence genomic regions.

**Bias-related filters (SOR, FS, MQRankSum, ReadPosRankSum)**: Identify systematic biases that could indicate technical artifacts rather than true genetic variation.

**Mapping-related filters (MQ)**: Focus on regions where alignment uncertainty could lead to false variant calls, particularly important in complex genomic regions.

#### Quality Control and Assessment
This filtering step enables evaluation of:
- **Dataset quality**: Overall proportion of variants passing quality thresholds
- **Technical issues**: Identification of systematic problems requiring investigation
- **Regional biases**: Areas of the genome with elevated filter failure rates
- **Sample consistency**: Uniform filtering performance across all tick samples

#### Expected Performance Metrics
Successful hard filtering should demonstrate:
- **Appropriate filter rates**: Neither too lenient (>20% failures) nor too stringent (<5% failures)
- **Quality improvement**: Better overall quality score distributions in passing variants
- **Biological consistency**: Filter failures concentrated in expected problematic regions
- **Sample uniformity**: Consistent filtering outcomes across all 61 tick samples

#### Integration with Downstream Analysis
The quality-annotated dataset serves as input for:
- **Step 14**: Selection of only variants passing all quality filters
- **Step 15**: Population-level filtering for missingness and allele frequency
- **Quality reporting**: Documentation of filtering outcomes for publication

#### Computational Considerations
- **Memory allocation**: 40GB accommodates large SNP datasets and annotation operations
- **Processing efficiency**: Single-threaded operation optimized for annotation tasks
- **I/O optimization**: Compressed input/output maintains storage efficiency
- **Quality preservation**: All original variant information maintained alongside filter annotations

This step represents the application of rigorously validated quality standards to ensure that subsequent population genetic analyses are based on high-confidence variant calls. By applying GATK best practices, it creates a quality-controlled dataset suitable for reliable biological interpretation and population structure analysis.

---

### Step 14: High-Quality Variant Selection
**Script**: `filter_after_genotypegvcfs_03_select_passing_only.sh`  
**Tools**: GATK4 v4.6 SelectVariants  
**Runtime**: 1 hour, single job, 20GB memory  

#### Purpose and Context
This step extracts only the variants that passed all quality filters from Step 13, creating a clean, high-confidence dataset for population genetic analysis. While the previous step flagged problematic variants with filter tags, this step definitively removes them, producing a streamlined dataset containing only reliable SNPs suitable for downstream population structure and diversity analyses.

#### Biological Rationale for Quality Selection
Creating a high-confidence variant subset is essential for reliable population genomics:

**Statistical validity**: Population genetic analyses assume that observed variants represent true biological differences rather than technical artifacts. Including low-quality variants can bias allele frequency estimates, population structure analyses, and demographic inferences.

**Reproducibility**: By using only variants that meet stringent quality criteria, analyses become more reproducible across different datasets and studies.

**Computational efficiency**: Smaller, high-quality datasets reduce computational requirements for intensive population genetic analyses while maintaining statistical power.

**Biological interpretation**: Reliable variants enable confident biological conclusions about population structure, gene flow, and evolutionary processes.

#### Technical Implementation and JEXL Expression
The script uses GATK's Java Expression Language (JEXL) to select variants:

**Selection criterion (`vc.isNotFiltered()`)**: This expression identifies variants where the FILTER field equals "PASS", meaning the variant passed all quality thresholds applied in Step 13.

**Exclusion logic**: Variants with any filter flags (QD2, QUAL30, SOR3, FS60, MQ40, MQRankSum-12.5, ReadPosRankSum-8) are excluded from the output.

**Preservation of quality**: All quality annotations and population genetic information are maintained for the selected high-confidence variants.

#### Expected Dataset Transformation
The selection process typically results in:

**Variant reduction**: 85-95% of original SNPs typically pass GATK hard filters in well-prepared datasets, meaning 5-15% of variants are removed.

**Quality improvement**: Dramatic improvement in overall dataset quality metrics, with remaining variants showing higher confidence scores and better statistical properties.

**Sample retention**: All 61 tick samples are maintained in the dataset, but some samples may have fewer called variants in regions where low-quality variants were removed.

**Genomic coverage**: High-quality variants should be well-distributed across the accessible tick genome, though some problematic regions may show reduced variant density.

#### Input Data Structure
Takes quality-annotated SNPs from Step 13:
- `cohort_ticks_june2025_snps_filtered.vcf.gz` (SNPs with filter annotations)

#### Output Structure
Generates high-confidence dataset:
- `cohort_ticks_june2025_snps_filtered_only.vcf.gz` (only variants passing all filters)

**Note**: The script filename suggests "filtered_only" but the JEXL expression `vc.isNotFiltered()` actually selects PASSING variants, not filtered variants. This creates a clean dataset of high-quality SNPs.

#### Quality Assessment and Validation
The high-confidence dataset should demonstrate:

**Improved quality metrics**:
- Higher average QUAL scores
- Better Ti/Tv ratios (closer to expected ~2.0-2.5)
- More uniform quality score distributions
- Reduced technical artifact signatures

**Population genetic validity**:
- Reasonable heterozygosity levels consistent with tick population genetics
- Allele frequency spectra consistent with demographic expectations
- Reduced singleton rates (variants found in only one sample)

#### Why This Step is Essential
1. **Analysis reliability**: Ensures population genetic analyses are based on high-confidence variants
2. **Bias reduction**: Removes technical artifacts that could skew population structure analyses
3. **Statistical power**: Concentrates analysis on variants with strongest statistical support
4. **Downstream compatibility**: Creates clean datasets suitable for population genetics software
5. **Quality assurance**: Provides final quality checkpoint before biological analysis

#### Expected Performance Metrics
Well-executed variant selection should show:
- **Retention rate**: 85-95% of original SNPs retained in high-quality dataset
- **Quality improvement**: Substantial improvement in average variant quality scores
- **Coverage maintenance**: Adequate variant density across the tick genome
- **Sample representation**: All 61 samples with reasonable variant counts

#### Integration with Population Analysis
The high-confidence SNP dataset serves as input for:
- **Step 15**: Population-level filtering for missingness and minor allele frequency
- **Step 16**: Format conversion to PLINK for population genetic software
- **Step 17**: Principal component analysis and population structure assessment

#### Computational Considerations
- **Reduced memory requirements**: 20GB sufficient for processing clean, filtered dataset
- **Processing efficiency**: Simple selection operation with minimal computational overhead
- **Storage optimization**: Smaller output files reduce storage requirements
- **Quality validation**: Enables assessment of filtering effectiveness

#### Quality Control Checkpoints
This step enables evaluation of:
- **Filter effectiveness**: Assessment of how many variants were removed and why
- **Dataset quality**: Overall improvement in variant quality metrics
- **Population coverage**: Verification that adequate variants remain for population analysis
- **Sample balance**: Confirmation that filtering didn't disproportionately affect specific samples

#### Expected Biological Outcomes
The final high-confidence dataset should provide:
- **Robust population structure**: Reliable data for PCA and population clustering
- **Accurate allele frequencies**: Unbiased estimates for population genetic analyses
- **Demographic insights**: Clean data for inferring population history and gene flow
- **Evolutionary analysis**: High-quality variants for selection and adaptation studies

This step completes the technical quality control phase of the pipeline, transforming the comprehensive variant dataset into a curated collection of high-confidence SNPs ready for biological analysis. The resulting dataset represents the optimal balance between stringent quality control and sufficient variant density for meaningful population genetic inferences about the 61 tick samples.

---

### Step 14b: Pre-Filter Dataset Characterization
**Script**: `downsample_to_100k_before_filtering.sh`  
**Tools**: bcftools, vcflib  
**Runtime**: 1 hours, single job, 50GB memory  

#### Purpose and Context
This step creates a representative subset of the raw, unfiltered variant dataset to establish baseline characteristics before applying quality filters. By downsampling the complete population VCF to approximately 100,000 sites, we can efficiently analyze the full spectrum of variant quality distributions, including low-quality variants that will be removed in subsequent filtering steps.

#### Technical Implementation
Uses `vcfrandomsample` with a sampling rate of 0.00037 to achieve approximately 100k variants from the complete joint-called dataset. This provides a statistically representative sample of all variants, including those that would fail quality filters, enabling comprehensive assessment of the unfiltered data landscape.

#### Input Data Structure
Takes the complete population VCF from Step 11:
- `cohort_ticks_june2025_final.vcf.gz` (unfiltered joint-called variants)

#### Output Structure
- `100Ksubset_before.cohort_ticks_june2025_final.vcf.gz` (pre-filter representative sample)

---

### Step 14b: Post-Filter Dataset Characterization  
**Script**: `downsample_to_100k_after_filtering.sh`  
**Tools**: bcftools, vcflib  
**Runtime**: 1 hour, single job, 50GB memory  

#### Purpose and Context
This step creates a representative subset of the high-quality, filtered variant dataset to assess the effectiveness of quality filtering and characterize the final data properties. By downsampling the PASS-only variants, we can evaluate how filtering has improved data quality and inform population-level filtering parameters.

#### Technical Implementation
Uses the same `vcfrandomsample` approach on the quality-filtered dataset, enabling direct comparison with the pre-filter analysis to quantify the impact of GATK hard filtering on variant quality distributions.

#### Input Data Structure
Takes high-quality variants from Step 14:
- `cohort_ticks_june2025_snps_passing_only.vcf.gz` (GATK-filtered variants)

#### Output Structure
- `100Ksubset_after.cohort_ticks_june2025_snps_passing_only.vcf.gz` (post-filter representative sample)

---

### Step 14c: Comparative Dataset Analysis
**Scripts**: `after_downsample_before_filtering_for_summary_VCF_R_plots.sh`, `after_downsample_after_filtering_for_summary_VCF_plots.sh`, `VCF_R_100k_before_filtering.sh`, `VCF_R_100k_after_filtering.sh`  
**Tools**: VCFtools, R (tidyverse)  
**Runtime**: 1 hour, statistical analysis and visualization  

#### Purpose and Context
This comprehensive analysis pipeline performs parallel statistical characterization of both pre-filter and post-filter datasets, generating comparative visualizations that demonstrate the effectiveness of quality filtering and inform population-level parameter optimization.

#### Before/After Analysis Workflow
The analysis includes two parallel tracks:

**Pre-Filter Analysis Track**:
1. **VCFtools statistics** (`after_downsample_before_filtering_for_summary_VCF_R_plots.sh`) on `100Ksubset_before.cohort_ticks_june2025_final.vcf.gz`
2. **R visualization** (`VCF_R_100k_before_filtering.sh`) generating blue-colored plots with "Before" labels

**Post-Filter Analysis Track**:
1. **VCFtools statistics** (`after_downsample_after_filtering_for_summary_VCF_plots.sh`) on `100Ksubset_after.cohort_ticks_june2025_snps_filtered_only.vcf.gz`
2. **R visualization** (`VCF_R_100k_after_filtering.sh`) generating orange-colored plots with "After" labels

#### Analytical Components Generated
Both analysis tracks produce comprehensive assessments:

**Variant Quality Metrics**: Distribution of variant confidence scores, demonstrating filtering effectiveness in removing low-confidence calls.

**Depth Analysis**: Site-level and individual-level coverage patterns, showing how filtering affects coverage uniformity and identifies optimal depth parameters.

**Missingness Assessment**: Before/after comparison of data completeness, crucial for determining appropriate missingness thresholds for population filtering.

**Allele Frequency Characterization**: MAF spectrum analysis comparing filtered vs. unfiltered datasets, informing population-level MAF filtering decisions.

**Individual-Level Quality**: Per-sample metrics including heterozygosity, inbreeding coefficients, and data completeness across the filtering process.

#### Comparative Visualization Outputs
- **Quality improvement plots**: Side-by-side comparison of variant quality distributions
- **Filtering effectiveness assessment**: Quantification of variants removed and quality gains achieved
- **Parameter optimization guidance**: Data-driven recommendations for population filtering thresholds
- **Publication-ready figures**: Professional comparative visualizations demonstrating pipeline effectiveness

#### Integration with Population Filtering
The comprehensive before/after analysis directly informs Step 15 parameter selection by:
- **Quantifying filtering impact**: Demonstrating improvement in data quality metrics
- **Optimizing MAF thresholds**: Showing how different MAF cutoffs would affect the filtered dataset
- **Informing missingness parameters**: Revealing optimal completeness requirements based on actual data patterns
- **Validating filtering effectiveness**: Confirming that GATK hard filters successfully improved dataset quality

This data-driven approach ensures that subsequent population-level filtering (Step 15) uses parameters optimized for the specific characteristics of the tick dataset rather than generic thresholds, maximizing both data quality and variant retention for robust population genetic analysis.

---

### Step 15: Population-Level Filtering
**Script**: `filter_after_genotypegvcfs_04_maf_and_missingness.sh`  
**Tools**: VCFtools  
**Runtime**: 1 hour 15 minutes, single job, 64GB memory  

#### Purpose and Context
This step applies population-level filters focusing on variant utility for population genetic analysis rather than technical quality. While previous steps addressed sequencing and alignment artifacts, this step removes variants that are either too rare to provide statistical power or have too much missing data across samples to enable reliable population inferences. These filters optimize the dataset for downstream population structure analysis and ensure robust statistical analyses.

#### Biological Rationale for Population-Level Filtering
Population genetic analyses require variants that meet specific criteria for statistical validity:

**Minor allele frequency requirements**: Very rare variants (MAF < 5%) provide little statistical power for population structure analysis and may represent sequencing errors, recent mutations, or population-specific variants that don't inform broad population relationships.

**Completeness requirements**: Variants with excessive missing data across samples can bias population genetic analyses, as missing genotypes are often non-random and can be confounded with population structure or technical factors.

**Statistical power considerations**: Population genetic methods require sufficient allele frequency and complete data to distinguish population structure from noise and technical artifacts.

#### Filter Parameters and Their Significance

**Minor Allele Frequency (MAF ≥ 0.05)**:
- Retains variants where the less common allele appears in at least 5% of chromosomes
- For 61 diploid samples (122 chromosomes), this means ≥6 copies of the minor allele
- Excludes singletons and very rare variants that provide minimal information for population structure
- Focuses analysis on variants with sufficient representation across the population

**Missingness Threshold (≥70% complete data)**:
- Requires variants to have successful genotype calls in at least 70% of samples
- For 61 samples, this means ≥43 samples must have called genotypes
- Excludes variants with excessive missing data that could bias population analyses
- Balances data completeness with variant retention

#### VCFtools Implementation Strategy
VCFtools provides efficient population genetic filtering:

**Frequency calculation**: Automatically calculates minor allele frequencies across all 61 samples
**Missing data assessment**: Evaluates completeness on a per-variant basis
**Combined filtering**: Applies both criteria simultaneously for efficiency
**Format preservation**: Maintains VCF format compatibility while applying population filters
**Information retention**: Preserves all variant annotations and quality information

#### Expected Dataset Transformation
Population-level filtering typically results in:

**Substantial variant reduction**: Often 50-80% reduction from technical filtering, as many variants are rare or have missing data
**Enhanced statistical power**: Remaining variants provide better resolution for population structure analysis
**Improved data completeness**: Higher average completion rates across samples
**Focus on informative variants**: Concentration on variants that distinguish populations rather than individuals

#### Input Data Structure
Takes high-confidence SNPs from Step 14:
- `cohort_ticks_june2025_snps_passing_only.vcf.gz` (technically validated variants)

#### Output Structure
Generates population-optimized dataset:
- `cohort_ticks_june2025_snps_passing_only.maf005.miss03.vcf.gz` (population-filtered variants)

#### Biological Implications of Filtering Choices
The filtering parameters reflect population genetic priorities:

**MAF 0.05 threshold**: 
- Excludes variants potentially representing sequencing errors or very recent mutations
- Focuses on variants likely to be segregating in multiple populations
- Provides sufficient statistical power for structure analyses
- Compatible with most population genetic software requirements

**70% completeness threshold**:
- Balances data quality with variant retention
- Removes variants in highly repetitive or difficult-to-sequence regions
- Ensures sufficient sample representation for population inferences
- Prevents bias from systematic missing data patterns

#### Why This Step is Essential
1. **Statistical optimization**: Creates dataset optimized for population genetic analysis methods
2. **Power enhancement**: Focuses on variants with sufficient frequency for meaningful analysis
3. **Bias reduction**: Removes variants with non-random missing data patterns
4. **Method compatibility**: Ensures compatibility with downstream population genetics tools
5. **Interpretability**: Concentrates on variants likely to reflect population processes

#### Expected Performance Metrics
Well-executed population filtering should show:
- **Variant reduction**: 50-80% of variants typically removed by MAF and missingness filters
- **Quality improvement**: Higher average completion rates and more uniform allele frequencies
- **Population representation**: Adequate variant density for population structure analysis
- **Statistical validity**: Sufficient variants for robust PCA and clustering analyses

#### Integration with Population Analysis
The population-optimized dataset serves as input for:
- **Step 16**: Format conversion to PLINK binary format for population genetics software
- **Step 17**: Principal component analysis and population structure assessment
- **Downstream analyses**: Phylogenetic analysis, admixture modeling, demographic inference

#### Computational Considerations
- **Memory allocation**: 64GB accommodates population-scale frequency calculations
- **Processing time**: 2-day allocation reflects computational intensity of population filtering
- **I/O efficiency**: Compressed input/output maintains storage efficiency
- **Statistical calculations**: VCFtools optimized for population genetic computations

#### Quality Control and Assessment
This filtering enables evaluation of:
- **Variant distribution**: Geographic and frequency distribution of retained variants
- **Sample balance**: Verification that filtering doesn't bias against specific samples
- **Population coverage**: Adequate variant density across the tick genome
- **Statistical readiness**: Dataset preparation for population genetic analyses

#### Expected Biological Outcomes
The population-filtered dataset should provide:
- **Clear population structure**: Enhanced resolution of genetic relationships between samples
- **Demographic signals**: Variants suitable for inferring population history
- **Geographic patterns**: Spatial genetic structure if present in the tick samples
- **Evolutionary insights**: Data for understanding tick population dynamics and gene flow

This step completes the transition from technical quality control to biological optimization, creating a dataset specifically tailored for population genetic discovery. The resulting variants represent the most informative subset for understanding population structure, demographic history, and evolutionary processes in the 61 tick samples.

---

### Step 16: Population Structure Analysis
**Script**: `analysis_01_plink.sh`  
**Tools**: PLINK v1.90, PLINK2  
**Runtime**: 1 hour, single job, 12 CPUs, 90GB memory  

#### Purpose and Context
This step converts the population-optimized VCF dataset to PLINK binary format and performs Principal Component Analysis (PCA) to reveal population structure and genetic relationships among the 61 tick samples. PCA is the gold standard method for population structure analysis, transforming the high-dimensional genetic data into interpretable components that capture the major axes of genetic variation. This analysis represents the culmination of the entire pipeline, revealing the underlying population genetic patterns in the tick dataset.

#### Biological Rationale for Population Structure Analysis
Understanding population structure is fundamental to population genomics:

**Evolutionary relationships**: PCA reveals genetic relatedness between samples, indicating shared evolutionary history and demographic processes.

**Geographic patterns**: Population structure often correlates with geographic distribution, revealing migration patterns, barriers to gene flow, and local adaptation.

**Demographic history**: The extent and pattern of population structure provides insights into historical population size changes, founder effects, and admixture events.

**Study design implications**: Understanding population structure is essential for interpreting subsequent analyses and controlling for confounding factors in association studies.

#### Technical Implementation Strategy
The analysis uses a two-tool approach optimized for different tasks:

**PLINK v1.90 for format conversion**:
- Efficient VCF to binary PLINK format conversion
- Biallelic variant filtering to ensure PCA compatibility
- Sample and variant quality control integration
- Optimized for handling large population datasets

**PLINK2 for PCA analysis**:
- Advanced PCA algorithms optimized for genomic data
- Efficient computation with large sample sizes and variant counts
- Standardized output formats for downstream visualization
- Superior performance for population genetic analyses

#### Processing Steps and Parameters

**LD Pruning and Variant ID Assignment**:
- `--vcf`: Specifies the input VCF file
- `--set-missing-var-ids @:#`: Corrects the critical error from the old pipeline. This command ensures that every variant, even those with a missing ID in the VCF, gets a unique, verifiable ID based on its genomic coordinates.
- `--indep-pairwise 50 10 0.1`: The core LD pruning parameters. It identifies and removes one of any two variants that are in high linkage disequilibrium (r^2 0.1) within a given sliding window, creating a list of independent variants for the next step.

**Variant Extraction, Format Conversion, and PCA Computation**:
- `--extract [prune.in file]`: Uses the list of independent variants generated in the first step
- `--make-bed`: Creates the binary PLINK format (.bed, .bim, .fam) for the pruned dataset
- `--pca 20`: Computes the first 20 principal components, providing comprehensive view of population structure


#### Expected PCA Outcomes and Interpretation
Principal Component Analysis reveals different aspects of population structure:

**PC1 and PC2**: Typically capture the major axes of population differentiation, often corresponding to geographic or demographic separation.

**Subsequent PCs (3-20)**: May reveal:
- Sub-population structure within major groups
- Recent admixture events
- Technical artifacts requiring further investigation
- Fine-scale geographic patterns

**Variance explained**: The proportion of genetic variance explained by each PC indicates the strength of population structure, with higher values suggesting stronger differentiation.

#### Input Data Structure
Takes population-optimized variants from Step 15:
- `cohort_ticks_june2025_snps_passing_only.maf005.miss03.vcf.gz` (filtered population dataset)

#### Output Structure
Generates population genetics files:
- `population_structure/tick_population.bed/.bim/.fam` (PLINK binary format)
- `population_structure/tick_pca.eigenvec` (principal component scores)
- `population_structure/tick_pca.eigenval` (variance explained by each PC)

#### Biological Interpretation Framework
PCA results can reveal various population genetic scenarios:

**No clear structure**: Random scatter suggests a single, panmictic population with high gene flow.

**Discrete clusters**: Separated groups indicate distinct populations with limited gene flow, possibly due to geographic barriers or demographic history.

**Continuous variation**: Gradual transitions suggest isolation-by-distance or continuous migration patterns.

**Outlier samples**: Individuals plotting separately may represent admixed individuals, recent migrants, or samples with technical issues.

#### Why This Step is Essential
1. **Population discovery**: Reveals the fundamental genetic structure of the tick population
2. **Quality assessment**: Identifies potential technical issues or sample problems
3. **Study design**: Informs subsequent analyses and interpretation frameworks
4. **Biological insight**: Provides immediate insights into tick population biology and evolution
5. **Visualization preparation**: Creates data for publication-quality population structure plots

#### Expected Performance Metrics
Successful population structure analysis should demonstrate:
- **Efficient conversion**: Complete VCF to PLINK format conversion without data loss
- **Adequate variants**: Sufficient SNPs (typically thousands to tens of thousands) for robust PCA
- **Clear patterns**: Interpretable population structure reflecting biological expectations
- **Technical validity**: PCA results consistent with known tick biology and geography

#### Computational Considerations
- **Memory allocation**: 90GB accommodates large variant matrices in memory for efficient PCA computation
- **Multi-threading**: 12 CPUs enable parallel processing for computationally intensive PCA calculations
- **I/O optimization**: Binary PLINK format provides efficient data access for iterative calculations
- **Algorithm efficiency**: PLINK2 implements optimized PCA algorithms for genomic data

#### Quality Control and Validation
This analysis enables assessment of:
- **Data integrity**: Verification that filtering steps preserved population genetic signal
- **Sample quality**: Identification of problematic samples through PCA outliers
- **Population expectations**: Comparison of observed structure with biological expectations
- **Technical artifacts**: Detection of batch effects or technical confounding

#### Integration with Visualization and Interpretation
The PCA results serve as input for:
- **Step 17**: Advanced visualization and population genetic analysis
- **Publication figures**: Population structure plots for manuscripts
- **Downstream analyses**: Population-aware association studies or selection analyses
- **Biological interpretation**: Understanding tick population biology and conservation genetics

#### Expected Biological Insights
The population structure analysis should reveal:
- **Genetic diversity**: Overall levels of genetic variation within and between tick populations
- **Population relationships**: Evolutionary relationships and demographic connections between samples
- **Geographic patterns**: Spatial structure reflecting tick ecology and dispersal patterns
- **Conservation implications**: Information relevant to tick population management and disease ecology

This step represents the analytical culmination of the entire pipeline, transforming millions of raw sequencing reads into interpretable population genetic insights. The PCA results provide the foundation for understanding the evolutionary biology, ecology, and population dynamics of the 61 tick samples, enabling both basic research discoveries and applied insights relevant to tick-borne disease management.

---

### Step 17: Advanced PCA Visualization and Population Analysis
**Script**: `pca_generation.py`  
**Tools**: Python (pandas, matplotlib, seaborn, numpy)  
**Runtime**: Interactive analysis and visualization  

#### Purpose and Context
This final step creates publication-quality visualizations of the population structure results and performs advanced interpretation of the PCA data within the context of tick geography and ecology. The analysis integrates the quantitative PCA results from Step 16 with geographical metadata to reveal spatial population structure patterns across Iowa, Kansas, and Nebraska. This visualization enables biological interpretation and provides the foundation for understanding tick population dynamics across the study region.

#### Biological and Geographic Framework
The analysis is designed to reveal population structure patterns related to:

**Geographic distribution**: The study includes tick samples from three states (Iowa, Kansas, Nebraska) with Nebraska further subdivided into northern (Thurston County) and southern regions (Dodge, Douglas, Sarpy Counties).

**Ecological considerations**: Different geographic regions may represent distinct ecological niches, host communities, or environmental conditions that could drive population differentiation.

**Dispersal patterns**: Tick population structure can reveal information about dispersal limitations, host movement patterns, and barriers to gene flow across the landscape.

**Disease ecology implications**: Population structure affects pathogen transmission dynamics and the spatial distribution of tick-borne disease risk.

#### Technical Implementation and Data Integration
The Python script implements a comprehensive analysis workflow:

**Data integration architecture**:
- Merges PCA eigenscores from PLINK2 output with geographic metadata
- Calculates percentage of variance explained by each principal component
- Assigns samples to meaningful geographic groupings for visualization

**Statistical processing**:
- Loads eigenvalues to calculate proportion of variance explained by each PC
- Standardizes sample identifiers across datasets for accurate merging
- Implements geographic grouping logic based on state and regional subdivisions

**Visualization optimization**:
- Creates publication-quality scatter plots with distinct colors for each geographic group
- Includes sample sizes in legend for statistical interpretation
- Applies consistent styling and formatting for professional presentation

#### Geographic Grouping Strategy
The analysis implements a sophisticated geographic classification:

**Iowa samples**: Represented as a single geographic unit reflecting the state-level sampling
**Kansas samples**: Treated as a distinct geographic population
**Nebraska North**: Thurston County samples representing northern Nebraska ecology
**Nebraska South**: Combined Dodge, Douglas, and Sarpy County samples representing the Omaha metropolitan region and southern Nebraska

This grouping strategy balances statistical power with meaningful ecological and geographic distinctions.

#### Expected Population Structure Patterns
Different scenarios could emerge from the PCA analysis:

**Geographic clustering**: Samples clustering by state or region would suggest limited gene flow and population structure driven by geographic distance or barriers.

**Continuous variation**: Gradual transitions across geographic space would indicate isolation-by-distance patterns typical of species with limited dispersal.

**Panmixia**: Random distribution regardless of geography would suggest high gene flow and a single, well-connected population across the study region.

**Complex structure**: Multiple clusters not clearly aligned with geography might indicate historical demographic events, adaptation to different environments, or host-associated differentiation.

#### Input Data Structure
Takes multiple data sources:
- PCA eigenscores from Step 16 (`tick_pca.eigenvec` converted to Excel format)
- Eigenvalues for variance calculation (`tick_pca.eigenval` converted to Excel format)
- Geographic metadata linking samples to collection locations and regional classifications

#### Visualization Features and Interpretation Tools
The analysis produces publication-ready visualizations with:

**Statistical information**: Percentage of variance explained by each principal component, enabling assessment of the strength of population structure
**Sample size documentation**: Legend includes sample sizes for each geographic group, important for interpreting statistical significance
**Professional formatting**: Publication-quality graphics suitable for manuscripts and presentations
**Interactive capability**: Python environment allows for easy modification of color schemes, groupings, or additional analyses

#### Why This Step is Essential
1. **Biological interpretation**: Transforms quantitative PCA results into interpretable biological patterns
2. **Geographic context**: Reveals spatial patterns relevant to tick ecology and disease transmission
3. **Publication preparation**: Creates publication-quality figures for manuscripts and presentations
4. **Hypothesis generation**: Identifies patterns that warrant further investigation or experimental validation
5. **Applied implications**: Provides insights relevant to tick-borne disease management and surveillance

#### Expected Analytical Outcomes
The visualization should reveal:

**Population genetic insights**:
- Degree of population structure across the study region
- Relationship between genetic differentiation and geographic distance
- Identification of potential barriers to gene flow or areas of high connectivity

**Ecological implications**:
- Patterns relevant to tick dispersal biology and host relationships
- Information about population connectivity that affects disease transmission
- Insights into tick population dynamics across different ecological regions

#### Computational and Technical Considerations
- **Interactive analysis**: Python environment enables rapid iteration and parameter adjustment
- **Data format flexibility**: Handles multiple input formats and enables easy data manipulation
- **Visualization quality**: Professional-grade graphics suitable for publication
- **Reproducibility**: Well-documented code enables reproducible analysis and modification

#### Integration with Broader Research Context
The population structure results inform:
- **Disease ecology studies**: Understanding how population structure affects pathogen transmission
- **Conservation genetics**: Information relevant to tick population management
- **Evolutionary biology**: Insights into tick population dynamics and demographic history
- **Public health**: Spatial patterns relevant to tick-borne disease surveillance and control

#### Expected Biological Discoveries
The final analysis should provide insights into:
- **Tick population connectivity**: How well-connected are tick populations across the study region?
- **Geographic barriers**: Are there identifiable barriers to gene flow between regions?
- **Ecological adaptation**: Do different geographic regions harbor genetically distinct tick populations?
- **Disease implications**: How might population structure affect tick-borne pathogen transmission patterns?

This final step completes the transformation from raw sequencing data to biological understanding, providing publication-ready results that contribute to our understanding of tick population biology, ecology, and disease transmission dynamics. The visualization enables immediate interpretation of population genetic patterns while providing the foundation for further research into tick-borne disease ecology and management.

---

## Pipeline Completion Summary

**Comprehensive Analysis Achieved:**
This 17-step pipeline has successfully transformed 61 tick samples from raw NovaSeq sequencing reads into detailed population genetic insights. The complete workflow encompasses:

- **Preprocessing (Steps 1-8)**: Quality control, alignment, and data preparation
- **Variant Calling (Steps 9-11)**: GATK-based population variant discovery  
- **Filtering (Steps 12-15)**: Rigorous quality control and population optimization
- **Population Analysis (Steps 16-17)**: Structure analysis and biological interpretation

**Key Achievements:**
- Analysis-ready BAM files with comprehensive quality metrics
- Population-scale variant calling with joint genotyping across all samples
- Rigorous quality filtering following GATK best practices
- Population-optimized SNP dataset for robust genetic analysis
- Publication-quality population structure visualization with geographic context

**Biological Insights Generated:**
The pipeline provides the foundation for understanding tick population dynamics, geographic structure, disease ecology implications, and evolutionary relationships across Iowa, Kansas, and Nebraska populations.

---
