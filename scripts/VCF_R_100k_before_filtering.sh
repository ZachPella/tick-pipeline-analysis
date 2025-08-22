#!/bin/bash
#SBATCH --job-name=vcf_qc_plots_before
#SBATCH --output=vcf_qc_plots_before_%j.out
#SBATCH --error=vcf_qc_plots_before_%j.err
#SBATCH --time=10:00:00
#SBATCH --mem=30G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=batch

# --- Set Working Directory and Variables ---
WORKDIR=/work/fauverlab/zachpella/scripts_ticksJune2025_10_scatter/genotyping
# Updated to use the output from your "before filtering" VCFtools script
OUTPUT_PREFIX="out.100Ksubset_before.cohort_ticks_june2025_final"
R_SCRIPT_FILE="vcf_plots_before_temp_${SLURM_JOB_ID}.R"

# Change to the working directory
cd "${WORKDIR}"

# Check for a required input file from the VCFtools script
if [ ! -f "${OUTPUT_PREFIX}.lqual" ]; then
    echo "Error: VCFtools output file not found: ${WORKDIR}/${OUTPUT_PREFIX}.lqual"
    echo "Please ensure the VCFtools script ran successfully first."
    exit 1
fi

# --- Load R module ---
echo "Loading R module..."
module purge
module load R

# --- Set and create personal R library path ---
R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.2
mkdir -p "${R_LIBS_USER}"
export R_LIBS_USER

# --- Check for and install tidyverse if needed ---
echo "Checking for 'tidyverse' package..."
Rscript -e '
if (!require("tidyverse", quietly = TRUE)) {
    message("tidyverse not found. Installing now...")
    install.packages("tidyverse", repos = "http://cran.us.r-project.org", lib = Sys.getenv("R_LIBS_USER"))
    if (!require("tidyverse", quietly = TRUE)) {
        stop("Failed to install tidyverse. Please check the output for errors.")
    }
    message("tidyverse installed successfully.")
} else {
    message("tidyverse is already installed. Skipping installation.")
}
'

# --- Write the R script code to a temporary file ---
echo "Creating temporary R script file: ${R_SCRIPT_FILE}"
cat <<EOF > "${R_SCRIPT_FILE}"
## load packages
library(tidyverse)

## set working directory
setwd("${WORKDIR}")

# Graph mean quality (BEFORE filtering)
var_qual <- read_delim(paste0("${OUTPUT_PREFIX}", ".lqual"), delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1, show_col_types = FALSE)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + ggtitle("Variant Quality (Before Filtering)")
ggsave("qual_density_plot_before.png", a)

# Graph & summarize mean depth (BEFORE filtering)
var_depth <- read_delim(paste0("${OUTPUT_PREFIX}", ".ldepth.mean"), delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1, show_col_types = FALSE)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + ggtitle("Mean Depth (Before Filtering)")
ggsave("depth_density_plot_before.png", a)
summary(var_depth\$mean_depth)

# Graph & summarize mean depth with upper x-axis limit of 35 (BEFORE filtering)
var_depth_limited <- var_depth %>% filter(mean_depth <= 35)
a <- ggplot(var_depth_limited, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + ggtitle("Mean Depth - Limited to 35x (Before Filtering)")
ggsave("depth_density_plot_limit35_before.png", a)

# Graph & summarize variant missingness (BEFORE filtering)
var_miss <- read_delim(paste0("${OUTPUT_PREFIX}", ".lmiss"), delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1, show_col_types = FALSE)
a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + ggtitle("Variant Missingness (Before Filtering)")
ggsave("missingness_density_plot_before.png", a)
summary(var_miss\$fmiss)

# Graph & summarize minor allele frequency (BEFORE filtering)
var_freq <- read_delim(paste0("${OUTPUT_PREFIX}", ".frq"), delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "freq1", "freq2"), skip = 1, show_col_types = FALSE)
var_freq\$maf <- pmin(var_freq\$freq1, var_freq\$freq2)
a <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + ggtitle("Minor Allele Frequency (Before Filtering)")
ggsave("maf_density_plot_before.png", a)
summary(var_freq\$maf)

# Plot mean depth per individual (BEFORE filtering)
ind_depth <- read_delim(paste0("${OUTPUT_PREFIX}", ".idepth"), delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1, show_col_types = FALSE)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + ggtitle("Individual Depth (Before Filtering)")
ggsave("ind_depth_histogram_before.png", a)

# Plot proportion of missing data per individual (BEFORE filtering)
ind_miss <- read_delim(paste0("${OUTPUT_PREFIX}", ".imiss"), delim = "\t",
                       col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1, show_col_types = FALSE)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + ggtitle("Individual Missingness (Before Filtering)")
ggsave("ind_missingness_histogram_before.png", a)

# Plot heterozygosity and inbreeding coefficients per individual (BEFORE filtering)
ind_het <- read_delim(paste0("${OUTPUT_PREFIX}", ".het"), delim = "\t",
                      col_names = c("ind", "ho", "he", "nsites", "f"), skip = 1, show_col_types = FALSE)
a <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + ggtitle("Individual Inbreeding Coefficient (Before Filtering)")
ggsave("ind_inbreeding_histogram_before.png", a)
EOF

# --- Execute the temporary R script ---
echo "Running R script..."
Rscript "${R_SCRIPT_FILE}"

# --- Cleanup ---
echo "Cleaning up temporary script file."
rm "${R_SCRIPT_FILE}"

echo "Script complete. Check the working directory for output files and plots."
