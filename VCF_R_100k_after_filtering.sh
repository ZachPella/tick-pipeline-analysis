#!/bin/bash
#SBATCH --job-name=vcf_qc_plots_after
#SBATCH --output=vcf_qc_plots_after_%j.out
#SBATCH --error=vcf_qc_plots_after_%j.err
#SBATCH --time=10:00:00
#SBATCH --mem=30G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=batch

# --- Set Working Directory and Variables ---
WORKDIR=/work/fauverlab/zachpella/scripts_ticksJune2025_10_scatter/genotyping
# Updated to use the output from your VCFtools script
OUTPUT_PREFIX="out.100Ksubset_after.cohort_ticks_june2025_snps_filtered_only"
R_SCRIPT_FILE="vcf_plots_after_temp_${SLURM_JOB_ID}.R"

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

# Graph mean quality (AFTER filtering)
var_qual <- read_delim(paste0("${OUTPUT_PREFIX}", ".lqual"), delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1, show_col_types = FALSE)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "orange", colour = "black", alpha = 0.3)
a + theme_light() + ggtitle("Variant Quality (After Filtering)")
ggsave("qual_density_plot_after.png", a)

# Graph & summarize mean depth (AFTER filtering)
var_depth <- read_delim(paste0("${OUTPUT_PREFIX}", ".ldepth.mean"), delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1, show_col_types = FALSE)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "orange", colour = "black", alpha = 0.3)
a + theme_light() + ggtitle("Mean Depth (After Filtering)")
ggsave("depth_density_plot_after.png", a)
summary(var_depth\$mean_depth)

# Graph & summarize mean depth with upper x-axis limit of 35 (AFTER filtering)
var_depth_limited <- var_depth %>% filter(mean_depth <= 35)
a <- ggplot(var_depth_limited, aes(mean_depth)) + geom_density(fill = "orange", colour = "black", alpha = 0.3)
a + theme_light() + ggtitle("Mean Depth - Limited to 35x (After Filtering)")
ggsave("depth_density_plot_limit35_after.png", a)

# Graph & summarize variant missingness (AFTER filtering)
var_miss <- read_delim(paste0("${OUTPUT_PREFIX}", ".lmiss"), delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1, show_col_types = FALSE)
a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "orange", colour = "black", alpha = 0.3)
a + theme_light() + ggtitle("Variant Missingness (After Filtering)")
ggsave("missingness_density_plot_after.png", a)
summary(var_miss\$fmiss)

# Graph & summarize minor allele frequency (AFTER filtering)
var_freq <- read_delim(paste0("${OUTPUT_PREFIX}", ".frq"), delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "freq1", "freq2"), skip = 1, show_col_types = FALSE)
var_freq\$maf <- pmin(var_freq\$freq1, var_freq\$freq2)
a <- ggplot(var_freq, aes(maf)) + geom_density(fill = "orange", colour = "black", alpha = 0.3)
a + theme_light() + ggtitle("Minor Allele Frequency (After Filtering)")
ggsave("maf_density_plot_after.png", a)
summary(var_freq\$maf)

# Plot mean depth per individual (AFTER filtering)
ind_depth <- read_delim(paste0("${OUTPUT_PREFIX}", ".idepth"), delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1, show_col_types = FALSE)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "orange", colour = "black", alpha = 0.3)
a + theme_light() + ggtitle("Individual Depth (After Filtering)")
ggsave("ind_depth_histogram_after.png", a)

# Plot proportion of missing data per individual (AFTER filtering)
ind_miss <- read_delim(paste0("${OUTPUT_PREFIX}", ".imiss"), delim = "\t",
                       col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1, show_col_types = FALSE)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "orange", colour = "black", alpha = 0.3)
a + theme_light() + ggtitle("Individual Missingness (After Filtering)")
ggsave("ind_missingness_histogram_after.png", a)

# Plot heterozygosity and inbreeding coefficients per individual (AFTER filtering)
ind_het <- read_delim(paste0("${OUTPUT_PREFIX}", ".het"), delim = "\t",
                      col_names = c("ind", "ho", "he", "nsites", "f"), skip = 1, show_col_types = FALSE)
a <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "orange", colour = "black", alpha = 0.3)
a + theme_light() + ggtitle("Individual Inbreeding Coefficient (After Filtering)")
ggsave("ind_inbreeding_histogram_after.png", a)
EOF

# --- Execute the temporary R script ---
echo "Running R script..."
Rscript "${R_SCRIPT_FILE}"

# --- Cleanup ---
echo "Cleaning up temporary script file."
rm "${R_SCRIPT_FILE}"

echo "Script complete. Check the working directory for output files and plots."
