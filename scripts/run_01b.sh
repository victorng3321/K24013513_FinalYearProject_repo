#!/bin/bash
#SBATCH -p cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH -J wgcna_uncorr
#SBATCH -o wgcna_uncorr_%j.out
# Task B: render the uncorrected-matrix WGCNA pipeline with module R 4.5.1 / Bioc 3.22
# (NOT the Singularity container). Submit with:  sbatch code/run_01b.sh
set -euo pipefail
cd /scratch/prj/ccc_dcis_ncrna/victor/data/code

# Load your WGCNA R module (swap for the exact name you normally use for R 4.5.1):
module load r/4.5.1 2>/dev/null || module load R/4.5.1 2>/dev/null || module load r 2>/dev/null || true
command -v Rscript >/dev/null || { echo "ERROR: no Rscript on PATH - load your R 4.5.1 module"; exit 1; }
echo "Using: $(Rscript --version 2>&1)"

# Bare module R has no pandoc, so we do NOT render to HTML. Instead tangle the Rmd
# to a plain R script and run it: every chunk executes and all CSV/PNG/RData outputs
# and console tables are produced, with no pandoc needed.
Rscript -e "knitr::purl('01b_WGCNA_uncorrected_HPC.Rmd', output='01b_tangled.R', documentation=0, quiet=TRUE); pdf(file.path('/scratch/prj/ccc_dcis_ncrna/victor/data/results/wgcna_uncorrected','01b_base_plots.pdf'), width=14, height=10); source('01b_tangled.R', echo=FALSE); dev.off()"
echo ">>> 01b execution complete"
