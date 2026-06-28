#!/usr/bin/env bash
# Run an R/Rmd script inside the DCIS isoform+network Singularity container.
# Usage: bash code/run_in_container.sh 04_isoform_switch_HPC.Rmd
set -euo pipefail
# Find the container runtime wherever it lives (PATH, module-loaded, or apptainer).
# Override by exporting SING=/full/path/to/singularity before running.
SING="${SING:-$(command -v singularity || command -v apptainer || true)}"
if [ -z "$SING" ]; then
  echo "ERROR: no 'singularity' or 'apptainer' on PATH." >&2
  echo "  Try: module load apptainer   (or singularity), then re-run." >&2
  echo "  And run on a compute node (srun ...), not the login node." >&2
  exit 127
fi
REAL=$(readlink -f /scratch/prj/ccc_dcis_ncrna/victor/data)
IMG="$REAL/bioconductor_3_19.sif"
export SINGULARITY_TMPDIR="/tmp/sing_$USER"; mkdir -p "$SINGULARITY_TMPDIR"
export SINGULARITYENV_R_LIBS_USER="$REAL/code/container_Rlib"
export SINGULARITYENV_DCIS_BASE="$REAL"          # config.R reads this; avoids the /scratch symlink
script="$1"
case "$script" in
  *.Rmd) cmd="rmarkdown::render('$script')" ;;
  *)     cmd="source('$script')" ;;
esac
"$SING" exec --bind /cephfs --pwd "$REAL/code" "$IMG" \
  Rscript -e ".libPaths(Sys.getenv('R_LIBS_USER')); $cmd"
