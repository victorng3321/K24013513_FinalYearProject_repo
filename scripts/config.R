# config.R — single source of paths and helpers for the DCIS project.
# Relocate everything by editing BASEDIR (or setting the $DCIS_BASE env var).

BASEDIR <- Sys.getenv("DCIS_BASE", "/scratch/prj/ccc_dcis_ncrna/victor/data")

# ---- Inputs (controlled-access; not for GitHub) ----
DIR_COUNTS <- file.path(BASEDIR, "counts")
DIR_CLIN   <- file.path(BASEDIR, "clinical")
DIR_SALMON <- file.path(BASEDIR, "salmon")
DIR_QC     <- file.path(BASEDIR, "qc")
DIR_REFS   <- file.path(BASEDIR, "refs")

# ---- Outputs ----
RESULTS     <- file.path(BASEDIR, "results")
DIR_WGCNA   <- file.path(RESULTS, "wgcna")
DIR_ENRICH  <- file.path(RESULTS, "enrichment")
DIR_CCDEG   <- file.path(RESULTS, "casecontrol_deg")
DIR_PAIRED  <- file.path(RESULTS, "paired_deg")
DIR_ISOFORM <- file.path(RESULTS, "isoform")
DIR_NETFIG  <- file.path(DIR_WGCNA, "network_figures")

# Curated final deliverables (named by thesis figure/table number)
DIR_FIG <- file.path(RESULTS, "figures")
DIR_TAB <- file.path(RESULTS, "tables")

# Large intermediate objects (git-ignored)
INTERMEDIATE <- file.path(BASEDIR, "intermediate")

# Create every output directory if missing
for (d in c(DIR_WGCNA, DIR_ENRICH, DIR_CCDEG, DIR_PAIRED, DIR_ISOFORM,
            DIR_NETFIG, DIR_FIG, DIR_TAB, INTERMEDIATE))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- Helpers: save a final figure / table by thesis ID ----
# e.g. save_fig(p_dendrogram, "Fig3.1a_dendrogram")
#      save_tab(recurrence_modules, "Table3.1_recurrence_modules")
save_fig <- function(p, id, w = 180, h = 120, dpi = 300)
  ggplot2::ggsave(file.path(DIR_FIG, paste0(id, ".png")), p,
                  width = w, height = h, units = "mm", dpi = dpi)

save_tab <- function(df, id)
  utils::write.csv(df, file.path(DIR_TAB, paste0(id, ".csv")), row.names = FALSE)

# ---- Reproducibility ----
set.seed(42)
