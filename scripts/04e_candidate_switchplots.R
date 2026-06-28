# 04e_candidate_switchplots.R
# Switch plots for the 4 pre-registered candidate genes (ESR1, PGR, ERBB2, AR)
# from the QC-98 analysed object. Rebuilt from switchPlot's component panels so
# the wasted right-hand legend column is removed and the lower panels fill the
# frame. Falls back to the standard switchPlot() if the custom layout errors.
source("config.R")
suppressMessages({
  library(IsoformSwitchAnalyzeR); library(ggplot2)
  library(gridExtra); library(grid)
})

OUT_DIR <- file.path(RESULTS, "isoform_qc98")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
sla <- readRDS(file.path(INTERMEDIATE, "switchListAnalyzed_qc98.rds"))

C1 <- "non_progressor"; C2 <- "progressor"
candidates <- c("ESR1", "PGR", "ERBB2", "AR")

make_plot <- function(g, f) {
  ok <- FALSE
  tryCatch({
    tx <- switchPlotTranscript(sla, gene = g) + theme(legend.position = "none")
    ge <- switchPlotGeneExp(sla, gene = g, condition1 = C1, condition2 = C2) +
            theme(legend.position = "none")
    ie <- switchPlotIsoExp(sla, gene = g, condition1 = C1, condition2 = C2) +
            theme(legend.position = "none")
    iu <- switchPlotIsoUsage(sla, gene = g, condition1 = C1, condition2 = C2) +
            theme(legend.position = "bottom")
    png(f, width = 10, height = 7.5, units = "in", res = 200)
    grid.arrange(tx, ge, ie, iu,
                 layout_matrix = rbind(c(1, 1, 1), c(2, 3, 4)),
                 heights = c(1.15, 1),
                 top = textGrob(paste0("Isoform usage in ", g,
                                       " (non_progressor vs progressor)"),
                                gp = gpar(fontface = "bold", fontsize = 14)))
    dev.off(); ok <- TRUE
    cat("plotted", g, "(custom layout)\n")
  }, error = function(e) {
    if (length(dev.list())) dev.off()
    cat(g, "custom layout failed:", conditionMessage(e), "\n")
  })
  if (!ok) {
    tryCatch({
      png(f, width = 9, height = 6.5, units = "in", res = 200)
      switchPlot(sla, gene = g)
      dev.off()
      cat("plotted", g, "(standard switchPlot fallback)\n")
    }, error = function(e) {
      if (length(dev.list())) dev.off()
      cat(g, "failed entirely:", conditionMessage(e), "\n")
    })
  }
}

for (g in candidates)
  make_plot(g, file.path(OUT_DIR, paste0("switch_candidate_", g, ".png")))
cat("done\n")
