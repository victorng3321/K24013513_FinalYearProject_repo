source("config.R")

OUT_DIR <- file.path(RESULTS, "isoform_qc98")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

rds <- file.path(INTERMEDIATE, "switchListAnalyzed_qc98.rds")
stopifnot(file.exists(rds))
sla <- readRDS(rds)
results <- sla$isoformFeatures

cat("QC-98 run — testable space:\n")
cat("  genes   :", length(unique(results$gene_id)), "\n")
cat("  isoforms:", length(unique(results$isoform_id)), "\n\n")

# Isoform-level significant switches (same thresholds as the 154-sample primary)
sig <- results[!is.na(results$isoform_switch_q_value) &
               results$isoform_switch_q_value < 0.05 &
               abs(results$dIF) > 0.10, ]
write.csv(results, file.path(OUT_DIR, "results_all_isoforms.csv"), row.names = FALSE)
write.csv(sig,     file.path(OUT_DIR, "significant_switches.csv"), row.names = FALSE)
cat("Significant switches (isoform-level, q<0.05 & |dIF|>0.10):", nrow(sig), "\n\n")

# Replication of the 7 primary switches
seven <- c("TOP2A","HLA-B","ACTR3","RTKN2","TMEM59","SLC4A1AP","ENSG00000264968")
chk <- results[results$gene_name %in% seven,
               c("gene_name","isoform_id","dIF","isoform_switch_q_value")]
chk <- chk[order(chk$gene_name, -abs(chk$dIF)), ]
cat("=== The 7 primary switches — their status in the QC-98 run ===\n")
print(chk, row.names = FALSE)

replicating <- sort(unique(sig$gene_name[sig$gene_name %in% seven]))
cat("\nReplicate at q<0.05 & |dIF|>0.10:",
    if (length(replicating)) paste(replicating, collapse = ", ") else "NONE", "\n")
cat("Replicating:", length(replicating), "of 7\n")
