# =============================================================================
# 06_hub_gene_barplots.R
# Hub-gene bar plots for the recurrence-associated WGCNA modules (Figure 3.6).
# Genes are ranked by intramodular connectivity (kIN = kWithin; Langfelder &
# Horvath, 2008). One plot per module, coloured by the module's WGCNA colour;
# a marker at each bar tip encodes biotype (circle = protein-coding, square = lncRNA).
#
# INPUT   <DIR_WGCNA>/hub_summary.csv, <DIR_WGCNA>/geneInfo_full.csv
# OUTPUT  <DIR_WGCNA>/hub_plots/hub_<module>.{pdf,png} + hub_recurrence_stacked.{pdf,png}
# =============================================================================

source("config.R")
library(readr); library(dplyr); library(ggplot2); library(forcats); library(patchwork)

data_dir     <- DIR_WGCNA
out_dir      <- file.path(DIR_WGCNA, "hub_plots")
n_top        <- 10
modules      <- c("paleturquoise", "violet", "darkorange")
label_values <- FALSE
pdf_device   <- cairo_pdf
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

hub <- read_csv(file.path(data_dir, "hub_summary.csv"), show_col_types = FALSE)
biotype <- read_csv(file.path(data_dir, "geneInfo_full.csv"), show_col_types = FALSE) |>
  select(gene_id = ensembl_id, gene_biotype)
hub <- hub |>
  left_join(biotype, by = "gene_id") |>
  mutate(Biotype     = factor(ifelse(gene_biotype == "lncRNA", "lncRNA", "Protein-coding"),
                              levels = c("Protein-coding", "lncRNA")),
         gene_symbol = ifelse(is.na(gene_symbol) | gene_symbol == "", gene_id, gene_symbol))

darken <- function(col, f = 0.55) { v <- col2rgb(col) * f; rgb(v[1], v[2], v[3], maxColorValue = 255) }
shape_map <- c("Protein-coding" = 21, "lncRNA" = 22)

plot_hub_module <- function(mod, n = n_top, show_values = label_values) {
  df <- hub |> filter(module == mod) |> slice_min(order_by = rank, n = n) |>
    mutate(gene_symbol = fct_reorder(gene_symbol, intramod_conn))
  ref <- data.frame(intramod_conn = df$intramod_conn[1], gene_symbol = df$gene_symbol[1],
                    Biotype = factor(c("Protein-coding", "lncRNA"),
                                     levels = c("Protein-coding", "lncRNA")))
  p <- ggplot(df, aes(intramod_conn, gene_symbol)) +
    geom_col(fill = mod, colour = "grey20", linewidth = 0.35, width = 0.72) +
    geom_point(data = ref, aes(shape = Biotype), size = 2.6, stroke = 0.6,
               fill = "white", colour = "black", alpha = 0) +
    geom_point(aes(shape = Biotype), size = 2.6, stroke = 0.6, fill = "white", colour = "black") +
    scale_shape_manual(values = shape_map, drop = FALSE, name = "Biotype") +
    guides(shape = guide_legend(override.aes = list(alpha = 1, fill = "white", colour = "black", size = 2.8))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = mod, x = "Intramodular connectivity (kIN)", y = NULL) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", colour = darken(mod), hjust = 0),
          axis.text.y = element_text(face = "italic", colour = "black"),
          axis.line = element_line(linewidth = 0.4),
          legend.position = "bottom", legend.key.size = unit(0.9, "lines"))
  if (show_values) p <- p + geom_text(aes(label = sprintf("%.2f", intramod_conn)),
                                      hjust = -0.55, size = 3, colour = "grey25")
  p
}

plots <- setNames(lapply(modules, plot_hub_module), modules)
for (m in modules) {
  ggsave(file.path(out_dir, paste0("hub_", m, ".pdf")), plots[[m]], width = 5.0, height = 3.6, device = pdf_device)
  ggsave(file.path(out_dir, paste0("hub_", m, ".png")), plots[[m]], width = 5.0, height = 3.6, dpi = 300, bg = "white")
}
if (length(modules) > 1) {
  stacked <- wrap_plots(plots, ncol = 1) + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold"), legend.position = "bottom")
  ggsave(file.path(out_dir, "hub_recurrence_stacked.pdf"), stacked,
         width = 6.5, height = 3.1 * length(modules), device = pdf_device)
  ggsave(file.path(out_dir, "hub_recurrence_stacked.png"), stacked,
         width = 6.5, height = 3.1 * length(modules), dpi = 300, bg = "white")
}
message("Done. Figures written to: ", normalizePath(out_dir))
