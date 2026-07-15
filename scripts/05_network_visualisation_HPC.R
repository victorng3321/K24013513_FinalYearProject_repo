source("config.R")

### !!!!! Note this part was conducted on local terminal !!!! ###

# plot_wgcna_networks 
#(biotype = shape, kIN = colour + size)


library(igraph); library(ggraph); library(ggrepel); library(scales); library(ggplot2)



######################## 1. PARAMETERS #######################
IN_DIR <- DIR_WGCNA
OUT_DIR <- DIR_NETFIG

MODULES <- "all"                 

TARGET_MEAN_DEGREE <- 8          
RESCUE_ISOLATES    <- TRUE

KIN_LOW  <- "#DEEBF7"             
KIN_HIGH <- "#084594"            

SIZE_RANGE   <- c(3, 12)          
NODE_BORDER  <- "grey25"          
BORDER_STROKE<- 0.5

SHAPE_VALUES <- c(protein_coding = 21, lncRNA = 22)   # 21 circle, 22 square

SKIP_ENSG_LABELS <- TRUE          # don't print unnamed ENSG... loci
LABEL_MIN_KIN_Q  <- 0             # 0 = label all named genes; e.g. 0.5 = top-half kIN only
FONT_SIZE        <- 2.7

PDF_W <- 9; PDF_H <- 8; PNG_DPI <- 300
SEED  <- 42


dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
hub_sum     <- read.csv(file.path(IN_DIR, "hub_summary.csv"),   stringsAsFactors = FALSE)
gene_info   <- read.csv(file.path(IN_DIR, "geneInfo_full.csv"), stringsAsFactors = FALSE)
biotype_map <- setNames(ifelse(gene_info$gene_biotype == "lncRNA", "lncRNA", "protein_coding"),
                        gene_info$geneSymbol)

prune_to_degree <- function(edges, node_ids, target_degree, rescue) {
  n <- length(node_ids); k <- max(1, round(n * target_degree / 2))
  kept <- edges[head(order(edges$weight, decreasing = TRUE), min(k, nrow(edges))), , drop = FALSE]
  if (rescue) {
    orphans <- setdiff(node_ids, unique(c(kept$from, kept$to)))
    if (length(orphans) > 0) {
      add <- do.call(rbind, lapply(orphans, function(g) {
        cand <- edges[edges$from == g | edges$to == g, , drop = FALSE]
        if (nrow(cand) == 0) NULL else cand[which.max(cand$weight), , drop = FALSE]
      }))
      kept <- unique(rbind(kept, add))
    }
  }
  kept
}

render_module <- function(mod) {
  ef <- file.path(IN_DIR, paste0("Cytoscape_edges_", mod, ".txt"))
  nf <- file.path(IN_DIR, paste0("Cytoscape_nodes_", mod, ".txt"))
  if (!file.exists(ef) || !file.exists(nf)) { cat("  [skip]", mod, "- missing pair\n"); return(invisible()) }
  
  edges <- read.delim(ef, stringsAsFactors = FALSE, check.names = FALSE)
  nodes <- read.delim(nf, stringsAsFactors = FALSE, check.names = FALSE)
  names(edges)[names(edges) == "fromNode"] <- "from"
  names(edges)[names(edges) == "toNode"]   <- "to"
  edges$weight <- as.numeric(edges$weight)
  names(nodes)[names(nodes) == "nodeName"] <- "name"
  kin_lookup <- setNames(hub_sum$intramod_conn[hub_sum$module == mod],
                         hub_sum$gene_symbol[hub_sum$module == mod])
  nodes$kIN  <- as.numeric(kin_lookup[nodes$name])
  if (all(is.na(nodes$kIN))) nodes$kIN <- 1
  nodes$kIN[is.na(nodes$kIN)] <- min(nodes$kIN, na.rm = TRUE)
  nodes$biotype <- biotype_map[nodes$name]
  nodes$biotype[is.na(nodes$biotype)] <- "protein_coding"
  if (nrow(nodes) < 2 || nrow(edges) < 1) { cat("  [skip]", mod, "- too small\n"); return(invisible()) }
  
  edges <- prune_to_degree(edges, nodes$name, TARGET_MEAN_DEGREE, RESCUE_ISOLATES)
  
  nodes$label <- nodes$name
  if (SKIP_ENSG_LABELS) nodes$label[grepl("^ENSG", nodes$name)] <- ""
  if (LABEL_MIN_KIN_Q > 0) {
    thr <- as.numeric(quantile(nodes$kIN, LABEL_MIN_KIN_Q, names = FALSE))
    nodes$label[nodes$kIN < thr] <- ""
  }
  
  nodes <- nodes[, c("name", setdiff(names(nodes), "name"))]
  g <- graph_from_data_frame(edges[, c("from", "to", "weight")], directed = FALSE, vertices = nodes)
  
  set.seed(SEED)
  lay <- layout_with_fr(g, weights = E(g)$weight)
  cat(sprintf("  %-16s %3d nodes | %4d edges (mean degree %.1f)\n",
              mod, vcount(g), ecount(g), 2 * ecount(g) / vcount(g)))
  
  p <- ggraph(g, layout = "manual", x = lay[, 1], y = lay[, 2]) +
    geom_edge_link(colour = "grey75", width = 0.18, alpha = 0.22) +
    geom_node_point(aes(fill = kIN, size = kIN, shape = biotype),
                    colour = NODE_BORDER, stroke = BORDER_STROKE) +
    geom_node_text(aes(label = label), repel = TRUE, size = FONT_SIZE,
                   max.overlaps = Inf, colour = "grey15",
                   bg.color = "white", bg.r = 0.15,
                   segment.size = 0.2, segment.colour = "grey70", min.segment.length = 0) +
    scale_fill_gradient(low = KIN_LOW, high = KIN_HIGH,
                        name = "Intramodular\nconnectivity (kIN)") +
    scale_size(range = SIZE_RANGE, guide = "none") +
    scale_shape_manual(values = SHAPE_VALUES, name = "Biotype",
                       labels = c(protein_coding = "protein-coding", lncRNA = "lncRNA")) +
    guides(fill = guide_colourbar(barheight = 6),
           shape = guide_legend(override.aes = list(fill = "grey70", size = 4,
                                                    colour = NODE_BORDER, stroke = BORDER_STROKE))) +
    labs(title = paste0(mod, " module")) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "right")
  
  base <- file.path(OUT_DIR, paste0("network_", mod))
  ggsave(paste0(base, ".pdf"), p, width = PDF_W, height = PDF_H, device = cairo_pdf)
  ggsave(paste0(base, ".png"), p, width = PDF_W, height = PDF_H, dpi = PNG_DPI, bg = "white")
}


####################### 3. DISCOVER + RENDER #######################
edge_files <- list.files(IN_DIR, pattern = "^Cytoscape_edges_.*\\.txt$", full.names = FALSE)
present <- sub("^Cytoscape_edges_(.*)\\.txt$", "\\1", edge_files)
if (!identical(MODULES, "all")) present <- intersect(present, MODULES)
if (length(present) == 0) stop("No Cytoscape_edges_*.txt files found in ", IN_DIR)
cat("Rendering", length(present), "module(s):\n  ", paste(present, collapse = ", "), "\n\n")

for (mod in present)
  tryCatch(render_module(mod),
           error = function(e) cat("  [ERROR]", mod, "-", conditionMessage(e), "\n"))

cat("\nDone. Figures in:", OUT_DIR, "\n")
