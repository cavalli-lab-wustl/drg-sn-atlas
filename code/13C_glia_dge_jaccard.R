#' ======================================================================
#' 13C_glia_dge_jaccard.R
#'
#' Direction-specific Jaccard similarity heatmaps of tissue-stratified
#' glial DGE, per glia Minor cluster. For each cluster and each tissue
#' (DRG, SN), builds three Jaccard matrices (UP / DOWN / ALL genes
#' passing p_adj<0.05 and |log2FC|>=logFC_threshold) across the 4 injury
#' conditions and saves each as a ComplexHeatmap PDF.
#'
#' A single combined DEG-count summary (both tissues) is written at the
#' end. Some clusters may be absent or too sparse in SN — wrapped in
#' tryCatch so missing pairs log a message and continue.
#'
#' Inputs:  ../r_objects/DRGSNMI_glia.RDS (only for cluster ordering)
#'          ../dge/glia_named/{DRG,SN}/DGE_<cluster>_<condition>.csv
#' Outputs: ../plots/glia_jaccard/{DRG,SN}/Jaccard_<cluster>_{UP,DOWN,ALL}_logFC<t>.pdf
#'          ../plots/glia_jaccard/DEG_summary_logFC<t>.csv
#'
#' Author: Michael B. Thomsen, Ph.D.
#' Created: 2024-07-01
#' ======================================================================


# Set working directory to script location (portable: Rscript or RStudio)
local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) setwd(dirname(normalizePath(sub("^--file=", "", f))))
  else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    setwd(dirname(rstudioapi::getSourceEditorContext()$path))
})

# Imports
source('utils/libraries.R')
source('utils/helper_functions.R')
source('utils/variables.R')
library(ComplexHeatmap)
library(circlize)
library(grid)


# ===== LOAD OBJECT (only for cluster ordering) =====
glia <- readRDS("../r_objects/DRGSNMI_glia.RDS")


# ===== CONFIG =====
injury_order    <- c("AIH", "DRC", "SCI", "SNC")
p_val_cutoff    <- 0.05
logFC_threshold <- 1   # stringent — robust changes only

# Canonical Minor ordering, intersected with what's actually in the object
clusters <- intersect(names(glia_cell_type_colors), unique(as.character(Idents(glia))))


# ===== HELPERS =====

# Extract a gene set from a DGE CSV with direction-specific filtering.
# Returns character(0) if the file is missing or no genes pass.
get_gene_set <- function(file, direction = c("up", "down", "both"),
                         logFC_col = "avg_log2FC",
                         padj_col  = "p_val_adj",
                         p_val_cutoff = 0.05,
                         logFC_threshold = 1) {
  direction <- match.arg(direction)
  if (!file.exists(file)) return(character(0))

  df <- read.csv(file, row.names = 1)
  df <- df[!is.na(df[[padj_col]]) & df[[padj_col]] < p_val_cutoff, , drop = FALSE]

  df <- switch(direction,
    up   = df[df[[logFC_col]] >=  logFC_threshold, , drop = FALSE],
    down = df[df[[logFC_col]] <= -logFC_threshold, , drop = FALSE],
    both = df[abs(df[[logFC_col]]) >= logFC_threshold, , drop = FALSE]
  )
  rownames(df)
}


# Build the Jaccard heatmap for one cluster + direction + tissue.
create_jaccard_heatmap <- function(cluster, gene_direction, dge_dir) {

  # Collect gene sets for each injury condition
  gene_sets <- lapply(injury_order, function(inj) {
    get_gene_set(
      file.path(dge_dir, paste0("DGE_", cluster, "_", inj, ".csv")),
      direction       = gene_direction,
      p_val_cutoff    = p_val_cutoff,
      logFC_threshold = logFC_threshold
    )
  })
  names(gene_sets) <- injury_order
  set_sizes <- vapply(gene_sets, length, integer(1))

  # If nothing to compare across conditions, bail
  if (sum(set_sizes) == 0) {
    return(list(ht = NULL, set_sizes = set_sizes))
  }

  # Build Jaccard matrix
  n_inj <- length(injury_order)
  jaccard_mat <- matrix(
    0, n_inj, n_inj,
    dimnames = list(injury_order, injury_order)
  )
  for (i in seq_len(n_inj)) {
    for (j in seq_len(n_inj)) {
      intersection <- length(intersect(gene_sets[[i]], gene_sets[[j]]))
      union_size   <- length(unique(c(gene_sets[[i]], gene_sets[[j]])))
      jaccard_mat[i, j] <- if (union_size == 0) 0 else intersection / union_size
    }
  }

  # Cluster-specific color, with grey fallback
  base_col <- glia_cell_type_colors[[cluster]]
  if (is.null(base_col)) base_col <- "#888888"

  # For the "down" direction only, shift hue to distinguish from "up"
  if (gene_direction == "down") {
    col_rgb <- col2rgb(base_col) / 255
    base_col <- rgb(
      min(col_rgb[1, ] + 0.3, 1),
      max(col_rgb[2, ] - 0.2, 0),
      max(col_rgb[3, ] - 0.2, 0)
    )
  }

  jaccard_col_fun <- colorRamp2(c(0, 1), c("white", base_col))

  cell_fun_labels <- function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", jaccard_mat[i, j]), x, y,
              gp = gpar(fontsize = 10, col = "black"))
  }

  ha_col <- HeatmapAnnotation(
    `DEG count` = anno_barplot(
      set_sizes, border = FALSE,
      gp     = gpar(fill = base_col, col = NA),
      height = unit(1.5, "cm")
    ),
    show_annotation_name = TRUE
  )

  direction_suffix <- c(up = "UP", down = "DOWN", both = "ALL")[[gene_direction]]

  ht <- Heatmap(
    jaccard_mat,
    name              = paste0("Jaccard (", direction_suffix, ")"),
    col               = jaccard_col_fun,
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    row_names_side    = "left",
    column_names_rot  = 45,
    top_annotation    = ha_col,
    heatmap_legend_param = list(
      title  = "Jaccard Index",
      at     = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "0.25", "0.5", "0.75", "1")
    ),
    cell_fun = cell_fun_labels
  )

  list(ht = ht, set_sizes = set_sizes)
}


# ===== MAIN LOOP: TISSUE x CLUSTER x DIRECTION =====

summary_stats <- data.frame()

for (tissue_label in c("DRG", "SN")) {
  cat("\n===== Tissue:", tissue_label, "=====\n")

  dge_dir <- file.path("../dge/glia_named", tissue_label)
  if (!dir.exists(dge_dir)) {
    message("  No DGE dir for tissue ", tissue_label, " — skipping (was phase 4 run?)")
    next
  }
  out_dir <- file.path("../plots/glia_jaccard", tissue_label)
  dir.create(out_dir, recursive = TRUE)

  for (cluster in clusters) {
    cat("  ", cluster, "\n", sep = "")

    for (direction in c("up", "down", "both")) {
      tryCatch({
        result <- create_jaccard_heatmap(cluster, direction, dge_dir = dge_dir)

        # Record DEG counts regardless of whether we plot
        summary_stats <- rbind(summary_stats, data.frame(
          Tissue    = tissue_label,
          Cluster   = cluster,
          Direction = direction,
          AIH = result$set_sizes[["AIH"]],
          DRC = result$set_sizes[["DRC"]],
          SCI = result$set_sizes[["SCI"]],
          SNC = result$set_sizes[["SNC"]],
          stringsAsFactors = FALSE
        ))

        if (is.null(result$ht)) {
          cat("    [", direction, "] no DEGs across any condition — no plot\n", sep = "")
          next
        }

        direction_suffix <- c(up = "UP", down = "DOWN", both = "ALL")[[direction]]
        pdf_path <- file.path(
          out_dir,
          paste0("Jaccard_", cluster, "_", direction_suffix, "_logFC", logFC_threshold, ".pdf")
        )
        pdf(pdf_path, width = 3, height = 3)
        draw(result$ht, annotation_legend_side = "bottom")
        dev.off()

      }, error = function(e) {
        cat("    [", direction, "] error — ", e$message, "\n", sep = "")
        if (dev.cur() != 1) dev.off()
      })
    }
  }
}


# ===== SAVE COMBINED SUMMARY =====
summary_dir <- "../plots/glia_jaccard"
dir.create(summary_dir, recursive = TRUE)
write.csv(
  summary_stats,
  file = file.path(summary_dir, paste0("DEG_summary_logFC", logFC_threshold, ".csv")),
  row.names = FALSE
)

cat("\nAnalysis parameters:\n")
cat("  p-value cutoff:  ", p_val_cutoff, "\n")
cat("  log2FC threshold:", logFC_threshold, "\n")
cat("  Output:          ", summary_dir, "\n")

message("\n13C_glia_dge_jaccard complete.")
