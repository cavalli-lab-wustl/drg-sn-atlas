#' ======================================================================
#' 13B_glia_dge_corrplots.R
#'
#' Per-cluster Pearson correlation of significant DGE gene sets across the
#' 4 injury conditions, at glia Minor resolution. For each glia cluster,
#' reads the injury-vs-control DGE CSVs produced by 13_glia_dge_go,
#' builds a presence/absence matrix of p_adj<0.05 genes across conditions,
#' computes Pearson correlation between condition columns, and saves a
#' heatmap colored by the cluster's canonical glia palette entry.
#'
#' Jaccard-based analogues live in 13C_glia_dge_jaccard.R.
#'
#' Inputs:  ../r_objects/DRGSNMI_glia.RDS (only for cluster ordering)
#'          ../dge/glia_named/DGE_<cluster>_<condition>.csv
#' Outputs: ../plots/glia_corrplots/corr_plot_<cluster>.pdf
#'
#' Author: Michael B. Thomsen, Ph.D.
#' Created: 2024-10-18
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


# ===== LOAD OBJECT =====
glia <- readRDS("../r_objects/DRGSNMI_glia.RDS")


# ===== OUTPUT DIRECTORY =====
output_directory <- "../plots/glia_corrplots"
dir.create(output_directory, recursive = TRUE)


# ===== CONDITIONS =====
conditions <- c("AIH", "DRC", "SCI", "SNC")


# ===== CELL TYPES (canonical glia Minor ordering) =====
# Intersect palette names with what's actually in the object to preserve
# the canonical order from variables.R while dropping absent labels.
clusters <- intersect(names(glia_cell_type_colors), unique(as.character(Idents(glia))))


# ===== CORRELATION LOOP =====
#
# Two-stage filter design (matches the original analysis):
#   - Stage 1 (strict, |log2FC|>1 & p_adj<0.05) identifies *robust* DEGs
#     in at least one condition; these form the universe of genes to study.
#   - Stage 2 (loose, p_adj<0.05 only) asks, for each of those robust
#     genes, whether it reaches *any* significance in each condition —
#     the binary vector used for the Pearson correlation.
# This captures "among robustly-changed genes, how consistently do they
# reach significance across conditions?" rather than the raw overlap of
# all-significant lists.

for (cluster in clusters) {

  # Stage 1: strict filter — collect robust DEGs per condition and their union
  condition_genes_strict <- list()
  condition_genes_loose  <- list()
  skip_cluster <- FALSE

  for (condition in conditions) {
    csv_path <- paste0("../dge/glia_named/DGE_", cluster, "_", condition, ".csv")
    if (!file.exists(csv_path)) {
      message("  Missing DGE for ", cluster, " x ", condition, " — skipping cluster")
      skip_cluster <- TRUE
      break
    }
    df <- read_csv(csv_path, show_col_types = FALSE)
    condition_genes_strict[[condition]] <- df %>%
      filter(p_val_adj < 0.05, abs(avg_log2FC) > 1) %>%
      pull(1)
    condition_genes_loose[[condition]] <- df %>%
      filter(p_val_adj < 0.05) %>%
      pull(1)
  }
  if (skip_cluster) next

  # Universe = union of strict-significant genes across all conditions
  unique_genes <- unique(unlist(condition_genes_strict))
  if (length(unique_genes) == 0) {
    message("  No robust DEGs (|log2FC|>1 & p_adj<0.05) for ", cluster, " — skipping")
    next
  }

  # Stage 2: presence/absence matrix — loose significance per condition
  corr_df <- data.frame(Genes = unique_genes, stringsAsFactors = FALSE)
  for (condition in conditions) {
    corr_df[[condition]] <- as.integer(corr_df$Genes %in% condition_genes_loose[[condition]])
  }

  # Pearson correlation between condition columns (phi coefficient on binary)
  corr_mat <- cor(corr_df[, -1], method = "pearson", use = "pairwise.complete.obs")
  melted_corr <- reshape2::melt(corr_mat)

  # Cluster-specific color; fall back to grey if missing from palette
  fill_hi <- glia_cell_type_colors[[cluster]]
  if (is.null(fill_hi)) fill_hi <- "grey40"

  p <- ggplot(melted_corr, aes(x = Var1, y = Var2)) +
    geom_tile(aes(fill = value), color = "white") +
    scale_fill_gradientn(colors = c("white", fill_hi), limits = c(0, 1)) +
    geom_text(aes(label = sprintf("%.2f", value)), vjust = 1, size = 3) +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(size = 10, color = "black"),
      axis.text.y  = element_text(size = 10, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )

  ggsave(
    filename = file.path(output_directory, paste0("corr_plot_", cluster, ".pdf")),
    plot = p, width = 3.3, height = 2.4
  )
}

message("13B_glia_dge_corrplots complete.")
