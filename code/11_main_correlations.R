#' ======================================================================
#' 11_main_correlations.R
#'
#' Per-celltype cross-condition correlation of significant DGE gene sets.
#' For each Major cell type, reads the 4 injury-vs-control DGE CSVs
#' produced by 10_main_dge_GO, builds a presence/absence matrix of
#' p_adj<0.05 genes across conditions, computes Pearson correlation
#' between condition columns, and saves a heatmap colored by the cell
#' type's canonical palette entry.
#'
#' Inputs:  ../r_objects/DRGSNMI_main_multilevel.RDS (only for Idents)
#'          ../dge/combined_named/DGE_<celltype>_<condition>.csv
#' Outputs: ../plots/corrplots/corr_plot_<celltype>.pdf
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


# ===== LOAD OBJECT (only to establish canonical cell-type ordering) =====
merged_named <- readRDS("../r_objects/DRGSNMI_main_multilevel.RDS")
Idents(merged_named) <- merged_named$Major


# ===== OUTPUT DIRECTORY =====
output_directory <- "../plots/corrplots"
dir.create(output_directory, recursive = TRUE)


# ===== CONDITIONS =====
conditions <- c(
  "cavalli_aih",
  "cavalli_drc",
  "cavalli_sci",
  "cavalli_snc"
)


# ===== CELL TYPES (canonical class-grouped order) =====
clusters <- intersect(class_ordered_majors, unique(as.character(Idents(merged_named))))


# ===== CLEAN LABEL HELPER =====
# Condition strings look like "cavalli_aih" → displays as "AIH"
clean_label <- function(x) {
  sapply(strsplit(x, "_"), function(v) toupper(v[2]))
}


# ===== CORRELATION LOOP =====
for (cluster in clusters) {

  # Read & filter all 4 condition DGEs for this cell type (skip missing)
  condition_genes <- list()
  skip_cluster <- FALSE
  for (condition in conditions) {
    csv_path <- paste0("../dge/combined_named/DGE_", cluster, "_", condition, ".csv")
    if (!file.exists(csv_path)) {
      message("  Missing DGE for ", cluster, " x ", condition, " — skipping cell type")
      skip_cluster <- TRUE
      break
    }
    df <- read_csv(csv_path, show_col_types = FALSE)
    condition_genes[[condition]] <- df %>%
      filter(p_val_adj < 0.05) %>%
      pull(1)
  }
  if (skip_cluster) next

  # Union of significant genes across all conditions
  unique_genes <- unique(unlist(condition_genes))
  if (length(unique_genes) == 0) {
    message("  No significant DEGs for ", cluster, " across any condition — skipping")
    next
  }

  # Presence/absence matrix (rows = genes, cols = conditions)
  corr_df <- data.frame(Genes = unique_genes, stringsAsFactors = FALSE)
  for (condition in conditions) {
    corr_df[[condition]] <- as.integer(corr_df$Genes %in% condition_genes[[condition]])
  }

  # Pearson correlation between condition columns (phi coefficient on binary)
  corr_mat <- cor(corr_df[, -1], method = "pearson", use = "pairwise.complete.obs")
  melted_corr <- melt(corr_mat)

  # Per-celltype color; fall back to grey if missing from palette
  fill_hi <- major_cell_type_colors[[cluster]]
  if (is.null(fill_hi)) fill_hi <- "grey40"

  p <- ggplot(melted_corr, aes(x = Var1, y = Var2)) +
    geom_tile(aes(fill = value), color = "white") +
    scale_fill_gradientn(colors = c("white", fill_hi), limits = c(0, 1)) +
    geom_text(aes(label = sprintf("%.2f", value)), vjust = 1, size = 3) +
    scale_x_discrete(labels = clean_label) +
    scale_y_discrete(labels = clean_label) +
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

message("11_main_correlations complete.")
