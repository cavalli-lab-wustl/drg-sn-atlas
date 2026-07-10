#' ======================================================================
#' 08_main_umaps_markers_scores.R
#'
#' Compute markers and GO enrichment at the Major cell type level for the
#' annotated atlas object, generate UMAP visualizations (overall +
#' tissue-split + per-cluster highlights), and produce cell-count tables
#' and tissue distribution plots.
#'
#' Inputs:  ../r_objects/DRGSNMI_main_multilevel.RDS
#'          ../r_objects/DRGSNMI_glia.RDS
#' Outputs: ../markers/combined_named_findmarkers/FindMarkers_*.csv + summary.xlsx
#'          ../markers/combined_named_findmarkers/GO/GOBP_*_up.csv
#'          ../plots/umaps/combined_named/*.pdf
#'          ../plots/umaps/combined_named_highlights/*.png
#'          ../plots/tissue_distribution/*.pdf
#'          ../metadata/combined_cellcount.csv
#'          ../metadata/glia_cellcount.csv
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


# ===== LOAD OBJECTS =====
merged_named <- readRDS("../r_objects/DRGSNMI_main_multilevel.RDS")
glia <- readRDS("../r_objects/DRGSNMI_glia.RDS")

cat("merged_named cells:", ncol(merged_named), "\n")
cat("glia cells:", ncol(glia), "\n")


# ===== SET MAJOR AS ACTIVE IDENT, ORDER BY CLASS =====
# class_ordered_majors and major_cell_type_colors come from utils/variables.R.
Idents(merged_named) <- merged_named$Major
major_levels <- intersect(class_ordered_majors, unique(as.character(Idents(merged_named))))
Idents(merged_named) <- factor(Idents(merged_named), levels = major_levels)

cat("Active Major levels (", length(major_levels), "):\n", sep = "")
print(major_levels)


# ===== FIND MARKERS =====
# min.pct = 0.1: gene detected in >=10% of cells in either group
# logfc.threshold = 0.5: require >=0.5 log2FC (stricter than Seurat default 0.25)
# skip_existing = TRUE: reuse CSVs from prior runs (summary.xlsx still regenerates).
# Delete the directory to force a fresh run.
output_directory <- "../markers/combined_named_findmarkers"
dir.create(output_directory, recursive = TRUE)
find_and_save_all_markers(merged_named, output_directory, min.pct = 0.1, logfc.threshold = 0.5, summary = TRUE, skip_existing = TRUE)


###########################################
###########################################
###                                     ###
###               UMAP PLOTS            ###
###                                     ###
###########################################
###########################################

dir.create("../plots/umaps/combined_named", recursive = TRUE)

p <- DimPlot(
  object = merged_named,
  reduction = "umap",
  cols = major_cell_type_colors,
  raster = FALSE
) +
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  ) +
  NoLegend()

psplit <- DimPlot(
  object = merged_named,
  reduction = "umap",
  cols = major_cell_type_colors,
  raster = FALSE,
  split.by = "tissue"
) +
  theme(
    plot.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none"
  ) +
  NoLegend()

p_labeled <- LabelClusters(
  plot = p,
  id = "ident",
  repel = TRUE,
  size = 2.7,
  fontface = "bold",
  colour = "black",
  box.padding = unit(0.35, "lines"),
  segment.color = "black",
  max.overlaps = Inf
)

psplit_labeled <- LabelClusters(
  plot = psplit,
  id = "ident",
  repel = TRUE,
  size = 2.7,
  fontface = "bold",
  colour = "black",
  box.padding = unit(0.35, "lines"),
  segment.color = "black",
  max.overlaps = Inf
)

print(p_labeled)
print(psplit_labeled)

ggsave("../plots/umaps/combined_named/combined_umap.pdf",
       p_labeled, width = 8, height = 8, device = cairo_pdf)
ggsave("../plots/umaps/combined_named/combined_umap_tissue_split.pdf",
       psplit_labeled, width = 12, height = 6, device = cairo_pdf)


# ===== CELL COUNT TABLES =====
cellcount_combined <- table(merged_named@meta.data$study, Idents(merged_named))
cellcount_glia     <- table(glia@meta.data$study,         Idents(glia))

write.csv(cellcount_combined, "../metadata/combined_cellcount.csv")
write.csv(cellcount_glia,     "../metadata/glia_cellcount.csv")


# ===== PER-CLUSTER HIGHLIGHT UMAPs =====
highlightAndSaveAllClusters(
  merged_named,
  custom_cluster_colors = major_cell_type_colors,
  save_directory = "../plots/umaps/combined_named_highlights/"
)


###########################################
###########################################
###                                     ###
###       TISSUE DISTRIBUTION PLOTS     ###
###                                     ###
###########################################
###########################################

dir.create("../plots/tissue_distribution", recursive = TRUE)

cell_types <- Idents(merged_named)
tissue     <- merged_named@meta.data$tissue

# Shared proportion table: per-CellType tissue share as % of that cell type's total
proportion_data_base <- data.frame(
  CellType = cell_types,
  Tissue   = tissue
) %>%
  group_by(CellType) %>%
  mutate(Total_CellType = n()) %>%
  group_by(CellType, Tissue) %>%
  summarise(
    Count      = n(),
    Total      = first(Total_CellType),
    Proportion = Count / Total * 100,
    .groups    = "drop"
  )

tissue_fill <- scale_fill_manual(values = tissue_colors)  # from utils/variables.R


# --- Plot 1: horizontal stacked bars ---
pd1 <- proportion_data_base %>%
  mutate(CellType = factor(CellType, levels = major_levels))

p1 <- ggplot(pd1, aes(x = CellType, y = Proportion, fill = Tissue)) +
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  tissue_fill +
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    minor_breaks = seq(0, 100, 10),
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold", color = "black"),
    legend.text = element_text(color = "black"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "black")
  ) +
  labs(
    title = "Relative Tissue Distribution by Cell Type",
    x = "Cell Type",
    y = "Relative Proportion (%)",
    fill = "Tissue"
  )

ggsave("../plots/tissue_distribution/tissue_distribution_horizontal.pdf",
       p1, width = 8, height = 4, device = cairo_pdf)


# --- Plot 2: vertical stacked bars (reversed order for top-to-bottom reading) ---
pd2 <- proportion_data_base %>%
  mutate(CellType = factor(CellType, levels = rev(major_levels)))

p2 <- ggplot(pd2, aes(x = Proportion, y = CellType, fill = Tissue)) +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  tissue_fill +
  scale_x_continuous(
    breaks = seq(0, 100, 20),
    minor_breaks = seq(0, 100, 10),
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold", color = "black"),
    legend.text = element_text(color = "black"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "black")
  ) +
  labs(
    title = "Relative Tissue Distribution by Cell Type",
    y = "Cell Type",
    x = "Relative Proportion (%)",
    fill = "Tissue"
  )

ggsave("../plots/tissue_distribution/tissue_distribution_vertical.pdf",
       p2, width = 2, height = 5, device = cairo_pdf)


# --- Plot 3: dodged bars with percentage labels ---
pd3 <- proportion_data_base %>%
  mutate(
    CellType = factor(CellType, levels = major_levels),
    CellType_Tissue = factor(
      paste(CellType, Tissue, sep = "_"),
      levels = unlist(lapply(major_levels, function(x) paste(x, c("DRG", "SN"), sep = "_")))
    )
  )

p3 <- ggplot(pd3, aes(x = CellType_Tissue, y = Proportion, fill = Tissue)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  tissue_fill +
  geom_text(
    aes(label = sprintf("%.1f%%", Proportion)),
    position = position_dodge(width = 0.8),
    vjust = -0.5, size = 2.5
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    minor_breaks = seq(0, 100, 10),
    expand = c(0, 0.1)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold", color = "black"),
    legend.text = element_text(color = "black"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "black")
  ) +
  labs(x = NULL, y = "Relative Proportion (%)", fill = "Tissue")

ggsave("../plots/tissue_distribution/tissue_distribution_dodge.pdf",
       p3, width = 12, height = 3, device = cairo_pdf)

# Print tissue distribution summary for log
cat("\nTissue distribution summary:\n")
pd1 %>%
  arrange(CellType, Tissue) %>%
  dplyr::select(CellType, Tissue, Count, Proportion) %>%
  print(n = Inf)


###########################################
###########################################
###                                     ###
###     GO ENRICHMENT ON UP-MARKERS     ###
###                                     ###
###########################################
###########################################

# Strict threshold (log2FC > 1) because inputs are already cluster-marker gene
# lists — GO on these highlights the strongest biology per cell type.
dir.create("../markers/combined_named_findmarkers/GO", recursive = TRUE)

unique_idents <- levels(Idents(merged_named))
for (ident in unique_idents) {
  marker_file <- paste0("../markers/combined_named_findmarkers/FindMarkers_", ident, ".csv")
  if (!file.exists(marker_file)) {
    cat("Skipping", ident, "(no marker file found)\n")
    next
  }

  marker_data <- read.csv(marker_file, stringsAsFactors = FALSE, row.names = 1)

  upregulated_genes <- row.names(marker_data[
    marker_data$p_val_adj < 0.05 & marker_data$avg_log2FC > 1, , drop = FALSE
  ])

  if (length(upregulated_genes) == 0) {
    cat("Skipping", ident, "(no up-regulated genes passing thresholds)\n")
    next
  }

  cat("Running enrichGO for", ident, "(", length(upregulated_genes), "genes)\n")
  go_enrich_up <- enrichGO(
    gene    = upregulated_genes,
    OrgDb   = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont     = "BP"
  )

  write.csv(
    go_enrich_up,
    file = paste0("../markers/combined_named_findmarkers/GO/GOBP_", ident, "_up.csv")
  )
}

message("08_main_umaps_markers_scores complete.")
