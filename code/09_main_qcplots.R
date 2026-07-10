#' ======================================================================
#' 09_main_qcplots.R
#'
#' Generate QC plots for the integrated atlas:
#'   - Absolute / relative cell type abundance bargraphs (by study_simple)
#'   - Violin plots of nFeature_RNA, nCount_RNA, percent.mt, percent.rpl
#'     split by study_simple across Major cell types
#'
#' Inputs:  ../r_objects/DRGSNMI_main_multilevel.RDS
#' Outputs: ../plots/qc/combined_abundances_bargraph_{absolute,relative}.pdf
#'          ../plots/qc/combined_{percentrpl,percentmt,nfeature,ncount}_vlnplot.pdf
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


# ===== LOAD OBJECT =====
combined_object <- readRDS("../r_objects/DRGSNMI_main_multilevel.RDS")
cat("combined_object cells:", ncol(combined_object), "\n")


# ===== SET MAJOR AS ACTIVE IDENT, ORDER BY CLASS =====
# class_ordered_majors and major_cell_type_colors come from utils/variables.R.
Idents(combined_object) <- combined_object$Major
major_levels <- intersect(class_ordered_majors, unique(as.character(Idents(combined_object))))
Idents(combined_object) <- factor(Idents(combined_object), levels = major_levels)


# ===== OUTPUT DIRECTORY =====
qc_dir <- "../plots/qc"
dir.create(qc_dir, recursive = TRUE)


###########################################
###########################################
###                                     ###
###       CELL TYPE ABUNDANCE BARS      ###
###                                     ###
###########################################
###########################################

# Build cell-type x study_simple count table
abundance_table <- table(Idents(combined_object), combined_object@meta.data$study_simple)
abundance_df <- as.data.frame.table(abundance_table)
names(abundance_df) <- c("Cell_Type", "Study", "Count")

plot_cell_type_abundance <- function(data) {
  # Order cell types by total count (descending) for visual readability
  data$Cell_Type <- with(data, reorder(Cell_Type, -Count))

  # Absolute abundance
  p1 <- ggplot(data, aes(x = Cell_Type, y = Count, fill = Study)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9),
             color = "black", linewidth = 0.5, width = 0.8) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Cell Type", y = "Count", fill = "Study",
         title = "Absolute Cell Type Abundances") +
    theme_classic() +
    theme(
      text = element_text(family = "sans", size = 12),
      title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
      legend.position = "inside",
      legend.position.inside = c(0.85, 0.95),
      legend.text = element_text(size = 10)
    ) +
    guides(fill = guide_legend(title = "Study"))

  # Relative abundance (within-cell-type stacked)
  data_relative <- data %>%
    group_by(Cell_Type) %>%
    mutate(Total = sum(Count)) %>%
    mutate(Relative_Count = Count / Total) %>%
    ungroup()

  p2 <- ggplot(data_relative, aes(x = Cell_Type, y = Relative_Count, fill = Study)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.5) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Cell Type", y = "Relative Abundance", fill = "Study",
         title = "Relative Cell Type Abundances") +
    theme_classic() +
    theme(
      text = element_text(family = "sans", size = 12),
      title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
      legend.position = "inside",
      legend.position.inside = c(0.85, 0.95),
      legend.text = element_text(size = 10)
    ) +
    guides(fill = guide_legend(title = "Study"))

  list(Absolute = p1, Relative = p2)
}

plots <- plot_cell_type_abundance(abundance_df)

ggsave(file.path(qc_dir, "combined_abundances_bargraph_absolute.pdf"),
       plot = plots$Absolute, width = 10, height = 8, dpi = 300)
ggsave(file.path(qc_dir, "combined_abundances_bargraph_relative.pdf"),
       plot = plots$Relative, width = 10, height = 8, dpi = 300)


###########################################
###########################################
###                                     ###
###           QC VIOLIN PLOTS           ###
###                                     ###
###########################################
###########################################

# Compute percent.rpl (ribosomal protein transcripts) — percent.mt already set upstream
combined_object[["percent.rpl"]] <- PercentageFeatureSet(object = combined_object, pattern = "^Rpl")

# Helper to produce a consistently-themed VlnPlot split by study
make_qc_vlnplot <- function(seurat_obj, feature, y.max = NULL) {
  VlnPlot(seurat_obj, features = feature, split.by = "study_simple",
          cols = major_cell_type_colors, pt.size = 0, y.max = y.max) +
    theme(
      axis.title.x = element_blank(),
      plot.margin = margin(5.5, 22, 5.5, 5.5, "points")
    )
}

rpl_vln_plot      <- make_qc_vlnplot(combined_object, "percent.rpl")
mt_vln_plot       <- make_qc_vlnplot(combined_object, "percent.mt")
nfeature_vln_plot <- make_qc_vlnplot(combined_object, "nFeature_RNA")
ncount_vln_plot   <- make_qc_vlnplot(combined_object, "nCount_RNA", y.max = 20000)

ggsave(file.path(qc_dir, "combined_percentrpl_vlnplot.pdf"),
       plot = rpl_vln_plot, width = 10, height = 3, dpi = 300)
ggsave(file.path(qc_dir, "combined_percentmt_vlnplot.pdf"),
       plot = mt_vln_plot, width = 10, height = 3, dpi = 300)
ggsave(file.path(qc_dir, "combined_nfeature_vlnplot.pdf"),
       plot = nfeature_vln_plot, width = 10, height = 3, dpi = 300)
ggsave(file.path(qc_dir, "combined_ncount_vlnplot.pdf"),
       plot = ncount_vln_plot, width = 10, height = 3, dpi = 300)

message("09_main_qcplots complete.")
