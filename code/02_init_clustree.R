#' ======================================================================
#' 02_init_clustree.R
#'
#' Test clustering resolutions (0.3-1.2) with clustree; selected res=0.9.
#'
#' Inputs:  ../r_objects/combined_init_object.RDS
#' Outputs: ../plots/clustree/combined_init/*.png
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


###########################################
###########################################
###                                     ###
###    CLUSTER RESOLUTION + CLUSTREE    ###
###                                     ###
###########################################
###########################################

combined_object <- readRDS("../r_objects/combined_init_object.RDS")
dir.create("../plots/clustree/combined_init", recursive = TRUE, showWarnings = FALSE)

# Do not set RNA assay or normalize - FindClusters requires SCT assay

# Array of resolutions
resolutions <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2)

for (res in resolutions) {
  # Find clusters and create plots
  combined_object <- FindClusters(combined_object, resolution = res)
  p <- DimPlot(combined_object, label = T, raster = F)
  ggplot2::ggsave(paste0("../plots/clustree/combined_init/clustree_combined_init_", res, "res.png"), plot = p, device = "png", width = 10, height = 10, dpi = 300)
}

# Run clustree to examine clustering layout
p <- clustree(combined_object, prefix = "SCT_snn_res.") + theme(text = element_text(size = 20))

# Saving the plot with specified parameters
ggsave(filename = "../plots/clustree/combined_init/clustree_combined_init_summary_plot.png", plot = p, device = "png", width = 18, height = 9, dpi = 300)

# Resolution 0.9 selected
message("02_init_clustree complete.")
