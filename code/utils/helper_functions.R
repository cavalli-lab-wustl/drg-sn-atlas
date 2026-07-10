#' ======================================================================
#' helper_functions.R
#'
#' Shared utility functions: markers export, plotting, RNA velocity helpers.
#'
#' Author: Michael B. Thomsen, Ph.D.
#' Created: 2024-07-01
#' ======================================================================


#' Set working directory to the calling script's location.
#' Works with Rscript (command line) and RStudio.
#' Note: Each script also has an inline bootstrap version of this logic
#' because setwd must happen before helper_functions.R can be sourced.
set_script_dir <- function() {
  # Rscript sets --file= in commandArgs
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) {
    setwd(dirname(normalizePath(sub("^--file=", "", f))))
    return(invisible(NULL))
  }
  # source()'d scripts have ofile in the call frame
  for (i in seq_len(sys.nframe())) {
    if (!is.null(sys.frame(i)$ofile)) {
      setwd(dirname(normalizePath(sys.frame(i)$ofile)))
      return(invisible(NULL))
    }
  }
  # RStudio interactive
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    setwd(dirname(rstudioapi::getSourceEditorContext()$path))
    return(invisible(NULL))
  }
  message("set_script_dir(): Could not determine script location; using current directory.")
}


# reverse_idents(seurat_obj)
# Reverse the order of the Idents in a Seurat object
reverse_idents <- function(seurat_obj) {
  seurat_obj <- seurat_obj %>% SetIdent(value = factor(Idents(seurat_obj), levels = rev(levels(Idents(seurat_obj)))))
  return(seurat_obj)
}

find_and_save_all_markers <- function(seurat_obj, output_directory, min.pct = 0.05, min.diff.pct = -Inf, logfc.threshold = -Inf, summary = FALSE, skip_existing = FALSE) {
  if(summary){
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      install.packages("openxlsx")
    }
    library(openxlsx)
  }

  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }

  clusters <- levels(seurat_obj)
  metadata_list <- data.frame(Cluster = character(), FindMarkers_Command = character(), stringsAsFactors = FALSE)
  skipped_clusters <- character()
  summary_list <- list()

  for (cluster in clusters) {
    cell_count <- sum(seurat_obj@meta.data$seurat_clusters == cluster)

    csv_filename <- paste0("FindMarkers_", cluster, ".csv")
    csv_filepath <- file.path(output_directory, csv_filename)

    find_markers_command <- paste0("FindMarkers(seurat_obj, ident.1 = '", cluster,
                                   "', min.pct = ", min.pct,
                                   ", min.diff.pct = ", min.diff.pct,
                                   ", logfc.threshold = ", logfc.threshold,")")

    if (skip_existing && file.exists(csv_filepath)) {
      cat("Skipping FindMarkers for cluster", cluster, "(existing CSV found at", csv_filepath, ")\n")
      if (summary) {
        markers <- read.csv(csv_filepath, row.names = 1, stringsAsFactors = FALSE)
      } else {
        markers <- NULL
      }
      find_markers_command <- paste0(find_markers_command, "  [skipped: existing CSV reused]")
    } else {
      cat("Running FindMarkers for cluster", cluster, "\n")
      markers <- eval(parse(text = find_markers_command))
      write.csv(markers, file = csv_filepath, row.names = TRUE)
    }

    metadata_list <- rbind(metadata_list, data.frame(Cluster = cluster, FindMarkers_Command = find_markers_command, stringsAsFactors = FALSE))

    if(summary){
      filtered_markers <- subset(markers, p_val_adj < 0.05)

      # Detect log fold change column name (varies by Seurat version)
      if("avg_log2FC" %in% colnames(filtered_markers)){
        logfc_col <- "avg_log2FC"
      } else if("avg_logFC" %in% colnames(filtered_markers)){
        logfc_col <- "avg_logFC"
      } else if("log2FoldChange" %in% colnames(filtered_markers)){
        logfc_col <- "log2FoldChange"
      } else {
        stop("No recognized log fold change column found in markers.")
      }

      sorted_markers <- filtered_markers[order(-filtered_markers[[logfc_col]]), ]
      top_markers <- head(sorted_markers, 100)
      summary_list[[cluster]] <- top_markers
    }
  }

  metadata <- data.frame(
    Seurat_Object_Name = deparse(substitute(seurat_obj)),
    Size = format(object.size(seurat_obj), units = "auto"),
    Cluster_Count = length(clusters),
    Cluster_Names = paste(clusters, collapse = ", "),
    Skipped_Clusters = paste(skipped_clusters, collapse = ", "),
    Date_Generated = Sys.Date(),
    Command_List = paste(metadata_list$FindMarkers_Command, collapse = "; "),
    stringsAsFactors = FALSE
  )

  write.csv(metadata, file = file.path(output_directory, "metadata.csv"), row.names = FALSE)

  if(summary && length(summary_list) > 0){
    summary_filepath <- file.path(output_directory, "summary.xlsx")
    wb <- createWorkbook()
    for(cluster in names(summary_list)){
      addWorksheet(wb, sheetName = cluster)
      writeData(wb, sheet = cluster, summary_list[[cluster]], rowNames = TRUE)
    }
    saveWorkbook(wb, file = summary_filepath, overwrite = TRUE)
  }

  return(metadata)
}

PlotGene <- function(obj, gene, dataset = "NA", output_directory = "../plots/") {
  p <- FeaturePlot(obj, features = gene, order = T, raster = F)
  p <- p + theme(
    plot.title = element_text(size = 20, face = "italic", hjust = 0.5),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 10, hjust = 0),
    axis.title.y = element_text(size = 10, vjust = 0, hjust = 0)
  )
  ggplot2::ggsave(paste0(output_directory, dataset, "_", gene, ".png"), plot = p, device = "png", width = 5, height = 5, dpi = 300)
}

PlotGeneZ1 <- function(obj, gene) {
    p <- FeaturePlot(obj, features = gene, order = T, raster = F)
    p <- p + theme(
        plot.title = element_text(size = 20, face = "italic", hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 10, hjust = 0),
        axis.title.y = element_text(size = 10, vjust = 0, hjust = 0)
    )
    ggplot2::ggsave(paste0("../plots/combined_featureplot_", gene, ".png"), plot = p, device = "png", width = 7, height = 7, dpi = 300)
}

PlotGeneZ2 <- function(obj, gene) {
  p <- FeaturePlot(obj, features = gene, order = T, raster = F)
  p <- p + theme(
    plot.title = element_text(size = 20, face = "italic", hjust = 0.5),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 10, hjust = 0),
    axis.title.y = element_text(size = 10, vjust = 0, hjust = 0)
  )
  ggplot2::ggsave(paste0("../plots/glia_featureplot_", gene, ".png"), plot = p, device = "png", width = 7, height = 7, dpi = 300)
}

highlightAndSaveAllClusters <- function(seurat_obj, save_plot = TRUE, custom_cluster_colors = NULL, save_directory = "../plots/") {
  # Create save_directory if it doesn't exist
  if (save_plot && !dir.exists(save_directory)) {
    dir.create(save_directory, recursive = TRUE)
  }

  # Get unique cluster ids in the order of factor levels
  cluster_ids <- levels(factor(Idents(seurat_obj)))

  # Set colors for clusters
  if (is.null(custom_cluster_colors)) {
    cluster_colors <- rainbow(length(cluster_ids))
  } else {
    # Ensure all clusters have a corresponding color
    if (!all(names(custom_cluster_colors) %in% cluster_ids)) {
      stop("Some clusters do not have a corresponding color.")
    }
    if (length(cluster_ids) != length(custom_cluster_colors)) {
      warning("Length of custom colors does not match number of clusters.")
    }
    # Match colors with cluster IDs based on the levels of Idents
    cluster_colors <- custom_cluster_colors[cluster_ids]
  }

  # For each cluster
  for (i in seq_along(cluster_ids)) {
    # Highlight the selected cluster
    highlight_colors <- setNames(rep("lightgrey", length(cluster_ids)), cluster_ids)
    highlight_colors[cluster_ids[i]] <- cluster_colors[i]

    # Generate the plot
    p <- DimPlot(seurat_obj, group.by = "ident", cells.highlight = WhichCells(seurat_obj, ident = cluster_ids[i]),
                 cols.highlight = cluster_colors[i], cols = "lightgrey", raster = FALSE) +
      ggtitle(paste(cluster_ids[i]))

    p <- p + theme(
      plot.title = element_text(size = 20, face = "italic", hjust = 0.5),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 10, hjust = 0),
      axis.title.y = element_text(size = 10, vjust = 0, hjust = 0)
    ) + NoLegend()

    if (save_plot) {
      filepath_png <- paste0(save_directory, cluster_ids[i], ".png")
      ggplot2::ggsave(filepath_png, plot = p, device = "png", width = 5, height = 5, dpi = 300)
      print(paste("PNG plot saved at", filepath_png))
    }
  }

  print("All cluster plots saved.")
}

###########################################
###########################################
###                                     ###
###           RNA VELOCITY              ###
###                                     ###
###########################################
###########################################

# Function to read and process a single loom file
process_loom <- function(file_path) {
  loom_data <- ReadVelocity(file = file_path)
  spliced <- loom_data$spliced
  unspliced <- loom_data$unspliced

  # Ensure unique gene names
  rownames(spliced) <- make.unique(rownames(spliced))
  rownames(unspliced) <- make.unique(rownames(unspliced))

  return(list(spliced = spliced, unspliced = unspliced))
}

# Function to modify barcode names
modify_barcodes <- function(colnames, sample_prefix, barcode_suffix) {
  colnames <- gsub(paste0(sample_prefix, ":"), "", colnames)
  colnames <- gsub("x$", paste0("-", barcode_suffix), colnames)
  return(colnames)
}

# Function to plot velocity of a single gene on UMAP
plot_gene_velocity <- function(seurat_obj, gene, vel_matrix) {
  # Extract UMAP coordinates
  umap_coords <- Embeddings(seurat_obj, reduction = "umap")

  # Extract velocity for the gene of interest
  gene_velocity <- vel_matrix[gene, ]

  # Create a data frame for plotting
  plot_data <- data.frame(
    UMAP_1 = umap_coords[, 1],
    UMAP_2 = umap_coords[, 2],
    Velocity = gene_velocity
  )

  # Create the plot
  p <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = Velocity)) +
    geom_point(size = 1, alpha = 0.7) +
    scale_color_gradient2(low = "blue", mid = "grey90", high = "red", midpoint = 0) +
    theme_bw() +
    labs(title = paste("Velocity of", gene),
         color = "Velocity") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  return(p)
}

# Function to create a grid of velocity plots
plot_velocity_grid <- function(seurat_obj, genes, vel_matrix, ncol = 3) {
  plots <- lapply(genes, function(gene) {
    plot_gene_velocity(seurat_obj, gene, vel_matrix) +
      theme(legend.position = "none")  # Remove individual legends
  })

  # Arrange plots in a grid
  plot_grid(plotlist = plots, ncol = ncol)
}


plot_total_expression_vs_velocity <- function(seurat_obj, gene, vel_matrix) {
  # Get expression data
  expression <- GetAssayData(seurat_obj, layer = "data")[gene, ]

  # Get velocity data
  velocity <- vel_matrix[gene, ]

  # Create a data frame for plotting
  plot_data <- data.frame(
    Expression = expression,
    Velocity = velocity
  )

  # Create the plot
  p <- ggplot(plot_data, aes(x = Expression, y = Velocity)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    theme_minimal() +
    labs(title = paste("Expression vs Velocity for", gene),
         x = "Expression", y = "Velocity")

  return(p)
}


plot_expression_vs_velocity <- function(seurat_obj, gene, vel_matrix, color_scale = NULL, point_size = 3, point_alpha = 0.5) {
  # Get expression data
  expression <- GetAssayData(seurat_obj, slot = "data")[gene, ]

  # Get velocity data
  velocity <- vel_matrix[gene, ]

  # Get cluster data
  clusters <- seurat_obj$CellType

  # Create a data frame for plotting
  plot_data <- data.frame(
    Expression = expression,
    Velocity = velocity,
    Cluster = clusters
  )

  # Create the plot
  p <- ggplot(plot_data, aes(x = Expression, y = Velocity, color = Cluster)) +
    geom_point(size = point_size, alpha = point_alpha) +
    theme_minimal() +
    labs(title = paste("Expression vs Velocity for", gene),
         x = "Expression", y = "Velocity", color = "Cell Type") +
    theme(legend.position = "right")

  # Add color scale if provided
  if (!is.null(color_scale)) {
    p <- p + scale_color_manual(values = color_scale)
  }

  # Add a density contour
  p <- p + geom_density_2d(color = "black", alpha = 0.5)

  return(p)
}


plot_spliced_unspliced_ratio <- function(seurat_obj, color_scale = NULL) {
  # Calculate the sum of spliced and unspliced counts for each cell
  spliced_sum <- colSums(GetAssayData(seurat_obj, assay = "spliced", slot = "counts"))
  unspliced_sum <- colSums(GetAssayData(seurat_obj, assay = "unspliced", slot = "counts"))

  # Calculate the ratio
  ratio <- spliced_sum / (spliced_sum + unspliced_sum)

  # Create a data frame for plotting
  plot_data <- data.frame(
    Ratio = ratio,
    Cluster = seurat_obj$CellType
  )

  # Create the plot
  p <- ggplot(plot_data, aes(x = Cluster, y = Ratio, fill = Cluster)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = "Spliced/Unspliced Ratio Across Clusters",
         x = "Cluster", y = "Spliced / (Spliced + Unspliced)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Add custom color scale if provided
  if (!is.null(color_scale)) {
    p <- p + scale_fill_manual(values = color_scale)
  }

  return(p)
}
