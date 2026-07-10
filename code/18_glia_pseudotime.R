#' ======================================================================
#' 18_glia_pseudotime.R
#'
#' Slingshot trajectory + tradeSeq on G_PROG population; gene-pseudotime correlations.
#'
#' Inputs:  ../r_objects/DRGSNMI_glia.RDS
#' Outputs: trajectory_genes_PROG2_fixed.csv, ../plots/pseudotime_*.pdf
#'
#' Author: Michael B. Thomsen, Ph.D.
#' Created: 2024-07-01
#' ======================================================================

# Fixed function to analyze gene expression along pseudotime
analyze_genes <- function(sds, expression_matrix, output_file = NULL) {
  # Get pseudotime values - fix the indexing issue
  pseudotime_matrix <- slingPseudotime(sds)

  # Handle case where there's only one trajectory
  if (is.vector(pseudotime_matrix)) {
    pseudotime_values <- pseudotime_matrix
  } else {
    # Use the first trajectory if multiple exist
    pseudotime_values <- pseudotime_matrix[, 1]
  }

  # Remove cells with NA pseudotime values
  valid_cells <- !is.na(pseudotime_values)
  pseudotime_clean <- pseudotime_values[valid_cells]
  expression_clean <- expression_matrix[, valid_cells]

  message("Analyzing ", nrow(expression_clean), " genes across ",
          length(pseudotime_clean), " cells with valid pseudotime")

  # Fixed correlation calculation function
  calc_correlation <- function(gene_expr) {
    # Remove any remaining NA values
    valid_indices <- !is.na(gene_expr) & !is.na(pseudotime_clean)

    if (sum(valid_indices) < 10) {
      # Not enough valid data points
      return(c(correlation = NA, pvalue = NA))
    }

    gene_clean <- gene_expr[valid_indices]
    pseudo_clean <- pseudotime_clean[valid_indices]

    # Check for constant values (no variation)
    if (var(gene_clean) == 0 || var(pseudo_clean) == 0) {
      return(c(correlation = 0, pvalue = 1))
    }

    # Calculate Spearman correlation (rank-based, robust to outliers)
    tryCatch({
      cor_result <- cor.test(gene_clean, pseudo_clean, method = "spearman", exact = FALSE)
      return(c(
        correlation = as.numeric(cor_result$estimate),
        pvalue = as.numeric(cor_result$p.value)
      ))
    }, error = function(e) {
      message("Error calculating correlation for gene: ", e$message)
      return(c(correlation = NA, pvalue = NA))
    })
  }

  # Calculate statistics for each gene
  message("Calculating correlations...")
  gene_stats <- t(apply(expression_clean, 1, function(x) {
    cor_stats <- calc_correlation(x)
    c(
      mean_expr = mean(x, na.rm = TRUE),
      var_expr = var(x, na.rm = TRUE),
      correlation = cor_stats["correlation"],
      pvalue = cor_stats["pvalue"]
    )
  }))

  # Convert to data frame and clean up
  gene_stats <- as.data.frame(gene_stats)
  gene_stats$gene <- rownames(gene_stats)

  # Calculate FDR only for non-NA p-values
  valid_pvals <- !is.na(gene_stats$pvalue)
  gene_stats$fdr <- NA
  if (sum(valid_pvals) > 0) {
    gene_stats$fdr[valid_pvals] <- p.adjust(gene_stats$pvalue[valid_pvals], method = "fdr")
  }

  # Calculate absolute correlation
  gene_stats$abs_correlation <- abs(gene_stats$correlation)

  # Sort by absolute correlation (NA values will go to end)
  gene_stats <- gene_stats[order(gene_stats$abs_correlation, decreasing = TRUE, na.last = TRUE), ]

  # Print summary statistics
  n_significant <- sum(gene_stats$fdr < 0.05, na.rm = TRUE)
  n_positive_cor <- sum(gene_stats$correlation > 0, na.rm = TRUE)
  n_negative_cor <- sum(gene_stats$correlation < 0, na.rm = TRUE)

  message("Summary:")
  message("- Total genes analyzed: ", nrow(gene_stats))
  message("- Genes with valid correlations: ", sum(!is.na(gene_stats$correlation.correlation)))
  message("- Significantly correlated genes (FDR < 0.05): ", n_significant)
  message("- Positively correlated: ", n_positive_cor)
  message("- Negatively correlated: ", n_negative_cor)

  # Save to file if specified
  if (!is.null(output_file)) {
    write.csv(gene_stats, output_file, row.names = FALSE)
    message("Gene statistics saved to: ", output_file)
  }

  return(gene_stats)
}

# Enhanced plotting function that handles the single trajectory case
plot_pseudotime_expression <- function(sds, expression_matrix, genes,
                                       max_genes = 9, ncol = 3) {
  # Get pseudotime values - handle single trajectory
  pseudotime_matrix <- slingPseudotime(sds)
  if (is.vector(pseudotime_matrix)) {
    pseudotime_values <- pseudotime_matrix
  } else {
    pseudotime_values <- pseudotime_matrix[, 1]
  }

  # Remove NA values
  valid_cells <- !is.na(pseudotime_values)
  pseudotime_clean <- pseudotime_values[valid_cells]
  expression_clean <- expression_matrix[, valid_cells]

  # Limit number of genes to plot
  genes_to_plot <- head(genes, max_genes)

  # Prepare data for plotting
  plot_data <- map_dfr(genes_to_plot, function(gene) {
    if (gene %in% rownames(expression_clean)) {
      data.frame(
        pseudotime = pseudotime_clean,
        expression = as.numeric(expression_clean[gene, ]),
        gene = gene,
        stringsAsFactors = FALSE
      )
    }
  })

  # Create plot
  if (nrow(plot_data) > 0) {
    ggplot(plot_data, aes(x = pseudotime, y = expression)) +
      geom_point(alpha = 0.4, size = 1, color = "steelblue") +
      geom_smooth(method = "loess", se = TRUE, color = "red") +
      facet_wrap(~gene, scales = "free_y", ncol = ncol) +
      theme_minimal() +
      labs(x = "Pseudotime", y = "Expression Level") +
      theme(
        strip.text = element_text(face = "bold", size = 10),
        axis.text = element_text(size = 8)
      )
  } else {
    message("No valid genes found for plotting")
    return(NULL)
  }
}

# Updated main analysis workflow
analyze_trajectories <- function(seurat_obj, prog_cluster = "PROG II",
                                 n_genes = 100, n_trajectories = 2,
                                 output_file = "trajectory_genes.csv") {
  # Prepare and clean data
  traj_data <- prepare_trajectory_data(seurat_obj, prog_cluster)

  # Run enhanced Slingshot
  sds <- run_slingshot(traj_data, prog_cluster, n_trajectories)

  # Create visualizations
  plots <- create_trajectory_plots(sds, traj_data)

  # Analyze genes with fixed correlation calculation
  gene_stats <- analyze_genes(sds, traj_data$expression, output_file)

  # Create expression plots for top correlated genes
  # Filter out genes with NA correlations first
  valid_genes <- gene_stats[!is.na(gene_stats$correlation), ]
  if (nrow(valid_genes) > 0) {
    top_genes <- head(valid_genes$gene, 9)
    gene_plots <- plot_pseudotime_expression(sds, traj_data$expression, top_genes)
  } else {
    gene_plots <- NULL
    message("No genes with valid correlations found for plotting")
  }

  return(list(
    slingshot_object = sds,
    gene_statistics = gene_stats,
    trajectory_plots = plots,
    gene_plots = gene_plots,
    data = traj_data
  ))
}

glia <- readRDS("../r_objects/DRGSNMI_glia.RDS")

results <- analyze_trajectories(glia, output_file = "trajectory_genes_PROG2_fixed.csv")

# Function to create volcano plot for pseudotime correlation analysis
plot_pseudotime_volcano <- function(gene_stats,
                                    genes_to_label = NULL,
                                    fdr_threshold = 0.05,
                                    correlation_threshold = 0.3,
                                    point_size = 1.5,
                                    label_size = 3,
                                    title = "Pseudotime Correlation Analysis") {

  require(ggplot2)
  require(ggrepel)
  require(dplyr)

  # Standardize column names (handle correlation.correlation from apply output)
  if ("correlation.correlation" %in% colnames(gene_stats)) {
    gene_stats$correlation <- gene_stats$correlation.correlation
  }
  if ("pvalue.pvalue" %in% colnames(gene_stats)) {
    gene_stats$pvalue <- gene_stats$pvalue.pvalue
  }

  # Clean the data - remove rows with NA values and invert correlation values
  plot_data <- gene_stats %>%
    filter(!is.na(correlation) & !is.na(fdr) & fdr > 0) %>%
    mutate(
      correlation = -correlation,  # Invert to match biological direction
      neg_log_fdr = -log10(fdr),
      significant = fdr < fdr_threshold,
      high_correlation = abs(correlation) > correlation_threshold,
      category = case_when(
        significant & high_correlation & correlation > 0 ~ "Significant Positive",
        significant & high_correlation & correlation < 0 ~ "Significant Negative",
        significant ~ "Significant",
        TRUE ~ "Not Significant"
      )
    )

  # Create the base plot
  p <- ggplot(plot_data, aes(x = correlation, y = neg_log_fdr)) +
    geom_point(aes(color = category), size = point_size, alpha = 0.7) +
    scale_color_manual(
      values = c(
        "Significant Positive" = "#e74c3c",
        "Significant Negative" = "#3498db",
        "Significant" = "#f39c12",
        "Not Significant" = "#95a5a6"
      ),
      name = "Category"
    ) +
    # Add threshold lines
    geom_hline(yintercept = -log10(fdr_threshold),
               linetype = "dashed", color = "gray50", alpha = 0.7) +
    geom_vline(xintercept = c(-correlation_threshold, correlation_threshold),
               linetype = "dashed", color = "gray50", alpha = 0.7) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    ) +
    labs(
      x = "Correlation with Pseudotime",
      y = "-log10(FDR)",
      title = title,
      subtitle = paste0("FDR threshold: ", fdr_threshold,
                        ", |Correlation| threshold: ", correlation_threshold)
    )

  # Add gene labels if provided
  if (!is.null(genes_to_label)) {
    label_data <- plot_data %>%
      filter(gene %in% genes_to_label)

    if (nrow(label_data) > 0) {
      p <- p +
        geom_point(data = label_data,
                   aes(x = correlation, y = neg_log_fdr),
                   color = "black", size = point_size + 0.5, shape = 21,
                   fill = "yellow", alpha = 0.8) +
        geom_text_repel(
          data = label_data,
          aes(label = gene),
          size = label_size,
          box.padding = 0.5,
          point.padding = 0.3,
          segment.color = "black",
          segment.size = 0.3,
          max.overlaps = Inf,
          force = 2
        )

      message("Labeled ", nrow(label_data), " genes out of ",
              length(genes_to_label), " requested")
    } else {
      message("None of the requested genes were found in the data")
    }
  }

  return(p)
}

# Function to get summary statistics for the volcano plot
summarize_volcano_data <- function(gene_stats,
                                   fdr_threshold = 0.05,
                                   correlation_threshold = 0.3) {

  # Standardize column names (handle correlation.correlation from apply output)
  if ("correlation.correlation" %in% colnames(gene_stats)) {
    gene_stats$correlation <- gene_stats$correlation.correlation
  }
  if ("pvalue.pvalue" %in% colnames(gene_stats)) {
    gene_stats$pvalue <- gene_stats$pvalue.pvalue
  }

  # Invert correlation values to match other figures
  gene_stats$correlation <- -gene_stats$correlation

  clean_data <- gene_stats %>%
    filter(!is.na(correlation) & !is.na(fdr))

  cat("Volcano Plot Summary:\n")
  cat("====================\n")
  cat("Total genes:", nrow(clean_data), "\n")
  cat("Significant genes (FDR <", fdr_threshold, "):",
      sum(clean_data$fdr < fdr_threshold, na.rm = TRUE), "\n")
  cat("High correlation genes (|r| >", correlation_threshold, "):",
      sum(abs(clean_data$correlation) > correlation_threshold, na.rm = TRUE), "\n")
  cat("Significant + High correlation:",
      sum(clean_data$fdr < fdr_threshold &
            abs(clean_data$correlation) > correlation_threshold, na.rm = TRUE), "\n")
  cat("Positive correlations:", sum(clean_data$correlation > 0, na.rm = TRUE), "\n")
  cat("Negative correlations:", sum(clean_data$correlation < 0, na.rm = TRUE), "\n")

  # Return top genes by category
  top_positive <- clean_data %>%
    filter(correlation > 0) %>%
    arrange(desc(correlation)) %>%
    head(5)

  top_negative <- clean_data %>%
    filter(correlation < 0) %>%
    arrange(correlation) %>%
    head(5)

  most_significant <- clean_data %>%
    arrange(fdr) %>%
    head(10)

  cat("\nTop 5 Positively Correlated Genes:\n")
  print(top_positive[, c("gene", "correlation", "fdr")])

  cat("\nTop 5 Negatively Correlated Genes:\n")
  print(top_negative[, c("gene", "correlation", "fdr")])

  cat("\nTop 10 Most Significant Genes:\n")
  print(most_significant[, c("gene", "correlation", "fdr")])

  return(list(
    top_positive = top_positive,
    top_negative = top_negative,
    most_significant = most_significant
  ))
}

# Define genes you want to label
my_genes_of_interest <- c("Mbp", "Csrp2", "Emid1", "Mpz", "Pmp22", "Prx", "Mag",
                          "Ednrb", "Ptn", "Dbi", "Csmd1", "Apoe",
                          "Scn7a", "Ncam1", "L1cam", "Tenm3")

# Create the volcano plot
volcano_plot <- plot_pseudotime_volcano(
 gene_stats = results$gene_statistics,
 genes_to_label = my_genes_of_interest,
 fdr_threshold = 0.05,
 correlation_threshold = 0,  # no correlation cutoff, show all significant
 title = "PROG II Pseudotime Trajectory Analysis"
)

# Display the plot
print(volcano_plot)

# Get summary statistics
volcano_summary <- summarize_volcano_data(results$gene_statistics)

# Save the plot
ggsave("../plots/pseudotime_volcano_plot.pdf", volcano_plot,
      width = 5, height = 5, dpi = 300)

# Enhanced function to overlay multiple genes on one plot with improved styling
plot_overlaid_genes <- function(sds, expression_matrix, genes,
                                smooth = TRUE,
                                scale_expression = TRUE,
                                point_alpha = 0.3,
                                point_size = 1,
                                line_size = 1.2,
                                title = "Gene Expression Along Pseudotime") {

  require(ggplot2)
  require(dplyr)
  require(purrr)
  require(scales)

  # Get pseudotime values and handle single trajectory case
  pseudotime_matrix <- slingPseudotime(sds)
  if (is.vector(pseudotime_matrix)) {
    pseudotime_values <- pseudotime_matrix
  } else {
    pseudotime_values <- pseudotime_matrix[, 1]
  }

  # INVERT pseudotime to match volcano plot correlation direction
  pseudotime_values <- -pseudotime_values

  # Remove cells with NA pseudotime values
  valid_cells <- !is.na(pseudotime_values)
  pseudotime_clean <- pseudotime_values[valid_cells]
  expression_clean <- expression_matrix[, valid_cells]

  # Check which genes are available - PRESERVE ORDER from input list
  available_genes <- genes[genes %in% rownames(expression_clean)]
  missing_genes <- genes[!genes %in% rownames(expression_clean)]

  if (length(missing_genes) > 0) {
    message("Warning: The following genes were not found: ", paste(missing_genes, collapse = ", "))
  }

  if (length(available_genes) == 0) {
    stop("None of the requested genes were found in the expression matrix")
  }

  # Prepare data for plotting - preserve gene order by using factor levels
  plot_data <- map_dfr(available_genes, function(gene) {
    expr <- as.numeric(expression_clean[gene, ])

    if (scale_expression) {
      # Scale expression to 0-1 range for better comparison across genes
      if (max(expr) != min(expr)) {
        scaled_expr <- (expr - min(expr)) / (max(expr) - min(expr))
      } else {
        scaled_expr <- rep(0.5, length(expr))  # Handle constant expression
      }
    } else {
      scaled_expr <- expr
    }

    data.frame(
      pseudotime = pseudotime_clean,
      expression = scaled_expr,
      raw_expression = expr,
      gene = gene,
      stringsAsFactors = FALSE
    )
  })

  # Set gene factor levels to preserve input order
  plot_data$gene <- factor(plot_data$gene, levels = available_genes)

  # Create a nice color palette - assign colors based on INPUT ORDER
  n_genes <- length(available_genes)
  if (n_genes <= 8) {
    # Colorblind-friendly palette for up to 8 genes
    colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
                "#9467bd", "#8c564b", "#e377c2", "#7f7f7f")[1:n_genes]
  } else {
    colors <- viridis_discrete(n_genes)
  }
  names(colors) <- available_genes

  # Create the base plot
  p <- ggplot(plot_data, aes(x = pseudotime, y = expression, color = gene)) +
    scale_color_manual(values = colors) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      legend.key.size = unit(1, "cm")
    ) +
    labs(
      x = "Inverted Pseudotime",
      y = ifelse(scale_expression, "Scaled Expression (0-1)", "Expression Level"),
      color = "Gene",
      title = title
    )

  # Add points and/or smooth lines
  if (smooth) {
    p <- p +
      geom_point(alpha = point_alpha, size = point_size) +
      geom_smooth(method = "loess", se = FALSE, size = line_size, span = 0.75)
  } else {
    p <- p + geom_point(alpha = point_alpha + 0.2, size = point_size + 0.5)
  }

  # Add subtle background
  p <- p + theme(panel.background = element_rect(fill = "white", color = NA))

  return(p)
}

# Enhanced function to create individual faceted plots for each gene
plot_individual_genes <- function(sds, expression_matrix, genes,
                                  smooth = TRUE,
                                  ncol = 3,
                                  point_alpha = 0.4,
                                  point_size = 1,
                                  title = "Individual Gene Expression Along Pseudotime") {

  require(ggplot2)
  require(dplyr)
  require(purrr)

  # Get pseudotime values and handle single trajectory case
  pseudotime_matrix <- slingPseudotime(sds)
  if (is.vector(pseudotime_matrix)) {
    pseudotime_values <- pseudotime_matrix
  } else {
    pseudotime_values <- pseudotime_matrix[, 1]
  }

  # INVERT pseudotime to match volcano plot correlation direction
  pseudotime_values <- -pseudotime_values

  # Remove cells with NA pseudotime values
  valid_cells <- !is.na(pseudotime_values)
  pseudotime_clean <- pseudotime_values[valid_cells]
  expression_clean <- expression_matrix[, valid_cells]

  # Check which genes are available - PRESERVE ORDER from input list
  available_genes <- genes[genes %in% rownames(expression_clean)]
  missing_genes <- genes[!genes %in% rownames(expression_clean)]

  if (length(missing_genes) > 0) {
    message("Warning: The following genes were not found: ", paste(missing_genes, collapse = ", "))
  }

  if (length(available_genes) == 0) {
    stop("None of the requested genes were found in the expression matrix")
  }

  # Prepare data for plotting - preserve gene order
  plot_data <- map_dfr(available_genes, function(gene) {
    data.frame(
      pseudotime = pseudotime_clean,
      expression = as.numeric(expression_clean[gene, ]),
      gene = gene,
      stringsAsFactors = FALSE
    )
  })

  # Set gene factor levels to preserve input order
  plot_data$gene <- factor(plot_data$gene, levels = available_genes)

  # Create the plot
  p <- ggplot(plot_data, aes(x = pseudotime, y = expression)) +
    geom_point(alpha = point_alpha, size = point_size, color = "steelblue") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 11, face = "bold")
    ) +
    labs(
      x = "Inverted Pseudotime",
      y = "Expression Level",
      title = title
    ) +
    facet_wrap(~gene, scales = "free_y", ncol = ncol)

  # Add smooth line if requested
  if (smooth) {
    p <- p + geom_smooth(method = "loess", se = TRUE, color = "red",
                         fill = "red", alpha = 0.2, size = 1, span = 0.75)
  }

  return(p)
}

# For overlaid plot (multiple genes on same plot)
genes_of_interest <- c("Ednrb", "Emid1", "Pou3f1")

overlaid_plot <- plot_overlaid_genes(
  sds = results$slingshot_object,
  expression_matrix = results$data$expression,
  genes = genes_of_interest,
  smooth = TRUE,
  scale_expression = TRUE,
  title = "Gene Expression Patterns - PROG II Trajectory"
)

print(overlaid_plot)

ggsave("../plots/pt_overlaid_genes_plot.pdf", overlaid_plot, width = 4, height = 3.5, dpi = 300)

# For overlaid plot (second gene set)
genes_of_interest <- c("Ptn", "Mbp", "Pou3f1")

overlaid_plot <- plot_overlaid_genes(
  sds = results$slingshot_object,
  expression_matrix = results$data$expression,
  genes = genes_of_interest,
  smooth = TRUE,
  scale_expression = TRUE,
  title = "Gene Expression Patterns - PROG II Trajectory"
)

print(overlaid_plot)

ggsave("../plots/pt_overlaid_genes_plot2.pdf", overlaid_plot, width = 4, height = 3.5, dpi = 300)

# For individual faceted plots
individual_plot <- plot_individual_genes(
  sds = results$slingshot_object,
  expression_matrix = results$data$expression,
  genes = genes_of_interest,
  smooth = TRUE,
  ncol = 2,
  title = "Individual Gene Expression - PROG II Trajectory"
)

print(individual_plot)

PlotGene <- function(obj, gene, dataset = "NA") {
  p <- FeaturePlot(
    obj,
    features  = gene,
    order     = F,
    raster    = FALSE
  ) &
    theme(
      plot.title  = element_text(size = 20, face = "italic", hjust = 0.5),
      axis.text   = element_blank(),
      axis.ticks  = element_blank(),
      axis.title.x = element_text(size = 10, hjust = 0),
      axis.title.y = element_text(size = 10, vjust = 0, hjust = 0)
    )

  print(p)
}

PlotGene(glia, "Ednrb")
