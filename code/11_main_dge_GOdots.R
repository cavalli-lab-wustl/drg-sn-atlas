#' ======================================================================
#' 11_main_dge_GOdots.R
#'
#' Per-celltype GO:BP dotplots across 4 injury conditions. For each Major
#' cell type, reads the GOBP_<celltype>_<condition>_up.csv outputs from
#' 10_main_dge_GO, applies per-celltype (Enrichment, -log10 pvalue)
#' thresholds, hierarchically clusters the passing terms, and generates
#' a dotplot showing enrichment across conditions. Exclude/custom GO
#' term lists (text files, one term per line) are honored per cell type.
#'
#' Thresholds and plot dimensions start as conservative defaults for every
#' cell type; tune per-cell-type via `custom_thresholds` / `custom_plot_dims`
#' after reviewing the first-pass figures.
#'
#' A final block produces a Macrophages-only "condition-unique GO terms"
#' dotplot — terms significant in exactly one condition.
#'
#' Inputs:  ../dge_GO/combined_named/GOBP_<celltype>_<condition>_up.csv
#'          ../dge_GO/combined_named/customization_lists/combined_<celltype>_{exclude,custom}.txt
#' Outputs: ../dge_GO/combined_named/customization_lists/<celltype>_injuries_goterms.csv
#'          ../plots/combined_godots/combined_<celltype>_injuries_godots.pdf
#'          ../plots/combined_godots/combined_MAC_injuries_godots_condition_unique.pdf
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


# ===== CONFIG =====
conditions <- c("cavalli_aih", "cavalli_drc", "cavalli_sci", "cavalli_snc")

# Default thresholds for every cell type: (Enrichment > x, -log10(pvalue) > y)
# Default plot dims (width, height). Override per-celltype below after
# reviewing first-run plots.
default_thresholds <- c(0.5, 5)
default_plot_dims  <- c(10, 14)

# Iterate over canonical Major cell-type labels (from utils/variables.R)
idents <- class_ordered_majors

# Starting point: every cell type uses defaults. Override individual entries
# once you've reviewed the first-run plots.
custom_thresholds <- setNames(replicate(length(idents), default_thresholds, simplify = FALSE), idents)
custom_plot_dims  <- setNames(replicate(length(idents), default_plot_dims,  simplify = FALSE), idents)

# Example overrides (uncomment / edit per your figure tuning):
# custom_thresholds[["SGC"]]   <- c(0.75, 7)
# custom_plot_dims[["NEU_A"]]  <- c(10, 20)


# ===== OUTPUT DIRECTORIES =====
customization_dir <- "../dge_GO/combined_named/customization_lists"
plots_dir         <- "../plots/combined_godots"
dir.create(customization_dir, recursive = TRUE)
dir.create(plots_dir,         recursive = TRUE)


# ===== ENSURE EXCLUDE / CUSTOM STUB FILES EXIST (one-time) =====
# Blank files are read as empty vectors and treated as no filtering.
for (ident in idents) {
  for (suffix in c("exclude", "custom")) {
    path <- file.path(customization_dir, paste0("combined_", ident, "_", suffix, ".txt"))
    if (!file.exists(path)) file.create(path)
  }
}


# ===== HELPER: GO term name → IDs (with unmatched reporting) =====
get_GO_IDs_from_names <- function(term_names, data) {
  term_names <- term_names[nzchar(trimws(term_names))]  # drop blank lines
  if (length(term_names) == 0) return(character())
  matched <- data$Description %in% term_names
  unmatched <- setdiff(term_names, data$Description)
  if (length(unmatched) > 0) {
    cat("  Unmatched exclude/custom terms: ", paste(unmatched, collapse = ", "), "\n")
  }
  unique(data$ID[matched])
}


###########################################
###########################################
###                                     ###
###       PER-CELLTYPE GO DOTPLOTS      ###
###                                     ###
###########################################
###########################################

for (ident in idents) {

  cat("Processing:", ident, "\n")

  # Read the 4 condition GO:BP-up results; skip cell type if any are missing
  all_data <- data.frame()
  skip_ident <- FALSE
  for (condition in conditions) {
    file_name <- paste0("../dge_GO/combined_named/GOBP_", ident, "_", condition, "_up.csv")
    if (!file.exists(file_name)) {
      cat("  Missing GO file for", ident, "x", condition, "— skipping cell type\n")
      skip_ident <- TRUE
      break
    }
    condition_data <- read_csv(file_name, show_col_types = FALSE)
    if (nrow(condition_data) == 0) next
    condition_data$Condition <- condition
    all_data <- rbind(all_data, condition_data)
  }
  if (skip_ident || nrow(all_data) == 0) next

  # Derive enrichment and log(p)
  all_data <- all_data %>%
    mutate(
      GeneNumerator   = as.numeric(str_extract(GeneRatio, "^[^/]+")),
      GeneDenominator = as.numeric(str_extract(GeneRatio, "[^/]+$")),
      BgNumerator     = as.numeric(str_extract(BgRatio,   "^[^/]+")),
      BgDenominator   = as.numeric(str_extract(BgRatio,   "[^/]+$")),
      Enrichment      = log10((GeneNumerator / GeneDenominator) / (BgNumerator / BgDenominator)),
      Log_pvalue      = -log10(pvalue)
    )

  # Threshold-filter to significant terms
  thresholds <- custom_thresholds[[ident]]
  significant_terms <- all_data %>%
    filter(Enrichment > thresholds[1], Log_pvalue > thresholds[2]) %>%
    dplyr::select(ID) %>%
    distinct(ID) %>%
    pull(ID)

  significant_data <- all_data %>% filter(ID %in% significant_terms)

  # Apply user exclude / custom lists
  excluded_terms <- readLines(file.path(customization_dir, paste0("combined_", ident, "_exclude.txt")))
  custom_terms   <- readLines(file.path(customization_dir, paste0("combined_", ident, "_custom.txt")))
  excluded_GO_IDs <- get_GO_IDs_from_names(excluded_terms, all_data)
  custom_GO_IDs   <- get_GO_IDs_from_names(custom_terms,   all_data)

  significant_data <- significant_data %>% filter(!(ID %in% excluded_GO_IDs))
  if (length(custom_GO_IDs) > 0) {
    custom_go_data <- all_data %>% filter(ID %in% custom_GO_IDs)
    if (nrow(custom_go_data) > 0) {
      significant_data <- rbind(significant_data, custom_go_data)
    }
  }

  if (nrow(significant_data) == 0) {
    cat("  No terms passing thresholds for", ident, "— skipping plot\n")
    next
  }

  # Hierarchical clustering of GO terms by binary presence across conditions
  descriptions <- unique(significant_data$Description)
  conditions_seen <- unique(significant_data$Condition)
  binary_matrix <- as.data.frame(matrix(
    0, nrow = length(descriptions), ncol = length(conditions_seen),
    dimnames = list(descriptions, conditions_seen)
  ))
  for (i in seq_len(nrow(significant_data))) {
    binary_matrix[as.character(significant_data$Description[i]),
                  as.character(significant_data$Condition[i])] <- 1
  }

  hc <- hclust(dist(binary_matrix, method = "binary"), method = "ward.D2")
  ordered_descriptions <- rownames(binary_matrix)[hc$order]

  # Save the clustered GO term order for reuse / customization
  write.csv(
    ordered_descriptions,
    file = file.path(customization_dir, paste0(ident, "_injuries_goterms.csv")),
    row.names = FALSE
  )

  # Reorder for plotting
  significant_data <- significant_data %>%
    mutate(Description = factor(Description, levels = ordered_descriptions)) %>%
    arrange(Condition, Description)

  # Dotplot
  p <- ggplot(significant_data, aes(x = as.factor(Condition), y = Description)) +
    geom_point(aes(size = Enrichment, color = Log_pvalue)) +
    scale_size_continuous(range = c(1, 6)) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(
      title = paste("GO:BP Enrichment —", ident),
      x = "Condition", y = "GO Term",
      size = "log10(Enrichment)", color = "-log10(pvalue)"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  dims <- custom_plot_dims[[ident]]
  ggsave(
    filename = file.path(plots_dir, paste0("combined_", ident, "_injuries_godots.pdf")),
    plot = p, width = dims[1], height = dims[2], dpi = 300
  )
}


###########################################
###########################################
###                                     ###
###    CONDITION-UNIQUE GO (MAC ONLY)   ###
###                                     ###
###########################################
###########################################

# Dotplot restricted to GO terms significant in exactly ONE of the 4 injury
# conditions for Macrophages — highlights condition-specific biology.
ident <- "MAC"

all_data <- data.frame()
for (condition in conditions) {
  file_name <- paste0("../dge_GO/combined_named/GOBP_", ident, "_", condition, "_up.csv")
  if (!file.exists(file_name)) {
    cat("Missing MAC GO file for", condition, "— skipping condition-unique analysis\n")
    all_data <- NULL
    break
  }
  condition_data <- read_csv(file_name, show_col_types = FALSE)
  if (nrow(condition_data) == 0) next
  condition_data$Condition <- condition
  all_data <- rbind(all_data, condition_data)
}

if (!is.null(all_data) && nrow(all_data) > 0) {
  all_data <- all_data %>%
    mutate(
      GeneNumerator   = as.numeric(str_extract(GeneRatio, "^[^/]+")),
      GeneDenominator = as.numeric(str_extract(GeneRatio, "[^/]+$")),
      BgNumerator     = as.numeric(str_extract(BgRatio,   "^[^/]+")),
      BgDenominator   = as.numeric(str_extract(BgRatio,   "[^/]+$")),
      Enrichment      = log10((GeneNumerator / GeneDenominator) / (BgNumerator / BgDenominator)),
      Log_pvalue      = -log10(pvalue)
    )

  # Looser threshold for this specialized plot
  significant_terms <- all_data %>%
    filter(Enrichment > 0.5, Log_pvalue > 3) %>%
    distinct(ID) %>%
    pull(ID)
  significant_data <- all_data %>% filter(ID %in% significant_terms)

  # Keep only terms appearing in exactly ONE condition
  term_counts <- significant_data %>%
    group_by(ID, Description) %>%
    summarise(Count = n_distinct(Condition), .groups = "drop") %>%
    filter(Count == 1)

  unique_significant_data <- significant_data %>%
    inner_join(term_counts, by = c("ID", "Description")) %>%
    arrange(Condition, Description) %>%
    mutate(Description = factor(Description, levels = unique(Description)))

  if (nrow(unique_significant_data) > 0) {
    p <- ggplot(unique_significant_data, aes(x = as.factor(Condition), y = Description)) +
      geom_point(aes(size = Enrichment, color = Log_pvalue)) +
      scale_size_continuous(range = c(1, 6)) +
      scale_color_gradient(low = "blue", high = "red") +
      labs(
        title = "MAC — Condition-Unique GO:BP Terms",
        x = "Condition", y = "GO Term",
        size = "log10(Enrichment)", color = "-log10(pvalue)"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(
      filename = file.path(plots_dir, "combined_MAC_injuries_godots_condition_unique.pdf"),
      plot = p, width = 10, height = 14, dpi = 300
    )
  } else {
    cat("MAC: no condition-unique GO terms passing thresholds.\n")
  }
}

message("11_main_dge_GOdots complete.")
