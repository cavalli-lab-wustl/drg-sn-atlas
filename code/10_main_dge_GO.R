#' ======================================================================
#' 10_main_dge_GO.R
#'
#' Injury-vs-control differential gene expression per Major cell type,
#' followed by GO:BP enrichment on all / up / down DGE hits.
#'
#' For each of 4 injury conditions (cavalli_aih/drc/sci/snc), cells matching
#' that condition's 'study' are compared against the pooled control set
#' (all non-injury studies). Per-celltype FindMarkers runs with no
#' thresholds so downstream consumers can re-filter. GO uses p_val_adj < 0.05
#' on each direction. This script is long (~30+ min DGE + ~20-30 min GO);
#' skip_existing toggles let you resume after interruption.
#'
#' Inputs:  ../r_objects/DRGSNMI_main_multilevel.RDS
#' Outputs: ../dge/combined_named/DGE_<celltype>_<condition>.csv
#'          ../dge_GO/combined_named/GOBP_<celltype>_<condition>_{all,up,down}.csv
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


# ===== TOGGLES =====
# Set to FALSE to force regeneration of that phase's outputs.
skip_existing_dge <- TRUE
skip_existing_go  <- TRUE


# ===== LOAD OBJECT =====
seurat_object <- readRDS("../r_objects/DRGSNMI_main_multilevel.RDS")
cat("Loaded DRGSNMI_main_multilevel:", ncol(seurat_object), "cells\n")

# Set Major as active ident in canonical class-grouped order
Idents(seurat_object) <- seurat_object$Major
all_labels <- unique(as.character(Idents(seurat_object)))
missing_from_order <- setdiff(all_labels, class_ordered_majors)
if (length(missing_from_order) > 0) {
  warning("Labels present in object but NOT in class_ordered_majors (will be skipped): ",
          paste(missing_from_order, collapse = ", "))
}
unique_idents <- intersect(class_ordered_majors, all_labels)
Idents(seurat_object) <- factor(Idents(seurat_object), levels = unique_idents)

DefaultAssay(seurat_object) <- "RNA"


# ===== OUTPUT DIRECTORIES =====
dge_dir <- "../dge/combined_named"
go_dir  <- "../dge_GO/combined_named"
dir.create(dge_dir, recursive = TRUE)
dir.create(go_dir,  recursive = TRUE)


# ===== CONDITIONS =====
injury_conditions <- c("cavalli_aih", "cavalli_drc", "cavalli_sci", "cavalli_snc")

control_conditions <- c(
  "berta_bbi_2023",
  "cavalli_contra",
  "cavalli_naive",
  "giger_elife_2022",
  "horste_pnas_2019",
  "kaminker_natcomm_2023",
  "kuruvilla_cellrep_2022",
  "milbrandt_natneuro_2022",
  "miller_fd_eneuro_2020",
  "miller_re_natcomm_2023",
  "suter_elife_2021"
)


# ===== PRE-COMPUTE CONTROL SUBSET (shared across all injury comparisons) =====
cat("\nBuilding pooled control subset...\n")
subset_control <- subset(seurat_object, subset = study %in% control_conditions)
subset_control$condition <- "control"
cat("Pooled control cells:", ncol(subset_control), "\n")


###########################################
###########################################
###                                     ###
###            DGE (PHASE 1)            ###
###                                     ###
###########################################
###########################################

for (condition in injury_conditions) {
  cat("\n===== DGE for", condition, "=====\n")

  subset_injury <- subset(seurat_object, subset = study == condition)
  subset_injury$condition <- "injury"
  cat(condition, "injury cells:", ncol(subset_injury), "\n")

  # Merge injury + shared control; join layers so FindMarkers can use RNA data
  comparison_object <- merge(
    subset_injury,
    y = subset_control,
    add.cell.ids = c("injury", "control"),
    project = "DGE_analysis"
  )
  comparison_object <- JoinLayers(comparison_object)
  DefaultAssay(comparison_object) <- "RNA"

  # Build celltype.condition ident for pairwise FindMarkers
  comparison_object$celltype.condition <- paste(
    Idents(comparison_object), comparison_object$condition, sep = "_"
  )
  Idents(comparison_object) <- "celltype.condition"

  for (celltype in unique_idents) {
    csv_path <- file.path(dge_dir, paste0("DGE_", celltype, "_", condition, ".csv"))

    if (skip_existing_dge && file.exists(csv_path)) {
      cat("  Skipping DGE:", celltype, "(existing CSV)\n")
      next
    }

    tryCatch({
      cat("  Running DGE:", celltype, "injury vs control\n")
      dge <- FindMarkers(
        comparison_object,
        ident.1         = paste(celltype, "injury",  sep = "_"),
        ident.2         = paste(celltype, "control", sep = "_"),
        min.pct         = -Inf,
        logfc.threshold = -Inf
      )
      write.csv(dge, csv_path)
    }, error = function(e) {
      cat(sprintf("  Error for %s x %s: %s\n", celltype, condition, e$message))
    })
  }

  rm(comparison_object, subset_injury)
  gc(verbose = FALSE)
}

rm(subset_control)
gc(verbose = FALSE)


###########################################
###########################################
###                                     ###
###       GO:BP ENRICHMENT (PHASE 2)    ###
###                                     ###
###########################################
###########################################

# Helper: run enrichGO on one gene list and save, with skip-existing + empty guards.
run_go <- function(genes, celltype, condition, suffix) {
  out_file <- file.path(
    go_dir,
    paste0("GOBP_", celltype, "_", condition, "_", suffix, ".csv")
  )

  if (skip_existing_go && file.exists(out_file)) {
    cat("  Skipping GO", suffix, "(existing):", celltype, "x", condition, "\n")
    return(invisible(NULL))
  }
  if (length(genes) == 0) {
    cat("  Skipping GO", suffix, "(no genes):  ", celltype, "x", condition, "\n")
    return(invisible(NULL))
  }

  cat("  Running GO", suffix, ":", celltype, "x", condition,
      "(", length(genes), "genes)\n")
  go <- enrichGO(
    gene    = genes,
    OrgDb   = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont     = "BP"
  )
  write.csv(go, out_file)
}


for (celltype in unique_idents) {
  for (condition in injury_conditions) {
    marker_file <- file.path(dge_dir, paste0("DGE_", celltype, "_", condition, ".csv"))

    if (!file.exists(marker_file)) {
      cat("Skipping GO (no DGE file):", celltype, "x", condition, "\n")
      next
    }

    marker_data <- read.csv(marker_file, stringsAsFactors = FALSE, row.names = 1)

    sig        <- !is.na(marker_data$p_val_adj) & marker_data$p_val_adj < 0.05
    up_genes   <- rownames(marker_data)[sig & marker_data$avg_log2FC > 0]
    down_genes <- rownames(marker_data)[sig & marker_data$avg_log2FC < 0]
    all_genes  <- rownames(marker_data)[sig]

    run_go(all_genes,  celltype, condition, "all")
    run_go(up_genes,   celltype, condition, "up")
    run_go(down_genes, celltype, condition, "down")
  }
}


message("10_main_dge_GO complete.")
