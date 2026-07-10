#' ======================================================================
#' 13_glia_dge_go.R
#'
#' Glia-specific differential expression, GO enrichment (two stringency
#' levels), tissue-stratified DGE, and SNC-specific KEGG enrichment.
#'
#' Phases:
#'   1. Main DGE: injury (AIH/DRC/SCI/SNC) vs control (WT) per minor glia
#'      cell type, using treatment_simple as the grouping.
#'   2. GO:BP (permissive): p_adj<0.05 on all / up / down DGE hits.
#'   3. GO:BP (stringent): |log2FC|>1, pct.1>5%, p_adj<0.05 on up / down.
#'   4. Tissue-stratified DGE: same as phase 1 but DRG-only and SN-only.
#'   5. KEGG enrichment on SNC only (relaxed threshold log2FC>0.25).
#'
#' Long-running. skip_existing toggles let you resume after interruption.
#'
#' Inputs:  ../r_objects/DRGSNMI_glia.RDS
#' Outputs: ../dge/glia_named/DGE_<celltype>_<condition>.csv
#'          ../dge/glia_named/{DRG,SN}/DGE_<celltype>_<condition>.csv
#'          ../dge_GO/glia_named/GOBP_<celltype>_<condition>_{all,up,down}.csv
#'          ../dge_GO/glia_named_stringent/GOBP_<celltype>_<condition>_{up,down}.csv
#'          ../dge_GO/glia_named_kegg/<direction>_genes_<celltype>_KEGG.csv
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
skip_existing_dge  <- TRUE
skip_existing_go   <- TRUE
skip_existing_kegg <- TRUE


# ===== CONFIG =====
injury_conditions  <- c("AIH", "DRC", "SCI", "SNC")
control_conditions <- c("WT")

# KEGG is run on SNC only (paper-focused, relaxed threshold)
kegg_condition <- "SNC"


# ===== LOAD GLIA OBJECT =====
glia_obj <- readRDS("../r_objects/DRGSNMI_glia.RDS")
DefaultAssay(glia_obj) <- "RNA"
cat("Loaded DRGSNMI_glia:", ncol(glia_obj), "cells\n")

# Safety: drop BATCH cells if any remain. They should have been removed by 07.
batch_mask <- as.character(Idents(glia_obj)) == "BATCH"
if (any(batch_mask)) {
  cat("Filtering BATCH cells from glia object (", sum(batch_mask), "cells)\n")
  glia_obj <- subset(glia_obj, idents = "BATCH", invert = TRUE)
}

# Canonical minor cell types present in the object
unique_celltypes <- as.character(unique(Idents(glia_obj)))
cat("Glia minor cell types (", length(unique_celltypes), "):\n", sep = "")
print(unique_celltypes)


# ===== OUTPUT DIRECTORIES =====
dge_dir        <- "../dge/glia_named"
go_dir         <- "../dge_GO/glia_named"
go_strict_dir  <- "../dge_GO/glia_named_stringent"
kegg_dir       <- "../dge_GO/glia_named_kegg"
dir.create(dge_dir,       recursive = TRUE)
dir.create(go_dir,        recursive = TRUE)
dir.create(go_strict_dir, recursive = TRUE)
dir.create(kegg_dir,      recursive = TRUE)


###########################################
###########################################
###                                     ###
###   HELPERS: DGE, GO, KEGG            ###
###                                     ###
###########################################
###########################################

# Run injury-vs-control DGE for every cell type in a given source object.
# Writes DGE_<celltype>_<condition>.csv into out_dir.
run_injury_vs_control_dge <- function(source_obj, out_dir,
                                      injury_conditions, control_conditions,
                                      celltypes, skip_existing = TRUE) {

  # Pre-compute the control subset once for this source object
  ctrl <- tryCatch(
    subset(source_obj, subset = treatment_simple %in% control_conditions),
    error = function(e) NULL
  )
  if (is.null(ctrl) || ncol(ctrl) == 0) {
    message("  No control cells found in this object — skipping DGE")
    return(invisible(NULL))
  }
  ctrl$condition <- "control"
  cat("  Control cells:", ncol(ctrl), "\n")

  for (condition in injury_conditions) {
    all_paths <- file.path(out_dir, paste0("DGE_", celltypes, "_", condition, ".csv"))
    if (skip_existing && all(file.exists(all_paths))) {
      cat("  Skipping", condition, "(all CSVs present)\n")
      next
    }

    inj <- tryCatch(
      subset(source_obj, subset = treatment_simple == condition),
      error = function(e) NULL
    )
    if (is.null(inj) || ncol(inj) == 0) {
      message("  No injury cells for condition ", condition, " — skipping")
      next
    }
    inj$condition <- "injury"
    cat("  ", condition, "injury cells:", ncol(inj), "\n")

    comparison_object <- merge(
      inj, y = ctrl, add.cell.ids = c("injury", "control"), project = "DGE_analysis"
    )
    comparison_object <- JoinLayers(comparison_object)
    DefaultAssay(comparison_object) <- "RNA"
    comparison_object$celltype.condition <- paste(
      Idents(comparison_object), comparison_object$condition, sep = "_"
    )
    Idents(comparison_object) <- "celltype.condition"

    for (celltype in celltypes) {
      out_path <- file.path(out_dir, paste0("DGE_", celltype, "_", condition, ".csv"))

      if (skip_existing && file.exists(out_path)) {
        cat("    Skipping DGE:", celltype, "(existing CSV)\n")
        next
      }

      tryCatch({
        cat("    Running DGE:", celltype, "injury vs control\n")
        dge <- FindMarkers(
          comparison_object,
          ident.1         = paste(celltype, "injury",  sep = "_"),
          ident.2         = paste(celltype, "control", sep = "_"),
          min.pct         = -Inf,
          logfc.threshold = -Inf
        )
        write.csv(dge, out_path)
      }, error = function(e) {
        cat(sprintf("    Error for %s x %s: %s\n", celltype, condition, e$message))
      })
    }

    rm(comparison_object, inj); gc(verbose = FALSE)
  }
}


# Run enrichGO on a single gene list and save to out_file.
run_go_for_gene_list <- function(genes, out_file, skip_existing = TRUE) {
  if (skip_existing && file.exists(out_file)) {
    cat("    Skipping GO (existing):", basename(out_file), "\n")
    return(invisible(NULL))
  }
  if (length(genes) == 0) {
    cat("    Skipping GO (no genes):  ", basename(out_file), "\n")
    return(invisible(NULL))
  }

  cat("    Running GO:", basename(out_file), "(", length(genes), "genes)\n")
  go <- tryCatch(
    enrichGO(gene = genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP"),
    error = function(e) { cat("      GO error:", e$message, "\n"); NULL }
  )
  if (!is.null(go)) write.csv(go, out_file)
}


# Loop GO over celltypes x conditions, applying provided up/down/all filters.
# up_filter and down_filter are functions taking a marker_data data frame and
# returning a logical vector of rows to keep.
run_go_phase <- function(dge_dir, go_dir, celltypes, conditions,
                         up_filter, down_filter, include_all = TRUE,
                         skip_existing = TRUE) {
  for (celltype in celltypes) {
    for (condition in conditions) {
      dge_file <- file.path(dge_dir, paste0("DGE_", celltype, "_", condition, ".csv"))
      if (!file.exists(dge_file)) {
        cat("  Skipping GO (no DGE):", celltype, "x", condition, "\n")
        next
      }
      marker_data <- read.csv(dge_file, row.names = 1, stringsAsFactors = FALSE)

      up_genes   <- rownames(marker_data)[up_filter(marker_data)]
      down_genes <- rownames(marker_data)[down_filter(marker_data)]

      run_go_for_gene_list(
        up_genes,
        file.path(go_dir, paste0("GOBP_", celltype, "_", condition, "_up.csv")),
        skip_existing
      )
      run_go_for_gene_list(
        down_genes,
        file.path(go_dir, paste0("GOBP_", celltype, "_", condition, "_down.csv")),
        skip_existing
      )

      if (include_all) {
        all_genes <- rownames(marker_data)[
          !is.na(marker_data$p_val_adj) & marker_data$p_val_adj < 0.05
        ]
        run_go_for_gene_list(
          all_genes,
          file.path(go_dir, paste0("GOBP_", celltype, "_", condition, "_all.csv")),
          skip_existing
        )
      }
    }
  }
}


# Convert mouse gene symbols to ENTREZ IDs for KEGG.
symbols_to_entrez <- function(gene_symbols) {
  mapIds(
    org.Mm.eg.db,
    keys    = gene_symbols,
    keytype = "SYMBOL",
    column  = "ENTREZID"
  )
}


run_kegg <- function(gene_symbols, celltype, direction, out_dir, skip_existing = TRUE) {
  out_file <- file.path(out_dir, paste0(direction, "_genes_", celltype, "_KEGG.csv"))

  if (skip_existing && file.exists(out_file)) {
    cat("  Skipping KEGG (existing):", basename(out_file), "\n")
    return(invisible(NULL))
  }
  if (length(gene_symbols) == 0) {
    cat("  Skipping KEGG (no genes):", basename(out_file), "\n")
    return(invisible(NULL))
  }

  cat("  Running KEGG:", celltype, direction, "(", length(gene_symbols), "genes)\n")
  gene_ids <- symbols_to_entrez(gene_symbols)

  kegg_res <- tryCatch({
    enrichKEGG(
      gene          = gene_ids,
      organism      = "mmu",
      keyType       = "ncbi-geneid",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05
    )
  }, error = function(e) {
    warning(paste0("KEGG enrichment failed for ", celltype, " (", direction, "): ", e$message))
    NULL
  })

  if (!is.null(kegg_res)) {
    write.csv(as.data.frame(kegg_res), file = out_file)
  }
}


# GO filters used below (kept as single-line lambdas for readability)
filter_up_permissive <- function(d)
  !is.na(d$p_val_adj) & d$p_val_adj < 0.05 & d$avg_log2FC > 0
filter_down_permissive <- function(d)
  !is.na(d$p_val_adj) & d$p_val_adj < 0.05 & d$avg_log2FC < 0

filter_up_stringent <- function(d)
  !is.na(d$p_val_adj) & d$p_val_adj < 0.05 & d$avg_log2FC > 1  & d$pct.1 > 0.05
filter_down_stringent <- function(d)
  !is.na(d$p_val_adj) & d$p_val_adj < 0.05 & d$avg_log2FC < -1 & d$pct.1 > 0.05


###########################################
###########################################
###                                     ###
###      PHASE 1 — MAIN GLIA DGE        ###
###                                     ###
###########################################
###########################################

message("\n===== Phase 1: main glia DGE =====")
run_injury_vs_control_dge(
  source_obj         = glia_obj,
  out_dir            = dge_dir,
  injury_conditions  = injury_conditions,
  control_conditions = control_conditions,
  celltypes          = unique_celltypes,
  skip_existing      = skip_existing_dge
)


###########################################
###########################################
###                                     ###
###   PHASE 2 — PERMISSIVE GO (BP)      ###
###                                     ###
###########################################
###########################################

message("\n===== Phase 2: permissive GO:BP =====")
run_go_phase(
  dge_dir       = dge_dir,
  go_dir        = go_dir,
  celltypes     = unique_celltypes,
  conditions    = injury_conditions,
  up_filter     = filter_up_permissive,
  down_filter   = filter_down_permissive,
  include_all   = TRUE,
  skip_existing = skip_existing_go
)


###########################################
###########################################
###                                     ###
###   PHASE 3 — STRINGENT GO (BP)       ###
###                                     ###
###########################################
###########################################

message("\n===== Phase 3: stringent GO:BP =====")
run_go_phase(
  dge_dir       = dge_dir,
  go_dir        = go_strict_dir,
  celltypes     = unique_celltypes,
  conditions    = injury_conditions,
  up_filter     = filter_up_stringent,
  down_filter   = filter_down_stringent,
  include_all   = FALSE,
  skip_existing = skip_existing_go
)


###########################################
###########################################
###                                     ###
###   PHASE 4 — TISSUE-STRATIFIED DGE   ###
###                                     ###
###########################################
###########################################

message("\n===== Phase 4: tissue-stratified DGE =====")
for (tissue_label in c("DRG", "SN")) {
  cat("\n--- Tissue:", tissue_label, "---\n")
  tissue_obj <- tryCatch(
    subset(glia_obj, subset = tissue == tissue_label),
    error = function(e) NULL
  )
  if (is.null(tissue_obj) || ncol(tissue_obj) == 0) {
    message("  No cells in tissue ", tissue_label, " — skipping")
    next
  }
  cat("  ", tissue_label, "cells:", ncol(tissue_obj), "\n")

  tissue_out_dir <- file.path(dge_dir, tissue_label)
  dir.create(tissue_out_dir, recursive = TRUE)

  run_injury_vs_control_dge(
    source_obj         = tissue_obj,
    out_dir            = tissue_out_dir,
    injury_conditions  = injury_conditions,
    control_conditions = control_conditions,
    celltypes          = unique_celltypes,
    skip_existing      = skip_existing_dge
  )

  rm(tissue_obj); gc(verbose = FALSE)
}


###########################################
###########################################
###                                     ###
###   PHASE 5 — SNC KEGG ENRICHMENT     ###
###                                     ###
###########################################
###########################################

message("\n===== Phase 5: SNC KEGG (relaxed log2FC>0.25) =====")
for (ct in unique_celltypes) {
  dge_file <- file.path(dge_dir, paste0("DGE_", ct, "_", kegg_condition, ".csv"))
  if (!file.exists(dge_file)) {
    cat("Skipping KEGG (no DGE file):", ct, "x", kegg_condition, "\n")
    next
  }

  cat("\n===== KEGG for", ct, "=====\n")
  dge <- read.csv(dge_file, row.names = 1, stringsAsFactors = FALSE)
  cat("  Total genes:", nrow(dge), "\n")

  up_genes <- rownames(dge)[
    !is.na(dge$p_val_adj) & dge$p_val_adj < 0.05 & dge$avg_log2FC > 0.25
  ]
  down_genes <- rownames(dge)[
    !is.na(dge$p_val_adj) & dge$p_val_adj < 0.05 & dge$avg_log2FC < -0.25
  ]

  run_kegg(up_genes,   celltype = ct, direction = "up",   out_dir = kegg_dir, skip_existing = skip_existing_kegg)
  run_kegg(down_genes, celltype = ct, direction = "down", out_dir = kegg_dir, skip_existing = skip_existing_kegg)
}


message("\n13_glia_dge_go complete.")
