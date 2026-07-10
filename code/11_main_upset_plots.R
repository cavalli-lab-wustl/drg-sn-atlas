#' ======================================================================
#' 11_main_upset_plots.R
#'
#' Per-celltype UpSet plots of DE gene overlap across the 4 injury
#' conditions. For each Major cell type, reads the DGE CSVs produced by
#' 10_main_dge_GO, extracts p_adj<0.05 gene sets per condition, and
#' generates an UpSet plot plus a tidy CSV of overlap groups → genes.
#'
#' Inputs:  ../dge/combined_named/DGE_<celltype>_<condition>.csv
#' Outputs: ../dge/combined_named/upset/<celltype>_upset.pdf
#'          ../dge/combined_named/upset/<celltype>_upset_groups_and_genes.csv
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
library(UpSetR)


# ===== CONFIG =====
conditions     <- c("cavalli_aih", "cavalli_drc", "cavalli_sci", "cavalli_snc")
condition_labs <- c("AIH",         "DRC",         "SCI",         "SNC")

dge_dir    <- "../dge/combined_named"
output_dir <- file.path(dge_dir, "upset")
dir.create(output_dir, recursive = TRUE)


# ===== HELPERS =====
# Read the 4 DGE CSVs for one cell type and return list of significant gene vectors.
load_dge_lists <- function(cell_type) {
  lists <- list()
  for (i in seq_along(conditions)) {
    f <- file.path(dge_dir, paste0("DGE_", cell_type, "_", conditions[i], ".csv"))
    if (!file.exists(f)) return(NULL)
    df <- read.csv(f)
    lists[[condition_labs[i]]] <- df[[1]][!is.na(df$p_val_adj) & df$p_val_adj < 0.05]
  }
  lists
}

# Enumerate all intersection groups and which genes belong to each.
process_overlap <- function(listInput, sort = TRUE) {
  listInputmat    <- fromList(listInput) == 1
  listInputunique <- unique(listInputmat)
  grouplist <- list()

  for (i in seq_len(nrow(listInputunique))) {
    currentRow <- listInputunique[i, ]
    myelements <- which(apply(listInputmat, 1, function(x) all(x == currentRow)))
    attr(myelements, "groups") <- currentRow
    grouplist[[paste(colnames(listInputunique)[currentRow], collapse = ":")]] <- myelements
  }

  if (sort) {
    grouplist <- grouplist[order(sapply(grouplist, length), decreasing = TRUE)]
  }
  attr(grouplist, "elements") <- unique(unlist(listInput))
  grouplist
}

# Flatten overlap groups to a long-form Group/Gene CSV.
save_overlap_csv <- function(li, filename) {
  df <- do.call(rbind, lapply(seq_along(li), function(i) {
    data.frame(Group = names(li)[i], Gene = attr(li, "elements")[li[[i]]])
  }))
  write.csv(df, filename, row.names = FALSE)
}


# ===== MAIN LOOP =====
# Iterate over canonical class-grouped cell-type labels (from utils/variables.R).
for (cell_type in class_ordered_majors) {
  degs_list <- tryCatch(load_dge_lists(cell_type), error = function(e) {
    cat("  Error loading DGE for", cell_type, ":", e$message, "\n"); NULL
  })

  if (is.null(degs_list)) {
    cat("Skipping", cell_type, "- missing DGE CSVs\n")
    next
  }
  if (length(unlist(degs_list)) == 0) {
    cat("Skipping", cell_type, "- no significant DEGs in any condition\n")
    next
  }

  cat("Processing", cell_type, "\n")

  tryCatch({
    pdf(file.path(output_dir, paste0(cell_type, "_upset.pdf")), width = 5, height = 4)
    print(upset(
      fromList(degs_list),
      sets          = condition_labs,
      order.by      = "freq",
      keep.order    = TRUE,
      matrix.color  = "#005C97",
      text.scale    = 1.5
    ))
    dev.off()

    li <- process_overlap(degs_list)
    save_overlap_csv(
      li,
      filename = file.path(output_dir, paste0(cell_type, "_upset_groups_and_genes.csv"))
    )
  }, error = function(e) {
    cat("  Error plotting / saving for", cell_type, ":", e$message, "\n")
    if (dev.cur() != 1) dev.off()
  })
}

message("11_main_upset_plots complete.")
