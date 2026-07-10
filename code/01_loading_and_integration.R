#' ======================================================================
#' 01_loading_and_integration.R
#'
#' Load 21 batch Seurat objects and integrate via SCTransform + RPCA.
#'
#' Inputs:  ../r_objects/batch_objects/*.RDS
#' Outputs: ../r_objects/combined_init_object.RDS
#'
#' Author: Michael B. Thomsen, Ph.D.
#' Created: 2024-07-01
#' ======================================================================


###########################################
###########################################
###                                     ###
###     OBJECT AND METADATA LOADING     ###
###                                     ###
###########################################
###########################################

obj_directory <- "../r_objects/batch_objects/"

# Paths to your .RDS files
# 21 batches included; commented-out batches excluded to avoid
# redundant conditions/genotypes (KO, cKO, NODICAM) or timepoints (1/7 DPI)
file_paths <- c(
  "berta_bbi_2023.RDS",
  "cavalli_aih_seurat.RDS",
  "cavalli_contra_seurat.RDS",
  "cavalli_drc_seurat.RDS",
  "cavalli_naive_seurat.RDS",
  "cavalli_sci_seurat.RDS",
  "cavalli_snc_seurat.RDS",
  "giger_elife_2022_SN_0DPI_seurat.RDS",
  #"giger_elife_2022_SN_1DPI_seurat.RDS",
  "giger_elife_2022_SN_3DPI_seurat.RDS",
  #"giger_elife_2022_SN_7DPI_seurat.RDS",
  "hackos_neuron_2023_NaV1_7_WT_seurat.RDS",
  #"hackos_neuron_2023_NaV1_7_KO_seurat.RDS",
  #"hackos_neuron_2023_NaV1_7_cKO_seurat.RDS",
  #"horste_pnas_2019_NODICAM_seurat.RDS",
  #"horste_pnas_2019_WT_NODCTRL_seurat.RDS",
  "horste_pnas_2019_WT_seurat.RDS",
  "kaminker_natcomm_2023_WT_seurat.RDS",
  "kuruvilla_cellrep_2022_DRG_seurat.RDS",
  "milbrandt_natneuro_2022_sciatic_nerve_seurat.RDS",
  "miller_re_natcomm_2023_seurat.RDS",
  "SN_CTRL_SN_DS_CELL_MM_Miller_FD_GSM4423509.RDS",
  "SN_INJ_SN_3d_DS_CELL_MM_Miller_FD_GSM4423506.RDS",
  "SN_JUV_SN_BEAD_DS_CELL_MM_Miller_FD_GSM4423508.RDS",
  "SN_JUV_SN_FACS_DS_CELL_MM_Miller_FD_GSM4423507.RDS",
  "suter_elife_2021_P1_seurat.RDS",
  "suter_elife_2021_P60_seurat.RDS"
)

# Metadata fields to extract
selected_metadata_fields <- c(
    "orig.ident",
    "study",
    "study_simple",
    "age",
    "study_sex",
    "tissue",
    "treatment",
    "treatment_simple",
    "tech",
    "tech_simple",
    "cellnuc"
)


# Function to load a Seurat object, extract counts and selected metadata, and create a new Seurat object
load_and_extract <- function(file_path, metadata_fields) {
  # Load Seurat object from file
  seurat_obj <- readRDS(paste0(obj_directory, file_path))


  DefaultAssay(seurat_obj) <- "RNA"
  tryCatch({
    seurat_obj <- JoinLayers(seurat_obj)
  }, error = function(e) {
    print(paste("Error in JoinLayers for", file_path))
  })
  # Extract raw counts
  counts <- seurat_obj[["RNA"]]$counts

  # Extract selected metadata fields
  metadata <- seurat_obj@meta.data[, metadata_fields, drop = FALSE]

  # Create a new Seurat object with only the extracted data
  # min.features = 3: exclude cells with fewer than 3 detected genes
  new_seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata, min.features = 3)
  rm(seurat_obj) # remove original object from memory
  return(new_seurat_obj)
}

# Load and extract data from each file
seurat_objects <- lapply(file_paths, load_and_extract, metadata_fields = selected_metadata_fields)

# Combine all objects into a single Seurat object
combined_object <- merge(seurat_objects[[1]], y = seurat_objects[-1])

combined_object[["percent.mt"]] <- PercentageFeatureSet(object = combined_object, pattern = "^mt-")


# Examine the combined object pre-trim
combined_object

# Filter cells with >10% mitochondrial reads (likely dead/damaged cells)
combined_object <- subset(combined_object, subset = percent.mt < 10)

# Examine the combined object post-trim
combined_object
# Examine the metadata
head(combined_object@meta.data)

###########################################
###########################################
###                                     ###
###           BEGIN INTEGRATION         ###
###                                     ###
###########################################
###########################################

# Integration, PCA, clustering, and UMAP
options(future.globals.maxSize = 3e+09)  # 3 GB — needed for large SCTransform jobs
combined_object <- SCTransform(combined_object)
combined_object <- RunPCA(combined_object, npcs = 50)
combined_object <- IntegrateLayers(
  object = combined_object,
  method = RPCAIntegration,
  normalization.method = "SCT"
)
# dims = 1:50: use all 50 PCs for neighbor graph and UMAP embedding
combined_object <- FindNeighbors(combined_object, dims = 1:50, reduction = "integrated.dr")
# resolution = 0.5: initial exploratory clustering (refined in 02_init_clustree.R)
combined_object <- FindClusters(combined_object, resolution = 0.5)
combined_object <- RunUMAP(combined_object, reduction = "integrated.dr", dims = 1:50)


# Plot UMAP, group by 'seurat_clusters'
p_init <- DimPlot(combined_object, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, raster = F) + NoLegend()

# Display plot (saved in next file/step)
#print(p_init)

# Save combined object
saveRDS(combined_object, "../r_objects/combined_init_object.RDS")

# LAST RUN 20241208
