#' ======================================================================
#' 16_glia_cellchat.R
#'
#' CellChat intercellular communication analysis (naive vs SNC conditions).
#'
#' Inputs:  ../r_objects/DRGSNMI_glia.RDS, ../r_objects/DRGSNMI_main_multilevel.RDS
#' Outputs: ../cellchat/combined/cellchat_subset_*.RDS, ../cellchat/combined/*_output_full.csv
#'
#' Author: Michael B. Thomsen, Ph.D.
#' Created: 2024-09-10
#' ======================================================================

###########################################
###########################################
###                                     ###
###    DATA LOADING AND PREP            ###
###                                     ###
###########################################
###########################################

# Load glia object
glia <- readRDS("../r_objects/DRGSNMI_glia.RDS")

table(glia@meta.data$orig.ident)

# Export idents (to be used in main object)
glia_labels <- Idents(glia)
glia_cells <- Cells(glia)

# Load main object
main <- readRDS("../r_objects/DRGSNMI_main_multilevel.RDS")

Idents(main, cells = glia_cells) <- glia_labels

# REMOVE background, technical artifacts
main <- subset(main, idents = c("Satellite Glial Cells", "Non-myelinating Schwann Cells", "Myelinating Schwann Cells"), invert = TRUE)

main <- subset(main, idents = c("Mast Cells", "Erythrocytes", "mSC III", "Lymphatic Endothelial Cells"), invert = TRUE)

naive_obj <- subset(main, subset = orig.ident %in% c(
    "DRG_naive_naive_unpb_10X_CELL_MM_Cavalli",
    "DRG_contra_contra_el21_10X_CELL_MM_Cavalli",
    "DRG_contra_contra_nc20_10X_CELL_MM_Cavalli",
    "DRG_contra_contra_unpb_10X_CELL_MM_Cavalli"
))

# load 3 subsets at once
snc_obj <- subset(main, subset = orig.ident %in% c(
    "DRG_snc_snc_nc20_10X_CELL_MM_Cavalli",
    "DRG_snc_snc_el21_10X_CELL_MM_Cavalli",
    "DRG_snc_snc_unpb_10X_CELL_MM_Cavalli"
))

naive_obj
snc_obj

DimPlot(naive_obj, label = T) + NoLegend()
DimPlot(snc_obj, label = T) + NoLegend()

table(snc_obj@meta.data$orig.ident)
table(naive_obj@meta.data$orig.ident)

table(Idents(snc_obj))
table(Idents(naive_obj))

###########################################
###########################################
###                                     ###
###    CELLCHAT                         ###
###                                     ###
###########################################
###########################################

library(CellChat)
packageVersion('CellChat')

output_directory <- "../cellchat/combined"
dir.create(output_directory, recursive = TRUE)

###########################################
###########################################
###                                     ###
###  NAIVE SUBSET                       ###
###                                     ###
###########################################
###########################################

# Subset for true control cells
subset_naive <- naive_obj

### 1. Prepare required input data for CellChat

# Select sample split metadata variable
sample_metadata_slot <- "orig.ident"

# Extract the CellChat input files from a Seurat V5 object
data.input <- GetAssayData(subset_naive, assay = "RNA", layer = "data") # normalized data matrix
labels <- Idents(subset_naive)
samples <- subset_naive@meta.data[[sample_metadata_slot]] # sample metadata
meta <- data.frame(group = labels, samples = samples, row.names = names(labels)) # create a dataframe with cell labels and samples

### 2. Create a CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group", assay = "RNA")

### 3. Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)
cellchat@DB <- CellChatDB.use

### 4. Identify over-represented ligands or receptors
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

### 5. Infer cell-cell communication at a ligand-receptor pair level
# NOTE: triMean is more stringent than truncatedMean (fewer false positives)
cellchat <- computeCommunProb(cellchat, type = "triMean")

### 6. Filter the cell-cell communication
cellchat <- filterCommunication(cellchat, min.cells = 10) # min 10 cells per group to be confident

df.net <- subsetCommunication(cellchat)
write.csv(df.net, "../cellchat/combined/subset_naive_cellchat_output_full.csv")

### 7. Infer cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

### 8. Calculate aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

### PAUSE AND SAVE OBJECT
saveRDS(cellchat, file = file.path(output_directory, "cellchat_subset_naive.RDS"))

###########################################
###########################################
###                                     ###
###  SNC                                ###
###                                     ###
###########################################
###########################################

# Subset for SNC condition cells
subset_snc <- snc_obj

### 1. Prepare required input data for CellChat

# Select sample split metadata variable
sample_metadata_slot <- "orig.ident"

# Extract the CellChat input files from a Seurat V5 object
data.input <- GetAssayData(subset_snc, assay = "RNA", layer = "data") # normalized data matrix
labels <- Idents(subset_snc)
samples <- subset_snc@meta.data[[sample_metadata_slot]] # sample metadata
meta <- data.frame(group = labels, samples = samples, row.names = names(labels)) # create a dataframe with cell labels and samples

### 2. Create a CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group", assay = "RNA")

### 3. Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)
cellchat@DB <- CellChatDB.use

### 4. Identify over-represented ligands or receptors
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

### 5. Infer cell-cell communication at a ligand-receptor pair level
# NOTE: triMean is more stringent than truncatedMean (fewer false positives)
cellchat <- computeCommunProb(cellchat, type = "triMean")

### 6. Filter the cell-cell communication
cellchat <- filterCommunication(cellchat, min.cells = 10) # min 10 cells per group to be confident

df.net <- subsetCommunication(cellchat)
write.csv(df.net, "../cellchat/combined/subset_snc_cellchat_output_full.csv")

### 7. Infer cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

### 8. Calculate aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

### PAUSE AND SAVE OBJECT
saveRDS(cellchat, file = file.path(output_directory, "cellchat_subset_snc.RDS"))

###########################################
###########################################
###                                     ###
###  SNC VS NAIVE                       ###
###                                     ###
###########################################
###########################################

cellchat.naive <- readRDS("../cellchat/combined/cellchat_subset_naive.RDS")
cellchat.snc <- readRDS("../cellchat/combined/cellchat_subset_snc.RDS")
object.list <- list(control = cellchat.naive, snc = cellchat.snc)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

save(cellchat, file = file.path(output_directory, "cellchat_snc_vs_naive.RData"))

# Load RData
load(file = file.path(output_directory, "cellchat_snc_vs_naive.RData"))

### 17. Compare the total number of interactions and interaction strength between groups

gg1 <- compareInteractions(cellchat, group = c(1,2))
gg2 <- compareInteractions(cellchat, group = c(1,2), measure = "weight")

# save png and pdf versions of gg1 + gg2 5"x5"
ggsave(file = file.path(output_directory, "cellchat_snc_vs_naive.png"), gg1 + gg2, width = 5, height = 5)
ggsave(file = file.path(output_directory, "cellchat_snc_vs_naive.pdf"), gg1 + gg2, width = 5, height = 5)

### 18. Compare the number of interactions and interaction strength among diffeent cell types

par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")

gg3 <- netVisual_heatmap(cellchat)
gg4 <- netVisual_heatmap(cellchat, measure = "weight")

pdf(file = file.path(output_directory, "cellchat_snc_vs_naive_heatmap.pdf"), width = 5, height = 4.5)
print(gg3)
dev.off()

pdf(file = file.path(output_directory, "cellchat_snc_vs_naive_heatmap_weight.pdf"), width = 5, height = 4.5)
print(gg4)
dev.off()

pdf(file = file.path(output_directory, "cellchat_snc_vs_naive_heatmap_combined.pdf"), width = 7.5, height = 4.5)
print(gg3+gg4)
dev.off()

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))

par(mfrow = c(1,2), xpd = TRUE)
for (i in 1:length(object.list)) {
    netVisual_circle(
        object.list[[i]]@net$count,
        weight.scale = TRUE,
        label.edge = FALSE,
        edge.weight.max = weight.max[2],
        edge.width.max = 12,
        title.name = paste0("Number of interactions - ", names(object.list[i]))
    )
}

# Here, CellChat categorize the cell populations into three celltypes, and then re-merge the list of CellChat objects.
group.cellType <- c(rep("Satellite Glial Cells", 4), rep("Fibroblasts", 4), rep("Macrophages", 4), rep("Non−Myelinating Schwann Cells", 4))
group.cellType <- factor(group.cellType, levels = c("Satellite Glial Cells", "Fibroblasts","Macrophages", "Non−Myelinating Schwann Cells"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Show the number of interactions or interaction strength between any two cell types in each dataset.
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))

pdf(file = file.path(output_directory, "cellchat_snc_vs_naive_circle.pdf"), width = 10, height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
    netVisual_circle(
        object.list[[i]]@net$count.merged,
        weight.scale = T,
        label.edge= T,
        edge.weight.max = weight.max[3],
        edge.width.max = 12,
        title.name = paste0("Number of interactions - ", names(object.list)[i])
    )
}
dev.off()

# Similarly, CellChat can also show the differential number of interactions or interaction strength between any two cell types using circle plot.
pdf(file = file.path(output_directory, "cellchat_snc_vs_naive_circle2.pdf"), width = 10, height = 10)
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", label.edge = T)
dev.off()

# Re-load CellChat objects for centrality and similarity analysis
cellchat.naive <- readRDS("../cellchat/combined/cellchat_subset_naive.RDS")
cellchat.snc <- readRDS("../cellchat/combined/cellchat_subset_snc.RDS")

# netAnalysis_computeCentrality
cellchat.naive <- netAnalysis_computeCentrality(cellchat.naive)
cellchat.snc <- netAnalysis_computeCentrality(cellchat.snc)

cellchat.naive <- computeNetSimilarity(cellchat.naive, type = "functional")
cellchat.naive <- netEmbedding(cellchat.naive, type = "functional")
cellchat.naive <- netClustering (cellchat.naive, type = "functional", do.parallel = FALSE)
cellchat.naive <- computeNetSimilarity(cellchat.naive, type = "structural")
cellchat.naive <- netEmbedding(cellchat.naive, type = "structural")
cellchat.naive <- netClustering(cellchat.naive, type = "structural",do.parallel = FALSE)

cellchat.snc <- computeNetSimilarity(cellchat.snc, type = "functional")
cellchat.snc <- netEmbedding(cellchat.snc, type = "functional")
cellchat.snc <- netClustering(cellchat.snc, type = "functional", do.parallel = FALSE)
cellchat.snc <- computeNetSimilarity(cellchat.snc, type = "structural")
cellchat.snc <- netEmbedding(cellchat.snc, type = "structural")
cellchat.snc <- netClustering(cellchat.snc, type = "structural", do.parallel = FALSE)

object.list <- list(naive = cellchat.naive, snc = cellchat.snc)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()

for (i in 1:length(object.list)) {
    gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

# save PDF
pdf(file = file.path(output_directory, "cellchat_snc_vs_naive_signalingRole.pdf"), width = 10, height = 5)
patchwork::wrap_plots(plots = gg)
dev.off()

# Create plots with fixed axes
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]],
                                               title = names(object.list)[i],
                                               weight.MinMax = weight.MinMax) +
    ggplot2::xlim(0, 31) +  # Set x-axis limit
    ggplot2::ylim(0, 31)    # Set y-axis limit
}

# Save PDF
pdf(file = file.path(output_directory, "cellchat_snc_vs_naive_signalingRole_scaled.pdf"), width = 10, height = 5)
patchwork::wrap_plots(plots = gg)
dev.off()

library(ggplot2)
library(ggrepel)

# Extract signaling role data for naive and snc conditions (v1: simple overlay, no color mapping)
extract_signaling_data <- function(cellchat_object, x.measure = "outdeg", y.measure = "indeg") {
  centr <- slot(cellchat_object, "netP")$centr
  outgoing <- rowSums(sapply(centr, function(x) x[[x.measure]]))
  incoming <- rowSums(sapply(centr, function(x) x[[y.measure]]))
  num.link <- rowSums(cellchat_object@net$count) + colSums(cellchat_object@net$count) - diag(cellchat_object@net$count)
  data.frame(cluster = names(outgoing), x = outgoing, y = incoming, size = num.link)
}

# Extract data for both naive and snc
df_naive <- extract_signaling_data(cellchat.naive)
df_snc <- extract_signaling_data(cellchat.snc)

# Merge the data frames by cluster
df_merged <- merge(df_naive, df_snc, by = "cluster", suffixes = c("_naive", "_snc"))

# Create the plot
p <- ggplot() +
  geom_point(data = df_merged, aes(x = x_naive, y = y_naive), color = "blue", size = 3, alpha = 0.6) +
  geom_point(data = df_merged, aes(x = x_snc, y = y_snc), color = "red", size = 3, alpha = 0.6) +
  geom_segment(data = df_merged, aes(x = x_naive, y = y_naive, xend = x_snc, yend = y_snc),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  xlim(0, 31) +
  ylim(0, 31) +
  labs(title = "Change in Signaling Role: Naive (blue) to SNC (red)",
       x = "Outgoing interaction strength",
       y = "Incoming interaction strength") +
  theme_minimal()

# Save the plot as PDF
pdf(file = file.path(output_directory, "cellchat_overlay_naive_snc.pdf"), width = 7, height = 7)
print(p)
dev.off()

# v2: Overlay with cluster-matched colors and size scaling
# Redefine extract_signaling_data with color parameter
extract_signaling_data <- function(cellchat_object, x.measure = "outdeg", y.measure = "indeg", color.use = NULL) {
  centr <- slot(cellchat_object, "netP")$centr
  outgoing <- rowSums(sapply(centr, function(x) x[[x.measure]]))
  incoming <- rowSums(sapply(centr, function(x) x[[y.measure]]))
  num.link <- rowSums(cellchat_object@net$count) + colSums(cellchat_object@net$count) - diag(cellchat_object@net$count)
  data.frame(cluster = names(outgoing), x = outgoing, y = incoming, size = num.link, color = color.use)
}

# Color mapping based on the number of clusters in cellchat.naive
color.use <- scPalette(nlevels(cellchat.naive@idents))

# Extract data for both naive and snc, keeping the same color scheme
df_naive <- extract_signaling_data(cellchat.naive, color.use = color.use)
df_snc <- extract_signaling_data(cellchat.snc, color.use = color.use)

# Merge the data frames by cluster
df_merged <- merge(df_naive, df_snc, by = "cluster", suffixes = c("_naive", "_snc"))

# Create the plot, using the same size and color as the original
p <- ggplot() +
  geom_point(data = df_merged, aes(x = x_naive, y = y_naive, size = size_naive, color = cluster), alpha = 0.6) +
  geom_point(data = df_merged, aes(x = x_snc, y = y_snc, size = size_snc, color = cluster), alpha = 0.6) +
  geom_segment(data = df_merged, aes(x = x_naive, y = y_naive, xend = x_snc, yend = y_snc),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  xlim(0, 31) +
  ylim(0, 31) +
  scale_size_continuous(range = c(2, 6)) +
  scale_color_manual(values = color.use) +
  labs(title = "Change in Signaling Role: Naive to SNC",
       x = "Outgoing interaction strength",
       y = "Incoming interaction strength") +
  theme_minimal()

# Save the plot as PDF
pdf(file = file.path(output_directory, "cellchat_overlay_naive_snc_same_size_color.pdf"), width = 12, height = 7)
print(p)
dev.off()

# v3: Overlay with positive-slope highlighting (solid vs dashed arrows)
df_naive <- extract_signaling_data(cellchat.naive, color.use = color.use)
df_snc <- extract_signaling_data(cellchat.snc, color.use = color.use)
df_merged <- merge(df_naive, df_snc, by = "cluster", suffixes = c("_naive", "_snc"))

# Flag clusters with increased outgoing or incoming signaling after SNC
df_merged$positive_slope <- (df_merged$x_snc - df_merged$x_naive > 0) | (df_merged$y_snc - df_merged$y_naive > 0)

p <- ggplot() +
  geom_point(data = df_merged, aes(x = x_naive, y = y_naive, size = size_naive, color = cluster), alpha = 0.6) +
  geom_point(data = df_merged, aes(x = x_snc, y = y_snc, size = size_snc, color = cluster), alpha = 0.6) +
  geom_segment(data = df_merged, aes(x = x_naive, y = y_naive, xend = x_snc, yend = y_snc,
                                     color = cluster, linetype = positive_slope),
               arrow = arrow(length = unit(0.2, "cm")), size = 1) +
  xlim(0, 31) +
  ylim(0, 31) +
  scale_size_continuous(range = c(2, 6)) +
  scale_color_manual(values = color.use) +
  scale_linetype_manual(values = c("dashed", "solid"), guide = FALSE) +  # solid = increased signaling
  labs(title = "Change in Signaling Role: Naive to SNC",
       x = "Outgoing interaction strength",
       y = "Incoming interaction strength") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

# Save the plot as PDF
pdf(file = file.path(output_directory, "cellchat_overlay_naive_snc_with_border_positive_slope.pdf"), width = 12, height = 7)
print(p)
dev.off()

# Loop over "idents" to generate plots for each cell type
for (cell_type in unique(Idents(naive_obj))) {
    gg <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = cell_type)
    pdf(file = file.path(output_directory, paste0("cellchat_snc_vs_naive_signalingChanges_", cell_type, ".pdf")), width = 7, height = 5)
    print(gg)
    dev.off()
}

### 21. Identify altered signaling pathways

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)

ggsave(file = file.path(output_directory, "cellchat_snc_vs_naive_rankNet_stacked.png"), gg1, width = 5, height = 10)
ggsave(file = file.path(output_directory, "cellchat_snc_vs_naive_rankNet.png"), gg2, width = 5, height = 10)

# PDFs
pdf(file = file.path(output_directory, "cellchat_snc_vs_naive_rankNet_stacked.pdf"), width = 5, height = 10)
print(gg1)
dev.off()

pdf(file = file.path(output_directory, "cellchat_snc_vs_naive_rankNet.pdf"), width = 5, height = 10)
print(gg2)
dev.off()

library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height= 24, font.size = 5)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 20, font.size = 5)

# save as PDF
pdf(file = file.path(output_directory, "cellchat_snc_vs_naive_signalingRole_heatmap.pdf"), width = 10, height = 15)
ht1 + ht2
dev.off()

# Compare the communication probabilities from certain cell groups to other cell groups
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), comparison= c(1, 2), angle.x = 45)
# Identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset.
gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in AIH",angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in IH",angle.x = 45, remove.isolate = T)

# save as PDFs
pdf(file = file.path(output_directory, "cellchat_snc_vs_naive_bubble_increased.pdf"), width = 7, height = 10)
print(gg1)
dev.off()

pdf(file = file.path(output_directory, "cellchat_snc_vs_naive_bubble_decreased.pdf"), width = 7, height = 10)
print(gg2)
dev.off()

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "snc"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "snc",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "snc",ligand.logFC = -0.1, receptor.logFC = -0.1)
# do further deconvolution to obtain the individual signaling genes
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
# Users can also find all the significant outgoing/incoming/both signaling according to the customized DEG features and cell groups of interest
df <- findEnrichedSignaling(object.list[[2]], features = c("JAM2", "JAM3"), idents = c("Neurons", "Satellite Glial Cells"), pattern ="outgoing")

### 24. Visualize the identified up-regulated and down-regulated signaling L-R pairs

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4,targets.use = c(5:11), comparison = c(1, 2), angle.x = 90, remove.isolate =T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use =4, targets.use = c(5:11), comparison = c(1, 2), angle.x = 90, remove.isolate= T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

# save PDFs
pdf(file = file.path(output_directory, "cellchat_snc_vs_naive_bubble_up.pdf"), width = 7, height = 10)
print(gg1)
dev.off()

pdf(file = file.path(output_directory, "cellchat_snc_vs_naive_bubble_down.pdf"), width = 7, height = 10)
print(gg2)
dev.off()

### 25. Visualize inferred cell-cell communication networks

pathways.show <- c("ApoE")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name =paste(pathways.show, names(object.list)[i]))
}

# save PDF
pdf(file = file.path(output_directory, "cellchat_COLLAGEN_snc_vs_naive_aggregate=.pdf"), width = 10, height = 5)
for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name =paste(pathways.show, names(object.list)[i]))
}
dev.off()

pathways.show <- c("COLLAGEN")
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
    ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# save PDF
pdf(file = file.path(output_directory, "cellchat_COLLAGEN_snc_vs_naive_heatmap.pdf"), width = 10, height = 5)
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
dev.off()

### 26. Visualize gene expression distribution

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("snc", "naive")) # set factor level
p <- plotGeneExpression(cellchat, signaling = "ApoE", split.by = "datasets", colors.ggplot = T)
ggsave(file = file.path(output_directory, "cellchat_apoe_snc_vs_naive_geneExpression.png"), p, width = 5, height = 7)
ggsave(file = file.path(output_directory, "cellchat_apoe_snc_vs_naive_geneExpression.pdf"), p, width = 5, height = 7)

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("snc", "naive")) # set factor level
p <- plotGeneExpression(cellchat, signaling = "MIF", split.by = "datasets")
ggsave(file = file.path(output_directory, "cellchat_mif_snc_vs_naive_geneExpression.png"), p, width = 5, height = 7)
ggsave(file = file.path(output_directory, "cellchat_mif_snc_vs_naive_geneExpression.pdf"), p, width = 5, height = 7)

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("snc", "naive")) # set factor level
p <- plotGeneExpression(cellchat, signaling = "ApoE", split.by = "datasets")
ggsave(file = file.path(output_directory, "cellchat_apoe_snc_vs_naive_geneExpression.png"), p, width = 5, height = 7)
ggsave(file = file.path(output_directory, "cellchat_apoe_snc_vs_naive_geneExpression.pdf"), p, width = 5, height = 7)

pathways <- c(
    "JAM",
    "MIF",
    "NCAM",
    "PSAP",
    "CADM",
    "ESAM",
    "CCL",
    "NRXN",
    "ADGRL",
    "COLLAGEN",
    "LAMININ",
    "PTN",
    "APP",
    "ApoE",
    "FN1"
)

for (pathway in pathways) {
    p <- plotGeneExpression(cellchat, signaling = pathway, split.by = "datasets", colors.ggplot = T)
    ggsave(file = file.path(output_directory, paste0("cellchat_", pathway, "_snc_vs_naive_geneExpression.png")), p, width = 5, height = 7)
    ggsave(file = file.path(output_directory, paste0("cellchat_", pathway, "_snc_vs_naive_geneExpression.pdf")), p, width = 5, height = 7)
}

###########################################
###########################################
###                                     ###
###    SARAH CUSTOM REQUEST             ###
###                                     ###
###########################################
###########################################

# Load required libraries
library(CellChat)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)

# Load RData
cellchat.naive <- readRDS("../cellchat/combined/cellchat_subset_naive.RDS")
cellchat.snc <- readRDS("../cellchat/combined/cellchat_subset_snc.RDS")

# netAnalysis_computeCentrality
cellchat.naive <- netAnalysis_computeCentrality(cellchat.naive)
cellchat.snc <- netAnalysis_computeCentrality(cellchat.snc)

# Compute network similarities and embeddings
cellchat.naive <- computeNetSimilarity(cellchat.naive, type = "functional")
cellchat.naive <- netEmbedding(cellchat.naive, type = "functional")
cellchat.naive <- netClustering(cellchat.naive, type = "functional", do.parallel = FALSE)
cellchat.naive <- computeNetSimilarity(cellchat.naive, type = "structural")
cellchat.naive <- netEmbedding(cellchat.naive, type = "structural")
cellchat.naive <- netClustering(cellchat.naive, type = "structural", do.parallel = FALSE)

cellchat.snc <- computeNetSimilarity(cellchat.snc, type = "functional")
cellchat.snc <- netEmbedding(cellchat.snc, type = "functional")
cellchat.snc <- netClustering(cellchat.snc, type = "functional", do.parallel = FALSE)
cellchat.snc <- computeNetSimilarity(cellchat.snc, type = "structural")
cellchat.snc <- netEmbedding(cellchat.snc, type = "structural")
cellchat.snc <- netClustering(cellchat.snc, type = "structural", do.parallel = FALSE)

# Merge the CellChat objects into a list
object.list <- list(naive = cellchat.naive, snc = cellchat.snc)

# Define the output directory
output_directory <- "../cellchat/combined/"

# Define cell types of interest
cell_types <- levels(object.list[[1]]@idents)

# Define separate color palettes for each condition
cell_type_colors <- list(
  naive = c(
    "SGC IMM" = "#E41A1C",        # Red
    "Macrophages" = "#377EB8",    # Blue
    "Neurons" = "#4DD091"         # Green
  ),
  snc = c(
    "SGC IMM" = "#FF7F00",        # Orange
    "Macrophages" = "#984EA3",    # Purple
    "Neurons" = "#FFEF59"         # Magenta
  )
)

# Ensure that the colors are assigned only to the cell types of interest
cell_type_colors$naive <- cell_type_colors$naive[names(cell_type_colors$naive) %in% cell_types]
cell_type_colors$snc <- cell_type_colors$snc[names(cell_type_colors$snc) %in% cell_types]

# Function to generate chord diagrams with condition-specific colors and updated legends
generate_chord_diagrams <- function(object.list, source_cell, target_cell, direction, output_filename, cell_type_colors) {
  # direction can be "forward" (source -> target) or "reverse" (target -> source)
  pdf(file = file.path(output_directory, output_filename), width = 12, height = 8)

  # Define layout: 2 chord diagrams on top, legends below
  layout_matrix <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
  layout(layout_matrix, heights = c(6, 1))  # Adjust heights as needed

  # Plot each chord diagram without legends
  for (i in seq_along(object.list)) {
    condition <- names(object.list)[i]
    colors <- cell_type_colors[[condition]]

    if (direction == "forward") {
      sources.use <- source_cell
      targets.use <- target_cell
      title_direction <- paste0(source_cell, " -> ", target_cell)
    } else {
      sources.use <- target_cell
      targets.use <- source_cell
      title_direction <- paste0(target_cell, " -> ", source_cell)
    }

    # Generate chord diagram without legend
    netVisual_chord_gene(
      object.list[[i]],
      sources.use = sources.use,
      targets.use = targets.use,
      slot.name = "net",  # Use "net" to visualize at the ligand-receptor level
      color.use = colors,
      title.name = paste0(title_direction, " - ", condition),
      legend.pos.x = -1,
      legend.pos.y = -1,
      lab.cex = .8,
      link.visible = TRUE,
      small.gap = 1,
      big.gap = 10,
      annotationTrackHeight = c(0.05),
      show.legend = FALSE
    )
  }

  # Create a combined legend with all 6 colors
  combined_labels <- c(
    "naive_SGC IMM",
    "snc_SGC IMM",
    "naive_Macrophages",
    "snc_Macrophages",
    "naive_Neurons",
    "snc_Neurons"
  )

  combined_colors <- c(
    cell_type_colors$naive["SGC IMM"],
    cell_type_colors$snc["SGC IMM"],
    cell_type_colors$naive["Macrophages"],
    cell_type_colors$snc["Macrophages"],
    cell_type_colors$naive["Neurons"],
    cell_type_colors$snc["Neurons"]
  )

  names(combined_colors) <- combined_labels

  # Create the legend using ComplexHeatmap
  lgd <- ComplexHeatmap::Legend(
    at = names(combined_colors),
    type = "grid",
    legend_gp = grid::gpar(fill = combined_colors),
    title = "Cell Type (Condition)",
    nrow = 2
  )

  # Draw the legend in the bottom layout area
  ComplexHeatmap::draw(lgd,
                       x = unit(0.5, "npc"),
                       y = unit(0.2, "npc"),
                       just = "center")

  dev.off()
}

# Generate chord diagrams for SGC IMM -> Macrophages
generate_chord_diagrams(
  object.list = object.list,
  source_cell = "SGC IMM",
  target_cell = "Macrophages",
  direction = "forward",
  output_filename = "SGC_IMM_to_Macrophages_chord_diagram.pdf",
  cell_type_colors = cell_type_colors
)

# Generate chord diagrams for Macrophages -> SGC IMM
generate_chord_diagrams(
  object.list = object.list,
  source_cell = "SGC IMM",
  target_cell = "Macrophages",
  direction = "reverse",
  output_filename = "Macrophages_to_SGC_IMM_chord_diagram.pdf",
  cell_type_colors = cell_type_colors
)

# Generate chord diagrams for SGC IMM -> Neurons
generate_chord_diagrams(
  object.list = object.list,
  source_cell = "SGC IMM",
  target_cell = "Neurons",
  direction = "forward",
  output_filename = "SGC_IMM_to_Neurons_chord_diagram.pdf",
  cell_type_colors = cell_type_colors
)

# Generate chord diagrams for Neurons -> SGC IMM
generate_chord_diagrams(
  object.list = object.list,
  source_cell = "SGC IMM",
  target_cell = "Neurons",
  direction = "reverse",
  output_filename = "Neurons_to_SGC_IMM_chord_diagram.pdf",
  cell_type_colors = cell_type_colors
)

# Function to generate chord diagrams for specific ligand-receptor pairs with updated legends
generate_chord_diagrams_lr <- function(object.list, source_cell, target_cell, output_filename, cell_type_colors, min_interaction_weight = 0) {
  pdf(file = file.path(output_directory, output_filename), width = 12, height = 8)

  # Define layout: 2 chord diagrams on top, legends below
  layout_matrix <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
  layout(layout_matrix, heights = c(6, 1))

  for (i in seq_along(object.list)) {
    condition <- names(object.list)[i]
    colors <- cell_type_colors[[condition]]

    # Extract interactions between source and target cells
    lr_network <- subsetCommunication(
      object.list[[i]],
      slot.name = "net",
      sources.use = source_cell,
      targets.use = target_cell
    )

    # Apply threshold
    lr_network <- lr_network[lr_network$prob > min_interaction_weight, ]

    if (nrow(lr_network) > 0) {
      netVisual_chord_gene(
        object.list[[i]],
        sources.use = source_cell,
        targets.use = target_cell,
        slot.name = "net",
        color.use = colors,
        title.name = paste0(source_cell, " -> ", target_cell, " - ", condition),
        legend.pos.x = -1,
        legend.pos.y = -1,
        lab.cex = 0.9,
        small.gap = 2,
        big.gap = 15,
        annotationTrackHeight = c(0.05),
        link.visible = TRUE,
        reduce = 0.02,         # Remove tiny sectors
        transparency = 0.4,
        scale = TRUE,
        thresh = 0.05,
        show.legend = FALSE
      )
    } else {
      print(paste("No interactions found between", source_cell, "and", target_cell, "in", condition))
    }
  }

  # Create a combined legend with all 6 colors
  combined_labels <- c(
    "naive_SGC IMM",
    "snc_SGC IMM",
    "naive_Macrophages",
    "snc_Macrophages",
    "naive_Neurons",
    "snc_Neurons"
  )

  combined_colors <- c(
    cell_type_colors$naive["SGC IMM"],
    cell_type_colors$snc["SGC IMM"],
    cell_type_colors$naive["Macrophages"],
    cell_type_colors$snc["Macrophages"],
    cell_type_colors$naive["Neurons"],
    cell_type_colors$snc["Neurons"]
  )

  names(combined_colors) <- combined_labels

  # Create the legend using ComplexHeatmap
  lgd <- ComplexHeatmap::Legend(
    at = names(combined_colors),
    type = "grid",
    legend_gp = grid::gpar(fill = combined_colors),
    title = "Cell Type (Condition)",
    nrow = 2
  )

  # Draw the legend in the bottom layout area
  ComplexHeatmap::draw(lgd,
                       x = unit(0.5, "npc"),
                       y = unit(0.2, "npc"),
                       just = "center")

  dev.off()
}

# Generate chord diagrams with evenly scaled interactions from SGC IMM to Macrophages
generate_chord_diagrams_lr(
  object.list = object.list,
  source_cell = "SGC IMM",
  target_cell = "Macrophages",
  output_filename = "SGC_IMM_to_Macrophages_even_chord_diagram.pdf",
  cell_type_colors = cell_type_colors
)

# Generate chord diagrams with evenly scaled interactions from Macrophages to SGC IMM
generate_chord_diagrams_lr(
  object.list = object.list,
  source_cell = "Macrophages",
  target_cell = "SGC IMM",
  output_filename = "Macrophages_to_SGC_IMM_even_chord_diagram.pdf",
  cell_type_colors = cell_type_colors
)

# Generate chord diagrams with thresholded interactions from SGC IMM to Neurons
generate_chord_diagrams_lr(
  object.list = object.list,
  source_cell = "SGC IMM",
  target_cell = "Neurons",
  output_filename = "SGC_IMM_to_Neurons_threshold_chord_diagram.pdf",
  cell_type_colors = cell_type_colors
)

# Generate chord diagrams with thresholded interactions from Neurons to SGC IMM
generate_chord_diagrams_lr(
  object.list = object.list,
  source_cell = "Neurons",
  target_cell = "SGC IMM",
  output_filename = "Neurons_to_SGC_IMM_threshold_chord_diagram.pdf",
  cell_type_colors = cell_type_colors
)

# List all unique signaling pathways in the first CellChat object
unique(cellchat.naive@netP$pathway_name)
# Define the pathway of interest
pathways.show <- c("ApoE")
# Get the maximum weight for the ApoE pathway across all CellChat objects
weight.max <- getMaxWeight(object.list, slot.name = "netP", attribute = pathways.show)

# Set up the plotting area to display two plots side by side
par(mfrow = c(1, 2), xpd = TRUE)

# Loop through each condition and generate a circle plot for the ApoE pathway
for (i in 1:length(object.list)) {
  condition.name <- names(object.list)[i]
  netVisual_aggregate(
    object = object.list[[i]],
    signaling = pathways.show,
    layout = "circle",
    edge.weight.max = weight.max[1],
    edge.width.max = 10,
    signaling.name = paste(pathways.show, condition.name)
  )
}

# Set up the plotting area to display two plots side by side
par(mfrow = c(1, 2), xpd = TRUE)

# Loop through each condition and generate a chord plot for the ApoE pathway
for (i in 1:length(object.list)) {
  condition.name <- names(object.list)[i]
  netVisual_aggregate(
    object = object.list[[i]],
    signaling = pathways.show,
    layout = "chord",
    signaling.name = paste(pathways.show, condition.name),
    small.gap = 0.5,
    big.gap = 5
  )
}

# Define cell types of interest
sources.use <- which(levels(object.list[[1]]@idents) == "SGC IMM")
targets.use <- which(levels(object.list[[1]]@idents) == "Macrophages")

# Define the signaling pathway of interest
pathways.show <- c("ApoE")

# Generate the chord diagram
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = sources.use, targets.use = targets.use, signaling = pathways.show, title.name = paste("ApoE signaling between SGC IMM and Macrophages -", names(object.list)[i]))
}
