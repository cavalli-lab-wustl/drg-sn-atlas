#' ======================================================================
#' 14_glia_wgcna.R
#'
#' hdWGCNA on glial cells: construct metacells, build network, identify 15 modules.
#'
#' Inputs:  ../r_objects/DRGSNMI_glia.RDS
#' Outputs: ../r_objects/glia_hdWGCNA_object.rds, ../wgcna_GO/glia_wgcna_hub_genes.csv
#'
#' Author: Michael B. Thomsen, Ph.D.
#' Created: 2024-07-01
#' ======================================================================

###########################################
###########################################
###                                     ###
###             GLIA WGCNA              ###
###                                     ###
###########################################
###########################################

# doParallel and doSNOW have problems loading from custom library directory

# be sure to load v5 Seurat and hdWGCNA from system libraries (IF CUSTOM LIBARIES LOADED RESTART SESSION)

glia <- readRDS("../r_objects/DRGSNMI_glia.RDS")

# set the seurat_clusters to the new Idents
glia$seurat_clusters <- Idents(glia)

# set random seed for reproducibility
set.seed(12345)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

DimPlot(glia, group.by='seurat_clusters', label=TRUE) +
   umap_theme() + ggtitle('Satellite Glia Subclusters') + NoLegend()


glia <- SetupForWGCNA(
  glia,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.03, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "glia" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
glia <- MetacellsByGroups(
  seurat_obj = glia,
  group.by = c("seurat_clusters"), # , "study" removed | specify the columns in seurat_obj@meta.data to group by
  reduction = 'integrated.dr', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'seurat_clusters' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
glia <- NormalizeMetacells(glia)

glia <- SetDatExpr(
  glia,
  group_name = c(
  "SGC",
  "SGC PRE",
  "SGC IMM",
  "nmSC",
  "nmSC PRE",
  "nmSC REP",
  "mSC",
  "mSC PRE",
  "mSC REP",
  "PROG I",
  "PROG II"
  ), # the name of the group of interest in the group.by column
  group.by='seurat_clusters', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data'
)

# Test different soft powers:
glia <- TestSoftPowers(
  glia,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(glia)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

# construct co-expression network:
glia <- ConstructNetwork(
  glia,
  #soft_power = 10, # the soft power to use, if not specified automatically selected
  tom_name = 'glia', # name of the topoligical overlap matrix written to disk
  overwrite_tom = TRUE
  )

PlotDendrogram(glia, main='Glia hdWGCNA Dendrogram')

TOM <- GetTOM(glia)

# need to run ScaleData first or else harmony throws an error:
glia <- ScaleData(glia, features=VariableFeatures(glia))

# compute all MEs in the full single-cell dataset
glia <- ModuleEigengenes(
 glia,
 group.by.vars="study"
)

# harmonized module eigengenes:
hMEs <- GetMEs(glia)

# compute eigengene-based connectivity (kME):
glia <- ModuleConnectivity(
  glia,
  group_name = c(
  "SGC",
  "SGC PRE",
  "SGC IMM",
  "nmSC",
  "nmSC PRE",
  "nmSC REP",
  "mSC",
  "mSC PRE",
  "mSC REP",
  "PROG I",
  "PROG II"
  ),
  group.by='seurat_clusters',
)


modules <- GetModules(glia)

# show the first 6 columns:
head(modules[,1:6])

# get hub genes
hub_df <- GetHubGenes(glia, n_hubs = Inf)

head(hub_df)

write.csv(hub_df, file='../wgcna_GO/glia_wgcna_hub_genes.csv')

saveRDS(glia, file='../r_objects/glia_hdWGCNA_object.rds')

glia <- readRDS('../r_objects/glia_hdWGCNA_object.rds')

# compute gene scoring for the top 30 hub genes by kME for each module
# with Seurat method
glia <- ModuleExprScore(
  glia,
  n_genes = 30,
  method='Seurat'
)

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  glia,
  features='hMEs', # plot the hMEs
  #order=TRUE # order so the points with highest hMEs are on top
)

plot_list_adjusted <- lapply(plot_list, function(p) {
  p + theme(
    aspect.ratio = 0.95,
    legend.background = element_rect(fill = "transparent", colour = "transparent"),
    legend.key = element_rect(fill = "transparent", colour = "transparent")
  )
})

for(i in seq_along(plot_list_adjusted)) {
  ggsave(
    filename = paste0("../plots/glia_hMEs_", i, ".png"),
    plot = plot_list_adjusted[[i]],
    device = 'png',
    width = 2,
    height = 2,
    dpi = 300,
    bg = "transparent"
  )
}

# Stitch together with patchwork
p <- wrap_plots(plot_list_adjusted, ncol=3)

# Save the plot
ggsave('../plots/glia_hMEs.pdf', plot=p, device='pdf', width=5, height=7, dpi=300)


# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  glia,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = FALSE # depending on Seurat vs UCell for gene scoring
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=6)

pdf(file = "../plots/glia_wgcna_correlations.pdf", width = 10, height = 10)
ModuleCorrelogram(glia, method = "square", outline = T, addgrid.col = "darkgray", order = "hclust", addrect = 4, mar = c(4,0,4,0), col = colorRampPalette(c("midnightblue","white","darkred"))(100))
dev.off()

# Open a pdf device to save the plot
pdf(file = "../plots/glia_eigengene_adjacency_heatmap.pdf", width = 10, height = 10)

# Create the correlogram with the desired style:
ModuleCorrelogram(
  glia,
  features = "hMEs",                         # Plot module eigengenes
  type = "full",                              # Show the full matrix (mirrored across the diagonal)
  method = "color",                           # Use a heatmap (color fill) rather than size/shape encoding
  order = "original",                         # Keep the original order (or use "hclust" for clustering)
  col = colorRampPalette(c("red", "white", "green"))(200),  # Set the palette: -1 = red, 0 = white, 1 = green
  tl.col = "black",                           # Text label color
  tl.srt = 45                                 # Text label rotation
)

# Close the pdf device
dev.off()

# get hMEs from seurat object
MEs <- GetMEs(glia, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
glia@meta.data <- cbind(glia@meta.data, MEs)

mods <- c(
  "brown",
  "yellow",
  "blue",
  "magenta",
  "black",
  "red",
  "salmon",
  "purple",
  "pink",
  "green",
  "greenyellow",
  "cyan",
  "turquoise",
  "tan"
)

# plot with Seurat's DotPlot function
p <- DotPlot(glia, features=mods, dot.scale = 7)

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  #coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
p

ggsave('../plots/glia_wgcna_dotplot.pdf', plot=p, device='pdf', width=6.5, height=4.5, dpi=300)
ggsave('../plots/glia_wgcna_dotplot.png', plot=p, device='png', width=6.5, height=4.5, dpi=300)

pdf(file = "../plots/glia_wgcna_correlations.pdf", width = 10, height = 10)
ModuleCorrelogram(
  glia,
  method = "square",
  outline = T,
  addgrid.col = "darkgray",
  order = "hclust",
  addrect = 4,
  mar = c(4,0,4,0),
  col = colorRampPalette(c("midnightblue","white","darkred"))(100)
  )
dev.off()

###########

# Load necessary libraries
library(Seurat)
library(hdWGCNA)
library(corrplot)

# Define all module eigengenes based on known module colors
moduleEigengenes <- c("pink", "yellow", "green", "greenyellow", "brown",
                      "grey", "blue", "black", "tan", "turquoise",
                      "red", "magenta", "purple")

# Define modules to exclude
excludeModules <- c("grey", "tan", "pink")  # Include all background modules

# Determine modules to include
modulesToInclude <- setdiff(moduleEigengenes, excludeModules)

# Verify the modules to include
print(modulesToInclude)
# [1] "yellow"      "green"       "greenyellow" "brown"       "blue"
# [6] "turquoise"   "red"         "magenta"     "purple"

# Extract the included module eigengenes from the Seurat object's metadata
MEs_filtered <- glia@meta.data[, modulesToInclude]

# Verify the structure
str(MEs_filtered)

# Ensure there are no duplicate column names
if(any(duplicated(colnames(MEs_filtered)))){
  stop("Duplicate module eigengene names found. Please ensure all module names are unique.")
}

# Ensure there are no missing values
if(any(is.na(MEs_filtered))){
  stop("Missing values detected in MEs_filtered. Please handle NA values before proceeding.")
}

pdf(file = "../plots/glia_wgcna_correlations_trimmed.pdf", width = 10, height = 10)

# Attempt to generate the Module Correlogram with adjusted parameters
tryCatch({
  ModuleCorrelogram(
    seurat_obj = glia,
    MEs2 = MEs_filtered,                  # Pass the filtered module eigengenes
    features = "MEs",                      # Use "MEs" instead of "hMEs"
    method = "square",                     # Correlation visualization method
    outline = TRUE,                        # Whether to draw outlines
    addgrid.col = "darkgray",              # Grid color
    order = "hclust",                      # Ordering method
    addrect = 4,                           # Number of rectangles in dendrogram
    mar = c(4, 0, 4, 0),                   # Plot margins
    col = colorRampPalette(c("midnightblue", "white", "darkred"))(100), # Color palette
    exclude_grey = FALSE,                  # Already excluded "grey" and "black"
    type = "upper",                        # Plot upper triangle
    tl.col = "black",                      # Text label color
    tl.srt = 45,                           # Text label rotation
    sig.level = c(1e-04, 0.001, 0.01, 0.05),# Significance levels
    pch.cex = 0.7,                         # Plot character expansion
    ncolors = 200                          # Number of colors in palette
  )
}, error = function(e){
  message("ModuleCorrelogram failed with error: ", e$message)
  message("Attempting to plot correlation matrix manually.")

  # Manually compute and plot the correlation matrix
  corr_matrix <- cor(MEs_filtered, method = "pearson")

  # Check if the correlation matrix is square
  if(nrow(corr_matrix) != ncol(corr_matrix)){
    stop("Correlation matrix is not square. Cannot proceed with plotting.")
  }

  # Plot using corrplot
  corrplot(
    corr_matrix,
    method = "square",
    type = "upper",
    order = "hclust",
    addrect = 4,
    col = colorRampPalette(c("midnightblue", "white", "darkred"))(100),
    tl.col = "black",
    tl.srt = 45,
    sig.level = c(1e-04, 0.001, 0.01, 0.05),
    pch.cex = 0.7
  )
})
dev.off()


################
set.seed(123)  # Fix the random seed

ModuleNetworkPlot(
  glia,
  outdir = '../plots/modules'
)

# hubgene network
HubGeneNetworkPlot(
  glia,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)


# get the list of modules:
modules <- GetModules(glia)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# hubgene network
HubGeneNetworkPlot(
  glia,
  n_hubs = 10, n_other=20,
  edge_prop = 0.75,
  return_graph = TRUE,
  mods = mods[1:5] # only select 5 modules
)

pdf(file = "../plots/combined_hub_network.pdf", width = 10, height = 10)
HubGeneNetworkPlot(
  glia,
  n_hubs = 3, n_other = 5,
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()

# Plot 1: Modules brown, yellow, blue, magenta
pdf(file = "../plots/combined_hub_network_group1.pdf", width = 10, height = 10)
HubGeneNetworkPlot(
  glia,
  n_hubs = 5,
  n_other = 5,
  edge_prop = 1,
  vertex.label.cex = 1,
  mods = c("brown", "yellow", "blue", "magenta")
)
dev.off()

# Plot 2: Modules black, red, salmon
pdf(file = "../plots/combined_hub_network_group2.pdf", width = 10, height = 10)
HubGeneNetworkPlot(
  glia,
  n_hubs = 5,
  n_other = 5,
  edge_prop = 1,
  vertex.label.cex = 1,
  mods = c("black", "red", "salmon")
)
dev.off()

# Plot 3: Modules purple, pink, green
pdf(file = "../plots/combined_hub_network_group3.pdf", width = 10, height = 10)
HubGeneNetworkPlot(
  glia,
  n_hubs = 5,
  n_other = 5,
  edge_prop = 0.75,
  mods = c("purple", "pink", "green")
)
dev.off()


#################

library(igraph)

# Generate the network graph for group 1 (brown, yellow, blue, magenta)
g <- HubGeneNetworkPlot(
  glia,
  n_hubs = 3,
  n_other = 5,
  edge_prop = 0.75,
  mods = c("brown", "yellow", "blue", "magenta"),
  return_graph = TRUE
)

# Remove vertices that are not connected (degree == 0)
isolated_nodes <- which(degree(g) == 0)
if(length(isolated_nodes) > 0){
  g <- delete.vertices(g, isolated_nodes)
}

# Optionally, you can re-calculate a layout if needed:
l <- layout_with_fr(g)

# Plot the modified graph with an increased font size for the labels
pdf(file = "../plots/combined_hub_network_group1_modified.pdf", width = 10, height = 10)
plot(g,
     layout = l,
     vertex.label.cex = 1.5,  # Increase font size of gene labels
     edge.arrow.mode = 0      # (Optional) Remove arrow heads for undirected networks
)
dev.off()
