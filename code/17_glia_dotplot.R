#' ======================================================================
#' 17_glia_dotplot.R
#'
#' Marker gene DotPlots for SGC subtypes, progenitors, and ErbB/NRG signaling.
#'
#' Inputs:  ../r_objects/DRGSNMI_glia.RDS
#' Outputs: ../plots/glia_dotplot_*.pdf
#'
#' Author: Michael B. Thomsen, Ph.D.
#' Created: 2024-09-10
#' ======================================================================

library(Seurat)
library(ggplot2)
library(pheatmap)

###########################################
###########################################
###                                     ###
###    DATA LOADING AND PREP            ###
###                                     ###
###########################################
###########################################

# Load glia object
glia <- readRDS("../r_objects/DRGSNMI_glia.RDS")

glia <- subset(glia, idents = c("SGC I", "SGC II", "SGC III", "SGC IEG", "SGC IMM"))

# This marker gene list is for the full size marker gene DotPlot of glia in the drgmi atlas
markergenes <- c(
  "Plp1",
  "Sox10",
  "Fabp7",
  "Ptprz1",
  "Pdpn",
  "Bcan",
  "Pou3f1",
  "Nkd1",
  "Cdkn1c",
  "Cxcl14",
  "Ngfr",
  "Hspa1a",
  "Fos",
  "Ifit1",
  "Ifit3",
  "Ifit3b",
  "Gbp2",
  "Stat1",
  "Cd274",
  "Irf7"
)

# Create the dot plot
dot_plot_collapsed <- DotPlot(glia, features = markergenes, dot.scale = 9, dot.min = 0.05) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x.top = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
  ) +
  labs(x = NULL, y = NULL) +
  guides(
    color = guide_colorbar(title = "Avg. Exp.", label = FALSE),
    size = guide_legend(title = "% Exp.")
  )

ggsave("../plots/sgc_marker_dotplot.pdf", plot = dot_plot_collapsed, width = 8.5, height = 3.5, dpi = 300)
ggsave("../plots/sgc_marker_dotplot.png", plot = dot_plot_collapsed, width = 8.5, height = 3.5, dpi = 300)

# Examining progenitor markers
markers <- c(
  "Sox2",    # Stemness and progenitor maintenance
  "Sox9",    # Astrocyte and oligodendrocyte lineage specification
  "Nes",  # Neural progenitors, including glial progenitors
  "Vim",# Immature glial progenitors
  "Aldh1l1", # Astrocyte progenitors
  "Fabp7",   # Radial glia and glial progenitors
  "Olig2",   # Oligodendrocyte progenitor marker
  "Pdgfra",  # OPC marker in CNS
  "Cspg4",     # OPC marker
  "Ascl1",   # Glial progenitor specification, especially oligodendrocytes
  "Cdh2",    # N-cadherin, progenitor cells in CNS
  "Gfap",    # Radial glia and astrocyte precursors
  "Sox10",   # PNS glial progenitors, Schwann cells
  "Egr2",  # Schwann cell progenitors to myelination transition
  "Pax3",    # Neural crest and Schwann cell precursors
  "Erbb3",   # Schwann cell progenitors
  "Pou3f1",    # Schwann cell progenitors
  "Mpz",     # Myelination marker in later Schwann cell progenitors
  "Dhh",     # Schwann cell progenitors
  "Fgfr3",   # Astrocyte progenitors
  "Id2",     # Inhibitors of differentiation, progenitor state
  "Id4",     # Inhibitors of differentiation, progenitor state
  "Cxcr4",   # Progenitor migration/proliferation
  "Hes1",    # Notch signaling, maintaining progenitor state
  "Efnb2",   # Cell signaling in progenitors
  "Mbp",     # Differentiating oligodendrocytes
  "Plp1",    # Oligodendrocyte maturation
  "S100b",   # Astrocyte and Schwann cell marker
  "Cnp",     # Oligodendrocyte and Schwann cell maturation
  "Nr2e1",     # Associated with glial stem cells and progenitor populations
  "Lif",     # Leukemia inhibitory factor, involved in progenitor plasticity
  "Tnfrsf19",# Associated with neural progenitor and plasticity
  "Bmp4",    # Promotes glial differentiation from progenitors
  "Notch1",  # Notch signaling important for maintaining glial progenitors
  "Prom1",   # Stemness marker in various progenitor cells
  "Gli1",    # Hedgehog signaling, involved in glial progenitor plasticity
  "Fgf2",    # Fibroblast growth factor 2, promotes proliferation of progenitors
  "Tcf7l2",  # Wnt signaling associated with glial progenitor state
  "Zbtb16",  # Plasticity-related gene in progenitors
  "Egr1",    # Associated with differentiation potential in progenitor cells
  "Stat3",    # Involved in signaling for glial progenitor maintenance
  "Sox4",
  "Sox8",
  "Foxd3",
  "Ednrb",
  "Mycn",
  "Hmga2",
  "Cux1"
)

# Refined progenitor/differentiation markers for SGC subtypes
markers <- c(
  "Fabp7",
  "Pdpn",
  "Ptprz1",
  "Sox10",
  "Sox2",
  "Mki67",
  "Top2a",
  "Ngfr",
  "Ncam1",
  "L1cam",
  "Pou3f1",
  "Cdkn1c",
  "Mpz",
  "Mbp",
  "Ncmap"
  )

# ErbB/NRG signaling pathway genes
  markers <- c(
    "Nrg1",
    "Nrg2",
    "Nrg3",
    "Nrg4",
    "Erbb2",
    "Erbb3",
    "Erbb4",
    "Egf",
    "Pik3ca",
    "Pik3cb",
    "Akt1",
    "Akt2",
    "Mapk1",
    "Mapk3",
    "Grb2",
    "Sos1",
    "Raf1",
    "Shc1",
    "Dock7",
    "Sos2",
    "Crk"
  )

# Create the dot plot
dot_plot <- DotPlot(glia, features = markers, dot.scale = 9, dot.min = 0.05) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x.top = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
  ) +
  labs(x = NULL, y = NULL) +
  guides(
    color = guide_colorbar(title = "Avg. Exp.", label = FALSE),
    size = guide_legend(title = "% Exp.")
)

print(dot_plot)

p <- FeaturePlot(glia, features = "Gfap", split.by = "Study", combine = F)
plot_grid(plotlist = p, ncol=4)

PlotGene <- function(obj, gene, output_directory = "../plots/") {
  p <- FeaturePlot(obj, features = gene, order = T, raster = F)
  p <- p + theme(
    plot.title = element_text(size = 20, face = "italic", hjust = 0.5),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 10, hjust = 0),
    axis.title.y = element_text(size = 10, vjust = 0, hjust = 0)
  )
  print(p)
  ggplot2::ggsave(paste0(output_directory, "glia_", gene, ".png"), plot = p, device = "png", width = 5, height = 5, dpi = 300)
}

# Define genes and conditions for condition-specific heatmap
genes_to_plot <- c(
  #SNC
  "Ctss",
  "H2-Q7",
  "C1qc",
  "Cd74",
  "C1qb",
  "BX470216.1",
  "Clba1",
  "Nanp",
  "Poglut2",
  "Ess2",
  "C1qa",
  "Lyz2",
  "H2-Eb1",
  "Baalc",
  "Rmc1",
  "Znrd2",
  "Ccn1",
  "Cfap410",
  "Crlf1",
  "Inka1",
  "Atpsckmt",
  "H2-Aa",
  "Ramac",

  #DRC
  "Gfap",
  "Timp1",
  "Crlf1",
  "Serpina3n",
  "Spp1",
  "Gm49980",
  "C4b",
  "Lyz2",
  "Hist2h2aa1",
  "Ciao2a",
  "Trf",
  "Gatd3a",
  "Cfap410",
  "Vps35l",
  "Fads2",
  "Ccnq",
  "Entr1",

  #AIH (duplicates with other conditions commented out)
  "Yju2",
  "Zfp983",
  "Pip5k1b",
  "I730030J21Rik",
  "Mtln",
  "Gm43305",
  "Dele1",
  "Ciao2b",
  "Pip4p1",

  #SCI (duplicates with other conditions commented out)
  "Aldoc",
  "Prxl2b",
  "Pip4p2",
  "Pip4p1",
  "Micos13",
  "Rab5if",
  "Cfap298"
)

conditions_of_interest <- c("cavalli-aih", "cavalli-snc", "cavalli-sci", "cavalli-drc", "cavalli-contra", "cavalli-naive")

# Subset the Seurat object to the cluster and conditions of interest
cluster_of_interest <- "SGC II"
subset_data <- subset(glia, idents = cluster_of_interest)

# Extract average expression for RNA assay, grouped by Study (conditions)
avg_expression <- AverageExpression(subset_data, features = genes_to_plot, group.by = "Study", assay = "RNA")

# Filter for conditions of interest
avg_expression_subset <- avg_expression$RNA[, conditions_of_interest]

# Manually specify the order of columns (conditions)
column_order <- c("cavalli-aih", "cavalli-snc", "cavalli-sci", "cavalli-drc", "cavalli-contra", "cavalli-naive")

# Ensure the data is in the custom column order
avg_expression_subset <- avg_expression_subset[, column_order]

# Generate Heatmap with row clustering and custom column order
pheatmap(avg_expression_subset,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",  # z-score normalization per gene
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Expression of Genes in Specified Cluster and Conditions")

# Generate Heatmap for all conditions
avg_expression_all <- avg_expression$RNA  # Includes all conditions

pheatmap(avg_expression_all, cluster_rows = TRUE, cluster_cols = TRUE,
         scale = "row",  # z-score normalization per gene
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Expression of Genes in Specified Cluster Across All Conditions")
