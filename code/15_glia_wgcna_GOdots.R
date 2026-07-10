#' ======================================================================
#' 15_glia_wgcna_GOdots.R
#'
#' GO enrichment analysis and dot plots for each WGCNA module.
#'
#' Inputs:  ../r_objects/glia_hdWGCNA_object.rds
#' Outputs: ../wgcna_GO/*_GOBP_results.csv
#'
#' Author: Michael B. Thomsen, Ph.D.
#' Created: 2024-07-01
#' ======================================================================

###########################################
###########################################
###                                     ###
###    WGCNA GODOTS                     ###
###                                     ###
###########################################
###########################################

library(dplyr)
library(readr)
library(clusterProfiler)
library(org.Mm.eg.db)

# Read the CSV file
wgcna_output <- read_csv("../glia_wgcna_hub_genes.csv")

# Group by module and create a named list of genes for each module
gene_lists <- wgcna_output %>%
  group_by(module) %>%
  summarise(genes = list(gene_name)) %>%
  deframe()

# Function for performing GO analysis and writing results to a file
perform_go_analysis <- function(gene_vector, module_name) {
  # Perform GO enrichment analysis
  go_result <- enrichGO(
    gene           = gene_vector,
    OrgDb          = org.Mm.eg.db,
    keyType        = "SYMBOL",
    ont            = "BP",  # Biological Process
    pAdjustMethod  = "BH",  # Adjust p-value using Benjamini & Hochberg method
    qvalueCutoff   = 0.05   # Cutoff for significant genes
  )

  # Write the results to a CSV file
  write.csv(as.data.frame(go_result), paste0("../wgcna_GO/", module_name, "_GOBP_results.csv"))
}

# Run GO analysis for each module and save the results
lapply(names(gene_lists), function(module_name) {
  perform_go_analysis(gene_lists[[module_name]], module_name)
})


# Initialize empty dataframe
all_data <- data.frame()

modules <- c(
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

# Loop through all cluster files to read and combine them
for (module in modules) {
  file_name <- paste0("../wgcna_GO/", module, "_GOBP_results.csv")
  cluster_data <- read_csv(file_name)
  cluster_data$Cluster <- module
  all_data <- rbind(all_data, cluster_data)
}

all_data <- all_data %>%
  mutate(
    GeneNumerator = as.numeric(str_extract(GeneRatio, "^[^/]+")),
    GeneDenominator = as.numeric(str_extract(GeneRatio, "[^/]+$")),
    BgNumerator = as.numeric(str_extract(BgRatio, "^[^/]+")),
    BgDenominator = as.numeric(str_extract(BgRatio, "[^/]+$")),
    Enrichment = (GeneNumerator / GeneDenominator) / (BgNumerator / BgDenominator)
  )

# Log-transform enrichment and p-value
all_data <- all_data %>%
  mutate(
    Enrichment = log10(Enrichment),
    Log_pvalue = -log10(pvalue)
  )

significant_data <- all_data #%>% filter(Log_pvalue > 1)

ggplot(significant_data, aes(x = as.factor(Cluster), y = Description)) +
  geom_point(aes(size = Enrichment, color = Log_pvalue)) +
  scale_size_continuous(range = c(1, 10)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "GO Term Enrichment Across Cell Types",
       x = "Cluster",
       y = "GO Term",
       size = "log10(Enrichment)",
       color = "-log10(pvalue)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Define N
N <- 5  # Number of top terms you want from each module

# Filter to get top N GO terms for each module
top_data <- significant_data %>%
  group_by(Cluster) %>%
  #top_n(N, Enrichment)  # Replace 'Enrichment' with the variable you want to sort by
  top_n(N, pvalue)

# Create an interaction term to group GO terms by their respective modules
top_data$Cluster_Term <- interaction(top_data$Cluster, top_data$Description, lex.order = TRUE)

# Reverse the levels of the Cluster_Term factor
top_data$Cluster_Term <- factor(top_data$Cluster_Term, levels = rev(levels(top_data$Cluster_Term)))

# Plotting
p <- ggplot(top_data, aes(x = as.factor(Cluster), y = Cluster_Term)) +
  geom_point(aes(size = Enrichment, color = Log_pvalue)) +
  scale_size_continuous(range = c(1, 7)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Top N GO Term Enrichment Across Cell Types",
       x = "Cluster",
       y = "GO Term",
       size = "log10(Enrichment)",
       color = "-log10(pvalue)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(labels = function(x) gsub(".*\\.", "", x))  # Remove module names from y-axis labels

p


#######################################################

# Initialize empty dataframe
all_data <- data.frame()

modules <- c(
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

# Loop through all cluster files to read and combine them
for (module in modules) {
  file_name <- paste0("../wgcna_GO/", module, "_GOBP_results.csv")
  cluster_data <- read_csv(file_name)
  cluster_data$Cluster <- module
  all_data <- rbind(all_data, cluster_data)
}

# Data transformations with log calculations
all_data <- all_data %>%
  mutate(
    GeneNumerator = as.numeric(str_extract(GeneRatio, "^[^/]+")),
    GeneDenominator = as.numeric(str_extract(GeneRatio, "[^/]+$")),
    BgNumerator = as.numeric(str_extract(BgRatio, "^[^/]+")),
    BgDenominator = as.numeric(str_extract(BgRatio, "[^/]+$")),
    Enrichment = log10((GeneNumerator / GeneDenominator) / (BgNumerator / BgDenominator)),
    Log_pvalue = -log10(pvalue)
  )

# Define N
N <- 5  # Number of top terms you want from each module

# Filter to get top N GO terms for each module
top_data <- all_data %>%
  group_by(Cluster) %>%
  top_n(N, wt = Log_pvalue) %>%
  ungroup() %>%
  arrange(Cluster, desc(Log_pvalue))  # Ensure data is ordered by Cluster and by top p-values

# Convert Cluster to a factor with the correct levels
top_data$Cluster <- factor(top_data$Cluster, levels = modules)

# Create a unique identifier for each term by combining Cluster and Description
top_data$Cluster_Term <- interaction(top_data$Cluster, top_data$Description, lex.order = TRUE)

# Plotting
p <- ggplot(top_data, aes(x = Cluster, y = fct_inorder(Cluster_Term))) +
  geom_point(aes(size = Enrichment, color = Log_pvalue)) +
  scale_size_continuous(range = c(1, 10)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Top N GO Term Enrichment Across Cell Types",
       x = "Cluster",
       y = "GO Term",
       size = "log10(Enrichment)",
       color = "-log10(pvalue)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(labels = function(x) gsub(".*\\.", "", x))  # Remove module names from y-axis labels

print(p)

ggsave("../plots/satglia_top_WGCNA_GO_terms.pdf", plot = p, device = "pdf", width = 14, height = 12, dpi = 300)
ggsave("../plots/satglia_top_WGCNA_GO_terms.png", plot = p, device = "png", width = 14, height = 12, dpi = 300)

####

library(tidyverse)
library(tidyr)

# Initialize empty dataframe
all_data <- data.frame()

modules <- c(
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

# Loop through all cluster files to read and combine them
for (module in modules) {
  file_name <- paste0("../wgcna_GO/", module, "_GOBP_results.csv")
  cluster_data <- read_csv(file_name)
  cluster_data$Cluster <- module
  all_data <- rbind(all_data, cluster_data)
}

# Ensure 'Description' is a character vector
all_data$Description <- as.character(all_data$Description)

# Data transformations with log calculations
all_data <- all_data %>%
  mutate(
    GeneNumerator = as.numeric(str_extract(GeneRatio, "^[^/]+")),
    GeneDenominator = as.numeric(str_extract(GeneRatio, "[^/]+$")),
    BgNumerator = as.numeric(str_extract(BgRatio, "^[^/]+")),
    BgDenominator = as.numeric(str_extract(BgRatio, "[^/]+$")),
    Enrichment = log10((GeneNumerator / GeneDenominator) / (BgNumerator / BgDenominator)),
    Log_pvalue = -log10(pvalue)
  )

# Define N
N <- 5  # Number of top terms you want from each module

# Get the top N GO terms for each module
top_terms_per_module <- all_data %>%
  group_by(Cluster) %>%
  top_n(N, wt = Log_pvalue) %>%
  ungroup()

# Get the unique GO terms from the top N terms across all modules
unique_terms <- unique(top_terms_per_module$Description)

# Ensure 'Description' is a character vector in 'unique_terms'
unique_terms <- as.character(unique_terms)

# Filter all_data to include only these unique terms across all modules
filtered_data <- all_data %>%
  filter(Description %in% unique_terms)

# Ensure 'Description' is a character vector in 'filtered_data'
filtered_data$Description <- as.character(filtered_data$Description)

# Create a data frame with all combinations of Clusters and unique terms
complete_data <- expand_grid(Cluster = modules, Description = unique_terms)

# Ensure 'Description' is a character vector in 'complete_data'
complete_data$Description <- as.character(complete_data$Description)

# Left join to add the data from filtered_data
complete_data <- left_join(complete_data, filtered_data, by = c("Cluster", "Description"))

# Replace NAs in Enrichment and Log_pvalue with zeros before further processing
complete_data$Enrichment[is.na(complete_data$Enrichment)] <- 0
complete_data$Log_pvalue[is.na(complete_data$Log_pvalue)] <- 0

# **Filter out rows where Enrichment or Log_pvalue is zero**
complete_data <- complete_data %>%
  filter(Enrichment != 0, Log_pvalue != 0)

# For ordering the terms on y-axis, order them based on the module where they have the highest Log_pvalue
term_module <- complete_data %>%
  group_by(Description) %>%
  arrange(desc(Log_pvalue)) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::select(Description, MaxCluster = Cluster, MaxLogP = Log_pvalue)

# Convert MaxCluster to factor with levels as modules
term_module$MaxCluster <- factor(term_module$MaxCluster, levels = modules)

# Order terms based on MaxCluster (module order) and MaxLogP
term_order <- term_module %>%
  arrange(MaxCluster, desc(MaxLogP))

# Set the levels of Description based on term_order
complete_data$Description <- factor(complete_data$Description, levels = term_order$Description)

# Convert Cluster to a factor with the correct levels
complete_data$Cluster <- factor(complete_data$Cluster, levels = modules)

# Plotting
p <- ggplot(complete_data, aes(x = Cluster, y = Description)) +
  geom_point(aes(size = Enrichment, color = Log_pvalue)) +
  scale_size_continuous(range = c(1, 10)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Top N GO Term Enrichment Across Cell Types",
       x = "Cluster",
       y = "GO Term",
       size = "log10(Enrichment)",
       color = "-log10(pvalue)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)


ggsave("../plots/glia_top_WGCNA_GO_terms.pdf", plot = p, device = "pdf", width = 11.5, height = 12, dpi = 300)
ggsave("../plots/glia_top_WGCNA_GO_terms.png", plot = p, device = "png", width = 11.5, height = 12, dpi = 300)

###########

library(tidyverse)
library(tidyr)

# Initialize empty dataframe
all_data <- data.frame()

modules <- c(
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

# Loop through all cluster files to read and combine them
for (module in modules) {
  file_name <- paste0("../wgcna_GO/", module, "_GOBP_results.csv")
  cluster_data <- read_csv(file_name)
  cluster_data$Cluster <- module
  all_data <- rbind(all_data, cluster_data)
}

# Ensure 'Description' is a character vector
all_data$Description <- as.character(all_data$Description)

# Data transformations with log calculations
all_data <- all_data %>%
  mutate(
    GeneNumerator = as.numeric(str_extract(GeneRatio, "^[^/]+")),
    GeneDenominator = as.numeric(str_extract(GeneRatio, "[^/]+$")),
    BgNumerator = as.numeric(str_extract(BgRatio, "^[^/]+")),
    BgDenominator = as.numeric(str_extract(BgRatio, "[^/]+$")),
    Enrichment = log10((GeneNumerator / GeneDenominator) / (BgNumerator / BgDenominator)),
    Log_pvalue = -log10(pvalue)
  )

# Define N
N <- 5  # Number of top terms you want from each module

# Get the top N GO terms for each module
top_terms_per_module <- all_data %>%
  group_by(Cluster) %>%
  top_n(N, wt = Log_pvalue) %>%
  ungroup()

# Identify duplicates (GO terms appearing in top N for multiple modules)
duplicate_terms <- top_terms_per_module %>%
  group_by(Description) %>%
  filter(n() > 1) %>%
  arrange(Description)

if(nrow(duplicate_terms) > 0) {
  cat("The following GO terms appear in the top N for multiple modules:\n")
  print(duplicate_terms %>% dplyr::select(Description, Cluster, Log_pvalue))
}

# For GO terms that are in multiple modules, select the module where they have the higher Log_pvalue
top_terms_per_module_unique <- top_terms_per_module %>%
  group_by(Description) %>%
  arrange(desc(Log_pvalue)) %>%
  dplyr::slice(1) %>%
  ungroup()

# Now, we have unique GO terms assigned to modules

# Collect all the unique GO terms
unique_terms <- top_terms_per_module_unique$Description

# Ensure 'Description' is a character vector in 'unique_terms'
unique_terms <- as.character(unique_terms)

# Filter all_data to include only these unique terms across all modules
filtered_data <- all_data %>%
  filter(Description %in% unique_terms)

# Ensure 'Description' is a character vector in 'filtered_data'
filtered_data$Description <- as.character(filtered_data$Description)

# Create a data frame with all combinations of Clusters and unique terms
complete_data <- expand_grid(Cluster = modules, Description = unique_terms)

# Ensure 'Description' is a character vector in 'complete_data'
complete_data$Description <- as.character(complete_data$Description)

# Left join to add the data from filtered_data
complete_data <- left_join(complete_data, filtered_data, by = c("Cluster", "Description"))

# Replace NAs in Enrichment and Log_pvalue with zeros before further processing
complete_data$Enrichment[is.na(complete_data$Enrichment)] <- 0
complete_data$Log_pvalue[is.na(complete_data$Log_pvalue)] <- 0

# Assign each GO term to its module (from top_terms_per_module_unique)
# We can create a mapping of Description to AssignedModule
term_module_assignment <- top_terms_per_module_unique %>%
  dplyr::select(Description, AssignedModule = Cluster)

# Merge this assignment back into complete_data
complete_data <- left_join(complete_data, term_module_assignment, by = "Description")

# Order the GO terms: group by AssignedModule, then by Log_pvalue within each group
term_order <- top_terms_per_module_unique %>%
  arrange(factor(Cluster, levels = modules), desc(Log_pvalue)) %>%
  dplyr::select(Description)

# Set the levels of Description based on term_order
complete_data$Description <- factor(complete_data$Description, levels = term_order$Description)

# Convert Cluster to a factor with the correct levels
complete_data$Cluster <- factor(complete_data$Cluster, levels = modules)

# Optional: Filter out rows where Enrichment and Log_pvalue are zero
# This will remove points where the GO term is not significant in a module
complete_data <- complete_data %>%
  filter(Enrichment != 0 & Log_pvalue != 0)

# Plotting
p <- ggplot(complete_data, aes(x = Cluster, y = Description)) +
  geom_point(aes(size = Enrichment, color = Log_pvalue)) +
  scale_size_continuous(range = c(1, 10)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Top N GO Term Enrichment Across Cell Types",
       x = "Cluster",
       y = "GO Term",
       size = "log10(Enrichment)",
       color = "-log10(pvalue)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)

ggsave("../plots/glia_top_WGCNA_GO_terms_ordered_v2.pdf", plot = p, device = "pdf", width = 11.5, height = 12, dpi = 300)
ggsave("../plots/glia_top_WGCNA_GO_terms_ordered_v2.png", plot = p, device = "png", width = 11.5, height = 12, dpi = 300)
