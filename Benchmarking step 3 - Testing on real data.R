# ------------------------------------------------------------------------------
# Author: Feargal Ryan
# GitHub: https://github.com/feargalr/
# Description: Script for downloading datasets from the curatedMetagenomicData 
#              package and performing differential expression and enrichment 
#              analysis using LINDA, FastANCOM, ALDEx2, and TaxSEA.
# ------------------------------------------------------------------------------

# Load required libraries for the analysis
library(ggplot2)                     # For data visualization
library(curatedMetagenomicData)       # To access curated metagenomics datasets
library(TaxSEA)                       # For taxon set enrichment analysis
library(MicrobiomeStat)               # For LINDA differential abundance testing
library(fastANCOM)                    # Fast implementation of ANCOM for DA testing
library(ALDEx2)                       # ALDEx2 for DA analysis with compositional data
library(ANCOMBC)                      # ANCOM-BC for differential abundance analysis
library(gridExtra)                    # For arranging multiple plots
library(TreeSummarizedExperiment)     # For tree-based data structures
library(ggpubr)                       # To create publication-ready plots
library(ggrepel)                      # For labeled scatterplots

# Set random seed for reproducibility
set.seed(123)

# Load metadata 
allmeta_data <- sampleMetadata

# -------------------------------------------------------------------------
# HMP_2019_ibdmdb Dataset Processing
# -------------------------------------------------------------------------

# Download the HMP 2019 IBD dataset
HMP_data <- curatedMetagenomicData(
  pattern = "2021-10-14.HMP_2019_ibdmdb.relative_abundance", 
  counts = TRUE, 
  dryrun = FALSE
)

# Extract the dataset and filter by abundance
HMP_counts <- assay(HMP_data$`2021-10-14.HMP_2019_ibdmdb.relative_abundance`)
HMP_counts <- HMP_counts[apply(HMP_counts > 1000, 1, sum) > 10, ]

# Match metadata with count data
HMP_metadata <- allmeta_data[allmeta_data$sample_id %in% colnames(HMP_counts), ]
rownames(HMP_metadata) <- HMP_metadata$sample_id

# Remove samples with antibiotic use
HMP_metadata <- HMP_metadata[HMP_metadata$antibiotics_current_use == "no", ]
valid_samples <- intersect(rownames(HMP_metadata), colnames(HMP_counts))
HMP_counts <- HMP_counts[, valid_samples]
HMP_metadata <- HMP_metadata[valid_samples, ]

# Handle missing disease subtypes and convert to factors
HMP_metadata$disease_subtype[is.na(HMP_metadata$disease_subtype)] <- "Control"
HMP_metadata$disease_subtype <- factor(HMP_metadata$disease_subtype, levels = c("Control", "UC", "CD"))

# Remove duplicate samples from the same participant
HMP_metadata <- HMP_metadata[order(HMP_metadata$visit_number, decreasing = TRUE), ]
HMP_metadata <- HMP_metadata[!duplicated(HMP_metadata$subject_id), ]
HMP_counts <- HMP_counts[, rownames(HMP_metadata)]

# Extract and clean species names from row names
species_names <- sapply(rownames(HMP_counts), function(x) strsplit(x, "\\|")[[1]][7])
rownames(HMP_counts) <- gsub("s__", "", species_names)

# -------------------------------------------------------------------------
# Differential Abundance Analysis Using LINDA, FastANCOM, and ALDEx2
# -------------------------------------------------------------------------

# Run LINDA on the HMP dataset
linda_res <- linda(HMP_counts, HMP_metadata, formula = '~study_condition')
linda_results <- linda_res$output$study_conditionIBD
write.csv(linda_results, "HMP_2019_ibdmdb_linda_results.csv")
table(linda_results$padj < 0.1)

# Extract LINDA ranks
linda_ranks <- linda_results$log2FoldChange
names(linda_ranks) <- rownames(linda_results)

# Run FastANCOM on the same dataset
fastancom_fit <- fastANCOM(Y = t(HMP_counts), x = HMP_metadata$study_condition)
fastancom_results <- fastancom_fit$results$final
table(fastancom_results$log2FC.qval < 0.1)
fastancom_ranks <- fastancom_results$log2FC
names(fastancom_ranks) = rownames(fastancom_results)

# Run ALDEx2
aldex_results <- aldex(HMP_counts, HMP_metadata$study_condition)
table(aldex_results$wi.eBH < 0.1)
aldex2_ranks <- aldex_results$effect
names(aldex2_ranks) = rownames(aldex_results)

# -------------------------------------------------------------------------
# Taxon Set Enrichment Analysis Using TaxSEA
# -------------------------------------------------------------------------

# Filter taxa that are present across all methods
common_taxa <- intersect(names(linda_ranks), names(fastancom_ranks))
common_taxa <- intersect(common_taxa, names(aldex2_ranks))

# Run TaxSEA for LINDA
linda_taxsea_results <- TaxSEA(taxon_ranks = linda_ranks)
linda_metabo_results <- linda_taxsea_results$Metabolite_producers
write.csv(linda_metabo_results, "HMP_2019_ibdmdb_TaxSEA_results.csv")

# Run TaxSEA for FastANCOM
fastancom_taxsea_results <- TaxSEA(taxon_ranks = fastancom_ranks)

# Run TaxSEA for ALDEx2
aldex2_taxsea_results <- TaxSEA(taxon_ranks = aldex2_ranks)

# -------------------------------------------------------------------------
# Plotting Correlations Between Methods
# -------------------------------------------------------------------------

# Combine results into a data frame
results_plot_data <- data.frame(
  FastANCOM = fastancom_taxsea_results$Metabolite_producers$PValue,
  LINDA = linda_metabo_results$PValue,
  ALDEx2 = aldex2_taxsea_results$Metabolite_producers$PValue
)

# Generate correlation plots
HMP_cor_plots <- grid.arrange(
  ggplot(results_plot_data, aes(x = LINDA, y = FastANCOM)) +
    geom_point() + geom_smooth(method = "lm") + theme_classic(),
  ggplot(results_plot_data, aes(x = LINDA, y = ALDEx2)) +
    geom_point() + geom_smooth(method = "lm") + theme_classic(),
  nrow = 1
)

# -------------------------------------------------------------------------
# Visualization of Results
# -------------------------------------------------------------------------

# Identify SCFA-producing taxa
scfa_taxa <- c(
  TaxSEA_db$GutMGene_producers_of_Butyrate, 
  TaxSEA_db$GutMGene_producers_of_Acetate, 
  TaxSEA_db$GutMGene_producers_of_Propionate
)

# Add SCFA annotation to LINDA results
ids_of_interest = rownames(linda_results)
ids_of_interest = get_ncbi_taxon_ids(ids_of_interest)
linda_results$Species = rownames(linda_results)
linda_results = linda_results[rownames(linda_results) %in% names(ids_of_interest),]
linda_results$Taxon = ids_of_interest[rownames(linda_results)]
linda_results = linda_results[!is.na(linda_results$Taxon),]
linda_results$SCFAs <- linda_results$Taxon %in% unlist(scfa_taxa)

# Plot LINDA results with SCFA annotation
ggplot(linda_results, aes(x = log2FoldChange, y = -log10(padj), fill = SCFAs)) +
  geom_point(shape = 21, size = 3, aes(alpha = SCFAs)) +
  scale_alpha_manual(values = c(0.6, 1)) +
  scale_fill_manual(values = c("grey9", "#4aba91")) +
  geom_text_repel(data = linda_results[linda_results$SCFAs, ], aes(label = Species), size = 2) +
  theme_classic() + 
  geom_hline(yintercept = -log10(0.1), linetype = 3)


#############################
#### QinJ_2012 Dataset ####
#############################

# ------------------------------------------------------------------------------
# Step 1: Download and Process the QinJ_2012 Dataset
# ------------------------------------------------------------------------------
# Download the QinJ_2012 dataset using the curatedMetagenomicData package
QinJ_2012_cmd_object <- curatedMetagenomicData(
  pattern = "QinJ_2012.relative_abundance", 
  counts = TRUE, 
  dryrun = FALSE
)

# Extract the abundance data as a matrix for further processing
QinJ_2012_counts <- assay(QinJ_2012_cmd_object$`2021-10-14.QinJ_2012.relative_abundance`)

# Filter: Keep taxa detected in at least 10 samples with >1000 reads
QinJ_2012_counts <- QinJ_2012_counts[apply(QinJ_2012_counts > 1000, 1, sum) > 10, ]

# ------------------------------------------------------------------------------
# Step 2: Match Metadata with the Count Data and Clean Data
# ------------------------------------------------------------------------------
# Match metadata based on sample IDs in the count data
QinJ_2012_metadata <- allmeta_data[allmeta_data$sample_id %in% colnames(QinJ_2012_counts), ]
rownames(QinJ_2012_metadata) <- QinJ_2012_metadata$sample_id

# Remove samples with antibiotic use
QinJ_2012_metadata <- QinJ_2012_metadata[QinJ_2012_metadata$antibiotics_current_use == "no", ]

# Ensure that sample IDs match between metadata and count data
valid_samples <- intersect(rownames(QinJ_2012_metadata), colnames(QinJ_2012_counts))
QinJ_2012_counts <- QinJ_2012_counts[, valid_samples]
QinJ_2012_metadata <- QinJ_2012_metadata[valid_samples, ]

# Remove duplicate entries from the same participant and any with missing conditions
QinJ_2012_metadata <- QinJ_2012_metadata[order(QinJ_2012_metadata$visit_number, decreasing = TRUE), ]
QinJ_2012_metadata <- QinJ_2012_metadata[!duplicated(QinJ_2012_metadata$subject_id), ]
QinJ_2012_metadata <- QinJ_2012_metadata[!is.na(QinJ_2012_metadata$study_condition), ]

# Ensure the count data matches the cleaned metadata
QinJ_2012_counts <- QinJ_2012_counts[, rownames(QinJ_2012_metadata)]

# ------------------------------------------------------------------------------
# Step 3: Assign Species Names to Taxa
# ------------------------------------------------------------------------------
# Extract and clean species names from the row names
species_names <- sapply(rownames(QinJ_2012_counts), function(x) strsplit(x, "\\|")[[1]][7])
species_names <- gsub("s__", "", species_names)  # Remove prefix 's__'
rownames(QinJ_2012_counts) <- species_names

# ------------------------------------------------------------------------------
# Step 4: Differential Abundance Analysis with LINDA, FastANCOM, and ALDEx2
# ------------------------------------------------------------------------------

## LINDA
# Run LINDA for differential abundance analysis
linda_res <- linda(QinJ_2012_counts, QinJ_2012_metadata, formula = '~study_condition')
linda_results <- linda_res$output$study_conditionT2D  # Extract relevant results
write.csv(linda_results, "QinJ_2012_linda_results.csv")  # Save results
table(linda_results$padj < 0.1)  # Check significant results

# Extract LINDA ranks
linda_ranks <- linda_results$log2FoldChange
names(linda_ranks) <- rownames(linda_results)

## FastANCOM
# Run FastANCOM for differential abundance analysis
fastancom_fit <- fastANCOM(Y = t(QinJ_2012_counts), x = QinJ_2012_metadata$study_condition)
fastancom_results <- fastancom_fit$results$final
table(fastancom_results$log2FC.qval < 0.1)  # Check significant results
fastancom_ranks <- fastancom_results$log2FC
names(fastancom_ranks) = rownames(fastancom_results)

## ALDEx2
# Run ALDEx2 for differential abundance analysis
aldex_results <- aldex(QinJ_2012_counts, QinJ_2012_metadata$study_condition)
table(aldex_results$wi.eBH < 0.1)  # Check significant results
aldex2_ranks <- aldex_results$effect
names(aldex2_ranks) = rownames(aldex_results)
# ------------------------------------------------------------------------------
# Step 5: Taxon Set Enrichment Analysis Using TaxSEA
# ------------------------------------------------------------------------------
# Filter common taxa across methods
common_taxa <- intersect(names(linda_ranks), names(fastancom_ranks))
common_taxa <- intersect(common_taxa, names(aldex2_ranks))

## LINDA + TaxSEA
# Run TaxSEA on LINDA results
linda_taxsea_results <- TaxSEA(taxon_ranks = linda_ranks)
linda_metabo_results <- linda_taxsea_results$Metabolite_producers

# Merge LINDA TaxSEA results and filter by FDR < 0.1
linda_taxsea_merged <- rbind(
  linda_taxsea_results$Metabolite_producers, 
  linda_taxsea_results$Health_associations
)
linda_taxsea_merged <- linda_taxsea_merged[linda_taxsea_merged$FDR < 0.1, ]
write.csv(linda_taxsea_merged, "QinJ_2012_Linda_TaxSEA_results_merged.csv")

## FastANCOM + TaxSEA
# Run TaxSEA on FastANCOM results
fastancom_taxsea_results <- TaxSEA(taxon_ranks = fastancom_ranks)

## ALDEx2 + TaxSEA
# Run TaxSEA on ALDEx2 results
aldex2_taxsea_results <- TaxSEA(taxon_ranks = aldex2_ranks)

# ------------------------------------------------------------------------------
# Step 6: Correlation Analysis Between Methods
# ------------------------------------------------------------------------------
# Perform Pearson correlation between LINDA, FastANCOM, and ALDEx2 results
cor.test(fastancom_taxsea_results$Metabolite_producers$PValue, 
         linda_metabo_results$PValue, method = "pearson")
cor.test(aldex2_taxsea_results$Metabolite_producers$PValue, 
         linda_metabo_results$PValue, method = "pearson")

# Create a data frame for plotting correlations
plot_data <- data.frame(
  FastANCOM = fastancom_taxsea_results$Metabolite_producers$PValue,
  LINDA = linda_metabo_results$PValue,
  ALDEx2 = aldex2_taxsea_results$Metabolite_producers$PValue
)

# Generate correlation plots
qin_cor_plots <- grid.arrange(
  ggplot(plot_data, aes(x = LINDA, y = FastANCOM)) +
    geom_point() + geom_smooth(method = "lm") + theme_classic(),
  ggplot(plot_data, aes(x = LINDA, y = ALDEx2)) +
    geom_point() + geom_smooth(method = "lm") + theme_classic(),
  nrow = 1
)

# Combine and display correlation plots with HMP dataset plots
grid.arrange(HMP_cor_plots, qin_cor_plots, nrow = 2)
