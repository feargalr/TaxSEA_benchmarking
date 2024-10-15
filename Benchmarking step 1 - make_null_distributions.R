# -----------------------------------------------------------------------------
# Author: Feargal Ryan
# GitHub: https://github.com/feargalr/
# Description: Script to generate null fold change distributions between random 
#              subsets of health controls as part of benchmarking for the TaxSEA 
#              R package. A signal implantation approach will be used in later 
#              steps for comparisons.
# -----------------------------------------------------------------------------

# Load required libraries
library(curatedMetagenomicData)
library(TaxSEA)
library(MicrobiomeStat)


# Load and prepare metadata
allmeta_data <- sampleMetadata
sort(table(allmeta_data$disease))
allmeta_data <- allmeta_data[allmeta_data$disease %in% c("healthy"),]
sort(table(allmeta_data$study_name))

# Load count data (pre-saved RDS)
LifeLifesDeep_2016counts.df <- readRDS("LifeLifesDeep_2016counts.RDS")

# Keep taxa with >1000 reads in more than 4 samples
LifeLifesDeep_2016counts.df <- LifeLifesDeep_2016counts.df[apply(LifeLifesDeep_2016counts.df > 1000, 1, sum) > 4,]

# Filter metadata to only relevant samples
LifeLifesDeep_2016md.df <- allmeta_data[allmeta_data$sample_id %in% colnames(LifeLifesDeep_2016counts.df),]
rownames(LifeLifesDeep_2016md.df) <- LifeLifesDeep_2016md.df$sample_id
table(LifeLifesDeep_2016md.df$subject_id)

# Ensure data matching between metadata and counts
LifeLifesDeep_2016counts.df <- LifeLifesDeep_2016counts.df[, rownames(LifeLifesDeep_2016md.df)]
LifeLifesDeep_2016md.df <- LifeLifesDeep_2016md.df[, c("study_name", "sample_id", "subject_id", "study_condition", "age_category", "DNA_extraction_kit")]

# Assign species names as rownames
spec_vec <- sapply(as.character(rownames(LifeLifesDeep_2016counts.df)), function(y) {strsplit(x = y, split = "\\|")[[1]][7]})
spec_vec <- gsub("s__", "", spec_vec)
rownames(LifeLifesDeep_2016counts.df) <- spec_vec

# Initialize vectors and lists for storing results
my_maxs <- numeric(3)
my_mins <- numeric(3)
my_means <- numeric(3)
my_medians <- numeric(3)
my_sd <- numeric(3)
shapiro_pvalues <- numeric(3)
my_distributions <- list()
dist_enr <- 0

# Main loop for generating null distributions
for (i in 1:1000) {
  sample_ids <- sample(colnames(LifeLifesDeep_2016counts.df), 50)
  test_data <- LifeLifesDeep_2016counts.df[, sample_ids]
  test_data <- test_data[apply(test_data > 1000, 1, sum) > 5, ]
  test_meta <- LifeLifesDeep_2016md.df[colnames(test_data), ]
  test_meta$Group <- as.factor(c(rep("Group1", 25), rep("Group2", 25)))
  
  # Run LINDA for differential abundance
  linda_results <- linda(test_data, test_meta, formula = '~Group', alpha = 0.1)
  test_ranks <- linda_results$output$GroupGroup2$log2FoldChange
  names(test_ranks) <- rownames(linda_results$output$GroupGroup2)
  
  # Run TaxSEA on the resulting ranks
  taxsea_res <- TaxSEA(taxon_ranks = test_ranks)
  meta_res <- taxsea_res$Metabolite_producers
  disease_res <- taxsea_res$Health_associations
  bugsigdb_res <- taxsea_res$BugSigdB
  meta_res <- meta_res[meta_res$FDR < 0.05, ]
  disease_res <- disease_res[disease_res$FDR < 0.05, ]
  bugsigdb_res <- bugsigdb_res[bugsigdb_res$FDR < 0.05, ]
  
  # Check for any enriched taxon sets
  if (nrow(meta_res) > 0 | nrow(disease_res) > 0 | nrow(bugsigdb_res) > 0) {
    next  # Skip to the next iteration
  }
  
  # Check for non-normality in fold change distribution
  if (shapiro.test(test_ranks)$p.value < 0.01) {
    next  # Skip to the next iteration
  }
  
  # Store results for the current iteration
  my_distributions <- c(my_distributions, list(test_ranks))
  my_maxs[i] <- max(test_ranks)
  my_mins[i] <- min(test_ranks)
  my_means[i] <- mean(test_ranks)
  my_medians[i] <- median(test_ranks)
  my_sd[i] <- sd(test_ranks)
  shapiro_pvalues[i] <- shapiro.test(test_ranks)$p.value
}

# Save the resulting null distributions to an RDS file
saveRDS(my_distributions, "Null_distributions.RDS")