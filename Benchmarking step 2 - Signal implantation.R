# ------------------------------------------------------------------------------
# Step 2: Signal Implementation for TaxSEA Benchmarking
# ------------------------------------------------------------------------------
# This script performs signal implantation by simulating taxon set enrichment 
# under various conditions. It generates fold change distributions and evaluates
# the performance of the TaxSEA package under different effect sizes and set sizes.
# ------------------------------------------------------------------------------

# Load necessary libraries and custom functions
library(TaxSEA)  # For TaxSEA benchmarking methods
source("TaxSEA_Benchmarking function.R")  # Custom function for signal implantation

# Load pre-computed null distributions from healthy controls
null_distributions <- readRDS("Null_distributions.RDS")

# Initialize an empty results data frame to store all results
benchmark_results <- data.frame(
  pvalue = c(1, 1, 1, 1),
  set_min_size = c(10, 10, 10, 10),
  set_max_size = c(10, 10, 10, 10),
  FC = c(10, 10, 10, 10),
  set_fraction = c(1, 1, 1, 1)
)

# Define parameters for benchmarking
min_set_sizes <- seq(5, 100, by = 10)  # Range of set sizes to evaluate
fold_changes <- c(-1, -0.42, 0.0, 1, 2, 3)  # Different fold changes to test
set_fractions <- c(1)  # Set fractions to use

# Perform benchmarking through nested loops
for (iteration in 1:1000) {  # Outer loop for multiple iterations
  for (set_size in min_set_sizes) {  # Loop through set sizes
    for (FC in fold_changes) {  # Loop through fold changes
      for (SF in set_fractions) {  # Loop through set fractions
        
        # Select a random null distribution from the list
        random_dist <- unlist(null_distributions[[sample(1:length(null_distributions), 1)]])
        
        # Run the TaxSEA benchmarking with the selected parameters
        benchmark_res <- benchmarking_TaxSEA(
          input_distribution = random_dist,
          set_fraction = SF,
          num_iterations = 1,
          set_min_size = set_size,
          set_max_size = set_size + 10,
          min_change = FC,
          max_change = FC + 2
        )
        
        # Store the results in a data frame
        res_df <- data.frame(
          pvalue = benchmark_res,
          set_min_size = rep(set_size, length(benchmark_res)),
          set_max_size = rep(set_size + 10, length(benchmark_res)),
          FC = rep(FC + 1, length(benchmark_res)),
          set_fraction = rep(SF, length(benchmark_res))
        )
        
        # Append the current results to the main results data frame
        benchmark_results <- rbind(benchmark_results, res_df)
      }
    }
  }
}

# ------------------------------------------------------------------------------
# Data Visualization: Plotting the Results
# ------------------------------------------------------------------------------

# Load ggplot2 for visualization
library(ggplot2)

# Prepare results for plotting
plot_data <- benchmark_results[benchmark_results$FC != 10, ]  # Filter out dummy data

# Define set size categories
plot_data$Set_size <- paste(plot_data$set_min_size, plot_data$set_max_size, sep = "-")
plot_data$Set_size[plot_data$set_min_size < 10] <- "1. Less than 10"
plot_data$Set_size[plot_data$set_min_size >= 10 & plot_data$set_min_size < 50] <- "2. 10 - 50"
plot_data$Set_size[plot_data$set_min_size >= 50] <- "3. 50+"

# Convert fold change values to readable factors
plot_data$Fold_change <- as.character(plot_data$FC)
plot_data$Fold_change[plot_data$Fold_change == "0"] <- "1x"
plot_data$Fold_change[plot_data$Fold_change == "0.58"] <- "1.5x"
plot_data$Fold_change[plot_data$Fold_change == "1"] <- "2x"
plot_data$Fold_change[plot_data$Fold_change == "2"] <- "4x"
plot_data$Fold_change[plot_data$Fold_change == "3"] <- "8x"
plot_data$Fold_change[plot_data$Fold_change == "4"] <- "16x"

# Order the factors for proper plotting
plot_data$Fold_change <- factor(plot_data$Fold_change, 
                                levels = c("1x", "1.5x", "2x", "4x", "8x", "16x"))

# Create a boxplot to visualize the results
ggplot(plot_data, aes(x = Fold_change, y = -log10(pvalue), color = Set_size)) +
  geom_boxplot(outlier.shape = NA) +  # Avoid plotting outliers for clarity
  geom_hline(yintercept = -log10(0.05), linetype = 3) +  # Add significance threshold
  theme_bw() +  # Use a clean theme
  xlab("Mean fold change of taxon set") +  # X-axis label
  scale_color_manual(values = c("goldenrod3", "deepskyblue3", "darkorchid4"))  # Set custom colors