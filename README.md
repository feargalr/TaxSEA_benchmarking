# TaxSEA_benchmarking
This repository contains code and workflows for benchmarking **[TaxSEA](https://github.com/feargalr/taxsea)**. The process consists of three key steps designed to evaluate its performance in detecting taxon set enrichment.

### Step 1: Null Distribution Creation
Generate log2 fold change distributions by comparing random subsets of healthy controls (Lifelines Deep cohort via `curatedMetagenomicData`). Ensure no significant enrichments by running TaxSEA on these distributions. Confirm they follow a normal distributions. These log2 fold change distributions represent a baseline for TaxSEA testing.

### Step 2: Signal Implantation
Implant artificial enrichment signals into the null distributions created in Step 1 to test TaxSEA’s ability to detect varying effect sizes.
Enrichment signals are added at different effect sizes across multiple sets in the TaxSEA database. TaxSEA is run on these distributions to measure how well it recovers the implanted signals.

### Step 3: Real Data Evaluation on real data
Apply TaxSEA to real datasets from `curatedMetagenomicData`. Perform differential abundance analysis between healthy controls and IBD samples.Use these fold changes as input for TaxSEA to assess its performance on real-world data.

## Additional steps
We also tested the ability of **[MSEA](https://msea.readthedocs.io/en/latest/)** to recover enrichments based upon the list of differential abundant taxa detected in Step 3. 
