"""
Author: Feargal Ryan
GitHub: https://github.com/feargalr/

Description:
This script uses the Microbe Set Enrichment Analysis (MSEA) library to test for 
the enrichment and depletion of taxon sets in various differential abundance outputs.
"""

# Import necessary libraries
from pprint import pprint
import msea
from msea import SetLibrary

# Define the file path to the microbe set library (GMT format)
gmt_filepath = 'https://bitbucket.org/wangz10/msea/raw/aee6dd184e9bde152b4d7c2f3c7245efc1b80d23/msea/data/human_genes_associated_microbes/set_library.gmt'

# Load the GMT file containing microbe sets
d_gmt = msea.read_gmt(gmt_filepath)

# Display the number of microbe sets available in the library
print('Number of microbe-sets:', len(d_gmt))

# Print a few reference sets to get a sense of the data structure
pprint(list(d_gmt.items())[:3])

# Define a sample input set of microbes
microbe_set_input = set([
    'Colwellia', 'Deinococcus', 'Idiomarina', 'Neisseria', 
    'Pseudidiomarina', 'Pseudoalteromonas'
])

# Perform enrichment analysis on the sample set
result = msea.enrich(microbe_set_input, d_gmt=d_gmt, universe=1000)
print(result)
print("---------------------------------------------------------------------------------------")

# HMP study: Test for depletion among decreased taxa
hmp_decreased_taxa = set([
    'Ruminococcus torques', 'Ruminococcus bromii', 
    'Firmicutes bacterium CAG 83', 'Roseburia sp CAG 182'
])
hmp_result = msea.enrich(hmp_decreased_taxa, d_gmt=d_gmt, universe=1000)
print("HMP Study - Depleted Sets:")
print(hmp_result)
print("---------------------------------------------------------------------------------------")

# Hall study: Test for depletion among decreased taxa
hall_decreased_taxa = set([
    'Roseburia hominis', 'Alistipes putredinis', 'Eubacterium rectale', 
    'Fusicatenibacter saccharivorans', 'Agathobaculum butyriciproducens', 
    'Anaerostipes hadrus', 'Collinsella aerofaciens', 'Alistipes finegoldii',
    # ... (other taxa omitted for brevity)
])
hall_result = msea.enrich(hall_decreased_taxa, d_gmt=d_gmt, universe=1000)
print("Hall Study - Depleted Sets:")
print(hall_result)
print("---------------------------------------------------------------------------------------")

# Ijaz study: Test for enrichment among increased taxa
ijaz_increased_taxa = set([
    'Eggerthella lenta', 'Klebsiella pneumoniae', 'Klebsiella variicola', 
    'Klebsiella quasipneumoniae', 'Eisenbergiella tayi',
    # ... (other taxa omitted for brevity)
])
ijaz_increased_result = msea.enrich(ijaz_increased_taxa, d_gmt=d_gmt, universe=1000)
print("Ijaz Study - Enriched Sets:")
print(ijaz_increased_result)
print("---------------------------------------------------------------------------------------")

# Ijaz study: Test for depletion among decreased taxa
ijaz_decreased_taxa = set([
    'Fusicatenibacter saccharivorans', 'Bifidobacterium longum', 
    'Collinsella aerofaciens', 'Ruminococcus torques', 
    'Anaerostipes hadrus', 'Agathobaculum butyriciproducens',
    # ... (other taxa omitted for brevity)
])
ijaz_decreased_result = msea.enrich(ijaz_decreased_taxa, d_gmt=d_gmt, universe=1000)
print("Ijaz Study - Depleted Sets:")
print(ijaz_decreased_result)
print("---------------------------------------------------------------------------------------")

# Nielson study: Test for enrichment among increased taxa
nielson_increased_taxa = set([
    'Bacteroides coprocola', 'Bacteroides stercoris', 'Collinsella intestinalis', 
    'Bifidobacterium bifidum', 'Olsenella scatoligenes',
    # ... (other taxa omitted for brevity)
])
nielson_increased_result = msea.enrich(nielson_increased_taxa, d_gmt=d_gmt, universe=1000)
print("Nielson Study - Enriched Sets:")
print(nielson_increased_result)
print("---------------------------------------------------------------------------------------")

# Nielson study: Test for depletion among decreased taxa
nielson_decreased_taxa = set([
    'Roseburia faecis', 'Butyrivibrio crossotus', 'Ruminococcus bromii', 
    'Faecalibacterium prausnitzii', 'Alistipes putredinis',
    # ... (other taxa omitted for brevity)
])
nielson_decreased_result = msea.enrich(nielson_decreased_taxa, d_gmt=d_gmt, universe=1000)
print("Nielson Study - Depleted Sets:")
print(nielson_decreased_result)
print("---------------------------------------------------------------------------------------")
