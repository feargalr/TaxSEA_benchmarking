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

# Perform enrichment analysis on the test taxa
result = msea.enrich(microbe_set_input, d_gmt=d_gmt, universe=1000)
print(result)
print("---------------------------------------------------------------------------------------")



# HMP study: Test for depletion among decreased taxa
hmp_decreased_taxa = set([
    'Ruminococcus torques', 'Ruminococcus bromii', 'Rosburia hominis', 
    'Firmicutes bacterium CAG 83', 'Roseburia sp CAG 182'
])
hmp_result = msea.enrich(hmp_decreased_taxa, d_gmt=d_gmt, universe=1000)
print("HMP Study - Depleted Sets:")
print(hmp_result)
print("---------------------------------------------------------------------------------------")

# Qin study: Test for depletion among decreased taxa
Qin_decreased_taxa = set([
    'Roseburia intestinalis', 'Haemophilus parainfluenzae', 'Anaerostipes hadrus', 
    'Prevotella stercorea', 'Clostridium sp CAG 299', 'Agathobaculum butyriciproducens',
    'Roseburia inulinivorans', 'Veillonella dispar', 'Ruminococcus callidus',
    'Faecalibacterium prausnitzii', 'Prevotella sp CAG 520', 'Citrobacter freundii',
    'Actinomyces sp ICM47', 'Citrobacter portucalensis', 'Alistipes finegoldii'
])
Qin_result = msea.enrich(Qin_decreased_taxa, d_gmt=d_gmt, universe=1000)
print("Qin Study - Depleted Sets:")
print(Qin_result)
print("---------------------------------------------------------------------------------------")

# Qin study: Test for depletion among increased taxa
Qin_increased_taxa = set([
    'Acidaminococcus_sp_CAG_542', 'Lactobacillus_amylovorus', 'Christensenella_minuta', 
    'Bifidobacterium_breve', 'Dielma_fastidiosa', 'Clostridium_aldenense', 
    'Olsenella_scatoligenes', 'Aeriscardovia_aeriphila', 'Anaerotruncus_colihominis', 
    'Clostridium_asparagiforme', 'Streptococcus_anginosus_group', 'Eisenbergiella_tayi', 
    'Enorma_massiliensis', 'Lawsonibacter_asaccharolyticus', 'Alistipes_indistinctus', 
    'Coprobacillus_cateniformis', 'Gordonibacter_pamelaeae', 'Clostridium_symbiosum', 
    'Clostridium_bolteae', 'Lactobacillus_mucosae', 'Clostridium_innocuum', 
    'Ruthenibacterium_lactatiformans', 'Intestinimonas_butyriciproducens', 
    'Clostridium_citroniae', 'Hungatella_hathewayi', 'Erysipelatoclostridium_ramosum'
])

Qin_result = msea.enrich(Qin_increased_taxa, d_gmt=d_gmt, universe=1000)
print("Qin Study - Enriched Sets:")
print(Qin_result)
print("---------------------------------------------------------------------------------------")












