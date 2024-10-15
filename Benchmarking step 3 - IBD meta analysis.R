# -----------------------------------------------------------------------------
# Author: Feargal Ryan
# GitHub: https://github.com/feargalr/
# Description: Script for downloading various datasets from the 
#              curatedMetagenomicData R package and performing differential 
#              expression and enrichment analysis using LINDA and TaxSEA.
# -----------------------------------------------------------------------------

library(ggplot2)
library(curatedMetagenomicData)
library(TaxSEA)
library(MicrobiomeStat)


allmeta_data = sampleMetadata

IBD_meta_data = allmeta_data[allmeta_data$study_condition=="IBD",]
IBD_meta_data = allmeta_data[allmeta_data$study_name %in% IBD_meta_data$study_name,]

#Identify the number of control samples per study
table(IBD_meta_data$study_condition,IBD_meta_data$study_name)["control",]



########################
 #### HMP_2019_ibdmdb ####
########################

# HMP_2019_ibdmdb dataset
HMP_2019_ibdmdb_cmd_object = curatedMetagenomicData(pattern = "2021-10-14.HMP_2019_ibdmdb.relative_abundance",counts = TRUE,dryrun = FALSE)
HMP_2019_ibdmdb_counts.df = HMP_2019_ibdmdb_cmd_object
HMP_2019_ibdmdb_counts.df = (HMP_2019_ibdmdb_counts.df$`2021-10-14.HMP_2019_ibdmdb.relative_abundance`)

HMP_2019_ibdmdb_counts.df = assay(HMP_2019_ibdmdb_counts.df)
HMP_2019_ibdmdb_counts.df = HMP_2019_ibdmdb_counts.df[apply(HMP_2019_ibdmdb_counts.df>1000,1,sum)>3,]
HMP_2019_ibdmdb_md.df = allmeta_data[allmeta_data$sample_id %in% colnames(HMP_2019_ibdmdb_counts.df),]
rownames(HMP_2019_ibdmdb_md.df) = HMP_2019_ibdmdb_md.df$sample_id

#Remove any samples with current antibiotics use
HMP_2019_ibdmdb_md.df = HMP_2019_ibdmdb_md.df[HMP_2019_ibdmdb_md.df$antibiotics_current_use == "no",]
s2k = intersect(rownames(HMP_2019_ibdmdb_md.df),colnames(HMP_2019_ibdmdb_counts.df))
HMP_2019_ibdmdb_counts.df = HMP_2019_ibdmdb_counts.df[,s2k]
HMP_2019_ibdmdb_md.df = HMP_2019_ibdmdb_md.df[s2k,]

## Ensure data is matching

HMP_2019_ibdmdb_md.df$disease_subtype[is.na(HMP_2019_ibdmdb_md.df$disease_subtype)] = "Control"
HMP_2019_ibdmdb_md.df$disease_subtype = factor(HMP_2019_ibdmdb_md.df$disease_subtype,levels=c("Control","UC","CD"))

## Remove duplicates from the same participant 
HMP_2019_ibdmdb_md.df = HMP_2019_ibdmdb_md.df[order(HMP_2019_ibdmdb_md.df$visit_number,decreasing = TRUE),]
HMP_2019_ibdmdb_md.df = HMP_2019_ibdmdb_md.df[!duplicated(HMP_2019_ibdmdb_md.df$subject_id),]
HMP_2019_ibdmdb_counts.df = HMP_2019_ibdmdb_counts.df[,rownames(HMP_2019_ibdmdb_md.df)]

## Assigning rows as species names ##
spec_vec = sapply(as.character(rownames(HMP_2019_ibdmdb_counts.df)),function(y) {strsplit(x = y,split="\\|")[[1]][7]})
names(spec_vec) = NULL
spec_vec = gsub("s__","",spec_vec)
rownames(HMP_2019_ibdmdb_counts.df) = spec_vec

HMP_2019_ibdmdb_linda_res = linda(HMP_2019_ibdmdb_counts.df, 
                              HMP_2019_ibdmdb_md.df, 
                              formula = '~age_category + study_condition')

HMP_2019_ibdmdb_ranks = HMP_2019_ibdmdb_linda_res$output$study_conditionIBD$log2FoldChange
write.csv(HMP_2019_ibdmdb_linda_res$output$study_conditionIBD,"HMP_2019_ibdmdb_LinDA_results.csv")
names(HMP_2019_ibdmdb_ranks) = rownames(HMP_2019_ibdmdb_linda_res$output$study_conditionIBD)
HMP_2019_ibdmdb_TaxSEA_results.df = TaxSEA(taxon_ranks = HMP_2019_ibdmdb_ranks)



########################
#### HallAB_2017 ####
########################
allmeta_data = sampleMetadata

# HallAB_2017 dataset
HallAB_2017_cmd_object = curatedMetagenomicData(pattern = "2021-10-14.HallAB_2017.relative_abundance",counts = TRUE,dryrun = FALSE)
HallAB_2017_counts.df = HallAB_2017_cmd_object
HallAB_2017_counts.df = (HallAB_2017_counts.df$`2021-10-14.HallAB_2017.relative_abundance`)

HallAB_2017_counts.df = assay(HallAB_2017_counts.df)
HallAB_2017_counts.df = HallAB_2017_counts.df[apply(HallAB_2017_counts.df>1000,1,sum)>4,]
HallAB_2017_md.df = allmeta_data[allmeta_data$sample_id %in% colnames(HallAB_2017_counts.df),]
rownames(HallAB_2017_md.df) = HallAB_2017_md.df$sample_id

## Ensure data is matching
colnames(HallAB_2017_counts.df) %in% rownames(HallAB_2017_md.df)
HallAB_2017_counts.df = HallAB_2017_counts.df[,rownames(HallAB_2017_md.df)]
HallAB_2017_md.df = HallAB_2017_md.df[,c("study_name","sample_id","subject_id","study_condition","age_category","DNA_extraction_kit","visit_number")]

## Remove duplicate samples 
HallAB_2017_md.df = HallAB_2017_md.df[order(HallAB_2017_md.df$visit_number,decreasing = FALSE),]
HallAB_2017_md.df = HallAB_2017_md.df[!duplicated(HallAB_2017_md.df$subject_id),]
HallAB_2017_counts.df = HallAB_2017_counts.df[,rownames(HallAB_2017_md.df)]

## Assigning rows as species names ##
spec_vec = sapply(as.character(rownames(HallAB_2017_counts.df)),function(y) {strsplit(x = y,split="\\|")[[1]][7]})
names(spec_vec) = NULL
spec_vec = gsub("s__","",spec_vec)
rownames(HallAB_2017_counts.df) = spec_vec

HallAB_2017_linda_res = linda(HallAB_2017_counts.df, 
                              HallAB_2017_md.df, 
                              formula = '~study_condition+DNA_extraction_kit')
HallAB_2017_ranks = HallAB_2017_linda_res$output$study_conditionIBD$log2FoldChange
names(HallAB_2017_ranks) = rownames(HallAB_2017_linda_res$output$study_conditionIBD)
HallAB_2017_TaxSEA_results.df = TaxSEA(taxon_ranks = HallAB_2017_ranks)
write.csv(HallAB_2017_linda_res$output$study_conditionIBD,"HallAB_2017_LinDA_results.csv")








##############################
#### NielsenHB_2014 ####
##############################
allmeta_data = sampleMetadata
# NielsenHB_2014 dataset
NielsenHB_2014_cmd_object = curatedMetagenomicData(pattern = "NielsenHB_2014.relative_abundance",counts = TRUE,dryrun = FALSE)
NielsenHB_2014_counts.df = NielsenHB_2014_cmd_object
NielsenHB_2014_counts.df = (NielsenHB_2014_counts.df$`2021-03-31.NielsenHB_2014.relative_abundance`)
#saveRDS(NielsenHB_2014_counts.df,"NielsenHB_2014_counts.df.RDS")
#NielsenHB_2014_counts.df = readRDS("NielsenHB_2014_counts.df.RDS")
NielsenHB_2014_counts.df = assay(NielsenHB_2014_counts.df)
NielsenHB_2014_counts.df = NielsenHB_2014_counts.df[apply(NielsenHB_2014_counts.df>1000,1,sum)>4,]
NielsenHB_2014_md.df = allmeta_data[allmeta_data$sample_id %in% colnames(NielsenHB_2014_counts.df),]
## Ensure data is matching
NielsenHB_2014_md.df = NielsenHB_2014_md.df[,c("study_name","sample_id","subject_id","study_condition","age_category","DNA_extraction_kit","visit_number")]
## Remove duplicate samples 
NielsenHB_2014_md.df = NielsenHB_2014_md.df[order(NielsenHB_2014_md.df$visit_number,decreasing = FALSE),]
NielsenHB_2014_md.df = NielsenHB_2014_md.df[!duplicated(NielsenHB_2014_md.df$subject_id),]
rownames(NielsenHB_2014_md.df) = NielsenHB_2014_md.df$sample_id
table(colnames(NielsenHB_2014_counts.df) %in% rownames(NielsenHB_2014_md.df))
NielsenHB_2014_md.df = NielsenHB_2014_md.df[rownames(NielsenHB_2014_md.df) %in% colnames(NielsenHB_2014_counts.df),]
NielsenHB_2014_counts.df = NielsenHB_2014_counts.df[,rownames(NielsenHB_2014_md.df)]
## Assigning rows as species names ##
spec_vec = sapply(as.character(rownames(NielsenHB_2014_counts.df)),function(y) {strsplit(x = y,split="\\|")[[1]][7]})
names(spec_vec) = NULL
spec_vec = gsub("s__","",spec_vec)
rownames(NielsenHB_2014_counts.df) = spec_vec
NielsenHB_2014_linda_res = linda(NielsenHB_2014_counts.df, 
                              NielsenHB_2014_md.df, 
                              formula = '~ study_condition + age_category')
NielsenHB_2014_ranks = NielsenHB_2014_linda_res$output$study_conditionIBD$log2FoldChange
names(NielsenHB_2014_ranks) = rownames(NielsenHB_2014_linda_res$output$study_conditionIBD)
NielsenHB_2014_TaxSEA_results.df = TaxSEA(taxon_ranks = NielsenHB_2014_ranks)
write.csv(NielsenHB_2014_linda_res$output$study_conditionIBD,"NielsenHB_2014_LinDA_results.csv")



##############################
#### IjazUZ_2017_2018 ####
##############################
allmeta_data = sampleMetadata
# IjazUZ_2017_2018 dataset
IjazUZ_2017_cmd_object = curatedMetagenomicData(pattern = "IjazUZ_2017.relative_abundance",counts = TRUE,dryrun = FALSE)
IjazUZ_2017_counts.df = IjazUZ_2017_cmd_object
IjazUZ_2017_counts.df = (IjazUZ_2017_counts.df$`2021-10-14.IjazUZ_2017.relative_abundance`)
#saveRDS(IjazUZ_2017_counts.df,"IjazUZ_2017_counts.df.RDS")
#IjazUZ_2017_counts.df = readRDS("IjazUZ_2017_counts.df.RDS")
IjazUZ_2017_counts.df = assay(IjazUZ_2017_counts.df)
IjazUZ_2017_counts.df = IjazUZ_2017_counts.df[apply(IjazUZ_2017_counts.df>1000,1,sum)>4,]
IjazUZ_2017_md.df = allmeta_data[allmeta_data$sample_id %in% colnames(IjazUZ_2017_counts.df),]
## Ensure data is matching
IjazUZ_2017_md.df = IjazUZ_2017_md.df[,c("study_name","sample_id","subject_id","study_condition","age_category","DNA_extraction_kit","visit_number")]
## Remove duplicate samples 
IjazUZ_2017_md.df = IjazUZ_2017_md.df[order(IjazUZ_2017_md.df$visit_number,decreasing = FALSE),]
IjazUZ_2017_md.df = IjazUZ_2017_md.df[!duplicated(IjazUZ_2017_md.df$subject_id),]
rownames(IjazUZ_2017_md.df) = IjazUZ_2017_md.df$sample_id
table(colnames(IjazUZ_2017_counts.df) %in% rownames(IjazUZ_2017_md.df))
IjazUZ_2017_md.df = IjazUZ_2017_md.df[rownames(IjazUZ_2017_md.df) %in% colnames(IjazUZ_2017_counts.df),]
IjazUZ_2017_counts.df = IjazUZ_2017_counts.df[,rownames(IjazUZ_2017_md.df)]
## Assigning rows as species names ##
spec_vec = sapply(as.character(rownames(IjazUZ_2017_counts.df)),function(y) {strsplit(x = y,split="\\|")[[1]][7]})
names(spec_vec) = NULL
spec_vec = gsub("s__","",spec_vec)
rownames(IjazUZ_2017_counts.df) = spec_vec
IjazUZ_2017_linda_res = linda(IjazUZ_2017_counts.df, 
                                      IjazUZ_2017_md.df, 
                                      formula = '~ study_condition + age_category')
IjazUZ_2017_ranks = IjazUZ_2017_linda_res$output$study_conditionIBD$log2FoldChange
names(IjazUZ_2017_ranks) = rownames(IjazUZ_2017_linda_res$output$study_conditionIBD)
IjazUZ_2017_TaxSEA_results.df = TaxSEA(taxon_ranks = IjazUZ_2017_ranks)
write.csv(IjazUZ_2017_linda_res$output$study_conditionIBD,"IjazUZ_2017_LinDA_results.csv")




IjazUZ_2017_metabolites.df = IjazUZ_2017_TaxSEA_results.df$Metabolite_producers
NielsenHB_2014_metabolites.df = NielsenHB_2014_TaxSEA_results.df$Metabolite_producers
HallAB_2017_metabolites.df = HallAB_2017_TaxSEA_results.df$Metabolite_producers
HMP_2019_ibdmdb_metabolites.df = HMP_2019_ibdmdb_TaxSEA_results.df$Metabolite_producers

IjazUZ_2017_metabolites.df$Study = "IjazUZ_2017"
NielsenHB_2014_metabolites.df$Study = "NielsenHB_2014"
HallAB_2017_metabolites.df$Study = "HallAB_2017"
HMP_2019_ibdmdb_metabolites.df$Study = "HMP_2019_ibdmdb"
all_metabolites.df = rbind(IjazUZ_2017_metabolites.df,NielsenHB_2014_metabolites.df,
                        HallAB_2017_metabolites.df,HMP_2019_ibdmdb_metabolites.df)

write.csv(all_metabolites.df,"TaxSEA_results_metabolites.csv")

all_metabolites.df$log10FDR = log10(all_metabolites.df$FDR)
all_metabolites.df$log10FDR[all_metabolites.df$NES > 0] = all_metabolites.df$log10FDR[all_metabolites.df$NES > 0]*-1
metabolites2plot = all_metabolites.df[all_metabolites.df$FDR < 0.1,"taxonSetName"]
metabolites2plot = table(metabolites2plot)
metabolites2plot = metabolites2plot[metabolites2plot>=2]
metabolites2plot = all_metabolites.df[all_metabolites.df$taxonSetName %in% names(metabolites2plot),]

ggplot(metabolites2plot,aes(x=log10FDR,y=taxonSetName,fill=Study))+
  geom_col(position = position_dodge())+theme_classic()+
  scale_fill_manual(values=c("#165241","grey44","#2d3a80","#bfa34c"))