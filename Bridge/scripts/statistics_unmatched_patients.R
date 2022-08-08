## Title: TNBC Bridge RNAseq Analysis-- statistic unmatched patients##
## Author: Qingwang Chen             ##
## Date: 2022-07-02                  ##  
## Version: V3.0                     ##
#######################################
# Set work directory and import libraries
rm(list = ls (all = TRUE))
source("./scripts/libraries.R")
source("./scripts/parameter.R")


Tes_Clusters_old <- read.csv("./results/tables/FUSCCTNBC_old_subtypes.csv")
colnames(Tes_Clusters_old) <- c("patient","subtype_old")

Tes_Clusters_new <- read.csv("./results/tables/FUSCCTNBC_new_subtypes.csv")
colnames(Tes_Clusters_new) <- c("patient","subtype_new")
Tes_Clusters_new$patient <- gsub("_rep","",Tes_Clusters_new$patient)

Tes_Clusters <- merge(Tes_Clusters_old,Tes_Clusters_new,by="patient")
MES_unmatched <- Tes_Clusters$patient[Tes_Clusters$subtype_old == "MES" &
                                      Tes_Clusters$subtype_new != "MES"  ]
saveRDS(MES_unmatched,"./Rdata/MES_unmatched.rds")

patients_unmatched <- Tes_Clusters$patient[Tes_Clusters$subtype_new != Tes_Clusters$subtype_old]
saveRDS(patients_unmatched,"./Rdata/patients_unmatched.rds")
