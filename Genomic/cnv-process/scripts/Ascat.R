#ASCAT
rm(list = ls (all = TRUE))
library(data.table)
devtools::install_github('VanLoo-lab/ascat/ASCAT')
library(ASCAT)

ascat.bc = ascat.loadData("./data/Tumor_LogR.txt","./data/Tumor_BAF.txt")
gg<-ascat.predictGermlineGenotypes(ascat.bc, platform = "AffySNP6")
ascat.bc = ascat.aspcf(ascat.bc, ascat.gg=gg)
ascat.output = ascat.runAscat(ascat.bc)

write.table(ascat.output$segments,"./data/segment_ascat.txt",col.names = T,row.names = F,sep = "\t")
write.csv(ascat.output$segments_raw,"./data/segment_ascat_raw20220107.csv",row.names = F)