## TNBC2019 RNAseq QC - 
## Author：Qingwang Chen
## Date：2022-6-6
## Version：1.0

# Set work directory and import libraries
rm(list = ls (all = TRUE))

source("./scripts/libraries.R")

######### meatadata import
qc_metadata <- read.csv("./data/TNBC2019_RNASeqQC_dir_metadata.csv",header=T)

######## batch info import
RNAseq_batch_info <- read.csv("./data/RNAseq_batch_info.csv",header=T)
RNAseq_batch_info$sample<- sapply(strsplit(as.character(RNAseq_batch_info$sample_id),"\\."),function(x){x[1]})
RNAseq_batch_info$type<- sapply(strsplit(as.character(RNAseq_batch_info$sample_id),"\\."),function(x){x[2]})
RNAseq_batch_info$type <- gsub("rep","TT",RNAseq_batch_info$type)
RNAseq_batch_info$type[which(is.na(RNAseq_batch_info$type))] <- "TT"
RNAseq_batch_info$sample_id <- paste0(RNAseq_batch_info$sample,"_",RNAseq_batch_info$type)
# batch2 here is actually batch1 in the origin
RNAseq_batch_info$batch <- gsub("batch1","x",RNAseq_batch_info$batch )
RNAseq_batch_info$batch <- gsub("batch2","batch1",RNAseq_batch_info$batch )
RNAseq_batch_info$batch <- gsub("x","batch2",RNAseq_batch_info$batch )

######### meatadata add batch info
qc_metadata_batch <- merge(qc_metadata,RNAseq_batch_info,by ="sample_id")

## export
write.csv(qc_metadata_batch,"./data/TNBC2019_RNASeqQC_batch_metadata.csv", row.names = F, quote = F)

