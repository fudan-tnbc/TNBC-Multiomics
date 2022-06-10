## TNBC2019 Genomic QC - batch information add
## Author：Qingwang Chen
## Date：2022-6-6
## Version：1.0

# Set work directory and import libraries
rm(list = ls (all = TRUE))

source("./scripts/libraries.R")
source("./scripts/parameter.R")

######### meatadata import
qc_metadata <- read.csv("./data/TNBC2019_WESQC_metadata.csv",header=T)

######## batch info import
WES_batch_info <- read.csv("./data/batch_info_TNBC_WES.csv",header=T)
WES_batch_info$sample <- sapply(strsplit(as.character(WES_batch_info$final_name),"\\."),function(x){x[1]})
WES_batch_info$type <- sapply(strsplit(as.character(WES_batch_info$final_name),"\\."),function(x){x[2]})
WES_batch_info$sample_id <- paste0(WES_batch_info$sample,"_",WES_batch_info$type)
sort(unique(WES_batch_info$batch))
WES_batch_info$batch <- Replace(data=WES_batch_info$batch,pattern = c("batch_2$:B01","batch_3$:B02","batch_4$:B03",
                                                                      "batch_5$:B04","batch_6$:B05","batch_7$:B06",
                                                                      "batch_8$:B07","batch_9$:B08","batch_10:B09",
                                                                      "batch_11:B10","batch_12:B11","batch_13:B12",
                                                                      "batch_14:B13","batch_17_2:B14","batch_18:B15",
                                                                      "batch_19:B16","batch_20:B17","batch_21:B18",
                                                                      "batch_22$:B19","batch_22_2:B20","batch_28:B21",
                                                                      "batch_29:B22","batch_30:B23","batch_31:B24"))


######### meatadata add batch info
qc_metadata_batch <- merge(qc_metadata,WES_batch_info,by ="sample_id")

## export
write.csv(qc_metadata_batch,"./data/TNBC2019_WESQC_batch_metadata.csv", row.names = F, quote = F)

