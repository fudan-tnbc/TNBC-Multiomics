## TNBC2019 RNAseq - exploration of insert size and RIN
## Author：Qingwang Chen
## Date：2022-5-21
## Version：1.0

# Set work directory and import libraries
rm(list = ls (all = TRUE))

source("./scripts/libraries.R")
source("./scripts/parameter.R")
library(do)

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

######## 5. Insert size
qc_ins_size <- qc_metadata[,c("sample_id", "median_insert_size", 
                              "tissue_type")] %>% distinct(sample_id, .keep_all = T)
qc_ins_size_batch <- merge(qc_ins_size,WES_batch_info,by="sample_id")

## plot
Ins_denplo_area <- ggdensity(qc_ins_size, x = "median_insert_size",
                             add = "median", rug = TRUE, 
                             color = "tissue_type", fill = "tissue_type", alpha = 0.5, 
                             palette = mypal[c(5,4)],xlab = "Insert Size", ylab = "Density", 
                             legend.title = "Type")+ mytheme;Ins_denplo_area

ggsave("./results/charts/Insert_sizeForSamples_WES.png",width=6,height=4,dpi=300)

aggregate(qc_ins_size$median_insert_size, by = list(qc_ins_size$tissue_type), FUN = median)

## batch
qc_ins_size_batch$batch <- as.factor(qc_ins_size_batch$batch)
Ins_vioplo <- ggviolin(qc_ins_size_batch, "batch", "median_insert_size",
                       fill = "tissue_type", alpha = 0.5,
                       palette = mypal[c(5,4)],
                       add = "boxplot",
                       ylab = "Insert Size", xlab = "Batch")+ mytheme; Ins_vioplo
ggsave("./results/charts/Insert_sizeForBatches_WES_violin.png",width=10,height=5,dpi=300)

Ins_denplo_area <- ggdensity(qc_ins_size_batch, x = "median_insert_size",
                             add = "median", rug = TRUE, 
                             color = "batch", fill = "batch", alpha = 0.5, 
                             xlab = "Insert Size", ylab = "Density", 
                             legend.title = "Batch")+ mytheme ;Ins_denplo_area
ggsave("./results/charts/Insert_sizeForBatches_WES.png",width=6,height=4,dpi=300)

table(qc_ins_size_batch$batch)
# batch_18   84; batch_20 44; batch_28  64

## inner batch
p1 <- ggdensity(qc_ins_size_batch[qc_ins_size_batch$batch=="B15",], x = "median_insert_size",
                             add = "median", rug = TRUE, 
                             color = "tissue_type", fill = "tissue_type", alpha = 0.5, 
                             palette = mypal[c(5,4)],xlab = "Insert Size", ylab = "Density", 
                             legend.title = "Type")+ mytheme ;p1 
ggsave("./results/charts/Insert_sizeForBatch18_WES.png",width=6,height=4,dpi=300)

p2 <- ggdensity(qc_ins_size_batch[qc_ins_size_batch$batch=="B08",], x = "median_insert_size",
                             add = "median", rug = TRUE, 
                             color = "tissue_type", fill = "tissue_type", alpha = 0.5, 
                             palette = mypal[c(5,4)],xlab = "Insert Size", ylab = "Density", 
                             legend.title = "Type")+ mytheme ;p2
ggsave("./results/charts/Insert_sizeForBatch9_WES.png",width=6,height=4,dpi=300)

p3 <- ggdensity(qc_ins_size_batch[qc_ins_size_batch$batch=="B21",], x = "median_insert_size",
                             add = "median", rug = TRUE, 
                             color = "tissue_type", fill = "tissue_type", alpha = 0.5, 
                             palette = mypal[c(5,4)],xlab = "Insert Size", ylab = "Density", 
                             legend.title = "Type")+ mytheme ;p3
ggsave("./results/charts/Insert_sizeForBatch3_RNAseq.png",width=6,height=4,dpi=300)


panel.ins_size_acrossbatch <- ggarrange(Ins_denplo_area,nrow = 1,labels = c("a"))+
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm")); panel.ins_size_acrossbatch

panel.ins_size_intrabatch <- ggarrange(p1,p2,p3,nrow = 1,labels = c("b","c","d"), 
                            widths = c(1,1,1))+theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm")); panel.ins_size_intrabatch
panel.ins_size <- ggarrange(panel.ins_size_acrossbatch,panel.ins_size_intrabatch,nrow = 2,
                           heights = c(1,1))+ theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm")); panel.ins_size
ggsave("./results/charts/panel.ins_size.pdf",panel.ins_size,width=8,height=8,dpi=300)



#------------------------------------------------------------------------------#
######## 6. CG content
qc_gc_content <- qc_metadata[,c("sample_id", "avg_gc", 
                                "tissue_type")] %>% distinct(sample_id, .keep_all = T)
qc_gc_content$tissue_type <- factor(qc_gc_content$tissue_type,levels = c("normal","tumor"))

GC_vioplo <- ggviolin(qc_gc_content, "tissue_type", "avg_gc",
                      fill = "tissue_type",add.params = list(fill = "white"), alpha = 0.5,
                      palette = mypal[c(5,4)],
                      add = "boxplot",
                      ylab = "% GC", xlab = "Tissue Type")+ mytheme + 
  theme(legend.position = "none"); GC_vioplo

ggsave("./results/charts/CG_content_vio_ForSamples_RNAseq.png",width=6,height=4,dpi=300)
aggregate(qc_gc_content$avg_gc, by = list(qc_gc_content$tissue_type), FUN = mean)

panel.def <- ggarrange(MR_histogram,Ins_denplo_area, GC_vioplo, nrow = 1, widths = c(1,1,1), 
                       labels = c("d","e","f"))+theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.def
ggsave("./results/charts/panel def.pdf",panel.def,width=16,height=6,dpi=300)