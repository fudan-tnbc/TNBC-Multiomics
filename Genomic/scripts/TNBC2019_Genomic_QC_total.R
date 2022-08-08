## TNBC2019 WES QC
## Author：Qingwang Chen; Yaqing Liu
## Date：2022-06-07
## Verision：3.0

# 设置工作目录并导入相关库
rm(list = ls (all = TRUE))
source("./scripts/libraries.R")
source("./scripts/parameter.R")
source("./scripts/theme_color.R")

######### meatadata import
qc_metadata <- read.csv("./data/TNBC2019_WESQC_batch_metadata.csv",header=T)

#------------------------------------------------------------------------------#
######### 1. Mapping ratio
qc_MR <- qc_metadata[,c("sample_id", "percentage_aligned",
                        "tissue_type")] %>% distinct(sample_id, .keep_all = T)
qc_MR$species <- "Human"
min(qc_MR$percentage_aligned) # 97.63026

#------------------------------------------------------------------------------#
######### 2. FastQ Screen (Contamination)
# colnames(qc_metadata)[grep("percentage",colnames(qc_metadata))]
fastqscreen <- qc_metadata[,c("fastqscreen_id", "Mouse",  "Yeast", "EColi", 
                              "Virus", "Vector")]

fastqscreen$fastqscreen_id <- gsub("_screen","",fastqscreen$fastqscreen_id,fix=T)
fastqscreen$sample_id <- sapply(strsplit(as.character(fastqscreen$fastqscreen_id),"_"),
                                function(x){paste(x[1],x[2],sep="_")})
fastqscreen$sample_id <- gsub("R1|R2|rep","TT",fastqscreen$sample_id)
fastqscreen[,2:6] <- fastqscreen[,2:6]/100
human_MR <- qc_MR[,c(1,2)]
MR <- merge(human_MR,fastqscreen,by="sample_id")  
rownames(MR)<-MR$fastqscreen_id 
MR$fastqscreen_id<-NULL; MR$sample_id <- NULL
colnames(MR)[1] <- "Human"

fastqscreen.forplot<-melt(as.matrix(MR))
colnames(fastqscreen.forplot)<-c("file","source","value")
fastqscreen.forplot$sample<-sapply(strsplit(as.character(fastqscreen.forplot$file),"_"),
                                   function(x){paste(x[1],x[2],sep="_")})
fastqscreen.forplot$type <- ifelse(fastqscreen.forplot$source == "Human","Homo sapiens","Other genomes")

###### plot
Con_boxplo <- ggboxplot(fastqscreen.forplot, x = "source", y = "value",
                        add = c("median","jitter"),
                        color = "source",
                        palette = mypal,
                        xlab = "", ylab = "% Mapping Ratio",
                        legend.title = "contamination") + mytheme+ 
  theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1))+
  ylim(0, 100)+facet_grid(.~type,scales = "free",space = "free"); Con_boxplo
ggsave("./results/charts/contaminForSamples_RNAseq.png",width=6,height=4,dpi=300)

panel.b <- ggarrange(Con_boxplo, nrow = 1, labels = c("b")) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.b
ggsave("./results/charts/panel b.pdf",panel.b,width=8,height=6,dpi=300)

#------------------------------------------------------------------------------#
######### 3. Coverage
# subset the coverage part of QC metadata
qc_coverage <- qc_metadata[,c("sample_id", "mean_coverage",
                              "tissue_type")]  %>% distinct(sample_id, .keep_all = T)

## plot
Cov_histogram <- gghistogram(qc_coverage, x = "mean_coverage",
                             add = "median", rug = TRUE, 
                             color = "tissue_type", fill = "tissue_type", alpha = 0.5,
                             palette = mypal[c(5,4)], xlab = "Coverage", ylab = "Sample counts", 
                             legend.title = "Type") + mytheme + 
  theme(legend.position = "none", legend.margin=unit(-.05,"cm")); Cov_histogram

## export
ggsave("./results/charts/coverageForSamples_WES.png",width=4,height=4,dpi=300)

aggregate(qc_coverage$mean_coverage, by = list(qc_coverage$tissue_type), FUN = max)
# rm(list = ls(pattern = "coverage|histogram"))

panel.c <- ggarrange(Cov_histogram, nrow = 1, labels = c("c")) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.c
ggsave("./results/charts/panel c.pdf",panel.c,width=8,height=6,dpi=300)

#------------------------------------------------------------------------------#
######## 4. Duplicate
qc_duplicate <- qc_metadata[,c("sample_id", "duplicate_picard", 
                               "tissue_type")] %>% distinct(sample_id, .keep_all = T)
## plot
qc_duplicate$tissue_type <- factor(qc_duplicate$tissue_type,levels = c("normal","tumor"))

Dup_vioplo <- ggviolin(qc_duplicate, "tissue_type", "duplicate_picard",
                       fill = "tissue_type", palette = mypal[c(5,4)], alpha = 0.5, 
                       add.params = list(fill = "white"),
                       add = "boxplot",
                       ylab = "% Duplication", xlab = "Tissue Type") + mytheme +
  theme(legend.position = "none"); Dup_vioplo

ggsave("./results/charts/duplicate_vio_ForSamples.png",width=4,height=4,dpi=300)

panel.d <- ggarrange(Dup_vioplo, nrow = 1, labels = c("d")) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.d
ggsave("./results/charts/panel d.pdf",panel.d,width=8,height=6,dpi=300)

panel.bcd <- ggarrange(panel.b, panel.c, panel.d, nrow = 1, widths = c(2, 1, 1)) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.bcd
ggsave("./results/charts/panel bcd.pdf",panel.bcd,width=16,height=4,dpi=300)

#------------------------------------------------------------------------------#
######## 5. Insert size
qc_ins_size <- qc_metadata[,c("sample_id", "median_insert_size", "batch",
                              "tissue_type")] %>% distinct(sample_id, .keep_all = T)
## plot
qc_ins_size$batch <- as.factor(qc_ins_size$batch)

Ins_boxplo <- ggboxplot(qc_ins_size, "batch", "median_insert_size",
                        fill = "tissue_type", alpha = 0.5,
                        palette = mypal[c(5,4)],
                        ylab = "Insert Size", xlab = "Batch") + mytheme +
  theme(legend.position = "none")+
  stat_summary(fun.y=median, geom="point", shape=18, size=3);Ins_boxplo
ggsave("./results/charts/Insert_sizeForBatches_WES_box.png",width=12,height=4,dpi=300)

# Ins_vioplo <- ggviolin(qc_ins_size, "batch", "median_insert_size",
#                        fill = "tissue_type", alpha = 0.5,
#                        palette = mypal[c(5,4)],
#                        add = "boxplot",
#                        ylab = "Insert Size", xlab = "")+ mytheme +
#                        theme(legend.position = "none"); Ins_vioplo
# ggsave("./results/charts/Insert_sizeForBatches_WES_violin.png",width=10,height=5,dpi=300)

# Ins_denplo_area <- ggdensity(qc_ins_size, x = "median_insert_size",
#                              add = "median", rug = TRUE, 
#                              color = "tissue_type", fill = "tissue_type", alpha = 0.5, 
#                              palette = mypal[c(5,4)],xlab = "Insert Size (bp)", ylab = "Density", 
#                              legend.title = "Type")+ mytheme + 
#   theme(legend.position = "none");Ins_denplo_area
# ggsave("./results/charts/Insert_sizeForSamples_WES.png",width=6,height=4,dpi=300)

aggregate(qc_ins_size$median_insert_size, by = list(qc_ins_size$tissue_type), FUN = mean)
# normal 195.4767; tumor 185.5484
aggregate(qc_ins_size$median_insert_size, by = list(qc_ins_size$tissue_type), FUN = sd)
# normal 15.66483; tumor 13.59844

panel.e <- ggarrange(Ins_boxplo, nrow = 1, labels = c("e")) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.e
ggsave("./results/charts/panel e.pdf",panel.e,width=12,height=4,dpi=300)

#------------------------------------------------------------------------------#
######## 6. CG content
qc_gc_content <- qc_metadata[,c("sample_id", "avg_gc", "batch", 
                                "tissue_type")] %>% distinct(sample_id, .keep_all = T)
qc_gc_content$tissue_type <- factor(qc_gc_content$tissue_type,levels = c("normal","tumor"))
qc_gc_content$batch <- as.factor(qc_gc_content$batch)

GC_boxplo <- ggboxplot(qc_gc_content, "batch", "avg_gc",
                       fill = "tissue_type", alpha = 0.5,
                       palette = mypal[c(5,4)],
                       ylab = "% GC", xlab = "Batch") + mytheme +
  theme(legend.position = "none")+
  stat_summary(fun.y=median, geom="point", shape=18, size=3);GC_boxplo
ggsave("./results/charts/gc_content_sizeForBatches_WES_box.png",width=5,height=4,dpi=300)

# GC_vioplo <- ggviolin(qc_gc_content, "batch", "avg_gc",
#                        fill = "tissue_type", alpha = 0.5,
#                        palette = mypal[c(5,4)],
#                        add = "boxplot",
#                        ylab = "% GC", xlab = "Batch")+ mytheme +
#   theme(legend.position = "none"); GC_vioplo
# ggsave("./results/charts/GC_ForBatches_WES_violin.png",width=10,height=5,dpi=300)

# GC_vioplo <- ggviolin(qc_gc_content, "tissue_type", "avg_gc",
#                       fill = "tissue_type",add.params = list(fill = "white"), alpha = 0.5,
#                       palette = mypal[c(5,4)],
#                       add = "boxplot",
#                       ylab = "% GC", xlab = "Tissue Type")+ mytheme + 
#   theme(legend.position = "none"); GC_vioplo
# ggsave("./results/charts/CG_content_vio_ForSamples_WES.png",width=6,height=4,dpi=300)

aggregate(qc_gc_content$avg_gc, by = list(qc_gc_content$tissue_type), FUN = mean)
# normal 51.22687; tumor 53.28124
aggregate(qc_gc_content$avg_gc, by = list(qc_gc_content$tissue_type), FUN = sd)
# normal 0.7648324; tumor 2.4243874

panel.f <- ggarrange(GC_boxplo, nrow = 1, labels = c("f")) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.f
ggsave("./results/charts/panel f.pdf",panel.f,width=12,height=4,dpi=300)

panel.ef <- ggarrange(Ins_boxplo,GC_boxplo, nrow = 2, heights = c(1,1), 
                     labels = c("e","f"))+theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.ef
ggsave("./results/charts/panel ef.pdf",panel.ef,width=12,height=6,dpi=300)

### FASTQC
img <- readPNG("./data/fastqc.png")
# plot with picture as layer
fastqc_p <- ggplot(mapping = aes(1:12, 1:12)) +
  annotation_raster(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
  theme_bw() + xlab("")+ylab("")+
  theme(panel.border = element_blank());fastqc_p

panel.a <- ggarrange(fastqc_p,nrow = 1,labels = c("a"))+
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm")); panel.a

#------------------------------------------------------------------------------#
######## 7. OncoScan QC
OncoScan_CNV_QC <- read.csv("./data/OncoScan_CNV_QC220524.csv",header = T)
OncoScan_CNV_QC$ndSNPQC_QC <- capitalize(paste0(OncoScan_CNV_QC$tissue_type,"_",OncoScan_CNV_QC$ndSNPQC_status))
OncoScan_CNV_QC$ndSNPQC_QC <- factor(OncoScan_CNV_QC$ndSNPQC_QC,levels = c("Tumor_Pass","Tumor_Fail","Normal_Pass","Normal_Fail"))
OncoScan_CNV_QC$MAPD_QC <- capitalize(paste0(OncoScan_CNV_QC$tissue_type,"_",OncoScan_CNV_QC$MAPD_status))
OncoScan_CNV_QC$MAPD_QC <- factor(OncoScan_CNV_QC$MAPD_QC,levels = c("Tumor_Pass","Tumor_Fail","Normal_Pass","Normal_Fail"))

ndSNPQC_histogram <- ggplot(OncoScan_CNV_QC)+
  geom_histogram(aes(x=ndSNPQC,fill=ndSNPQC_QC,color=ndSNPQC_QC),alpha=0.5,bins = 50)+
  scale_fill_manual(values = c("#001E6C","#035397",'#A64B2A','#AC7D88'))+
  scale_color_manual(values = c("#001E6C","#035397",'#A64B2A','#AC7D88'))+
  theme_classic() + geom_vline(aes(xintercept=26), colour="#990000", linetype="dashed") + 
                          theme(plot.background = element_rect(colour = "white"),
                          plot.title = element_text(hjust=0.5,size=12,face='bold'),
                          axis.title = element_text(size=12),
                          axis.text = element_text(colour= "black",size=11),
                          strip.text=element_text(size=11),
                          legend.title = element_text(size=10),
                          legend.position = "none");ndSNPQC_histogram

MAPD_histogram <- ggplot(OncoScan_CNV_QC)+
  geom_histogram(aes(x=MAPD,fill=MAPD_QC,color=MAPD_QC),alpha=0.5,bins = 50)+
  scale_fill_manual(values = c("#001E6C","#035397",'#A64B2A','#AC7D88'))+
  scale_color_manual(values = c("#001E6C","#035397",'#A64B2A','#AC7D88'))+
  theme_classic() + geom_vline(aes(xintercept=0.3), colour="#990000", linetype="dashed")+ 
                          theme(plot.background = element_rect(colour = "white"),
                          plot.title = element_text(hjust=0.5,size=12,face='bold'),
                          axis.title = element_text(size=12),
                          axis.text = element_text(colour= "black",size=11),
                          strip.text=element_text(size=11),
                          legend.title = element_text(size=10),
                          legend.position = "none");MAPD_histogram

panel.gh <- ggarrange(ndSNPQC_histogram, MAPD_histogram, nrow = 1,labels = c("g","h"))+
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm")); panel.gh
ggsave("./results/charts/panel gh.pdf",panel.gh,width=6,height=4,dpi=300)

fig3 <- ggarrange(panel.a,panel.bcd,panel.e,panel.f,panel.gh,nrow = 5, widths = c(2,2,2,2,2))+
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));fig3
ggsave("./results/charts/fig3_0616.pdf",fig3,width=12,height=20,dpi=300)

######## Delete unrelated variables
rm(list = ls())

