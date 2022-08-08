## TNBC2019 RNAseq - exploration of insert size and RIN (GC Content)
## Author：Qingwang Chen
## Date：2022-5-18
## Version：1.0

# Set work directory and import libraries
rm(list = ls (all = TRUE))

source("./scripts/libraries.R")
source("./scripts/parameter.R")

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

######## 5. Insert size
qc_ins_size <- qc_metadata[,c("sample_id", "median_insert_size", 
                              "tissue_type")] %>% distinct(sample_id, .keep_all = T)
qc_ins_size_batch <- merge(qc_ins_size,RNAseq_batch_info,by="sample_id")

## plot
Ins_denplo_area <- ggdensity(qc_ins_size, x = "median_insert_size",
                             add = "median", rug = TRUE, 
                             color = "tissue_type", fill = "tissue_type", alpha = 0.5, 
                             palette = mypal[c(5,4)],xlab = "Insert Size", ylab = "Density", 
                             legend.title = "Type")+ mytheme;Ins_denplo_area

ggsave("./results/charts/Insert_sizeForSamples_RNAseq.png",width=6,height=4,dpi=300)

aggregate(qc_ins_size$median_insert_size, by = list(qc_ins_size$tissue_type), FUN = median)

## batch
qc_ins_size_batch$batch <- as.factor(qc_ins_size_batch$batch)
Ins_vioplo <- ggviolin(qc_ins_size_batch, "batch", "median_insert_size",
                      fill = "tissue_type", alpha = 0.5,
                      palette = mypal[c(5,4)],
                      add = "boxplot",
                      ylab = "Insert Size", xlab = "Batch")+ mytheme; Ins_vioplo
ggsave("./results/charts/Insert_sizeForBatches_RNAseq_violin.png",width=5,height=4,dpi=300)

Ins_denplo_area <- ggdensity(qc_ins_size_batch, x = "median_insert_size",
                             add = "median", rug = TRUE, 
                             color = "batch", fill = "batch", alpha = 0.5, 
                             palette = mypal[c(1,2,3)],xlab = "Insert Size", ylab = "Density", 
                             legend.title = "Batch")+ mytheme ;Ins_denplo_area
ggsave("./results/charts/Insert_sizeForBatches_RNAseq.png",width=6,height=4,dpi=300)


## inner batch
p1 <- ggdensity(qc_ins_size_batch[qc_ins_size_batch$batch=="batch1",], x = "median_insert_size",
                             add = "median", rug = TRUE, 
                             color = "tissue_type", fill = "tissue_type", alpha = 0.5, 
                             palette = mypal[c(5,4)],xlab = "Insert Size", ylab = "Density", 
                             legend.title = "Type")+ mytheme ;p1
ggsave("./results/charts/Insert_sizeForBatch1_RNAseq.png",width=6,height=4,dpi=300)

p2 <- ggdensity(qc_ins_size_batch[qc_ins_size_batch$batch=="batch2",], x = "median_insert_size",
                             add = "median", rug = TRUE, 
                             color = "tissue_type", fill = "tissue_type", alpha = 0.5, 
                             palette = mypal[c(5,4)],xlab = "Insert Size", ylab = "Density", 
                             legend.title = "Type")+ mytheme ;p2
ggsave("./results/charts/Insert_sizeForBatch2_RNAseq.png",width=6,height=4,dpi=300)

p3 <- ggdensity(qc_ins_size_batch[qc_ins_size_batch$batch=="batch3",], x = "median_insert_size",
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

qc_gc_content_batch <- merge(qc_gc_content,RNAseq_batch_info,by="sample_id")


GC_vioplo <- ggviolin(qc_gc_content, "tissue_type", "avg_gc",
                      fill = "tissue_type",add.params = list(fill = "white"), alpha = 0.5,
                      palette = mypal[c(5,4)],
                      add = "boxplot",
                      ylab = "% GC", xlab = "Tissue Type")+ mytheme + 
  theme(legend.position = "none"); GC_vioplo

ggsave("./results/charts/CG_content_vio_ForSamples_RNAseq.png",width=6,height=4,dpi=300)
aggregate(qc_gc_content$avg_gc, by = list(qc_gc_content$tissue_type), FUN = mean)

## acrossbatch
qc_gc_content_batch$batch <- factor(qc_gc_content_batch$batch)

GC_vioplo <- ggviolin(qc_gc_content_batch, "batch", "avg_gc",
                      fill = "tissue_type", alpha = 0.5,
                      palette = mypal[c(5,4)],
                      add = "boxplot",
                      ylab = "% GC", xlab = "Batch")+ mytheme; GC_vioplo
ggsave("./results/charts/gc_content_ForBatches_RNAseq.png",width=6,height=4,dpi=300)

## import RIN value 
RNAseq_RIN <- read.csv("./data/TNBC_RNASeq_RIN.csv",header=T)
qc_gc_content_batch_RIN <- merge(qc_gc_content_batch,RNAseq_RIN,by="sample_id")

qc_gc_content_batch1_RIN <- qc_gc_content_batch_RIN[qc_gc_content_batch_RIN$batch=="batch1",]

RIN_vioplo <- ggviolin(qc_gc_content_batch1_RIN, "type", "RIN",
                      fill = "tissue_type", alpha = 0.5,
                      palette = mypal[c(5,4)],
                      add = "boxplot",
                      ylab = "RIN", xlab = "Tissue Type")+ mytheme; RIN_vioplo
ggsave("./results/charts/RIN_ForSamples_RNAseq.png",width=4,height=4,dpi=300)

## further exploration
fpkm <- fread("./data/tnbc2019_448samples_dir_fpkm.csv") %>% as.data.frame()
row.names(fpkm) <- fpkm$gene_id; fpkm$gene_id <- NULL
count <- readRDS("./Rdata/FUSCCTNBC_RNAseqShi_Counts_Unfiltered_448x56858_200521.rds")
count_col <- colnames(count) %>% as.data.frame(); colnames(count_col) <- "name"
count_col$sample <- sapply(strsplit(as.character(count_col$name),"_"),function(x){x[1]})
count_col$type <- sapply(strsplit(as.character(count_col$name),"_"),function(x){x[2]})
count_col$type <- gsub("rep","TT",count_col$type)
count_col$type[which(is.na(count_col$type))] <- "TT"
count_col$sample_id <- paste0(count_col$sample,"_",count_col$type)
colnames(count) <- count_col$sample_id

count_batch1 <- count[which( colnames(count) %in% qc_gc_content_batch1_RIN$sample_id),]

expd <- ifelse(count_batch1>3,1,0)
expdN <- data.frame(colSums(expd))
expdratio <- apply(count,2,function(x){x/sum(x)})

nn2 <- apply(expdratio,2,function(x){sum(x[which(x>0.02)])})

qc_gc_content_batch1_RIN[qc_gc_content_batch1_RIN$sample_id %in% names(which(nn2>0.15)),]

qc_gc_content_batch1_RIN$highexpdo <- nn2[match(qc_gc_content_batch1_RIN$sample_id,names(nn2))]
plot(qc_gc_content_batch1_RIN$highexpdo,qc_gc_content_batch1_RIN$RIN)
cor.test(qc_gc_content_batch1_RIN$highexpdo,qc_gc_content_batch1_RIN$RIN)

qc_gc_content_batch1_RIN$highexpdo <- qc_gc_content_batch1_RIN$highexpdo*100
highexpdo_vioplo <- ggviolin(qc_gc_content_batch1_RIN, "type", "highexpdo",
                       fill = "tissue_type", alpha = 0.5,
                       palette = mypal[c(5,4)],
                       add = "boxplot",
                       ylab = "Top Gene (>2%)/All Gene (%)", xlab = "Tissue Type")+ mytheme; highexpdo_vioplo
ggsave("./results/charts/highexpdo_vioplo_RNAseq.png",width=4,height=4,dpi=300)



#### gene GC content
library("biomaRt")
ensembl=useMart("ensembl")

datasets <- listDatasets(ensembl)
ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)

attributes <-listAttributes(ensembl)
head(attributes)

ssids = getBM(attributes=c('ensembl_gene_id', 'description', 'external_gene_name', 
                           "percentage_gene_gc_content", "transcript_gencode_basic", 
                           "gene_biotype", "transcript_biotype"), 
              filters='ensembl_gene_id', values=rownames(fpkm), mart=ensembl)
# saveRDS(ssids,"./Rdata/GC_content_ensembl_gene_id.rds")

## filter diffrent samples
normal_samples <- as.character(qc_gc_content_batch1_RIN$sample_id[qc_gc_content_batch1_RIN$type=="PT"])
tumor_samples <- as.character(qc_gc_content_batch1_RIN$sample_id[qc_gc_content_batch1_RIN$type=="TT"])

normal_genes<-c()
for ( i in 1:length(normal_samples)){
  normal_genes <- c(normal_genes,rownames(expdratio)[expdratio[,colnames(expdratio) %in% normal_samples[i]]>0.02])
}
normal_genes<-unique(normal_genes)
plot(density(ssids[ssids$external_gene_name %in% normal_genes,"percentage_gene_gc_content"]))

tumor_genes<-c()
for ( i in 1:length(tumor_samples)){
  tumor_genes <- c(tumor_genes,rownames(expdratio)[expdratio[,colnames(expdratio) %in% tumor_samples[i]]>0.02])
}
tumor_genes<-unique(tumor_genes)
plot(density(ssids[ssids$external_gene_name %in% tumor_genes,"percentage_gene_gc_content"]))

dd_normal <-ssids[ssids$external_gene_name %in% normal_genes,c("ensembl_gene_id","percentage_gene_gc_content")]
dd_normal$type<-"top genes normal"

dd_tumor <-ssids[ssids$external_gene_name %in% tumor_genes,c("ensembl_gene_id","percentage_gene_gc_content")]
dd_tumor$type<-"top genes tumor"

dd_all<-ssids[,c("ensembl_gene_id","percentage_gene_gc_content")]
dd_all$type<-"all genes"

dd_merge <- rbind(dd_normal,dd_tumor,dd_all)
dd_merge$type <- factor(dd_merge$type,levels=c("top genes normal","top genes tumor","all genes"),ordered=T)

Top_gene_GC_Content <- ggplot(dd_merge,aes(x=percentage_gene_gc_content,fill=type,color=type))+
  geom_density(alpha=0.3)+
  theme_bw()+
  mytheme+
  scale_fill_npg(name="")+
  scale_color_npg(name="")+
  labs(x="Gene GC Conetent",y="Density")+
  theme(legend.text = element_text(size=11),
        legend.position = c(0.20,0.7));Top_gene_GC_Content
ggsave("./results/charts/Top_gene_GC_Content_RNAseq.png",width=6,height=6,dpi=300)
