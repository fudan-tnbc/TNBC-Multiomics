## Title: TNBC Bridge RNAseq Analysis##
## Author: Qingwang Chen             ##
## Date: 2022-05-20                  ##  
## Version: V3.0                     ##
#######################################

# Set work directory and import libraries
rm(list = ls (all = TRUE))
source("./scripts/libraries.R")
source("./scripts/parameter.R")

# ----------------------Subtyping---------------------- #
# 1. Data process
load("./data/old_data_old_pipeline_rnaseq/FUSCCTNBC_RNAseqShi.Complete_log2_448x45308_V15_190209.Rdata")
fpkm_old_pipe <- FUSCCTNBC_RNAseqShi.Complete_log2
fpkm_old_pipe <- 2^fpkm_old_pipe-1
filtered_fpkm_old_pipe <- fpkm_old_pipe[apply(fpkm_old_pipe, MARGIN=1, 
                                              FUN=function(x) {(sum(x==0)/length(x))<=0.3}),]
filtered_fpkm_old_pipe_PT <- filtered_fpkm_old_pipe[,grep("PT",colnames(filtered_fpkm_old_pipe))]
filtered_fpkm_old_pipe_TT <- filtered_fpkm_old_pipe[,-grep("PT",colnames(filtered_fpkm_old_pipe))]

fpkm_new_pipeline <- readRDS("./data/FUSCCTNBC_new_fpkm_448x56858_220302.rds")
## filer the genes whose tpm=0 in 30% samples
# filtered_fpkm_new_pipe <- fpkm_new_pipeline[apply(fpkm_new_pipeline, MARGIN=1, 
#                                                 FUN=function(x) {(sum(x==0)/length(x))<=0.3}),]
## save the filtered tpm witn annotation
# saveRDS(filtered_tpm_new_pipe_log2,"./data/FUSCCTNBC_new_filtered_tpm_log2_448x25589_220520.rds")
filtered_fpkm_new_pipe_PT <- filtered_fpkm_new_pipe[,grep("PT",colnames(filtered_fpkm_new_pipe))]
filtered_fpkm_new_pipe_TT <- filtered_fpkm_new_pipe[,-grep("PT",colnames(filtered_fpkm_new_pipe))]

## method from Dr. Ma (method before)
### filter the protein coding genes
Gene_annotation_TNBC <- read.csv("./data/RNA_seq_Gene_annotation_FUSCCTNBC.csv", header = T)
row.names(Gene_annotation_TNBC) <- Gene_annotation_TNBC$GeneSymbol
Gene_annotation_TNBC <- Gene_annotation_TNBC[which(!is.na(Gene_annotation_TNBC$Category)), ]
ProteinCoding <- Gene_annotation_TNBC$GeneSymbol[Gene_annotation_TNBC$Category == "protein_coding"]

### Expression dataset for subsequent test
# filtered_fpkm_old_pipe_TT 
logexpr_old_TT <- log2(filtered_fpkm_old_pipe_TT +1)
Tes_Expr <- logexpr_old_TT[ProteinCoding, ]
dim(Tes_Expr)

### Select genes for clustering (sd top 2000)
SDs <- NULL
for(i in 1:nrow(Tes_Expr))
{
  SDs <- c(SDs, sd(as.numeric(Tes_Expr[i,])))
  # if(i/1000 == round(i/1000)) {print(i/nrow(Tes_Expr))}
}
names(SDs) <- row.names(Tes_Expr)
SDs        <- sort(SDs, decreasing = T)
Clustering_Genes <- names(SDs[1:2000])
### 2000 genes

### Expression dataset for clustering
Tes_Expr <- Tes_Expr[Clustering_Genes, ]
dim(Tes_Expr)

### clustering by k-means
set.seed(1234)
km1 <- kmeans(t(Tes_Expr),4,iter.max = 1000,nstart = 100)
table(km1$cluster)
Tes_Clusters <- km1$cluster

#### Rename clusters
Tes_Clusters[Tes_Clusters == "1"] <- "BLIS"  ; Tes_Clusters[Tes_Clusters == "2"] <- "LAR"
Tes_Clusters[Tes_Clusters == "3"] <- "IM" ; Tes_Clusters[Tes_Clusters == "4"] <- "MES"

#### Outputs
write.csv(Tes_Clusters, "./results/tables/FUSCCTNBC_old_subtypes.csv")

###################--------------new pipeline-----------------##################
# filtered_fpkm_new_pipe_TT
logexpr_new_TT <- log2(filtered_fpkm_new_pipe_TT+1) %>% as.data.frame()
Tes_Expr <- logexpr_new_TT[ProteinCoding, ]
dim(Tes_Expr)

### Expression dataset for clustering
Tes_Expr <- Tes_Expr[Clustering_Genes, ] %>% na.omit()
dim(Tes_Expr)
### 1854 genes(consistent with old pipeline)

### clustering by k-means
set.seed(1234)
km1 <- kmeans(t(Tes_Expr),4,iter.max = 1000,nstart = 100)
table(km1$cluster)
Tes_Clusters <- km1$cluster

#### Rename clusters
Tes_Clusters[Tes_Clusters == "1"] <- "BLIS"  ; Tes_Clusters[Tes_Clusters == "2"] <- "LAR"
Tes_Clusters[Tes_Clusters == "3"] <- "IM" ; Tes_Clusters[Tes_Clusters == "4"] <- "MES"

#### Outputs
write.csv(Tes_Clusters, "./results/tables/FUSCCTNBC_new_subtypes.csv")
rm(list=ls())

# ----------------------Subtyping Column Table---------------------- #
#1. data import
FUSCCTNBC_old_subtypes <- read.csv("./results/tables/FUSCCTNBC_old_subtypes.csv")
colnames(FUSCCTNBC_old_subtypes) <- c("patient","subtype_old")
table(FUSCCTNBC_old_subtypes$subtype_old)
FUSCCTNBC_new_subtypes <- read.csv("./results/tables/FUSCCTNBC_new_subtypes.csv")
colnames(FUSCCTNBC_new_subtypes) <- c("patient","subtype_new")
table(FUSCCTNBC_new_subtypes$subtype_new)
FUSCCTNBC_new_subtypes$patient <- gsub("_rep|TPM.|_TT|_TP","",FUSCCTNBC_new_subtypes$patient)
FUSCCTNBC_subtypes <- merge(FUSCCTNBC_old_subtypes,FUSCCTNBC_new_subtypes,by="patient")
rm(FUSCCTNBC_old_subtypes,FUSCCTNBC_new_subtypes)
write.csv(FUSCCTNBC_subtypes,"./results/tables/FUSCCTNBC_subtypes.csv")

# 2. Compare the results
## check the subtype status
result = 0
patients = c()
for(i in 1:length(FUSCCTNBC_subtypes$patient)){
  temp <- sum(FUSCCTNBC_subtypes$subtype_old[i] %in% FUSCCTNBC_subtypes$subtype_new[i])
  if(temp == 1){
    check <- FUSCCTNBC_subtypes$patient[i]
    patients <- append(patients,check)
  }
  result = result+temp
}
### 340 patients matched; 20 patients unmatched

df <- xtabs(~ subtype_old + subtype_new, data = FUSCCTNBC_subtypes)
tabel1 <- as.data.frame(df) %>% melt()
library(plyr)
tabel1 <- ddply(tabel1, "subtype_old", transform,
                percent_old = value / sum(value) * 100)
tabel1 <- ddply(tabel1, "subtype_new", transform,
                percent_new = value / sum(value) * 100)

## plot 
library(ggplot2)
lev <- c("BLIS","IM","MES","LAR")
tabel1$subtype_old <- factor(tabel1$subtype_old,levels = rev(lev))
tabel1$subtype_new <- factor(tabel1$subtype_new,levels = lev)
p1 <- ggplot(tabel1, aes(x=subtype_old, y=percent_old/100, fill=subtype_new)) +
  geom_bar(stat="identity", colour="black")+
  guides(fill=guide_legend(reverse=TRUE))+
  scale_fill_manual(values = c("#ef4922","#7bc242","#3baade","#8f80ba"))+theme_classic()+ 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size=12, face="bold"),
        axis.title.x = element_text(angle = 0, hjust = 0.5, size=14, face="bold"),
        panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        legend.position="none",
        plot.margin=unit(rep(1,4),'lines')) + 
  scale_y_continuous(expand = c(0,0),position = "right")+
  ylab("Proportion")+coord_flip();p1
# +scale_y_reverse()
ggsave("./results/charts/subtype_old_new.pdf",p1,width = 4,height = 4,dpi = 300,bg = "transparent")


tabel1$subtype_old <- factor(tabel1$subtype_old,levels = rev(lev))
tabel1$subtype_new <- factor(tabel1$subtype_new,levels = rev(lev))
p2 <- ggplot(tabel1, aes(x=subtype_new, y=percent_new/100, fill=subtype_old)) +
  geom_bar(stat="identity", colour="black")+
  scale_fill_manual(values = c("#ef4922","#7bc242","#3baade","#8f80ba"))+theme_classic()+ 
  ylab("Proportion")+
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, size=12, face="bold"),
        axis.title.y = element_text(angle = 90, hjust = 0.5, size=14, face="bold"),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        legend.position="none")+
  scale_y_continuous(expand = c(0,0),position = "left");p2

ggsave("./results/charts/subtype_new_old.pdf",p2,width = 4,height = 4,dpi = 300,bg = "transparent")

## remove all
rm(list=ls())


