## Title: TNBC Bridge RNAseq Analysis-- MES subtypig further analysis##
## Author: Qingwang Chen             ##
## Date: 2022-06-07                  ##  
## Version: V3.0                     ##
#######################################

# Set work directory and import libraries
rm(list = ls (all = TRUE))
source("./scripts/libraries.R")
source("./scripts/parameter.R")
library(scatterplot3d)
library(ggalt)

# ----------------------Subtyping---------------------- #
Tes_Clusters <- read.csv("./results/tables/FUSCCTNBC_old_subtypes.csv")
colnames(Tes_Clusters) <- c("patient","subtype")
### filter the protein coding genes
Gene_annotation_TNBC <- read.csv("./data/RNA_seq_Gene_annotation_FUSCCTNBC.csv", header = T)
row.names(Gene_annotation_TNBC) <- Gene_annotation_TNBC$GeneSymbol
Gene_annotation_TNBC <- Gene_annotation_TNBC[which(!is.na(Gene_annotation_TNBC$Category)), ]
ProteinCoding <- Gene_annotation_TNBC$GeneSymbol[Gene_annotation_TNBC$Category == "protein_coding"]

### Expression dataset for subsequent test
#### old expression dataset
load("./data/old_data_old_pipeline_rnaseq/FUSCCTNBC_RNAseqShi.Complete_log2_448x45308_V15_190209.Rdata")
fpkm_old_pipe <- FUSCCTNBC_RNAseqShi.Complete_log2
fpkm_old_pipe <- 2^fpkm_old_pipe-1
filtered_fpkm_old_pipe <- fpkm_old_pipe[apply(fpkm_old_pipe, MARGIN=1, 
                                              FUN=function(x) {(sum(x==0)/length(x))<=0.3}),]
filtered_fpkm_old_pipe_PT <- filtered_fpkm_old_pipe[,grep("PT",colnames(filtered_fpkm_old_pipe))]
filtered_fpkm_old_pipe_TT <- filtered_fpkm_old_pipe[,-grep("PT",colnames(filtered_fpkm_old_pipe))]

# filtered_fpkm_old_pipe_TT 
# logexpr_old_TT <- log2(filtered_fpkm_old_pipe_TT +1)
logexpr_old <- log2(filtered_fpkm_old_pipe +1)
# Tes_Expr <- logexpr_old_TT[ProteinCoding, ]
Tes_Expr <- logexpr_old[ProteinCoding, ]

dim(Tes_Expr)

### Select genes for clustering (sd top 2000)
# SDs <- NULL
# for(i in 1:nrow(Tes_Expr))
# {
#   SDs <- c(SDs, sd(as.numeric(Tes_Expr[i,])))
#   # if(i/1000 == round(i/1000)) {print(i/nrow(Tes_Expr))}
# }
# names(SDs) <- row.names(Tes_Expr)
# SDs        <- sort(SDs, decreasing = T)
# Clustering_Genes <- names(SDs[1:2000])
### 2000 genes
Clustering_Genes <- readRDS("./Rdata/Clustering_Genes_old.rds")

### Expression genes for clustering
Tes_Expr <- Tes_Expr[Clustering_Genes, ] %>% na.omit()
dim(Tes_Expr)
### 2000 genes(old pipeline)

### clustering by k-means
# set.seed(1234)
# km1 <- kmeans(t(Tes_Expr),4,iter.max = 1000,nstart = 100)
# table(km1$cluster)
# Tes_Clusters <- km1$cluster

### PCA
pca.all <- prcomp(t(Tes_Expr),retx=T)

### Create horizontal and vertical labels
# PC3comp<-summary(pca.all)$importance[2,3]*100
lev <- c("BLIS","IM","MES","LAR")
subtype_pal <- c("#ef4922","#7bc242","#3baade","#8f80ba")
Tes_Clusters$subtype <- factor(Tes_Clusters$subtype,levels = lev)
pca_cluster <- fviz_pca_ind(pca.all, 
                            axes = c(1, 2),
                            label = "none",
                            addEllipses = T,
                            geom.ind = c("point", "text"),
                            repel = T,
                            habillage = Tes_Clusters$subtype,
                            ellipse.level=0.95,
                            palette = subtype_pal)+
  theme(axis.text = element_text(size=16,face='bold',color='gray40'),
        axis.title = element_text(size=16,face='bold'),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16,color='gray40'),
        plot.subtitle = element_text(hjust=0.5,size=16),
        plot.title = element_text(hjust=0.5,size=16))+
  scale_shape_manual(values = c(15:19))+
  labs(x = sprintf("PC1 (%.2f%%)", summary(pca.all)$importance[2,1]*100),
       y = sprintf("PC2 (%.2f%%)", summary(pca.all)$importance[2,2]*100),
       title = "")+
  theme(aspect.ratio = 1/1);pca_cluster
#,title=paste('PCA-All ',ncol(log_filtered.expr.mat.ft),'samples',' N=',nrow(log_filtered.expr.mat.ft))



ggsave("./results/charts/PCA_all_subtype_old_448samples_2000genes.pdf",pca_cluster,width=6,height=4,dpi=300)

### total samples
# MES_unmatched <- readRDS("./Rdata/MES_unmatched.rds")
patients_unmatched <- readRDS("./Rdata/patients_unmatched.rds")

pcs <- pca.all$x %>% as.data.frame()
pcs <- pcs[,1:3]
pcs$patient <- rownames(pcs)

tnbc_mes_tumor_ratio <- read.csv("./results/tables/TNBC_mes_tumor_ratio.txt")
table(tnbc_mes_tumor_ratio$tumor_purity)
low_tp_mes_patients <- tnbc_mes_tumor_ratio$sample_id[which(tnbc_mes_tumor_ratio$tumor_purity<75)]
high_tp_mes_patients <- tnbc_mes_tumor_ratio$sample_id[which(tnbc_mes_tumor_ratio$tumor_purity>=75)]


for_label <- pcs[pcs$patient %in% tnbc_mes_tumor_ratio$sample_id,]
for_label$subtype <- "MES"
colnames(for_label)[4] <- "sample_id"

for_label <- merge(for_label,tnbc_mes_tumor_ratio,by="sample_id")
saveRDS(for_label,"./Rdata/pca_all_samples_360patients_MES_annotated_for_label.rds")


pcs$patient <- sapply(strsplit(as.character(rownames(pcs)),"\\."),function(x){x[1]})
pcs$type <- sapply(strsplit(as.character(rownames(pcs)),"\\."),function(x){x[2]})
pcs$type[is.na(pcs$type)] <- "TT"
pcs$type <- gsub("TT","tumor",pcs$type); pcs$type <- gsub("PT","normal",pcs$type)
pcs$compare <- ifelse(pcs$patient %in% patients_unmatched, "inconsistent", "consistent")
pcs$class <- paste0(pcs$type,"-subtype_",pcs$compare)

pcs_an <- merge(pcs,Tes_Clusters,by="patient")
pcs_an$type <- factor(pcs_an$type,levels = c("tumor","normal"))

write.csv(pcs_an,"./results/tables/pca_all_samples_360patients_MES_annotated.csv",quote = F,row.names = F)
saveRDS(pca.all,"./Rdata/pca_all_samples_360patients_MES_annotated.rds")


pcs_an$subtype <- factor(pcs_an$subtype,levels = lev)
pca_cluster_total_old <- ggplot(pcs_an,aes(x=PC1,y=PC2,color=subtype))+
  geom_point(aes(shape=class),size=2.5)+
  # geom_encircle(aes(color=type,fill=type), expand=0.02, alpha=.2)+
  theme_few()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.position = "right",
        plot.subtitle = element_text(hjust=0.5,size=14),
        plot.title = element_text(hjust=0.5,size=16))+
  scale_color_manual(values = subtype_pal)+
  scale_fill_manual(values = subtype_pal)+
  scale_shape_manual(values = c(0,15,2,17))+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  guides(shape=guide_legend(override.aes = list(size=3)))+
  labs(x = sprintf("PC1 (%.2f%%)", summary(pca.all)$importance[2,1]*100),
       y = sprintf("PC2 (%.2f%%)", summary(pca.all)$importance[2,2]*100),
       title = "");pca_cluster_total_old
ggplotly(pca_cluster_total_old)

pca_cluster_total_old+geom_text_repel(data = for_label, aes(PC1,PC2,
                  color=factor(subtype),
                  label=tumor_purity),
                  arrow = arrow(length=unit(0.01, "npc")),
                  max.overlaps = 30)
ggsave("./results/charts/PCA_all_subtype_old_448samples(all)_2000genes_TP.pdf",width=10,height=7,dpi=300)

ggsave("./results/charts/PCA_all_subtype_old_448samples(all)_2000genes.pdf",pca_cluster_total,width=7,height=4,dpi=300)
ggsave("./results/charts/PCA_all_subtype_old_448samples(all)_2000genes_nolegend.pdf",pca_cluster_total,width=4,height=4,dpi=300)

###################--------------new pipeline-----------------##################
Tes_Clusters <- read.csv("./results/tables/FUSCCTNBC_new_subtypes.csv")
colnames(Tes_Clusters) <- c("patient","subtype")

#### new expression dataset
fpkm_new_pipe <- readRDS("./Rdata/FUSCCTNBC_new_fpkm_448x56858_220302.rds")
#### unfiltered
fpkm_new_pipe_PT <- fpkm_new_pipe[,grep("PT",colnames(fpkm_new_pipe))]
fpkm_new_pipe_TT <- fpkm_new_pipe[,-grep("PT",colnames(fpkm_new_pipe))]

# filtered_fpkm_new_pipe_TT
# logexpr_new_TT <- log2(fpkm_new_pipe_TT+1)
logexpr_new <- log2(fpkm_new_pipe+1)

# Tes_Expr <- logexpr_new_TT[ProteinCoding, ]
Tes_Expr <- logexpr_new[ProteinCoding, ]

dim(Tes_Expr)

### Expression genes for clustering
# Tes_Expr <- logexpr_new_TT[which(rownames(Tes_Expr) %in% Clustering_Genes), ] %>% na.omit()
Tes_Expr <- Tes_Expr[Clustering_Genes, ] %>% na.omit()
dim(Tes_Expr)
### 1889 genes(consistent with old pipeline)

### PCA
pca.all <- prcomp(t(Tes_Expr),retx=T)

### Create horizontal and vertical labels
# PC3comp<-summary(pca.all)$importance[2,3]*100
lev <- c("BLIS","IM","MES","LAR")
subtype_pal <- c("#ef4922","#7bc242","#3baade","#8f80ba")
Tes_Clusters$subtype <- factor(Tes_Clusters$subtype,levels = lev)
pca_cluster <- fviz_pca_ind(pca.all, 
                            axes = c(1, 2),
                            label = "none",
                            addEllipses = T,
                            geom.ind = c("point", "text"),
                            repel = T,
                            habillage = Tes_Clusters$subtype,
                            ellipse.level=0.95,
                            palette = subtype_pal)+
  theme(axis.text = element_text(size=16,face='bold',color='gray40'),
        axis.title = element_text(size=16,face='bold'),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16,color='gray40'),
        plot.subtitle = element_text(hjust=0.5,size=16),
        plot.title = element_text(hjust=0.5,size=16))+
  scale_shape_manual(values = c(15:19))+
  labs(x = sprintf("PC1 (%.2f%%)", summary(pca.all)$importance[2,1]*100),
       y = sprintf("PC2 (%.2f%%)", summary(pca.all)$importance[2,2]*100),
       title = "")+
  theme(aspect.ratio = 1/1);pca_cluster
#,title=paste('PCA-All ',ncol(log_filtered.expr.mat.ft),'samples',' N=',nrow(log_filtered.expr.mat.ft))
ggsave("./results/charts/PCA_all_subtype_new_448samples_1889genes.pdf",pca_cluster,width=6,height=4,dpi=300)

### total samples
patients_unmatched <- readRDS("./Rdata/patients_unmatched.rds")

pcs <- pca.all$x %>% as.data.frame()
pcs <- pcs[,1:3]
pcs$patient <- sapply(strsplit(as.character(rownames(pcs)),"_"),function(x){x[1]})
pcs$type <- sapply(strsplit(as.character(rownames(pcs)),"_"),function(x){x[2]})
pcs$type[pcs$type =="rep"] <- "TT"; pcs$type[is.na(pcs$type)] <- "TT"
pcs$type <- gsub("TT","tumor",pcs$type); pcs$type <- gsub("PT","normal",pcs$type)
pcs$compare <- ifelse(pcs$patient %in% patients_unmatched, "inconsistent", "consistent")
pcs$class <- paste0(pcs$type,"-subtype_",pcs$compare)

pcs_an <- merge(pcs,Tes_Clusters,by="patient")
pcs_an$type <- factor(pcs_an$type,levels = c("TT","PT"))
pcs_an$subtype <- factor(pcs_an$subtype,levels = lev)

# MES_unmatched <- readRDS("./Rdata/MES_unmatched.rds")
# for_label <- pcs_an[pcs_an$patient %in% MES_unmatched,]

pca_cluster_total_new <- ggplot(pcs_an,aes(x=PC1,y=-PC2,color=subtype))+
  geom_point(aes(shape=class),size=2.5)+
  # geom_encircle(aes(color=type,fill=type), expand=0.02, alpha=.2)+
  theme_few()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.position = "none",
        plot.subtitle = element_text(hjust=0.5,size=14),
        plot.title = element_text(hjust=0.5,size=16))+
  scale_color_manual(values = subtype_pal)+
  scale_fill_manual(values = subtype_pal)+
  scale_shape_manual(values = c(0,15,2,17))+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  guides(shape=guide_legend(override.aes = list(size=3)))+
  labs(x = sprintf("PC1 (%.2f%%)", summary(pca.all)$importance[2,1]*100),
       y = sprintf("PC2 (%.2f%%)", summary(pca.all)$importance[2,2]*100),
       title = ""); pca_cluster_total_new
  # geom_text_repel(data = for_label, aes(PC1,PC2,
  #                                       # color=factor(subtype),
  #                                       label=patient),
  #                 arrow = arrow(length=unit(0.01, "npc")),
  #                 max.overlaps = 30);pca_cluster_total
ggsave("./results/charts/PCA_all_subtype_new_448samples(all)_1889genes_labeled.pdf",pca_cluster_total,width=7,height=4,dpi=300)
ggsave("./results/charts/PCA_all_subtype_new_448samples(all)_1889genes_labeled_nolegend.pdf",pca_cluster_total,width=4,height=4,dpi=300)


###################----------111 genes not in new pipeline----------##################
genes_111 <- setdiff(Clustering_Genes,rownames(Tes_Expr))
# test <- as.data.frame(Clustering_Genes)
# test1 <- rownames(Tes_Expr) %>% as.data.frame()
genes_111 <- genes_111[-which(genes_111 %in% c("IGJ","KIAA175"))] 

### Expression genes for clustering
# Tes_Expr <- logexpr_old_TT[ProteinCoding, ]
Tes_Expr <- logexpr_old[ProteinCoding, ]

dim(Tes_Expr)
Tes_Expr <- Tes_Expr[genes_111, ] %>% na.omit()

dim(Tes_Expr)
### 111 genes(old pipeline)

### PCA
pca.all <- prcomp(t(Tes_Expr),retx=T)

### Create horizontal and vertical labels
# PC3comp<-summary(pca.all)$importance[2,3]*100
lev <- c("BLIS","IM","MES","LAR")
subtype_pal <- c("#ef4922","#7bc242","#3baade","#8f80ba")
Tes_Clusters$subtype <- factor(Tes_Clusters$subtype,levels = lev)

# pca.result<-as.data.frame(pca.all$x)
# pca.result$patient<-rownames(pca.result)
# pca.result.total <- merge(pca.result,Tes_Clusters,by="patient")
# pca.result.total = pca.result.total %>% mutate(colour = case_when(
#   pca.result.total$subtype == "BLIS" ~ "#ef4922",
#   pca.result.total$subtype == "IM" ~ "#7bc242",
#   pca.result.total$subtype == "MES" ~ "#3baade",
#   pca.result.total$subtype == "LAR" ~ "#8f80ba"
# ))
# 
# scatterplot3d(pca.result.total[,2:4],
#               color = pca.result.total$colour,
#               pch = 16,angle=30,
#               cex.symbols = 2,
#               box=T,type="p",
#               main = "3D PCA Plot",
#               lty.hide=2,lty.grid = 2)

pca_cluster <- fviz_pca_ind(pca.all,
                            axes = c(2, 3),
                            label = "none",
                            addEllipses = T,
                            geom.ind = c("point", "text"),
                            repel = T,
                            habillage = Tes_Clusters$subtype,
                            ellipse.level=0.95,
                            palette = subtype_pal)+
  theme(axis.text = element_text(size=16,face='bold',color='gray40'),
        axis.title = element_text(size=16,face='bold'),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16,color='gray40'),
        plot.subtitle = element_text(hjust=0.5,size=16),
        plot.title = element_text(hjust=0.5,size=16))+
  scale_shape_manual(values = c(15:19))+
  labs(x = sprintf("PC2 (%.2f%%)", summary(pca.all)$importance[2,2]*100),
       y = sprintf("PC3 (%.2f%%)", summary(pca.all)$importance[2,3]*100),
       title = "")+
  theme(aspect.ratio = 1/1);pca_cluster
#,title=paste('PCA-All ',ncol(log_filtered.expr.mat.ft),'samples',' N=',nrow(log_filtered.expr.mat.ft))
ggsave("./results/charts/PCA_all_subtype_old_448samples_111genes.pdf",pca_cluster,width=6,height=4,dpi=300)

### total samples
Tes_Clusters <- read.csv("./results/tables/FUSCCTNBC_old_subtypes.csv")
colnames(Tes_Clusters) <- c("patient","subtype")

patients_unmatched <- readRDS("./Rdata/patients_unmatched.rds")

pcs <- pca.all$x %>% as.data.frame()
pcs <- pcs[,1:3]
pcs$patient <- sapply(strsplit(as.character(rownames(pcs)),"\\."),function(x){x[1]})
pcs$type <- sapply(strsplit(as.character(rownames(pcs)),"\\."),function(x){x[2]})
pcs$type[is.na(pcs$type)] <- "TT"
pcs$type <- gsub("TT","tumor",pcs$type); pcs$type <- gsub("PT","normal",pcs$type)
pcs$compare <- ifelse(pcs$patient %in% patients_unmatched, "inconsistent", "consistent")
pcs$class <- paste0(pcs$type,"-subtype_",pcs$compare)

pcs_an <- merge(pcs,Tes_Clusters,by="patient")
pcs_an$type <- factor(pcs_an$type,levels = c("tumor","normal"))
pcs_an$subtype <- factor(pcs_an$subtype,levels = lev)

pca_cluster_total_111 <- ggplot(pcs_an,aes(x=PC1,y=-PC2,color=subtype))+
  geom_point(aes(shape=class),size=2.5)+
  # geom_encircle(aes(color=type,fill=type), expand=0.02, alpha=.2)+
  theme_few()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        # legend.position = "none",
        plot.subtitle = element_text(hjust=0.5,size=14),
        plot.title = element_text(hjust=0.5,size=16))+
  scale_color_manual(values = subtype_pal)+
  scale_fill_manual(values = subtype_pal)+
  scale_shape_manual(values = c(0,15,2,17))+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  guides(shape=guide_legend(override.aes = list(size=3)))+
  labs(x = sprintf("PC1 (%.2f%%)", summary(pca.all)$importance[2,1]*100),
       y = sprintf("PC2 (%.2f%%)", summary(pca.all)$importance[2,2]*100),
       title = "");pca_cluster_total_111
ggsave("./results/charts/PCA_all_subtype_old_448samples(all)_111genes.pdf",pca_cluster_total,width=7,height=4,dpi=300)
ggsave("./results/charts/PCA_all_subtype_old_448samples(all)_111genes_nolegend.pdf",pca_cluster_total,width=4,height=4,dpi=300)

panel.all <- ggarrange(pca_cluster_total_old, pca_cluster_total_new,
                      pca_cluster_total_111,nrow = 1,labels = c("a","b","c"),
                      widths = c(1,1,1.7)) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.all
ggsave("./results/charts/PCA_all_subtype_old_new_111genes_448samples(all).pdf",panel.all,width=15,height=4,dpi=300)


## remove all
rm(list=ls())
