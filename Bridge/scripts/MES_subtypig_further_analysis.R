## Title: TNBC Bridge RNAseq Analysis-- MES subtypig further analysis##
## Author: Qingwang Chen             ##
## Date: 2022-06-07                  ##  
## Version: V3.0                     ##
#######################################

# Set work directory and import libraries
rm(list = ls (all = TRUE))
source("./scripts/libraries.R")
source("./scripts/parameter.R")

# ----------------------Subtyping---------------------- #
Tes_Clusters <- read.csv("./results/tables/FUSCCTNBC_new_subtypes.csv")
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

###################--------------new pipeline-----------------##################
#### new expression dataset
fpkm_new_pipe <- readRDS("./Rdata/FUSCCTNBC_new_fpkm_448x56858_220302.rds")
#### unfiltered
fpkm_new_pipe_PT <- fpkm_new_pipe[,grep("PT",colnames(fpkm_new_pipe))]
fpkm_new_pipe_TT <- fpkm_new_pipe[,-grep("PT",colnames(fpkm_new_pipe))]

# filtered_fpkm_new_pipe_TT
logexpr_new_TT <- log2(fpkm_new_pipe_TT+1)
Tes_Expr <- logexpr_new_TT[ProteinCoding, ] %>% na.omit()
dim(Tes_Expr)

### Expression genes for clustering
Tes_Expr <- logexpr_new_TT[Clustering_Genes, ] %>% na.omit()
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
ggsave("./results/charts/PCA_all_subtype_448samples_1889genes.pdf",pca_cluster,width=6,height=4,dpi=300)

## remove all
rm(list=ls())
