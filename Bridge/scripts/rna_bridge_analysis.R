## Title: TNBC Bridge RNAseq Analysis##
## Author: Qingwang Chen             ##
## Date: 2022-06-07                  ##  
## Version: V3.0                     ##
#######################################

# Set work directory and import libraries
rm(list = ls (all = TRUE))
source("./scripts/libraries.R")
source("./scripts/parameter.R")

# 1. Data process
# counts from old analysis
load("./data/FUSCCTNBC_RNAseqShi_Counts_Unfiltered_448x45308_V15_190223.Rdata")
counts_old_pipe <- FUSCCTNBC_RNAseqShi_Counts
sample_list_old <- colnames(counts_old_pipe) %>% as.data.frame()
colnames(sample_list_old) <- "sample"
sample_list_old$type <- sapply(strsplit(as.character(sample_list_old$sample),"\\."),function(x){x[2]})
sample_list_old$type[which(is.na(sample_list_old$type))] <- "TT"

# counts from new analysis
counts_new_pipe_max <- readRDS("./Rdata/counts_RUVg_r56858c448_20220604.rds") %>% as.data.frame()
# fwrite(counts_new_pipe_max,"./Rdata/FUSCCTNBC_Expression_RNAseqCount_RUVg_hg38.tsv",sep = "\t",row.names = T)
# counts_new_pipe <- fread("./data/FUSCCTNBC_Gene_Counts_new_20220521.csv") %>% as.data.frame()
# counts_new_pipe$gene_symbol <- sapply(strsplit(as.character(counts_new_pipe$gene_id),"\\|"),function(x){x[2]})
# counts_new_pipe$gene_id <- NULL
# ## keep the max counts if the same gene symbol
# counts_new_pipe_max <- aggregate(counts_new_pipe[,-449], list(counts_new_pipe$gene_symbol), FUN=max)
# rownames(counts_new_pipe_max) <- counts_new_pipe_max$Group.1
# counts_new_pipe_max$Group.1 <- NULL
# dim(counts_new_pipe_max)
# colnames(counts_new_pipe_max)
# saveRDS(counts_new_pipe_max,"./Rdata/FUSCCTNBC_RNAseqShi_Counts_Unfiltered_448x56858_200521.rds")
sample_list_new <- colnames(counts_new_pipe_max) %>% as.data.frame()
colnames(sample_list_new) <- "sample"
sample_list_new$type <- sapply(strsplit(as.character(sample_list_new$sample),"_"),function(x){x[2]})
sample_list_new$type[which(is.na(sample_list_new$type))] <- "TT"
sample_list_new$type <- gsub("rep","TT",sample_list_new$type)

## limma-voom
p_load(edgeR) 
DEG_cal <- function(exprMat,group,compare,p.thr=0.05,FC.thr=2,n.label=10,use.p.adjust=F){
  ## Parameter explanation.
  ## exprMat: the counts expression matrix as input.
  ## group: a vector that represents the grouping information of the samples, which needs to correspond to the column names of counts expression matrix. 
  ## compare: a vector for comparing differences, e.g., "groupTumor-groupNormal" (must be enclosed in quotes and preceded by group), which is needed to perform contrast.fit in limma.
  ## p.thr:cut-off of p-value, default 0.05.
  ## FC.thr:cut-off of fold change, default 2.
  ## n.label: the number of genes you want to label top in the volcano map, default is top10 for each up and down-regulated.
  ## use.p.adjust: if or not to use corrected p-value for up and down-regulated, default FALSE.
  group <- factor(group,ordered = F)
  design <- model.matrix(~0+group)
  ## Filter the counts expression matrix
  layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))
  dge <- DGEList(counts = exprMat)
  keep <- filterByExpr(dge, design)
  dge <- dge[keep,keep.lib.sizes=F]
  ## Calculate the normalization factor
  dge <- calcNormFactors(dge)
  ## View the grouping factor variance, dim1 is the farthest away with the largest potential variance.
  ## mds analysis can be alternative.
  col.group <- group
  levels(col.group) <- brewer.pal(nlevels(col.group), "Paired") 
  col.group <- as.character(col.group)
  plotMDS(dge, labels=group, col=col.group)
  title("sample group")
  ## voom conversion and variance weight calculation
  ## First the original counts are converted to log2 CPM (counts per million reads)
  ## where per million reads was specified based on the norm.factors previously calculated by calcNormFactors.
  ## A linear model was then produced based on the log2CPM for each gene, and residuals were calculated; ##
  ## The sqrt(residual standard deviation) was fitted using the mean expression (red line).
  ## The final smoothing curve obtained can be used to obtain the weights for each gene and sample
  ## The data have been filtered out for low expression genes, the smooth curve indicates good quality, 
  ## and the curve has an upward bump indicating that further filtering is needed.
  v <- voom(dge, design, plot=T) #plot可以等于F
  ## Generalized Linear Fitting
  fit1 <- lmFit(v, design)
  ## Difference Comparison Matrix
  contr <- makeContrasts(contrasts=compare,levels=design)
  ## Analysis of Variance Fitting
  fit2 <- contrasts.fit(fit1,contr)
  ## Empirical Bayesian fit
  fit3 <- eBayes(fit2)
  ## Plot log2 residual standard deviation versus log-CPM mean
  plotSA(fit3, main="Final model: Mean-variance trend")
  ## Results
  deg.data <- topTable(fit3,sort.by = "P",number = Inf)
  ## Annotate the expression up-down and top.gene, the following parts do not affect the calculation results
  logFC.thr <- log2(FC.thr)
  deg.data$Symbol <- rownames(deg.data)
  if(use.p.adjust){
    deg.data$logP.adj <- -log10(deg.data$adj.P.Val)
    deg.data$group <-rep("not-significant",nrow(deg.data));
    deg.data$group[which((deg.data$adj.P.Val<p.thr)&(deg.data$logFC> logFC.thr))] <- "up-regulated" ;
    deg.data$group[which((deg.data$adj.P.Val<p.thr)&(deg.data$logFC< -logFC.thr))] <- "down-regulated" ;
    deg.data$label <- rep("",nrow(deg.data))
    deg.data <- deg.data[order(deg.data$adj.P.Val),]
    up.genes <- head(deg.data$Symbol[which(deg.data$group=="up-regulated")],n.label)
    down.genes <- head(deg.data$Symbol[which(deg.data$group=="down-regulated")],n.label)
    deg.top.genes <- c(as.character(up.genes),
                       as.character(down.genes))
    deg.data$label[match(deg.top.genes,deg.data$Symbol)] <- deg.top.genes
  }else{
    deg.data$logP <- -log10(deg.data$P.Value)
    deg.data$group <- rep("not-significant",nrow(deg.data));
    deg.data$group[which((deg.data$P.Value<p.thr)&(deg.data$logFC> logFC.thr))] <- "up-regulated" ;
    deg.data$group[which((deg.data$P.Value<p.thr)&(deg.data$logFC< -logFC.thr))] <- "down-regulated" ;
    deg.data$label <- rep("",nrow(deg.data))
    deg.data <- deg.data[order(deg.data$logFC,decreasing = T),]
    up.genes <- head(deg.data$Symbol[which(deg.data$group=="up-regulated")],n.label)
    deg.data <- deg.data[order(deg.data$logFC),]
    down.genes <- head(deg.data$Symbol[which(deg.data$group=="down-regulated")],n.label)
    deg.top.genes <- c(as.character(up.genes),
                       as.character(down.genes))
    deg.data$label[match(deg.top.genes,deg.data$Symbol)] <- deg.top.genes
  }
  return(deg.data)
}

deg_old <- DEG_cal(exprMat=counts_old_pipe,sample_list_old$type,('groupTT-groupPT'))
deg_new <-  DEG_cal(exprMat=counts_new_pipe_max,sample_list_new$type,('groupTT-groupPT'))

## Filter DEG and Log2FC Comparison
deg_old_pipe <- deg_old[,c("logFC","P.Value","Symbol")]
filter_deg_old_pipe <- deg_old_pipe[abs(deg_old_pipe$logFC)>1 & deg_old_pipe$P.Value <0.05,]

deg_new_pipe <- deg_new[,c("logFC","P.Value","Symbol")]
filter_deg_new_pipe <- deg_new_pipe[abs(deg_new_pipe$logFC)>1 & deg_new_pipe$P.Value <0.05,]

compare_FC <- merge(filter_deg_new_pipe,filter_deg_old_pipe,by="Symbol")
plot(compare_FC$logFC.x,as.numeric(compare_FC$logFC.y))  
cor.test(compare_FC$logFC.x,as.numeric(compare_FC$logFC.y))
#### cor=0.9777226 
write.csv(compare_FC,"./results/tables/compare_FC.csv")

## plot
compare_FC_plot <- ggscatter(compare_FC, x = "logFC.x", y = "logFC.y",
                             add = "reg.line", conf.int = TRUE, alpha = 0.1,
                             add.params = list(linetype=2, color = "lightblue", fill = "lightgray"))+
  stat_cor(method = "pearson",alternative = "two.sided",r.digits = 4,
           p.digits = 4,label.x = 0, label.y = -2.5)+
  labs(title="", x="New pipeline", y="Old pipeline")+
  theme(plot.title = element_text(hjust = 0.5));compare_FC_plot

# for_label <- compare_FC[(-2 < compare_FC$logFC.x & compare_FC$logFC.x < 0) & 
#                           (0 < compare_FC$logFC.y & compare_FC$logFC.y < 4),]
# 
# compare_FC_plot + geom_label_repel(data = for_label, aes(label = Symbol), 
#                                    alpha = 0.7, max.overlaps = 15) 

ggsave('./results/charts/compare_FC_plot_scatter.png',compare_FC_plot,width=5,height=5)
ggsave('./results/charts/compare_FC_plot_scatter.pdf',compare_FC_plot,width=5,height=5)


# ----------------------Subtyping---------------------- #
## method from Dr. Ma (method before)
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
# saveRDS(Clustering_Genes,"./Rdata/Clustering_Genes_old.rds")
Clustering_Genes <- readRDS("./Rdata/Clustering_Genes_old.rds")

### Expression profile for clustering
Tes_Expr <- Tes_Expr[which(rownames(Tes_Expr) %in% Clustering_Genes), ]
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
### filter the protein coding genes
# Gene_annotation_TNBC <- readRDS("./Rdata/gene_anotation_ensembl_gene_id.rds") %>% distinct(external_gene_name, .keep_all = T)
# row.names(Gene_annotation_TNBC) <- Gene_annotation_TNBC$external_gene_name
# Gene_annotation_TNBC <- Gene_annotation_TNBC[which(!is.na(Gene_annotation_TNBC$gene_biotype)), ]
# ProteinCoding <- Gene_annotation_TNBC$external_gene_name[Gene_annotation_TNBC$gene_biotype == "protein_coding"]

#### new expression dataset
fpkm_new_pipe <- readRDS("./Rdata/FUSCCTNBC_new_fpkm_448x56858_220302.rds")
#### unfiltered
fpkm_new_pipe_PT <- fpkm_new_pipe[,grep("PT",colnames(fpkm_new_pipe))]
fpkm_new_pipe_TT <- fpkm_new_pipe[,-grep("PT",colnames(fpkm_new_pipe))]

# filtered_fpkm_new_pipe_TT
logexpr_new_TT <- log2(fpkm_new_pipe_TT+1)
Tes_Expr <- logexpr_new_TT[ProteinCoding, ] %>% na.omit()
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

### Expression genes for clustering
Tes_Expr <- logexpr_new_TT[Clustering_Genes, ] %>% na.omit()
# Tes_Expr <- logexpr_new_TT[which(rownames(Tes_Expr) %in% Clustering_Genes), ] %>% na.omit()
dim(Tes_Expr)
### 1889 genes(consistent with old pipeline)

### clustering by k-means
set.seed(1234)
km1 <- kmeans(t(Tes_Expr),4,iter.max = 1000,nstart = 100)
table(km1$cluster)
Tes_Clusters <- km1$cluster

#### Rename clusters
Tes_Clusters[Tes_Clusters == "1"] <- "IM"  ; Tes_Clusters[Tes_Clusters == "2"] <- "LAR"
Tes_Clusters[Tes_Clusters == "3"] <- "BLIS" ; Tes_Clusters[Tes_Clusters == "4"] <- "MES"

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
FUSCCTNBC_new_subtypes$patient <- gsub("_rep","",FUSCCTNBC_new_subtypes$patient)
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

