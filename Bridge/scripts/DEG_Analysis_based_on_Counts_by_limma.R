## TNBC2019 RNAseq - DEG Analysis based on Counts by limma
## Author：Qingwang Chen
## Date：2022-5-21
## Version：1.0

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
counts_new_pipe <- fread("./data/FUSCCTNBC_Gene_Counts_new_20220521.csv") %>% as.data.frame()
counts_new_pipe$gene_symbol <- sapply(strsplit(as.character(counts_new_pipe$gene_id),"\\|"),function(x){x[2]})
counts_new_pipe$gene_id <- NULL
## keep the max counts if the same gene symbol
counts_new_pipe_max <- aggregate(counts_new_pipe[,-449], list(counts_new_pipe$gene_symbol), FUN=max)
rownames(counts_new_pipe_max) <- counts_new_pipe_max$Group.1
counts_new_pipe_max$Group.1 <- NULL
dim(counts_new_pipe_max)
colnames(counts_new_pipe_max)
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
#### cor=0.974329 
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
