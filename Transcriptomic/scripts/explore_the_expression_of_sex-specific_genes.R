## TNBC2019 RNAseq - explore the expression of sex-specific genes
## Author：Qingwang Chen
## Date：2022-5-18
## Version：1.0

# Set work directory and import libraries
rm(list = ls (all = TRUE))

source("./scripts/libraries.R")
source("./scripts/parameter.R")

######## Batch Information Process
######## batch info import
RNAseq_batch_info <- read.csv("./data/RNAseq_batch_info.csv",header=T)
RNAseq_batch_info$sample<- sapply(strsplit(as.character(RNAseq_batch_info$sample_id),"\\."),function(x){x[1]})
RNAseq_batch_info$type<- sapply(strsplit(as.character(RNAseq_batch_info$sample_id),"\\."),function(x){x[2]})
RNAseq_batch_info$type <- gsub("rep","TT",RNAseq_batch_info$type)
RNAseq_batch_info$type[which(is.na(RNAseq_batch_info$type))] <- "TT"
RNAseq_batch_info$sample_id <- paste0("TPM.",RNAseq_batch_info$sample,"_",RNAseq_batch_info$type)
RNAseq_batch_info$sample_name <- paste0(RNAseq_batch_info$sample,"_",RNAseq_batch_info$type)
# batch2 here is actually batch1 in the origin
RNAseq_batch_info$batch <- gsub("batch1","x",RNAseq_batch_info$batch )
RNAseq_batch_info$batch <- gsub("batch2","batch1",RNAseq_batch_info$batch )
RNAseq_batch_info$batch <- gsub("x","batch2",RNAseq_batch_info$batch )


######## Exp Data Process
## 1. Gender-checked sex cluster
### Import of gender-checked related genes and exp files
sexgene <- read.table("./data/sexgenelist.txt",header=T,sep="\t",stringsAsFactors = FALSE)
gene_tpm <- readRDS("./Rdata/gene_tpm_r56858c448_genesymbol_TNBC2019_20220518.rds")
dim(gene_tpm)
### Calculate log(expr+0.01)
log_gene_expr <- apply(gene_tpm,2,function(x){log2(x+0.01)})

### Screening sample for lines associated with sex-checked genes
logexpr_sex <- log_gene_expr[rownames(log_gene_expr) %in% sexgene$GeneSymbol,]

annotation_col <- RNAseq_batch_info[,c("sample_id","batch")]
row.names(annotation_col) <- annotation_col$sample_id
annotation_col$sample_id <- NULL
colnames(annotation_col) <- "Batch"

# plot
anocolor <- list(Batch = c(batch1="#2ea15c",batch2="#6c90c7",batch3="#7e25cf"))
sex_heatmap_withbatch <- pheatmap(logexpr_sex,scale = "column",distance_rows = "euclidean",
                        clustering_distance_cols = "euclidean",
                        annotation_col=annotation_col, annotation_colors = anocolor,
                        colorRampPalette(c('#436eee','white','#EE0000'))(100),
                        clustering_method = "ward.D",
                        show_rownames = TRUE,show_colnames=FALSE);sex_heatmap_withbatch

# pheatmap to ggplot
sex_heatmap_withbatch_g = as.ggplot(sex_heatmap_withbatch)

ggsave("./results/charts/sex check total with batch information.pdf",sex_heatmap_withbatch_g,width=6,height=4,dpi=300)

### Filter out samples with abnormal expression of RPS4Y1 gene
sample_gene_sex <- t(logexpr_sex) %>% as.data.frame()
max(sample_gene_sex$RPS4Y1)
# 1.064635 TPM.FUSCCTNBC475 (TT)
# rm(list=ls())

## Remove batch2
logexpr_sex_batch13 <- logexpr_sex[,-which(colnames(logexpr_sex) %in% RNAseq_batch_info$sample_id[RNAseq_batch_info$batch == "batch2"])]
sex_heatmap_withbatch <- pheatmap(logexpr_sex_batch13,scale = "column",distance_rows = "euclidean",
                                  clustering_distance_cols = "euclidean",
                                  annotation_col=annotation_col, annotation_colors = anocolor,
                                  colorRampPalette(c('#436eee','white','#EE0000'))(100),
                                  clustering_method = "ward.D",
                                  show_rownames = TRUE,show_colnames=FALSE);sex_heatmap_withbatch


######## Count
count <- readRDS("./Rdata/FUSCCTNBC_RNAseqShi_Counts_Unfiltered_448x56858_200521.rds")
count_col <- colnames(count) %>% as.data.frame(); colnames(count_col) <- "name"
count_col$sample <- sapply(strsplit(as.character(count_col$name),"_"),function(x){x[1]})
count_col$type <- sapply(strsplit(as.character(count_col$name),"_"),function(x){x[2]})
count_col$type <- gsub("rep","TT",count_col$type)
count_col$type[which(is.na(count_col$type))] <- "TT"
count_col$sample_id <- paste0(count_col$sample,"_",count_col$type)
colnames(count) <- count_col$sample_id

### Calculate log(count+1)
log_count <- apply(count,2,function(x){log2(x+1)})

### Screening sample for lines associated with sex-checked genes
logcount_sex <- log_count[rownames(log_count) %in% sexgene$GeneSymbol,]

annotation_col <- RNAseq_batch_info[,c("sample_name","batch")]
row.names(annotation_col) <- annotation_col$sample_name
annotation_col$sample_name <- NULL
colnames(annotation_col) <- "Batch"

# plot
anocolor <- list(Batch = c(batch1="#2ea15c",batch2="#6c90c7",batch3="#7e25cf"))
sex_heatmap_withbatch <- pheatmap(logcount_sex,scale = "column",distance_rows = "euclidean",
                                  clustering_distance_cols = "euclidean",
                                  annotation_col=annotation_col, annotation_colors = anocolor,
                                  colorRampPalette(c('#436eee','white','#EE0000'))(100),
                                  clustering_method = "ward.D",
                                  show_rownames = TRUE,show_colnames=FALSE);sex_heatmap_withbatch

# pheatmap to ggplot
sex_heatmap_withbatch_g = as.ggplot(sex_heatmap_withbatch)

ggsave("./results/charts/sex check total with batch information (Count).pdf",sex_heatmap_withbatch_g,width=6,height=4,dpi=300)

######## RUVseq
sample <- as.character(count_col$type)
x <- as.factor(sample)
set <- newSeqExpressionSet(as.matrix(count),
                           phenoData = data.frame(x, row.names=colnames(count)))

## ----empirical--------------------------------------------
design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design) # consume long time

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)

top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:round(nrow(set)*1/2)]))]

## ----emp_ruvg, fig.cap="RUVg normalization based on empirical controls.", fig.subcap=c("RLE plot","PCA plot")----
ubatch <- unique(RNAseq_batch_info$batch)
set2 <- RUVg(set, empirical, k=length(ubatch))

#plotPCA(set2, col=colors.sample.Quartet.fill[x], cex=0.7)

RUVg_normed_Counts<-normCounts(set2)
logRUVg_counts <-apply(RUVg_normed_Counts,2,function(x){log2(x+1)})
dim(logRUVg_counts)
saveRDS(logRUVg_counts,"./Rdata/log_counts_RUVg_r56858c448_20220604.rds")

### Screening sample for lines associated with sex-checked genes
logRUVgcount_sex <- logRUVg_counts[rownames(logRUVg_counts) %in% sexgene$GeneSymbol,]

# annotation_col <- RNAseq_batch_info[,c("sample_name","batch","type")]
# row.names(annotation_col) <- annotation_col$sample_name
# annotation_col$sample_name <- NULL
# colnames(annotation_col) <- c("Batch","Type")

# plot
anocolor <- list(Batch = c(batch1="#2ea15c",batch2="#6c90c7",batch3="#7e25cf"))
sex_heatmap_withbatch <- pheatmap(logRUVgcount_sex,scale = "column",distance_rows = "euclidean",
                                  clustering_distance_cols = "euclidean",
                                  annotation_col=annotation_col, annotation_colors = anocolor,
                                  colorRampPalette(c('#436eee','white','#EE0000'))(100),
                                  clustering_method = "ward.D",
                                  show_rownames = TRUE,show_colnames=FALSE);sex_heatmap_withbatch

# pheatmap to ggplot
sex_heatmap_withbatch_g = as.ggplot(sex_heatmap_withbatch)

ggsave("./results/charts/sex check total with batch information (RUVg normalized Count).pdf",sex_heatmap_withbatch_g,width=6,height=4,dpi=300)

rm(list = ls())


