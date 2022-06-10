## TNBC2019 RNAseq - Consistency of gene mutation frequency in different analyses for WES 
## Author：Qingwang Chen
## Date：2022-5-21
## Version：1.0


################-------------------------Preperation---------------------################
# Set work directory and import libraries
rm(list = ls (all = TRUE))
source("./scripts/libraries.R")
source("./scripts/parameter.R")
library(maftools)


# 1. Data process
# maf from old analysis
maf_old_txt <- read.table("./data/FUSCCTNBC_Mutations_incIntronIGR_V15.1.txt",header = T)
maf_old <- read.maf(maf_old_txt)
write.mafSummary(maf = maf_old, basename = './results/tables/maf_old')

# maf from new analysis
maf_new <- read.maf("./data/data_mutations_extended_somatic_hg38.maf")
write.mafSummary(maf = maf_new, basename = './results/tables/maf_new')

# oncoplot(maf = test1, top = 20,writeMatrix = T)

# gene mutation info
gene_mutation_old <- read.table("./results/tables/maf_old_geneSummary.txt",header = T)
gene_mutation_old$GeneMF <- round(gene_mutation_old$MutatedSamples/279,2)

gene_mutation_new <- read.table("./results/tables/maf_new_geneSummary.txt",header = T)
gene_mutation_new$GeneMF <- round(gene_mutation_new$MutatedSamples/279,2)

compare_gene_mutation <- merge(gene_mutation_new,gene_mutation_old,by="Hugo_Symbol")

# simple test
plot(compare_gene_mutation$GeneMF.x,compare_gene_mutation$GeneMF.y)
cor.test(compare_gene_mutation$GeneMF.x,compare_gene_mutation$GeneMF.y)
write.csv(compare_gene_mutation,"./results/tables/compare_genes_info.csv")

p <- ggscatter(compare_gene_mutation, x = "GeneMF.x", y = "GeneMF.y",
                add = "reg.line", conf.int = TRUE, 
                add.params = list(linetype=2, color = "lightblue", fill = "lightgray"))+
  stat_cor(method = "pearson",alternative = "two.sided",r.digits = 4,
           p.digits = 4,label.x = 0.4, label.y = 0.2)+
  labs(title="", x="New pipeline", y="Old pipeline")+
  theme(plot.title = element_text(hjust = 0.5))+mytheme;p
# geom_label_repel(aes(percentile.x,percentile.y,label=Var1));p

ggsave('./results/charts/all_genes_mutation_frequency_scatter.png',p,width=4,height=4)

