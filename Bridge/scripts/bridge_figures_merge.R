## Title: TNBC Bridge Figures Merge  ##
## Author: Qingwang Chen             ##
## Date: 2022-03-02                  ##  
## Version: V1.0                     ##
#######################################

# Set work directory and import libraries
rm(list = ls (all = TRUE))
source("./scripts/libraries.R")
source("./scripts/parameter.R")

#### --------------------------------RNAseq--------------------------------#### 
## a. FC Correlation & b. Subtype consistency
compare_FC <- read.csv("./results/tables/compare_FC.csv")

compare_FC_plot <- ggscatter(compare_FC, x = "logFC.x", y = "logFC.y",
                             add = "reg.line", conf.int = TRUE, alpha = 0.1,
                             add.params = list(linetype=2, color = "lightblue", fill = "lightgray"))+
  stat_cor(method = "pearson",alternative = "two.sided",r.digits = 4,
           p.digits = 4,label.x = 0, label.y = -2.5)+
  labs(title="", x="Log2FC (New)", y="Log2FC (Old)")+
  theme(plot.title = element_text(hjust = 0.5));compare_FC_plot

panel.a <- ggarrange(compare_FC_plot, nrow = 1, labels = c("a")) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.a
ggsave("./results/charts/panela.pdf",panel.a,width=6,height=6,dpi=300)

panel.ab <- ggarrange(panel.a,NULL,widths = c(1,2)) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.ab
ggsave("./results/charts/panel ab.pdf",panel.ab,width=12,height=4,dpi=300)

## c. tumor purity statistics by four subtypes
tnbc_tumor_ratio <- fread("./data/Supplementary Table 5.csv") %>% as.data.frame()

tnbc_tumor_ratio_sub <- drop_na(tnbc_tumor_ratio)

lev <- c("BLIS","IM","MES","LAR")
subtype_pal <- c("#ef4922","#7bc242","#3baade","#8f80ba")
tnbc_tumor_ratio_sub$subtype_old <- factor(tnbc_tumor_ratio_sub$subtype_old,levels = lev)
# 315 samples
aggregate(tnbc_tumor_ratio_sub$tumor_purity, by = list(tnbc_tumor_ratio_sub$subtype_old), FUN = sd)

#------------------------plot---------------------------------#
my_comparisons <- list(c("BLIS", "IM"),
                       c("BLIS", "MES"),
                       c("BLIS", "LAR"),
                       c("IM", "MES"),
                       c("IM", "LAR"),
                       c("MES", "LAR"))

#  tumor purity by four subtypes(Boxplot)
tp4sub_plot <- ggboxplot(tnbc_tumor_ratio_sub,x = "subtype_old",y = "tumor_purity",
          color = "subtype_old", palette = subtype_pal,add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,label="p.signif",
                     label.y=100)+   # Pairwise comparison against all
  stat_compare_means(method = "anova", label.y=110)+ # Add global p-value
  xlab("")+ylab("Tumor Purity (%)")+
  theme(legend.position = "none"); tp4sub_plot

panel.c <- ggarrange(tp4sub_plot, nrow = 1, labels = c("c")) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.c
ggsave("./results/charts/panelc.pdf",panel.c,width=6,height=6,dpi=300)

## d. PCA for the all samples annotated with the tumor purity of MES
pcs_an <- read.csv("./results/tables/pca_all_samples_360patients_MES_annotated.csv")
pca.all <- readRDS("./Rdata/pca_all_samples_360patients_MES_annotated.rds")
for_label <- readRDS("./Rdata/pca_all_samples_360patients_MES_annotated_for_label.rds")

lev <- c("BLIS","IM","MES","LAR")
subtype_pal <- c("#ef4922","#7bc242","#3baade","#8f80ba")
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

pca_cluster_total_old_anotated <- pca_cluster_total_old+
                        geom_text_repel(data = for_label, aes(PC1,PC2,
                                        color=factor(subtype),
                                        label=tumor_purity),
                                      arrow = arrow(length=unit(0.01, "npc")),
                                      max.overlaps = 30);pca_cluster_total_old_anotated
ggsave("./results/charts/PCA_all_subtype_old_448samples(all)_2000genes_TP.pdf",pca_cluster_total_old_anotated,width=10,height=7,dpi=300)

panel.d <- ggarrange(pca_cluster_total_old_anotated, nrow = 1, labels = c("d")) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.d
ggsave("./results/charts/paneld.pdf",panel.d,width=7,height=4,dpi=300)

panel.cd <- ggarrange(panel.c,panel.d,widths = c(1,2)) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.cd
ggsave("./results/charts/panel cd.pdf",panel.cd,width=12,height=4,dpi=300)


#### ----------------------------------WES----------------------------------#### 
## e. the consistency of variant genes
### whole genes
FUSCCTNBC_JI <- read.csv("./results/tables/FUSCCTNBC_JI.csv",row.names = 1)
FUSCCTNBC_JI_forplot <- melt(FUSCCTNBC_JI)
FUSCCTNBC_JI_forplot$variable <- factor(FUSCCTNBC_JI_forplot$variable, 
                                        levels = c("Final.Mutation","Filter.Mutation",
                                                   "VarScan","TNscope","TNseq"))

FUSCCTNBC_JI_plo <- ggplot(FUSCCTNBC_JI_forplot,aes(x=variable, y=value,
                                                    color=variable,fill=variable))+
  geom_jitter(alpha = 0.5,size = 0.8)+
  geom_boxplot(alpha = 0.5,width = 0.8,
               outlier.shape = NA)+
  scale_fill_npg()+
  scale_color_npg()+
  theme_bw()+mytheme+
  labs(title="", x="", y="Jaccard Index")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylim(0,1);FUSCCTNBC_JI_plo
ggsave('./results/charts/FUSCCTNBC_JI_boxplot.png',FUSCCTNBC_JI_plo,width=4,height=3)

panel.e <- ggarrange(FUSCCTNBC_JI_plo, nrow = 1, labels = c("e")) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.e
ggsave("./results/charts/panel e.pdf",panel.e,width=4,height=4,dpi=300)

## f. the consistency of all mutated genes
## Known hotspot cancer-related genes from the article
cancer_genes <- c("TP53","PIK3CA","KMT2C","PTEN","RB1","KMT2D","PIK3R1","FAT3","FOXA1",
                  "ATR","NF1","ATK1","ATRX","NOTCH2","APC","NCOR1","HRAS","ERBB2","CBFB")
# percentile <- c(74,18,7,6,4,4,4,4,3,3,3,3,3,3,3,3,2,2,1)
# cancer_genes_info <- data.frame(cancer_genes = cancer_genes, percentile=percentile)
# write.csv(cancer_genes_info,"./results/tables/cancer_genes_info.csv")

compare_genes_info <- read.csv("./results/tables/compare_genes_info.csv",row.names = 1) %>% as.data.frame()
compare_genes_info <- compare_genes_info[-which(compare_genes_info$GeneMF.x <= 0.01 | compare_genes_info$GeneMF.y <= 0.01),]
for_label <- compare_genes_info[compare_genes_info$Hugo_Symbol %in% cancer_genes,]

p1 <- ggscatter(compare_genes_info, x = "GeneMF.x", y = "GeneMF.y",
                add = "reg.line", conf.int = TRUE, 
                add.params = list(linetype=2, color = "lightblue", fill = "lightgray"))+
  # stat_cor(method = "pearson",alternative = "two.sided",r.digits = 4,
  #          p.digits = 4,label.x = 0.4, label.y = 0.2)+
  labs(title="", x="Frequency (New)", y="Frequency (Old)")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")+mytheme;p1
# ggsave('./results/charts/genes_MF_scatter.pdf',p1,width=4,height=4)
# R = 0.9316, p = 2.2e-16
p2 <- ggpar(p1,xscale = "log10",yscale = "log10");p2
p3 <- p2 + geom_label_repel(data = for_label, aes(color=factor(Hugo_Symbol),
                                                  label = Hugo_Symbol),
                            alpha = 0.7, 
                            arrow = arrow(length=unit(0.01, "npc")),
                            max.overlaps = 30);p3

ggsave('./results/charts/total_genes_with_hotspotgens_scatter.png',p3,width=4,height=4)

panel.f <- ggarrange(p3, nrow = 1, labels = c("f")) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.f
ggsave("./results/charts/panel f.pdf",panel.f,width=4,height=4,dpi=300)

panel.ef<- ggarrange(FUSCCTNBC_JI_plo,p3,ncol = 2,labels = c("e","f")) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.ef
ggsave("./results/charts/panel ef.pdf",panel.ef,width=12,height=4,dpi=300)

#### ----------------------------------CNV----------------------------------#### 
# the consistency of CNV Peaks(venn.diagram)
amp_inter_venn <- read.csv("./data/amp_inter_venn.csv",header = T)
del_inter_venn <- read.csv("./data/del_inter_venn.csv",header = T)

data1 <- list()
data1$Amp_peak_new <- amp_inter_venn[,4]
data1$Amp_peak_old <- amp_inter_venn[,2]

data2 <- list()
data2$Del_peak_new <- del_inter_venn[,4]
data2$Del_peak_old <- del_inter_venn[,2]

venn.plot <- venn.diagram(
  x = list(Amp_peak_new = data1$Amp_peak_new,
           Amp_peak_old = data1$Amp_peak_old),
  filename = NULL,
  fill = c("#be482c", "#be1b32"),
  alpha = 0.6,
  label.col = "white",
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("#be482c", "#be1b32"),
  cat.cex = 2,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  margin = 0.05,
  cat.dist = c(0.03, 0.03),
  cat.pos = c(360, 180),
  rotation.degree=90
)
p8 <- as.ggplot(plot_grid(grobTree(venn.plot)));p8

venn.plot <- venn.diagram(
  x = list(Del_peak_new = data2$Del_peak_new,
           Del_peak_old = data2$Del_peak_old),
  filename = NULL,
  fill = c("#64999f", "#93b18d"),
  alpha = 0.6,
  label.col = "white",
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("#64999f", "#93b18d"),
  cat.cex = 2,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  margin = 0.05,
  cat.dist = c(0.03, 0.03),
  cat.pos = c(360, 180),
  rotation.degree=90
)
p9 <- as.ggplot(plot_grid(grobTree(venn.plot)));p9

panel.g <- ggarrange(p8,p9, nrow = 1, labels = c("g"),widths = c(1,1)) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.g
ggsave("./results/charts/panel g.pdf",panel.g,width=12,height=6,dpi=300)

## the consistency of CNV genes (venn.diagram)
# new_gene_cnv <- read.csv("./results/tables/new_gene_cnv.csv",header = T,row.names = 1)
# old_gene_cnv <- read.csv("./results/tables/old_gene_cnv.csv",header = T,row.names = 1)

# data <- list()
# data$newgene <- rownames(new_gene_cnv)
# data$oldgene <- rownames(old_gene_cnv)
# 
# venn.plot <- venn.diagram(
#   x = list(new_gene = data$newgene,
#            old_gene = data$oldgene),
#   filename = NULL,
#   fill = c("#f39b7f", "#3c5488"),
#   alpha = 0.6,
#   label.col = "white",
#   cex = 1.5,
#   fontfamily = "serif",
#   fontface = "bold",
#   cat.col = c("#f39b7f", "#3c5488"),
#   cat.cex = 2,
#   cat.fontfamily = "serif",
#   cat.fontface = "bold",
#   margin = 0.05,
#   cat.dist = c(0.03, 0.03),
#   cat.pos = c(0, 180),
#   rotation.degree=270
# );grid.draw(venn.plot)
# p10 <- as.ggplot(plot_grid(grobTree(venn.plot)));p10
# ggsave("./results/charts/p10.pdf",p10,width=6,height=6,dpi=300)

## the consistency of CNV genes (Jaccard Index)
CNV_ji <- read.csv("./results/tables/CNV_ji.csv",header = T,row.names = 1,)
CNV_ji_forplot <-melt(CNV_ji)
CNV_ji_forplot$variable <- factor(CNV_ji_forplot$variable, 
                                        levels = c("Total.Genes","Netural","Gain",
                                                   "Loss"))
## Comparison of old and new processes
FUSCCTNBC_CNV_JI_plo <- ggplot(CNV_ji_forplot,aes(x=variable, y=value,
                                              color=variable,fill=variable))+
  # geom_jitter(alpha = 0.5,size = 1)+
  geom_boxplot(alpha = 0.5,width = 0.8,
               outlier.shape = NA)+
  scale_fill_aaas()+
  scale_color_aaas()+
  theme_bw()+mytheme+
  labs(title="", x="", y="Jaccard Index")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")+theme(axis.text.x = element_text(angle = 45, hjust = 1,size =12))+theme(axis.text.y = element_text(size =12))+
  ylim(0,1);FUSCCTNBC_CNV_JI_plo


panel.h <- ggarrange(FUSCCTNBC_CNV_JI_plo, nrow = 1, labels = c("h"),widths = c(2,4)) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.h
ggsave("./results/charts/panel.h.png",panel.h,width=8,height=6,dpi=300)
ggsave("./results/charts/panel.h.pdf",panel.h,width=8,height=6,dpi=300)

panel.gh <- ggarrange(panel.g, panel.h,nrow = 1, widths = c(2,3)) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));panel.gh

fig4 <- ggarrange(panel.ab,panel.cd,panel.ef,panel.gh,nrow = 4,heights = c(1,1,1,1)) +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"));fig4
ggsave("./results/charts/fig4_0807.pdf",fig4,width=16,height=12,dpi=300)

## remove all
rm(list=ls())
