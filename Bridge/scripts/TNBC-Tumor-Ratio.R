#' Title: TNBC Tumor Ratio
#' Author: Qingwang Chen
#' Date: 2022-07-23
#' Version: V1.0
#' 

# Set work directory and import libraries
rm(list = ls (all = TRUE))
source("./scripts/libraries.R")

# data import
tnbc_tumor_ratio <- read.table("./data/TNBC_tumor_ratio.txt")
colnames(tnbc_tumor_ratio) <- tnbc_tumor_ratio[1,]
tnbc_tumor_ratio <- tnbc_tumor_ratio[-1,]
tnbc_tumor_ratio$tumor_purity <- as.numeric(tnbc_tumor_ratio$tumor_purity)

tnbc_tumor_ratio$num <- sapply(strsplit(as.character(tnbc_tumor_ratio$sample_id),"V"),function(x){x[2]})
tnbc_tumor_ratio$num <- str_pad(tnbc_tumor_ratio$num,3,side = "left",0)
tnbc_tumor_ratio$sample_id <- paste0("FUSCCTNBC",tnbc_tumor_ratio$num)

write.csv(tnbc_tumor_ratio,"./results/tables/TNBC_tumor_ratio.txt",quote = F,row.names = F)


Tes_Clusters <- read.csv("./results/tables/FUSCCTNBC_old_subtypes.csv")
colnames(Tes_Clusters) <- c("patient","subtype")

mes_patients <- Tes_Clusters$patient[which(Tes_Clusters$subtype=="MES")]
setdiff(mes_patients,tnbc_tumor_ratio$sample_id)
tnbc_mes_tumor_ratio <- tnbc_tumor_ratio[which(tnbc_tumor_ratio$sample_id %in% mes_patients),]
write.csv(tnbc_mes_tumor_ratio,"./results/tables/TNBC_mes_tumor_ratio.txt",quote = F,row.names = F)

table(tnbc_mes_tumor_ratio$tumor_purity)
low_tp_mes_patients <- tnbc_mes_tumor_ratio$sample_id[which(tnbc_mes_tumor_ratio$tumor_purity<75)]

### 
# tumor purity statistics by four subtypes
tnbc_tumor_ratio <- fread("./results/tables/TNBC_tumor_ratio.txt") %>% as.data.frame()
Tes_Clusters <- read.csv("./results/tables/FUSCCTNBC_subtypes.csv")
Tes_Clusters$X <- NULL
colnames(Tes_Clusters)[1] <- "sample_id"

tnbc_tumor_ratio_sub <- left_join(Tes_Clusters,tnbc_tumor_ratio) 
tnbc_tumor_ratio_sub <- tnbc_tumor_ratio_sub[,-c(4,6)]
write.csv(tnbc_tumor_ratio_sub,"./results/tables/TNBC_315samples_tumor_ratio.csv",quote = F,row.names = F)

tnbc_tumor_ratio_sub <- drop_na(tnbc_tumor_ratio_sub)

lev <- c("BLIS","IM","MES","LAR")
subtype_pal <- c("#ef4922","#7bc242","#3baade","#8f80ba")
tnbc_tumor_ratio_sub$subtype <- factor(tnbc_tumor_ratio_sub$subtype,levels = lev)
# 315 samples
aggregate(tnbc_tumor_ratio_sub$tumor_purity, by = list(tnbc_tumor_ratio_sub$subtype), FUN = sd)

#------------------------plot---------------------------------#
my_comparisons <- list(c("BLIS", "IM"),
                       c("BLIS", "MES"),
                       c("BLIS", "LAR"),
                       c("IM", "MES"),
                       c("IM", "LAR"),
                       c("MES", "LAR"))

#  tumor purity by four subtypes(Boxplot)
ggboxplot(tnbc_tumor_ratio_sub,x = "subtype",y = "tumor_purity",
          color = "subtype", palette = subtype_pal,add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,label="p.signif",
                     label.y=100)+   # Pairwise comparison against all
  stat_compare_means(method = "anova", label.y=110)+ # Add global p-value
  xlab("Subtype")+ylab("Tumor Purity (%)")

ggsave('./results/charts/TumorPurityByFourSubtypes (Boxplot).png',width=4,height=3)
ggsave('./results/charts/TumorPurityByFourSubtypes (Boxplot).pdf',width=5,height=4)
