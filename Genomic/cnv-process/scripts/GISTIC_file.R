#Eliminate false positive segments in Ascat segment result files
#Correction of results to GISTIC2.0 input format

segment_ascat <- read.csv("./data/segment_ascat_raw20220107.csv",header = T)
chas_segment_falseposit <- false_posit 
segment_ascat_cut <- read.table("./data/segment_ascat.txt",header = T)

chas_segment_falseposit <- na.omit(chas_segment_falseposit)
chas_segment_falsepositive <- chas_segment_falseposit[which(chas_segment_falseposit$count >= 12),]
chr <- unique(segment_ascat$chr)

temp <- rbind(segment_ascat_cut[,1:3],segment_ascat[,1:3])
intersect_result <- temp[duplicated(temp), ]
test <- as.matrix(duplicated(temp))

segment_ascat$index <- c(1:nrow(segment_ascat))
n <- 1
list_index <- list()

for (i in chr)
{
  ascati <- segment_ascat[which(segment_ascat$chr == i),]
  false_positi <- chas_segment_falsepositive[which(chas_segment_falsepositive$Chromosome == i),]
  cat("the chromosome id is",i,'\n')
  for (k in (1:nrow(ascati)))
  {
    cat("segment id",k,'\n')
    select_seg <- ascati[k,c(3,4,9)]
    flagseccess <- 0
    l <- 1
    while(l<nrow(false_positi)&&flagseccess == 0)
    {
      if(false_positi[l,8] >= select_seg[1,2]&&false_positi[l,7] <= select_seg[1,1])
      {list_index[n] <- select_seg[1,3];n <- n+1;flagseccess <- flagseccess+1;break}
      else
      {l <- l+1;next}
    }
  }
}

segment_ascat_true <- segment_ascat[-match(as.matrix(list_index),segment_ascat$index),]
segment_ascat_true$cpnb <- segment_ascat_true$nAraw+segment_ascat_true$nBraw

segment_ascat_true[which(segment_ascat_true$cpnb == 0),10] <- 0.05
segment_ascat_true$count_probe <- NA
chrom <- unique(segment_ascat_true[,2])
SNP <- probeset[,1:3]

temp_result <- matrix(NA,nrow(segment_ascat_true),1)
m <- 1
for(i in chrom)
{
  segi <- segment_ascat_true[which(segment_ascat_true$chr == i),]
  snpi <- SNP[which(SNP$Chromosome == i),]
  
  for(k in (1:nrow(segi)))
  {
    temp_result[m,] <- length(which(snpi[,3] >= segi[k,3]&snpi[,3] <= segi[k,4]))
    m <- m+1
  }
}
segment_ascat_true$count_probe <- temp_result
segment_ascat_true$cpnb1 <- log2(segment_ascat_true$cpnb)-1


SNP[which(SNP$Chromosome == 25),2] <- "Y"
SNP[which(SNP$Chromosome == 24),2] <- "X"
write.table(SNP,"./data/markersfile220226.txt",col.names = F,row.names = F,quote = F,sep = "\t")
write.table(seg_ascat_ture[,c(1:4,11,12)],"./data/segmentationfile220307.txt",col.names = F,row.names = F,quote = F,sep = "\t")