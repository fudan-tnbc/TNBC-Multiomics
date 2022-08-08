#After obtaining CHAS results, using log2ratio of segments to remove false positive probesets

library(data.table)
probeset <- as.data.frame(fread("./data/result.probeset.txt"))
segment <- as.data.frame(fread("./data/Result.segment.txt"))

probset_wbc <- probeset[,c(1:3,grep("WBC",colnames(probeset)))]
probset_wbc_log2 <- probset_wbc[,c(1:3,grep("Log2Ratio",colnames(probset_wbc)))]
Probset_wbc_log2 <- probset_wbc_log2[,-grep("WeightedLog2Ratio ",colnames(probset_wbc_log2))]

segment_WBC <- segment[grep("WBC",segment[,12]),]
segment_WBC_ttCN <- segment_WBC[which(segment_WBC$Type == "TotalCN"),]
segment_WBC_ttCN_copy <- segment_WBC_ttCN
segment_WBC_ttCN_copy[,8] <- as.numeric(segment_WBC_ttCN_copy[,8])
datasetid <- as.matrix(as.numeric(substr(colnames(Probset_wbc_log2)[4:26],21,23)))
chrom <- c(1:22,24,25)
Probset_wbc_log2_cpseg <- Probset_wbc_log2

output_Result <- Probset_wbc_log2[1,]
for (i in chrom)
{
  output <- Probset_wbc_log2[which(Probset_wbc_log2$Chromosome == i),]    
  segmenti <- segment_WBC_ttCN_copy[which(segment_WBC_ttCN_copy$Chromosome == i),]
  for (k in (1:23))
  {
    segmentk <- segmenti[grep(datasetid[k],segmenti[,12]),]
    for (j in (1:nrow(segmentk)))
    {output[which(output[,3] >= segmentk[j,3]&output[,3] <= segmentk[j,4]),(k+3)] <- as.numeric(segmentk[j,8])}
  }
  output_Result <- rbind(output_Result,output)
}
output_Result <- output_Result[-1,]

cutoff_count <- as.matrix(rowSums(abs(output_Result[,-c(1:3)]) > 0.1))
false_posi <- output_Result[which(cutoff_count >= 12),]
probset_select <- probeset[-match(false_posi[,1],probeset[,1]),]

##准备输入ascat的文件
prob_bafdata <- probset_select[,grep("BAF",colnames(probset_select))]
all_log2 <- probset_select[,grep("Log2Ratio ",colnames(probset_select))]
prob_log2 <- all_log2[,-grep("WeightedLog2Ratio",colnames(all_log2))]

tumor_baf <- cbind(probset_select[,1:3],prob_bafdata[,-grep("WBC",colnames(prob_bafdata))])
tumor_log2 <- cbind(probset_select[,1:3],prob_log2[,-grep("WBC",colnames(prob_log2))])
datasetidtumor <- as.matrix(substr(colnames(tumor_log2)[-(1:3)],12,23))
colnames(tumor_baf) <- c("ProbeSetName","Chromosome","Position",datasetidtumor)
colnames(tumor_log2) <- c("ProbeSetName","Chromosome","Position",datasetidtumor)

write.table(tumor_log2,"./data/Tumor_LogR.txt",row.names = F,col.names = T,sep = "\t")
write.table(tumor_baf,"./data/Tumor_BAF.txt",row.names = F,col.names = T,sep = "\t")