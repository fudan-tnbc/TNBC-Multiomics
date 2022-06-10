#Counting false positive segments using Result.segment profile from CHAS results
segment <- as.data.frame(fread("./data/Result.segment.txt"))

segment_WBC <- segment[grep("WBC",segment[,12]),]
segment_WBC_ttCN <- segment_WBC[which(segment_WBC$Type == "TotalCN"),]
segment_WBC_ttCN$sampleid <- as.numeric(substr(segment_WBC_ttCN[,12],10,12))
segment_select_orig <- segment_WBC_ttCN[,c(2:4,8,13)]
segment_select_orig[,4] <- as.numeric(segment_select_orig[,4])
segment_select_orig_log2cut <- segment_select_orig[which(abs(segment_select_orig[,4]) >= 0.1),]

segment_select_orig_log2cut$segment_length <- segment_select_orig_log2cut[,3]-segment_select_orig_log2cut[,2]+1
segment_select_orig_log2cut$start2 <- segment_select_orig_log2cut[,2]-floor(0.5*segment_select_orig_log2cut$segment_length)
segment_select_orig_log2cut$stop2 <- segment_select_orig_log2cut[,3]+floor(0.5*segment_select_orig_log2cut$segment_length)
segment_select_orig_log2cut[which(segment_select_orig_log2cut[,7] <= 0),7] <- 1

result_mat <- as.data.frame(matrix(NA,7725,10))
colnames(result_mat) <- c(colnames(segment_select_orig_log2cut),"count")
m <- 1
segment_select_orig_log2cut$index <- 0

for (i in chrom)
{
  segmenti <- segment_select_orig_log2cut[which(segment_select_orig_log2cut$Chromosome == i),]
  iter_mat <- segmenti
  cat("the Chromosome id is",i,'\n')
  while(nrow(iter_mat) > 1)
  {
    select_seg <- iter_mat[1,c(7,8,6,9)]
    result_mat[m,1:9] <- iter_mat[1,]
    datasetid <- as.matrix(unique(iter_mat[,5]))
    dataset_other <- setdiff(datasetid ,as.matrix(iter_mat[1,5]))
    flagseccess <- 0
    count <- 0
    list_index <- list()
    n <- 1
    list_index[n] <- select_seg[1,4];n <- n+1;
    for(j in dataset_other)
    {
      cat("other datasets",j,'\n')
      l <- 1
      set <- iter_mat[which(iter_mat[,5] == j),]
      flagseccess <- 0
      while (l<nrow(set)&&flagseccess == 0)
      {
        cat("other segments",l,'\n')
        if(set[l,2] >= select_seg[1,1]&&set[l,3] <= select_seg[1,2])
        {count <- count+1;
        list_index[n] <- set[l,9];n <- n+1;
        flagseccess <- flagseccess+1;
        break}
        else
        {l <- l+1;next}
      }
    }
    
    result_mat[m,10] <- count
    m <- m+1
    flagseccess <- 0
    iter_mat <- iter_mat[-(match(as.matrix(list_index),iter_mat[,9])),]
  }
}

false_posit <- result_mat