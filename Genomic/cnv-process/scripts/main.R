#After obtaining CHAS results, using log2ratio of segments to remove false positive probesets
source("./scripts/FP_prob.R")


#ASCAT
source("./scripts/Ascat.R")


#Counting false positive segments using Result.segment profile from CHAS results
source("./scripts/FP_seg.R")

#Eliminate false positive segments in Ascat segment result files
#Correction of results to GISTIC2.0 input format
source("./scripts/GISTIC_file.R")