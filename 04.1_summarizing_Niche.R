
# rm(list=ls())
setwd("~/Desktop/WCVP_special_issue/Eve_MyrtalesPAFTOL/myrtales")
library(matrixStats)
library(ape)

all_niche_files <- list.files("datasets/Niche",".csv$", full.names = T)
all_niches <- lapply(all_niche_files, read.csv)
all_niches <- do.call(rbind, all_niches)

# Minimum volume for imputation
min_vol <- all_niches
min_vol <- subset(min_vol, min_vol$Vol!=0)
min_vol <- subset(min_vol, !is.na(min_vol$Vol))
min_vol <- min(min_vol$Vol) # for imputation
all_niches$Vol[which(all_niches$Vol==0)] <- min_vol
all_niches$Vol[which(is.na(all_niches$Vol))] <- min_vol
#-------------------------------

genera_per_row <- gsub("\\ .*","", all_niches$X)

# Matching with genera in tree and summarizing
tree <- read.tree("tree/myrtales_pruned.tre")
genera <- gsub("^[^_]*_","", tree$tip.label)

all_results <- matrix(ncol=ncol(all_niches), nrow=0)
for(i in 1:length(genera)){
  genus <- genera[i]
  one_subset <- all_niches[grep(genus, genera_per_row),]
  one_subset <- subset(one_subset, !is.na(one_subset$CHELSA_bio10_01))
  if(nrow(one_subset)>0){
  n <- sum(one_subset$N)
  all_medians <- as.data.frame(t(colMedians(as.matrix(one_subset[,-c(1,2,3)]))))
  vol <- median(subset(one_subset$Vol, !is.na(one_subset$Vol)))
  all_results <- rbind(all_results, cbind(genus,n,vol, all_medians))
  }
}
colnames(all_results) <- colnames(all_niches)
all_results <- all_results[,-c(4:6)]
colnames(all_results)[c(1,2)] <- c("genus","n_points")

write.csv(all_results, file="myrtales_niche_summaries.csv",row.names=F)
#View(all_results)

