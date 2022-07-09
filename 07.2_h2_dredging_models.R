# Regression analyses 2: Dredging models
# rm(list=ls())
setwd("~/2022_myrtales")
#################################################################################################
library(phylolm)
library(ape)
library(phytools)
library(MuMIn)

# Functions to get r squares and to organize results
get.rqrs <- function(organized_table, full_dataset, phy, dep.var="div_rate_eps0.9") {
  vars <- setdiff(colnames(organized_table),c("df","logLik","AICc","delta","weight"))
  only_var <- organized_table[,vars]
  colnames(full_dataset)[which(colnames(full_dataset)==dep.var)] <- "focal_var"
  rqsrs <- c()
  for(i in 1:nrow(organized_table)) {
    to_include_in_model <- colnames(only_var)[which(only_var[i,] != "")]
    tmp_dataset <- full_dataset[c("focal_var", to_include_in_model)]
    model <- phylolm(focal_var~.,data=tmp_dataset, phy=tree)
    rqsrs[i] <- round(model$r.squared,3)
    cat(i,"\r")
  }
  organized_table$rsqs <- rqsrs
  return(organized_table)
}

organize.table <- function(table_results, thrsh=T) {
  if(thrsh) {
    subset_best_fit <- subset(table_results, table_results$delta < 2)
  } else {
    subset_best_fit <- table_results
  }
  subset_best_fit <- subset_best_fit[,-1]
  vars <- setdiff(colnames(subset_best_fit),c("df","logLik","AICc","delta","weight"))
  only_var <- subset_best_fit[,vars]
  stats <- subset_best_fit[,c("df","logLik","AICc","delta","weight")]
  for(i in 1:nrow(only_var)) {
    only_var[i,][which(!is.na(only_var[i,]))] <- "x"
    only_var[i,][which(is.na(only_var[i,]))] <- ""
  }
  tmp_df <- cbind(only_var, stats)
  tmp_df$logLik <- round(tmp_df$logLik, 1)
  tmp_df$AICc <- round(tmp_df$AICc, 1)
  tmp_df$delta <- round(tmp_df$delta, 2)
  tmp_df$weight <- round(tmp_df$weight, 2)
  rownames(tmp_df) <- paste0("model ", 1:nrow(tmp_df))
  cat(i,"\r")
  return(tmp_df)
}

# Tree
tree <- read.tree("tree/myrtales_pruned.tre")

# Master table
master_table <- readRDS("datasets/Myrtales_full_dataset.Rdata") 

############################################
############################################
# Key traits and environment are correlated and determine niche breadth and diversification rates of extant taxa at a global scale. 

############################################
# Building global model for diversification (div rates as dependent variable)
# here note that we're keeping only CHELSA_bio10_17 and CHELSA_bio10_11 for 
# precipitation and temperature
master_table <- subset(master_table, !is.na(master_table$most_common_life_form))
master_table <- subset(master_table, !is.na(master_table$fm_scoring_fruit))
master_table <- subset(master_table, !is.na(master_table$seed.length.mean))
master_table <- subset(master_table, !is.na(master_table$fm_scoring_seed_number))
master_table <- subset(master_table, !is.na(master_table$fm_scoring_corolla_diam))

model_div_full <- phylolm(div_rate_eps0.9~
                            most_common_life_form+
                            fm_scoring_fruit+
                            seed.length.mean+
                            fm_scoring_seed_number+ 
                            fm_scoring_corolla_diam+
                            CHELSA_bio10_11+
                            CHELSA_bio10_02+
                            CHELSA_bio10_17+
                            GLOBAL_SLOPE_10MIN+
                            depthtobedrock2+
                            meanwatercap+
                            meancarbon+
                            meanpH, data=master_table, phy=tree)

# model_div_full <- phylolm(div_rate_eps0.9~
#                             most_common_life_form,
#                           data=master_table, phy=tree)
# plot(model_div_full)
# Dredging full model for "best" combinations
# dredge_div <- dredge(model_div_full)
# save(dredge_div, file = "results/h2/dredge_div.Rsave")
# coefTableList <- lapply(dredge_div, coefTable)
# write.csv(dredge_div, file="results/h2/dredged_divrate_full.csv", row.names=F)
#----
# dredge_div <- read.csv("results/h2/dredged_divrate_full.csv") 
load("results/h2/dredge_div.Rsave")
subset(dredge_div, delta < 4)
model.avg(dredge_div, subset = delta < 4)
# model.avg(dredge_div, subset = cumsum(weight) <= .95) # get averaged coefficients
summary(get.models(dredge_div, 1)[[1]])


# pdf(file = "h2-results.pdf", height = 20, width = 20)
# plot(dredge_div)
# dev.off()

# dredge_div <- organize.table(dredge_div, thrsh=F)
# dredge_div <- get.rqrs(organized_table=dredge_div, full_dataset=master_table, phy=tree, dep.var="div_rate_eps0.9")
# write.csv(dredge_div, file="results/h2/dredged_divrate_organized_table.csv")

dredge_div[order(dredge_div$rsqs,decreasing = T),]


############################################
############################################
# Building global model for niche breadth (niche breadth  as dependent variable)
model_vol_full <- phylolm(niche_through_time~
                            most_common_life_form+
                            fm_scoring_fruit+
                            seed.length.mean+
                            fm_scoring_seed_number+ 
                            fm_scoring_corolla_diam+
                            CHELSA_bio10_11+
                            CHELSA_bio10_02+
                            CHELSA_bio10_17+
                            GLOBAL_SLOPE_10MIN+
                            depthtobedrock2+
                            meanwatercap+
                            meancarbon+
                            meanpH
                          , data=master_table, phy=tree)

# Dredging full model for "best" combinations
dredge_vol <- dredge(model_vol_full)
write.csv(dredge_vol, file="results/h2/dredged_niche_full.csv", row.names=F)

#----
dredge_vol <- read.csv("results/h2/dredged_niche_full.csv") 
dredge_vol <- organize.table(dredge_vol)
dredge_vol <- get.rqrs(organized_table=dredge_vol, full_dataset=master_table, phy=tree, dep.var="niche_through_time")
write.csv(dredge_vol, file="results/h2/dredge_niche_organized_table.csv")



