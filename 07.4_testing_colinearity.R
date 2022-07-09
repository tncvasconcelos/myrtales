# Regression analyses 2: Dredging models
# rm(list=ls())
setwd("~/Desktop/WCVP_special_issue/Eve_MyrtalesPAFTOL/myrtales")
#################################################################################################
library(phylolm)
library(ape)
library(phytools)
library(MuMIn)
library(car)


# Tree
tree <- read.tree("tree/myrtales_pruned.tre")

# Master table
master_table <- readRDS("datasets/Myrtales_full_dataset.Rdata") 
#master_table_by_clade <- readRDS("datasets/Myrtales_by_clade_dataset.Rdata") 

############################################
############################################
# Key traits and environment are correlated and determine niche breadth and diversification rates of extant taxa at a global scale. 

############################################
# Preparing table
subset_master_table <- master_table
subset_master_table$seed.length.mean <- NA
# Getting a mean value for seed length
for(i in 1:nrow(subset_master_table)) {
  seed_min <- subset_master_table$seed.length.min[i]
  seed_max <- subset_master_table$seed.length.max[i]
  subset_master_table$seed.length.mean[i] <- mean(c(seed_min, seed_max))
}

subset_master_table$most_common_life_form[which(subset_master_table$most_common_life_form=="annual")] <- "herbaceous"
subset_master_table$most_common_life_form[which(subset_master_table$most_common_life_form=="woody perennial")] <- "woody"
subset_master_table$most_common_life_form[which(subset_master_table$most_common_life_form=="herbaceous perennial")] <- "herbaceous"
subset_master_table$most_common_life_form[which(subset_master_table$most_common_life_form=="epiphyte")] <- "herbaceous"

############################################
# Excluding variables with colinearity problems
# Selecting traits we would keep in the largest model
div_model_var_to_keep <- c("most_common_life_form","fm_scoring_fruit","seed.length.mean",
                           "fm_scoring_seed_number","fm_scoring_corolla_diam","CHELSA_bio10_01",
                           "CHELSA_bio10_02","CHELSA_bio10_10","CHELSA_bio10_11","CHELSA_bio10_12",
                           "CHELSA_bio10_15","CHELSA_bio10_16","CHELSA_bio10_17","GLOBAL_SLOPE_10MIN",
                           "depthtobedrock2","meancarbon","meanpH","meanwatercap")

subset_master_table_test <- subset_master_table[,c("div_rate_eps0.9",div_model_var_to_keep)]
test_cor <- lm(div_rate_eps0.9~., data=subset_master_table_test)
vif(test_cor) # vif test

# testing removing vars
div_model_vars_to_exclude <- names(vif(test_cor))[which(vif(test_cor)>10)]
div_model_vars_to_include <- names(vif(test_cor))[which(vif(test_cor)<10)]