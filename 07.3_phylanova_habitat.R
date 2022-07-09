# Regression analyses 2: Dredging models
# rm(list=ls())
setwd("~/Desktop/WCVP_special_issue/Eve_MyrtalesPAFTOL/myrtales")
#################################################################################################
library(phylolm)
library(ape)
library(phytools)
library(MuMIn)
library(corHMM)

# Tree
tree <- read.tree("tree/myrtales_pruned.tre")
#tree$tip.label <- gsub("^[^_]*_","", tree$tip.label)
# plot(tree)
# axisPhylo()

# Master table
master_table <- readRDS("datasets/Myrtales_full_dataset.Rdata") 
master_table_by_clade <- readRDS("datasets/Myrtales_by_clade_dataset.Rdata") 
############################################
# Hypothesis 3: 
# The occupation of habitats depend on intrinsic traits related on both survival and reproduction. For example, wind pollinated species and/or those with larger, more conspicuous flowers are more species-rich in open habitats. Species with fleshy fruits and/or species with fruits with larger, fewer seeds are more species-rich in closed habitats. 

subset_master_table <- master_table
subset_master_table$seed.length.mean <- NA

# Getting a mean value for seed length
for(i in 1:nrow(subset_master_table)) {
  seed_min <- subset_master_table$seed.length.min[i]
  seed_max <- subset_master_table$seed.length.max[i]
  subset_master_table$seed.length.mean[i] <- mean(c(seed_min, seed_max))
}
#subset_master_table$seed.trade.off <- subset_master_table$seed.length.mean/subset_master_table$fm_scoring_seed_number
#x <- lm(log(subset_master_table$seed.length.mean)~log(subset_master_table$fm_scoring_seed_number))
#plot(log(subset_master_table$fm_scoring_seed_number)~ log(subset_master_table$seed.length.mean))

# colnames(master_table)
subset_master_table <- subset(subset_master_table, subset_master_table$div_rate_eps0.9!=0) #removing 0s before loging
# transforming temperature in Kelvin before logging
subset_master_table$CHELSA_bio10_01 <- subset_master_table$CHELSA_bio10_01+273.15
subset_master_table$CHELSA_bio10_10 <- subset_master_table$CHELSA_bio10_10+273.15
subset_master_table$CHELSA_bio10_11 <- subset_master_table$CHELSA_bio10_11+273.15

# Logging all continuous vars
subset_master_table$div_rate_eps0.9 <- log(subset_master_table$div_rate_eps0.9)
subset_master_table$Vol <- log(subset_master_table$Vol)

# abiotic:
subset_master_table$CHELSA_bio10_01 <- log(subset_master_table$CHELSA_bio10_01)
subset_master_table$CHELSA_bio10_10 <- log(subset_master_table$CHELSA_bio10_10)
subset_master_table$CHELSA_bio10_11 <- log(subset_master_table$CHELSA_bio10_11)
subset_master_table$CHELSA_bio10_12 <- log(subset_master_table$CHELSA_bio10_12)
subset_master_table$CHELSA_bio10_15 <- log(subset_master_table$CHELSA_bio10_15)
subset_master_table$CHELSA_bio10_16 <- log(subset_master_table$CHELSA_bio10_16)
subset_master_table$CHELSA_bio10_17 <- log(subset_master_table$CHELSA_bio10_17)
subset_master_table$GLOBAL_SLOPE_10MIN <- log(subset_master_table$GLOBAL_SLOPE_10MIN)
subset_master_table$depthtobedrock2 <- log(subset_master_table$depthtobedrock2)
subset_master_table$meancarbon <- log(subset_master_table$meancarbon)
subset_master_table$meanpH <- log(subset_master_table$meanpH)
subset_master_table$meanwatercap <- log(subset_master_table$meanwatercap)

# traits:
subset_master_table$fm_scoring_seed_number <- log(subset_master_table$fm_scoring_seed_number)
subset_master_table$seed.length.mean <- log(subset_master_table$seed.length.mean)
subset_master_table$fm_scoring_corolla_diam <- log(subset_master_table$fm_scoring_corolla_diam)

# Building phyloanova?
test1 <- subset_master_table[c("fm_scoring_corolla_diam","main_habitat")]
test1 <- subset(test1, !is.na(test1$fm_scoring_corolla_diam))
test1 <- subset(test1, !is.na(test1$main_habitat))

corolla_diam <- test1$fm_scoring_corolla_diam
habitat <- test1$main_habitat
names(corolla_diam) <- names(habitat) <- rownames(test1)

tree_pruned <- keep.tip(tree, rownames(test1))
phylANOVA(tree_pruned, habitat, corolla_diam)
boxplot(corolla_diam~habitat)


###
test2 <- subset_master_table[c("seed.length.mean","main_habitat")]
test2 <- subset(test2, !is.na(test2$seed.length.mean))
test2 <- subset(test2, !is.na(test2$main_habitat))

seed_length <- test2$seed.length.mean
habitat <- test2$main_habitat
names(seed_length) <- names(habitat) <- rownames(test2)

tree_pruned <- keep.tip(tree, rownames(test2))
phylANOVA(tree_pruned, habitat, seed_length)
boxplot(seed_length~habitat)

###############################
colnames(subset_master_table)
test3 <- subset_master_table[c("fm_scoring_fruit","main_habitat")]
test3 <- subset(test3, !is.na(test3$fm_scoring_fruit))
test3 <- subset(test3, !is.na(test3$main_habitat))
test3$species <- rownames(test3)

dataset_traits <- test3[,c("species","fm_scoring_fruit","main_habitat")]
#habitat <- test3[,c("species","main_habitat")]

tree_pruned <- keep.tip(tree, rownames(test3))

dependent_model_matrix <- getStateMat4Dat(dataset_traits, collapse = FALSE)$rate.mat 
independent_model_matrix <- getStateMat4Dat(dataset_traits, collapse = FALSE, indep = TRUE)$rate.mat 

# data legend
getStateMat4Dat(dataset_traits, collapse = FALSE)$legend

dependent_model_fit <- corHMM(tree_pruned, dataset_traits, rate.cat = 1, rate.mat = dependent_model_matrix, node.states = "none", collapse = FALSE)

corHMM:::fitCorrelationTest(tree_pruned, dataset_traits) 
corHMM:::fitCorrelationTest
