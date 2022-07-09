# Regression analyses 2: Dredging models
# rm(list=ls())
setwd("~/Desktop/WCVP_special_issue/Eve_MyrtalesPAFTOL/myrtales")
setwd("~/2022_myrtales/")
#################################################################################################
library(phylolm)
library(ape)
library(phytools)
library(MuMIn)
library(corHMM)

# Tree
tree <- read.tree("tree/myrtales_pruned.tre")

# Master table
master_table <- readRDS("datasets/Myrtales_full_dataset.Rdata") 
#master_table_by_clade <- readRDS("datasets/Myrtales_by_clade_dataset.Rdata") 
############################################
# Hypothesis 3: 
# The occupation of habitats depend on intrinsic traits related on both survival and reproduction. For example, wind pollinated species and/or those with larger, more conspicuous flowers are more species-rich in open habitats. Species with fleshy fruits and/or species with fruits with larger, fewer seeds are more species-rich in closed habitats. 
#plot(sort(subset_master_table$open_canopy))
subset_master_table <- master_table

###############################
### Habitat vs. corolla diameter
test1 <- subset_master_table[c("fm_scoring_corolla_diam","main_habitat")]
test1 <- subset(test1, !is.na(test1$fm_scoring_corolla_diam))
test1 <- subset(test1, !is.na(test1$main_habitat))

corolla_diam <- test1$fm_scoring_corolla_diam
habitat <- test1$main_habitat
names(corolla_diam) <- names(habitat) <- rownames(test1)

tree_pruned <- keep.tip(tree, rownames(test1))
phylANOVA(tree_pruned, habitat, corolla_diam)
boxplot(corolla_diam~habitat)

# 40% cutoff
test1 <- subset_master_table[c("fm_scoring_corolla_diam","cutoff_40")]
test1 <- subset(test1, !is.na(test1$fm_scoring_corolla_diam))
test1 <- subset(test1, !is.na(test1$cutoff_40))

corolla_diam <- test1$fm_scoring_corolla_diam
habitat <- test1$cutoff_40
names(corolla_diam) <- names(habitat) <- rownames(test1)

tree_pruned <- keep.tip(tree, rownames(test1))
phylANOVA(tree_pruned, habitat, corolla_diam)
boxplot(corolla_diam~habitat)

# 60% cutoff
test1 <- subset_master_table[c("fm_scoring_corolla_diam","cutoff_60")]
test1 <- subset(test1, !is.na(test1$fm_scoring_corolla_diam))
test1 <- subset(test1, !is.na(test1$cutoff_60))

corolla_diam <- test1$fm_scoring_corolla_diam
habitat <- test1$cutoff_60
names(corolla_diam) <- names(habitat) <- rownames(test1)

tree_pruned <- keep.tip(tree, rownames(test1))
phylANOVA(tree_pruned, habitat, corolla_diam)
boxplot(corolla_diam~habitat)

###############################
### Habitat vs. seed length
test2 <- subset_master_table[c("seed.length.mean","main_habitat")]
test2 <- subset(test2, !is.na(test2$seed.length.mean))
test2 <- subset(test2, !is.na(test2$main_habitat))

seed_length <- test2$seed.length.mean
habitat <- test2$main_habitat
names(seed_length) <- names(habitat) <- rownames(test2)

tree_pruned <- keep.tip(tree, rownames(test2))
phylANOVA(tree_pruned, habitat, seed_length)
boxplot(seed_length~habitat)

# 40% cutoff
test2 <- subset_master_table[c("seed.length.mean","cutoff_40")]
test2 <- subset(test2, !is.na(test2$seed.length.mean))
test2 <- subset(test2, !is.na(test2$cutoff_40))

seed_length <- test2$seed.length.mean
habitat <- test2$cutoff_40
names(seed_length) <- names(habitat) <- rownames(test2)

tree_pruned <- keep.tip(tree, rownames(test2))
phylANOVA(tree_pruned, habitat, seed_length)
boxplot(seed_length~habitat)

# 60% cutoff
test2 <- subset_master_table[c("seed.length.mean","cutoff_60")]
test2 <- subset(test2, !is.na(test2$seed.length.mean))
test2 <- subset(test2, !is.na(test2$cutoff_60))

seed_length <- test2$seed.length.mean
habitat <- test2$cutoff_60
names(seed_length) <- names(habitat) <- rownames(test2)

tree_pruned <- keep.tip(tree, rownames(test2))
phylANOVA(tree_pruned, habitat, seed_length)
boxplot(seed_length~habitat)

###############################
### Habitat vs. seed number
test3 <- subset_master_table[c("fm_scoring_seed_number","main_habitat")]
test3 <- subset(test3, !is.na(test3$fm_scoring_seed_number))
test3 <- subset(test3, !is.na(test3$main_habitat))

seed_number <- test3$fm_scoring_seed_number
habitat <- test3$main_habitat
names(seed_number) <- names(habitat) <- rownames(test3)

tree_pruned <- keep.tip(tree, rownames(test3))
phylANOVA(tree_pruned, habitat, seed_number)
boxplot(seed_number~habitat)

# 40% cutoff
test3 <- subset_master_table[c("fm_scoring_seed_number","cutoff_40")]
test3 <- subset(test3, !is.na(test3$fm_scoring_seed_number))
test3 <- subset(test3, !is.na(test3$cutoff_40))

seed_number <- test3$fm_scoring_seed_number
habitat <- test3$cutoff_40
names(seed_number) <- names(habitat) <- rownames(test3)

tree_pruned <- keep.tip(tree, rownames(test3))
phylANOVA(tree_pruned, habitat, seed_number)
boxplot(seed_number~habitat)

# 60% cutoff
test3 <- subset_master_table[c("fm_scoring_seed_number","cutoff_60")]
test3 <- subset(test3, !is.na(test3$fm_scoring_seed_number))
test3 <- subset(test3, !is.na(test3$cutoff_60))

seed_length <- test3$fm_scoring_seed_number
habitat <- test3$cutoff_60
names(seed_length) <- names(habitat) <- rownames(test3)

tree_pruned <- keep.tip(tree, rownames(test3))
phylANOVA(tree_pruned, habitat, seed_length)
boxplot(seed_length~habitat)


###############################
### Habitat vs. fruit type
test4 <- subset_master_table[c("fm_scoring_fruit","main_habitat")]
test4 <- subset(test4, !is.na(test4$fm_scoring_fruit))
test4 <- subset(test4, !is.na(test4$main_habitat))
test4$species <- rownames(test4)

dataset_traits <- test4[,c("species","fm_scoring_fruit","main_habitat")]
tree_pruned <- keep.tip(tree, rownames(test4))

corHMM:::fitCorrelationTest(tree_pruned, dataset_traits) 


# library
library(ggplot2)

# create a dataset
specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(specie,condition,value)

# Stacked + percent
ggplot(dataset_traits, aes(fill=fm_scoring_fruit, x=main_habitat)) + 
  geom_bar(stat = "count",
           position = "stack")

?geom_bar
