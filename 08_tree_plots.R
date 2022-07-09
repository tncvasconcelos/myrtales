# tree plots

# Regression analyses 2: Dredging models
# rm(list=ls())
setwd("~/Desktop/WCVP_special_issue/Eve_MyrtalesPAFTOL/myrtales")
#################################################################################################
library(phylolm)
library(ape)
library(phytools)
library(MuMIn)

# Tree
tree <- read.tree("tree/myrtales_pruned.tre")
#tree$tip.label <- gsub("^[^_]*_","", tree$tip.label)
# plot(tree)
# axisPhylo()
