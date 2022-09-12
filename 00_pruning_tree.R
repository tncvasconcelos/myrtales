
tidyng.Myrtales.tree <- function(tree, tips_to_drop) {
  # Remove tips from "tips to drop" file
  tree_pruned <- ape::drop.tip(tree, tips_to_drop)
  # Remove family name from tips
  tmp <- strsplit(tree_pruned$tip.label, "_")
  tips <- c()
  for(i in 1:length(tmp)){
   tip0 <- tmp[[i]]
    family <- tip0[1]
    genus <- tip0[2]
    if(i == grep("Pimenta_pseudocaryophyllus",tree_pruned$tip.label)) {
      tips[i] <- paste(family, "Pseudocaryophyllus", sep="_")
    } else {
    if(!is.na(tip0[3])){
      species <- tip0[3:length(tip0)]
      tips[i] <- paste(family, genus, paste(species, collapse="_"), sep="_")
    } else {
      species <- ""
      tips[i] <- paste0(family,"_", genus, species)
    }
    }
  }
  tips <- sapply(strsplit(tips, "_"), function(x) {
    g <- seq_along(x)
    g[g < 1] <- 1
    g[g > 2 ] <- 2
    paste(tapply(x, g, paste, collapse = "-"), collapse = "_")
  })
  tree_pruned$tip.label <- tips
  # Removing epiphet
  tree_pruned$tip.label <- gsub("\\-.*", "",tree_pruned$tip.label)
  return(tree_pruned)
}


#-------------------------------------------------
#-------------------------------------------------
#-------------------------------------------------
# Set wd as the repo
setwd("~/Desktop/Pubs_inprep/WCVP_special_issue/Eve_MyrtalesPAFTOL/myrtales/")
# rm(list=ls())

# Load packages
library(ape)
library(phytools)
library(picante)

# Load tree
tree <- read.nexus("./tree/Myrtales_dated_Aug2022.tre")
tree <- ladderize(tree)
tips_to_drop_table <- read.csv("./tree/tips_to_drop_final.csv", h=T)
tips_to_drop <- as.character(tips_to_drop_table[tips_to_drop_table$to_drop=="x",][,"tip"])
tree_pruned <- tidyng.Myrtales.tree(tree, tips_to_drop)

# fixing typos
tree_pruned$tip.label[which(tree_pruned$tip.label=="Myrtaceae_Astereomyrtus")] <- "Myrtaceae_Asteromyrtus"
tree_pruned$tip.label[which(tree_pruned$tip.label=="Melastomataceae_Potheranthera")] <- "Melastomataceae_Poteranthera"

#
write.tree(tree_pruned, file="./tree/myrtales_pruned.tre")



