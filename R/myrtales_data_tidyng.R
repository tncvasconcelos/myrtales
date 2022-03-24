
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
setwd("~/Desktop/WCVP_special_issue/Eve_MyrtalesPAFTOL/myrtales/")
# rm(list=ls())

# Load packages
library(ape)
library(phytools)
library(picante)

# Load tree
tree <- read.tree("./tree/myrtales_dated_final.tre")
tips_to_drop_table <- read.csv("./tree/tips_to_drop_final.csv", h=T)
tips_to_drop <- as.character(tips_to_drop_table[tips_to_drop_table$to_drop=="x",][,"tip"])
tree_pruned <- tidyng.Myrtales.tree(tree, tips_to_drop)
write.tree(tree_pruned, file="./tree/myrtales_pruned.tre")

# Load most up to date datasets
traits <- read.csv("./datasets/fruit_example.csv")

# Making datasets talk to each other
rownames(traits) <- traits[,1]
phy0 <- picante::match.phylo.data(tree_pruned, traits)$phy
trait0 <- picante::match.phylo.data(tree_pruned, traits)$data
mode <- trait0[,2]
names(mode) <- trait0[,1]
colors_states <- c("midnightblue", "goldenrod", "green")


# Plot states at tips
pdf("./plots/trait_fruit.pdf", width= 4, height= 12)
plot(phy0, show.tip.label=T, edge.width=0.2, adj=1, cex=0.08)
par(fg="transparent")
tiplabels(pie=to.matrix(mode, sort(unique(mode))),piecol=colors_states,cex= 0.3,lwd=0.2, frame = "n")
par(fg="black")
#tiplabels(pch=1, bg=tip.cols, adj=1, cex=0.1, width = 0.1)
legend("topleft", legend=sort(unique(mode)), pt.bg = colors_states, pch=21, cex=0.8)
axisPhylo()
dev.off()



