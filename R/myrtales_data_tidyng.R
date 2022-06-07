
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
tree <- ladderize(tree)
tips_to_drop_table <- read.csv("./tree/tips_to_drop_final.csv", h=T)
tips_to_drop <- as.character(tips_to_drop_table[tips_to_drop_table$to_drop=="x",][,"tip"])
tree_pruned <- tidyng.Myrtales.tree(tree, tips_to_drop)
write.tree(tree_pruned, file="./tree/myrtales_pruned.tre")

# Load most up to date datasets
# traits <- read.csv("./datasets/fruit_example.csv")
traits <- read.csv("./datasets/Myrtales_master_table.csv")
match_tables <- read.csv("match_trait_tree_datasets.csv")

# adding missing tips:
#traits$trait_species <- NA
#for(i in 1:nrow(traits)) {
#  traits$trait_species[i] <- paste(c(traits$Family[i], traits$Genus[i], traits$Species[i]), collapse="_")
#}
#
#match_tables$trait_species <- NA
#for(i in 1:nrow(match_tables)) {
#  match_tables$trait_species[i] <- paste(c(match_tables$Family[i], match_tables$Genus[i], match_tables$Species[i]), collapse="_")
#}
#write.csv(merged_tables, file="merged_tables.csv", row.names=F)

trait_dataset_curated <- read.csv("adjusting_datasets/merged_tables_curated_final.csv")
trait_dataset_curated <- subset(trait_dataset_curated, trait_dataset_curated$match_to_tree %in% tips_to_drop_table$tip)
to_keep <- subset(tips_to_drop_table, tips_to_drop_table$to_drop=="keep")
trait_dataset_curated <- subset(trait_dataset_curated, trait_dataset_curated$match_to_tree %in% to_keep$tip)

add_table <- as.data.frame(matrix(nrow=6, ncol=ncol(trait_dataset_curated)))
colnames(add_table) <- colnames(trait_dataset_curated)
add_table$match_to_tree <- to_keep$tip[!(to_keep$tip %in% trait_dataset_curated$match_to_tree)]
trait_dataset_curated <- rbind(trait_dataset_curated, add_table)

write.csv(trait_dataset_curated, file="2022-25-04_Myrtales-master-table.csv", row.names = F)
# remove duplicated as in 
# trait_dataset_curated$match_to_tree[duplicated(trait_dataset_curated$match_to_tree)]
# Load back and adjust names according to tree

final_dataset <- read.csv("2022-25-04_Myrtales-master-table.csv")
tips <- final_dataset$match_to_tree
tips[grep("Pimenta_pseudocaryophyllus",tips)] <- "Myrtaceae_Pseudocaryophyllus"
tips <- sapply(strsplit(tips, "_"), function(x) {
  g <- seq_along(x)
  g[g < 1] <- 1
  g[g > 2 ] <- 2
  paste(tapply(x, g, paste, collapse = "-"), collapse = "_")
})
tips <- gsub("\\-.*", "",tips)
final_dataset$match_to_tree <- tips
final_dataset$Genus.b <- unlist(lapply(strsplit(tips, "_"), "[[", 2))

write.csv(final_dataset, file="2022-25-04_Myrtales-master-table.csv", row.names = F)

# including life_form
life_form <- read.csv("most_common_life_form.csv")
life_form <- subset(life_form, life_form$genera %in% final_dataset$Genus)
write.csv(life_form, "2022-04-25_most_common_life_form.csv")

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



