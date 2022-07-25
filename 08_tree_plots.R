# tree plots

# Regression analyses 2: Dredging models
# rm(list=ls())
# setwd("~/Desktop/WCVP_special_issue/Eve_MyrtalesPAFTOL/myrtales")
setwd("2022_myrtales/")
#################################################################################################
library(phylolm)
library(ape)
library(phytools)
library(MuMIn)
require(reshape2)
require(ggplot2)
require(ggplotify)
require(gridExtra)
require(ggtree)
require(aplot)
require(cowplot)

# Tree
tree <- read.tree("tree/myrtales_pruned.tre")
data <- readRDS("datasets/Myrtales_full_dataset.Rdata")
tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% rownames(data)])
data <- data[match(tree$tip.label, rownames(data)),]
#tree$tip.label <- gsub("^[^_]*_","", tree$tip.label)
# plot(tree)
# axisPhylo()

family_names <- unlist(lapply(strsplit(tree$tip.label, "_"), function(x) x[[1]]))
unique(family_names)
pal <- hcl.colors(6, palette = "Viridis", alpha = 0.7)
names(pal) <- c("Myrtaceae", "Onagraceae", "Combretaceae", "Melastomataceae_CAPclade", "Vochysiaceae", "Lythraceae")
to_rename <- unique(family_names)[!unique(family_names) %in% names(pal)]

family_colors <- viridis::viridis(length(unique(family_names)))
names(family_colors) <- unique(family_names)

a <- ggtree(tree) +
  ggtitle("a) Myrtales Phylogeny")

plot_data <- data.frame(id = tree$tip.label, div_rate = exp(data$div_rate_eps0.9), niche = exp(data$niche_through_time), clade = family_names)
plot_data$clade[plot_data$clade %in% to_rename] <- "Melastomataceae_CAPclade"
plot_data$clade <- as.character(plot_data$clade)

cols <- pal[match(plot_data$clade, names(pal))]

b <- ggplot(plot_data, aes(x = id, y = div_rate, color = clade)) + 
  geom_col(width = .00001) +
  scale_color_manual(values = cols) +
  coord_flip() + 
  theme_tree2() + 
  ylab("Div rate") +
  ggtitle("b) Diversification rates") +
  theme(legend.position='none')

c <- ggplot(plot_data, aes(x = id, y = niche, color = clade)) + 
  geom_col(width = .00001) +
  scale_color_manual(values = cols) +
  coord_flip() + 
  theme_tree2() + 
  ylab("Niche thingy") +
  ggtitle("c) Niche Bredth") +
  theme(legend.position='none')

# bc <- grid.arrange(b, c, nrow = 1)
# bc <- plot_grid(b, c, nrow = 1)

ab <- as.grob(b %>% insert_left(a, width = 3))
plot(ab)

ac <- as.grob(c %>% insert_left(a, width = 3))
plot(ac)

