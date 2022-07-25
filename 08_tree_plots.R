# tree plots
# rm(list=ls())
# setwd("~/Desktop/WCVP_special_issue/Eve_MyrtalesPAFTOL/myrtales")
setwd("2022_myrtales/")
#################################################################################################
library(ape)
<<<<<<< HEAD

=======
library(phytools)
library(MuMIn)
require(reshape2)
require(ggplot2)
require(ggplotify)
require(gridExtra)
require(ggtree)
require(aplot)
require(cowplot)
>>>>>>> 5defe376c44f4ef289d0a9e623b1bef59b007752

# Tree
tree <- read.tree("tree/myrtales_pruned.tre")
data <- readRDS("datasets/Myrtales_full_dataset.Rdata")
tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% rownames(data)])
data <- data[match(tree$tip.label, rownames(data)),]
#tree$tip.label <- gsub("^[^_]*_","", tree$tip.label)
# plot(tree)
# axisPhylo()
<<<<<<< HEAD
master_table <- readRDS("datasets/Myrtales_full_dataset.Rdata") 


#pdf("mu_plot.pdf", width=4, height=5)

  one_tree <- keep.tip(tree, row.names(master_table))
  one_label <- "Myrtales tree"

  #if(nrow(one_set)>0) {
    
    tree_pruned <- one_tree
    # Small corrections so that tips appear in the right order in the plot
    is_tip <- tree_pruned$edge[,2] <= length(tree_pruned$tip.label)
    ordered_tips <- tree_pruned$edge[is_tip, 2]
    right_order <- as.character(tree_pruned$tip.label[ordered_tips])
    
    # Organizing so tip rates are in the same order as tips of the tree
    cleaned_table <- master_table[match(as.character(right_order), row.names(master_table)),]

      # Plotting
      color_breaks = 4 
      layout.matrix <- matrix(c(1,2), nrow = 1, ncol = 2)
      layout(mat = layout.matrix,widths = c(2,1))
      par(mar=c(5,0.5,3,0.5))
      # Tree
      plot(tree_pruned, show.tip.label=T, edge.width=0.2, adj=1, cex=0.1)
      title(main=one_label)
      ape::axisPhylo()
      
      palette_name1 = "Inferno"
      
      divrates <- as.numeric(cleaned_table$div_rate_eps0.9)
      divrates <- exp(divrates)
      rounded_rates <- round(divrates, color_breaks)
      pal <- hcl.colors(length(levels(as.factor(rounded_rates))), palette = palette_name1, alpha = 0.75)
      pal <- pal[match(rounded_rates, as.numeric(levels(as.factor(rounded_rates))))] 
      x <- 1:length(divrates)
      plot(divrates, x,  lwd = 0.2, xlim=range(c(min(divrates), max(divrates))),
           pch=19, yaxt = "n", xlab="div_rates", ylab="", frame.plot=T, cex=0.75, col=pal)
      segments(min(divrates), 1:length(divrates), divrates[1:length(divrates)], 1:length(divrates), col= pal,lwd = 0.2)
      

dev.off()


# Figure 3: results from global model relative variable importance

div_global <- read.csv("results/h2/param_table_div.csv")

div_global$X[which(div_global$X=="CHELSA_bio10_17")]
div_global$X[which(div_global$X=="most_common_life_formwoody")]
div_global$X[which(div_global$X=="CHELSA_bio10_02")]
div_global$X[which(div_global$X=="meanpH")]
div_global$X[which(div_global$X=="depthtobedrock2")]
div_global$X[which(div_global$X=="GLOBAL_SLOPE_10MIN")]
div_global$X[which(div_global$X=="meancarbon")]
div_global$X[which(div_global$X=="fm_scoring_fruitFleshy")]
div_global$X[which(div_global$X=="meanwatercap")]
div_global$X[which(div_global$X=="fm_scoring_seed_number")]
div_global$X[which(div_global$X=="seed.length.mean")]
div_global$X[which(div_global$X=="CHELSA_bio10_11")]
div_global$X[which(div_global$X=="fm_scoring_corolla_diam")]




# Figure 4: phyloanovas and corhmm




=======

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
>>>>>>> 5defe376c44f4ef289d0a9e623b1bef59b007752

