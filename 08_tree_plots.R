# tree plots
# rm(list=ls())
setwd("~/Desktop/Pubs_inprep/WCVP_special_issue/Eve_MyrtalesPAFTOL/myrtales")
# setwd("2022_myrtales/")
#################################################################################################
library(ape)
library(phytools)
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
#master_table <- readRDS("datasets/Myrtales_full_dataset.Rdata") 


# #pdf("mu_plot.pdf", width=4, height=5)
# 
#   one_tree <- keep.tip(tree, row.names(master_table))
#   one_label <- "Myrtales tree"
# 
#   #if(nrow(one_set)>0) {
#     
#     tree_pruned <- ladderize(one_tree)
#     
#     # Small corrections so that tips appear in the right order in the plot
#     is_tip <- tree_pruned$edge[,2] <= length(tree_pruned$tip.label)
#     ordered_tips <- tree_pruned$edge[is_tip, 2]
#     right_order <- as.character(tree_pruned$tip.label[ordered_tips])
#     
#     # Organizing so tip rates are in the same order as tips of the tree
#     cleaned_table <- master_table[match(as.character(right_order), row.names(master_table)),]
#     tip_names <- rownames(cleaned_table)
#     tip_names <- gsub("_",") ", tip_names)
#     tip_names <- paste0("(",tip_names)
#     write.csv(tip_names, file="tip_names.csv", row.names=F)
#       
#       # Plotting
#       color_breaks = 4 
#       layout.matrix <- matrix(c(1,2), nrow = 1, ncol = 2)
#       layout(mat = layout.matrix,widths = c(2,1))
#       par(mar=c(5,0.5,3,0.5))
#       # Tree
#       plot(tree_pruned, show.tip.label=T, edge.width=0.2, adj=1, cex=0.01)
#       
#      
#       title(main=one_label)
#       ape::axisPhylo()
#       
#       palette_name1 = "Inferno"
#       
#       divrates <- as.numeric(cleaned_table$div_rate_eps0.9)
#       divrates <- exp(divrates)
#       rounded_rates <- round(divrates, color_breaks)
#       pal <- hcl.colors(length(levels(as.factor(rounded_rates))), palette = palette_name1, alpha = 0.75)
#       pal <- pal[match(rounded_rates, as.numeric(levels(as.factor(rounded_rates))))] 
#       x <- 1:length(divrates)
#       plot(divrates, x,  lwd = 0.2, xlim=range(c(min(divrates), max(divrates))),
#            pch=19, yaxt = "n", xlab="div_rates", ylab="", frame.plot=T, cex=0.75, col=pal)
#       segments(min(divrates), 1:length(divrates), divrates[1:length(divrates)], 1:length(divrates), col= pal,lwd = 0.2)
#       
# 
# dev.off()


# Figure 3: results from global model relative variable importance
# 
# div_global <- read.csv("results/h2/param_table_div.csv")
# 
# div_global$X[which(div_global$X=="CHELSA_bio10_17")]
# div_global$X[which(div_global$X=="most_common_life_formwoody")]
# div_global$X[which(div_global$X=="CHELSA_bio10_02")]
# div_global$X[which(div_global$X=="meanpH")]
# div_global$X[which(div_global$X=="depthtobedrock2")]
# div_global$X[which(div_global$X=="GLOBAL_SLOPE_10MIN")]
# div_global$X[which(div_global$X=="meancarbon")]
# div_global$X[which(div_global$X=="fm_scoring_fruitFleshy")]
# div_global$X[which(div_global$X=="meanwatercap")]
# div_global$X[which(div_global$X=="fm_scoring_seed_number")]
# div_global$X[which(div_global$X=="seed.length.mean")]
# div_global$X[which(div_global$X=="CHELSA_bio10_11")]
# div_global$X[which(div_global$X=="fm_scoring_corolla_diam")]
# 
# 
# 
# 
# # Figure 4: phyloanovas and corhmm
# 


family_names <- unlist(lapply(strsplit(tree$tip.label, "_"), function(x) x[[1]]))
unique(family_names)
pal <- hcl.colors(6, palette = "Viridis", alpha = 1)
names(pal) <- c("Myrtaceae", "Onagraceae", "Combretaceae", "Melastomataceae_CAPclade", "Vochysiaceae", "Lythraceae")
to_rename <- unique(family_names)[!unique(family_names) %in% names(pal)]

family_colors <- viridis::viridis(length(unique(family_names)))
names(family_colors) <- unique(family_names)

#tip_names <- tree$tip.label
#tip_names <- gsub("_",") ", tip_names)
#tip_names <- paste0("(",tip_names)
#tree$tip.label <- tip_names

tree$tip.label <- gsub(".*_", "",tree$tip.label)


a <- ggtree(tree) +
  geom_tiplab(geom="text", size=0.3) +
  ggtitle("a) Myrtales Phylogeny")


plot_data <- data.frame(id = tree$tip.label, div_rate = data$div_rate_eps0.9, niche_expansion = exp(data$niche_through_time),most_common_life_form=data$Most.Common.Life.Form, fruit_type=data$Dry.or.Fleshy.Fruit, seed_length=exp(data$Mean.Seed.Length), seed_number=exp(data$Mean.Seed.Number.per.Fruit), corolla_diam=exp(data$Mean.Corolla.Diameter), niche=exp(data$Vol), main_habitat=data$main_habitat, clade = family_names)

plot_data$clade[plot_data$clade %in% to_rename] <- "Melastomataceae_CAPclade"
plot_data$clade <- as.character(plot_data$clade)

cols <- pal[match(plot_data$clade, names(pal))]

b <- ggplot(plot_data, aes(x = id, y = div_rate, color = clade)) + 
  geom_col(width = .0001) +
  scale_color_manual(values = cols) +
  coord_flip() + 
  theme_tree2() + 
  ylab("Div rate") +
  ggtitle("b) Diversification rates") +
  theme(legend.position='none')

c <- ggplot(plot_data, aes(x = id, y = niche, color = clade)) + 
  geom_col(width = .0001) +
  scale_color_manual(values = cols) +
  coord_flip() + 
  theme_tree2() + 
  ylab("Niche thingy") +
  ggtitle("c) Niche Bredth") +
  theme(legend.position='none')

d <- ggplot(plot_data, aes(x = id, y = niche_expansion, color = clade)) + 
  geom_col(width = .0001) +
  scale_color_manual(values = cols) +
  coord_flip() + 
  theme_tree2() + 
  ylab("Niche expansion") +
  ggtitle("d) Niche Bredth/time") +
  theme(legend.position='none')

e <- ggplot(plot_data, aes(x = id, y = seed_length, color = clade)) + 
  geom_col(width = .0001) +
  scale_color_manual(values = cols) +
  coord_flip() + 
  theme_tree2() + 
  ylab("length (mm)") +
  ggtitle("e) Seed length") +
  theme(legend.position='none')

f <- ggplot(plot_data, aes(x = id, y = seed_number, color = clade)) + 
  geom_col(width = .0001) +
  scale_color_manual(values = cols) +
  coord_flip() + 
  theme_tree2() + 
  ylab("number") +
  ggtitle("f) Seed number") +
  theme(legend.position='none')

g <- ggplot(plot_data, aes(x = id, y = corolla_diam, color = clade)) + 
  geom_col(width = .0001) +
  scale_color_manual(values = cols) +
  coord_flip() + 
  theme_tree2() + 
  ylab("diam mm") +
  ggtitle("g) Corolla diam") +
  theme(legend.position='none')

h <- ggplot(plot_data, aes(x = id, y = 1, fill = factor(fruit_type),color = fruit_type)) + 
  geom_point(aes(shape = fruit_type, color = clade, stroke=0.001), size = 0.8) + 
  scale_shape_manual(values = c(0, 15)) + 
  scale_color_manual(values = cols) +
  coord_flip() + 
  theme_tree2() + 
  ggtitle("h) Fruit type") +
  theme(legend.position='none')

i <- ggplot(plot_data, aes(x = id, y = 1, fill = factor(most_common_life_form),color = most_common_life_form)) + 
  geom_point(aes(shape = most_common_life_form, color = clade, stroke=0.001), size = 0.8) + 
  scale_shape_manual(values = c(0, 15)) + 
  scale_color_manual(values = cols) +
  coord_flip() + 
  theme_tree2() + 
  ggtitle("i) Life form") +
  theme(legend.position='none')

j <- ggplot(plot_data, aes(x = id, y = 1, fill = factor(main_habitat),color = main_habitat)) + 
  geom_point(aes(shape = main_habitat, color = clade, stroke=0.001), size = 0.8) + 
  scale_shape_manual(values = c(15, 0)) + 
  scale_color_manual(values = cols) +
  coord_flip() + 
  theme_tree2() + 
  ggtitle("j) Habitat") +
  theme(legend.position='none')

# PLOTS
pdf("plots/figure1div.pdf", height=8, width=5)
ab <- as.grob(b %>% insert_left(a, width = 2))
plot(ab)
dev.off()

pdf("plots/figure1niche1.pdf", height=8, width=5)
ac <- as.grob(c %>% insert_left(a, width = 2))
plot(ac)
dev.off()

pdf("plots/figure1niche2.pdf", height=8, width=5)
ad <- as.grob(d %>% insert_left(a, width = 2))
plot(ad)
dev.off()

pdf("plots/figure1seedlength.pdf", height=8, width=5)
ae <- as.grob(e %>% insert_left(a, width = 2))
plot(ae)
dev.off()

pdf("plots/figure1seednumber.pdf", height=8, width=5)
af <- as.grob(f %>% insert_left(a, width = 2))
plot(af)
dev.off()

pdf("plots/figure1corolla.pdf", height=8, width=5)
ag <- as.grob(g %>% insert_left(a, width = 2))
plot(ag)
dev.off()

pdf("plots/figure1fruittype.pdf", height=8, width=5)
ah <- as.grob(h %>% insert_left(a, width = 2))
plot(ah)
dev.off()

pdf("plots/figure1lifeform.pdf", height=8, width=5)
ai <- as.grob(i %>% insert_left(a, width = 2))
plot(ai)
dev.off()

pdf("plots/figure1habitat.pdf", height=8, width=5)
aj <- as.grob(j %>% insert_left(a, width = 2))
plot(aj)
dev.off()

