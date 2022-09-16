# Regression analyses 2: Dredging models
# rm(list=ls())
setwd("~/Desktop/Pubs_inprep/WCVP_special_issue/Eve_MyrtalesPAFTOL/myrtales")
#setwd("~/2022_myrtales/")
#################################################################################################
library(phylolm)
library(ape)
library(phytools)
library(MuMIn)
library(corHMM)
library(ggplot2)
library(gridExtra)

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
phylanova_results <- phylANOVA(tree_pruned, habitat, corolla_diam)

plot_corolla50 <- ggplot(test1, aes(x=main_habitat, y=exp(fm_scoring_corolla_diam), fill=main_habitat)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.25, alpha=0.8) + 
  scale_fill_manual( values = c("#66a182","#edae49")) +
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=0.5, y=max(exp(corolla_diam)), size=4, hjust=0, label=paste0("p=",phylanova_results$Pf)) + 
  xlab("Habitat") +
  ylab("Mean corolla diameter (mm)") +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.x = element_text(colour = 'black', size = 10)) 


# 40% cutoff
test1 <- subset_master_table[c("fm_scoring_corolla_diam","cutoff_40")]
test1 <- subset(test1, !is.na(test1$fm_scoring_corolla_diam))
test1 <- subset(test1, !is.na(test1$cutoff_40))

corolla_diam <- test1$fm_scoring_corolla_diam
habitat <- test1$cutoff_40
names(corolla_diam) <- names(habitat) <- rownames(test1)

tree_pruned <- keep.tip(tree, rownames(test1))
phylanova_results <- phylANOVA(tree_pruned, habitat, corolla_diam)

plot_corolla40 <- ggplot(test1, aes(x=cutoff_40, y=exp(fm_scoring_corolla_diam), fill=cutoff_40)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.25, alpha=0.8) + 
  #stat_compare_means(comparisons = 3) + # Add pairwise comparisons p-value
  #geom_jitter(colour = 2, alpha=0.3, size=0.9) +
  #coord_cartesian(ylim = c(0, 1)) +
  #coord_flip()  + 
  scale_fill_manual( values = c("#66a182","#edae49")) +
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=0.5, y=max(exp(corolla_diam)), size=4, hjust=0, label=paste0("p=",phylanova_results$Pf)) + 
  xlab("Habitat") +
  ylab("Mean corolla diameter (mm)") +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.x = element_text(colour = 'black', size = 10)) 


# 60% cutoff
test1 <- subset_master_table[c("fm_scoring_corolla_diam","cutoff_60")]
test1 <- subset(test1, !is.na(test1$fm_scoring_corolla_diam))
test1 <- subset(test1, !is.na(test1$cutoff_60))

corolla_diam <- test1$fm_scoring_corolla_diam
habitat <- test1$cutoff_60
names(corolla_diam) <- names(habitat) <- rownames(test1)

tree_pruned <- keep.tip(tree, rownames(test1))
phylanova_results <- phylANOVA(tree_pruned, habitat, corolla_diam)

plot_corolla60 <- ggplot(test1, aes(x=cutoff_60, y=exp(fm_scoring_corolla_diam), fill=cutoff_60)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.25, alpha=0.8) + 
  #stat_compare_means(comparisons = 3) + # Add pairwise comparisons p-value
  #geom_jitter(colour = 2, alpha=0.3, size=0.9) +
  #coord_cartesian(ylim = c(0, 1)) +
  #coord_flip()  + 
  scale_fill_manual( values = c("#66a182","#edae49")) +
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=0.5, y=max(exp(corolla_diam)), size=4, hjust=0, label=paste0("p=",phylanova_results$Pf)) + 
  xlab("Habitat") +
  ylab("Mean corolla diameter (mm)") +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.x = element_text(colour = 'black', size = 10)) 


pdf("plots/corolla_diam_habitat.pdf",height=3.5,width=5)
grid.arrange(plot_corolla40, plot_corolla50, plot_corolla60, ncol=3, nrow = 1)
dev.off()

##############################################################
##############################################################
### Habitat vs. seed length
test2 <- subset_master_table[c("seed.length.mean","main_habitat")]
test2 <- subset(test2, !is.na(test2$seed.length.mean))
test2 <- subset(test2, !is.na(test2$main_habitat))

seed_length <- test2$seed.length.mean
habitat <- test2$main_habitat
names(seed_length) <- names(habitat) <- rownames(test2)

tree_pruned <- keep.tip(tree, rownames(test2))
phylanova_results <- phylANOVA(tree_pruned, habitat, seed_length)

plot_seedlength50 <- ggplot(test2, aes(x=main_habitat, y=exp(seed_length), fill=main_habitat)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.25, alpha=0.8) + 
  scale_fill_manual( values = c("#66a182","#edae49")) +
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=0.5, y=max(exp(seed_length)), size=4, hjust=0, label=paste0("p=",phylanova_results$Pf)) + 
  xlab("Habitat") +
  ylab("Mean seed length (mm)") +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.x = element_text(colour = 'black', size = 10)) 

# 40% cutoff
test2 <- subset_master_table[c("seed.length.mean","cutoff_40")]
test2 <- subset(test2, !is.na(test2$seed.length.mean))
test2 <- subset(test2, !is.na(test2$cutoff_40))

seed_length <- test2$seed.length.mean
habitat <- test2$cutoff_40
names(seed_length) <- names(habitat) <- rownames(test2)

tree_pruned <- keep.tip(tree, rownames(test2))
phylanova_results <- phylANOVA(tree_pruned, habitat, seed_length)

plot_seedlength40 <- ggplot(test2, aes(x=cutoff_40, y=exp(seed_length), fill=cutoff_40)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.25, alpha=0.8) + 
  scale_fill_manual( values = c("#66a182","#edae49")) +
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=0.5, y=max(exp(seed_length)), size=4, hjust=0, label=paste0("p=",phylanova_results$Pf)) + 
  xlab("Habitat") +
  ylab("Mean seed length (mm)") +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.x = element_text(colour = 'black', size = 10)) 

# 60% cutoff
test2 <- subset_master_table[c("seed.length.mean","cutoff_60")]
test2 <- subset(test2, !is.na(test2$seed.length.mean))
test2 <- subset(test2, !is.na(test2$cutoff_60))

seed_length <- test2$seed.length.mean
habitat <- test2$cutoff_60
names(seed_length) <- names(habitat) <- rownames(test2)

tree_pruned <- keep.tip(tree, rownames(test2))
phylanova_results <- phylANOVA(tree_pruned, habitat, seed_length)

plot_seedlength60 <- ggplot(test2, aes(x=cutoff_60, y=exp(seed_length), fill=cutoff_60)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.25, alpha=0.8) + 
  scale_fill_manual( values = c("#66a182","#edae49")) +
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=0.5, y=max(exp(seed_length)), size=4, hjust=0, label=paste0("p=",phylanova_results$Pf)) + 
  xlab("Habitat") +
  ylab("Mean seed length (mm)") +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.x = element_text(colour = 'black', size = 10)) 

pdf("plots/seed_length_habitat.pdf",height=3.5,width=5)
grid.arrange(plot_seedlength40, plot_seedlength50, plot_seedlength60, ncol=3, nrow = 1)
dev.off()

###############################
### Habitat vs. seed number
test3 <- subset_master_table[c("fm_scoring_seed_number","main_habitat")]
test3 <- subset(test3, !is.na(test3$fm_scoring_seed_number))
test3 <- subset(test3, !is.na(test3$main_habitat))

seed_number <- test3$fm_scoring_seed_number
habitat <- test3$main_habitat
names(seed_number) <- names(habitat) <- rownames(test3)

tree_pruned <- keep.tip(tree, rownames(test3))
phylanova_results<-phylANOVA(tree_pruned, habitat, seed_number)

plot_seednumber50 <- ggplot(test3, aes(x=main_habitat, y=exp(fm_scoring_seed_number), fill=main_habitat)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.25, alpha=0.8) + 
  scale_fill_manual( values = c("#66a182","#edae49")) +
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=0.5, y=max(exp(seed_number)), size=4, hjust=0, label=paste0("p=",phylanova_results$Pf)) + 
  xlab("Habitat") +
  ylab("Seed number") +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.x = element_text(colour = 'black', size = 10)) 


# 40% cutoff
test3 <- subset_master_table[c("fm_scoring_seed_number","cutoff_40")]
test3 <- subset(test3, !is.na(test3$fm_scoring_seed_number))
test3 <- subset(test3, !is.na(test3$cutoff_40))

seed_number <- test3$fm_scoring_seed_number
habitat <- test3$cutoff_40
names(seed_number) <- names(habitat) <- rownames(test3)

tree_pruned <- keep.tip(tree, rownames(test3))
phylanova_results<-phylANOVA(tree_pruned, habitat, seed_number)

plot_seednumber40 <- ggplot(test3, aes(x=cutoff_40, y=exp(fm_scoring_seed_number), fill=cutoff_40)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.25, alpha=0.8) + 
  scale_fill_manual( values = c("#66a182","#edae49")) +
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=0.5, y=max(exp(seed_number)), size=4, hjust=0, label=paste0("p=",phylanova_results$Pf)) + 
  xlab("Habitat") +
  ylab("Seed number") +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.x = element_text(colour = 'black', size = 10)) 

# 60% cutoff
test3 <- subset_master_table[c("fm_scoring_seed_number","cutoff_60")]
test3 <- subset(test3, !is.na(test3$fm_scoring_seed_number))
test3 <- subset(test3, !is.na(test3$cutoff_60))

seed_length <- test3$fm_scoring_seed_number
habitat <- test3$cutoff_60
names(seed_length) <- names(habitat) <- rownames(test3)

tree_pruned <- keep.tip(tree, rownames(test3))
phylanova_results <- phylANOVA(tree_pruned, habitat, seed_length)

plot_seednumber60 <- ggplot(test3, aes(x=cutoff_60, y=exp(fm_scoring_seed_number), fill=cutoff_60)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.25, alpha=0.8) + 
  scale_fill_manual( values = c("#66a182","#edae49")) +
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=0.5, y=max(exp(seed_number)), size=4, hjust=0, label=paste0("p=",phylanova_results$Pf)) + 
  xlab("Habitat") +
  ylab("Seed number") +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.x = element_text(colour = 'black', size = 10)) 

pdf("plots/seed_number_habitat.pdf",height=3.5,width=5)
grid.arrange(plot_seednumber40, plot_seednumber50, plot_seednumber60, ncol=3, nrow = 1)
dev.off()

###############################
### Habitat vs. fruit type
pal <- hcl.colors(5, palette = "Viridis", alpha = 0.7)

test4 <- subset_master_table[c("fm_scoring_fruit","main_habitat")]
test4 <- subset(test4, !is.na(test4$fm_scoring_fruit))
test4 <- subset(test4, !is.na(test4$main_habitat))
test4$species <- rownames(test4)

dataset_traits <- test4[,c("species","fm_scoring_fruit","main_habitat")]
tree_pruned <- keep.tip(tree, rownames(test4))

t0 <- ggplot(dataset_traits, aes(x=fm_scoring_fruit, fill=main_habitat)) + 
  geom_bar(position = "fill", alpha=0.8)+
  theme_bw(base_size = 8) +
  theme(legend.position = "none") + 
  scale_fill_manual(values = c("#66a182","#edae49")) +
  #annotate(geom="text", x=0.5, y=max(corolla_diam), size=4, hjust=0, label=paste0("p=",phylanova_results$Pf)) + 
  xlab("Fruit type") +
  ylab("Proportion in habitat") +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.x = element_text(colour = 'black', size = 10)) 


corhmm_fits <- corHMM:::fitCorrelationTest(tree_pruned, dataset_traits) 
save(corhmm_fits, file = "results/h3/corhmm_fits.Rsave")
load("results/h3/corhmm_fits.Rsave")
corhmm_tbl <- corHMM:::getModelTable(corhmm_fits)

# conducting a lrt
teststat <- -2 * (corhmm_tbl$lnLik[1] - corhmm_tbl$lnLik[2])
p.val50 <- pchisq(teststat, df = 1, lower.tail = FALSE)
print(p.val)

# 40% cutoff
test4 <- subset_master_table[c("fm_scoring_fruit","cutoff_40")]
test4 <- subset(test4, !is.na(test4$fm_scoring_fruit))
test4 <- subset(test4, !is.na(test4$cutoff_40))
test4$species <- rownames(test4)

dataset_traits <- test4[,c("species","fm_scoring_fruit","cutoff_40")]
tree_pruned <- keep.tip(tree, rownames(test4))

t_40cutoff <- ggplot(dataset_traits, aes(x=fm_scoring_fruit, fill=cutoff_40)) + 
  geom_bar(position = "fill", alpha=0.8)+
  theme_bw(base_size = 8) +
  theme(legend.position = "none") + 
  scale_fill_manual(values = c("#66a182","#edae49")) +
  #annotate(geom="text", x=0.5, y=max(corolla_diam), size=4, hjust=0, label=paste0("p=",phylanova_results$Pf)) + 
  xlab("Fruit type") +
  ylab("Proportion in habitat") +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.x = element_text(colour = 'black', size = 10)) 


corhmm_fits <- corHMM:::fitCorrelationTest(tree_pruned, dataset_traits) 
save(corhmm_fits, file = "results/h3/corhmm_fits_cutoff40.Rsave")

load("results/h3/corhmm_fits_cutoff40.Rsave")
corhmm_tbl <- corHMM:::getModelTable(corhmm_fits)

# conducting a lrt
teststat <- -2 * (corhmm_tbl$lnLik[1] - corhmm_tbl$lnLik[2])
p.val40 <- pchisq(teststat, df = 1, lower.tail = FALSE)
print(p.val40)

# 60% cutoff
test4 <- subset_master_table[c("fm_scoring_fruit","cutoff_60")]
test4 <- subset(test4, !is.na(test4$fm_scoring_fruit))
test4 <- subset(test4, !is.na(test4$cutoff_60))
test4$species <- rownames(test4)

dataset_traits <- test4[,c("species","fm_scoring_fruit","cutoff_60")]
tree_pruned <- keep.tip(tree, rownames(test4))

t_60cutoff <- ggplot(dataset_traits, aes(x=fm_scoring_fruit, fill=cutoff_60)) + 
  geom_bar(position = "fill", alpha=0.8)+
  theme_bw(base_size = 8) +
  theme(legend.position = "none") + 
  scale_fill_manual(values = c("#66a182","#edae49")) +
  #annotate(geom="text", x=0.5, y=max(corolla_diam), size=4, hjust=0, label=paste0("p=",phylanova_results$Pf)) + 
  xlab("Fruit type") +
  ylab("Proportion in habitat") +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.x = element_text(colour = 'black', size = 10)) 


corhmm_fits <- corHMM:::fitCorrelationTest(tree_pruned, dataset_traits) 
save(corhmm_fits, file = "results/h3/corhmm_fits_cutoff60.Rsave")

load("results/h3/corhmm_fits_cutoff60.Rsave")
corhmm_tbl <- corHMM:::getModelTable(corhmm_fits)

# conducting a lrt
teststat <- -2 * (corhmm_tbl$lnLik[1] - corhmm_tbl$lnLik[2])
p.val60 <- pchisq(teststat, df = 1, lower.tail = FALSE)

pdf("plots/fruittype_habitat.pdf" ,height=3.5,width=5)
grid.arrange(t_40cutoff, t0, t_60cutoff, ncol=3, nrow = 1)
dev.off()

sink("results/h3/p_values.txt")
cat("p value 40% cutoff")
print(round(p.val40,4))
cat("p value 50% cutoff")
print(round(p.val50,4))
cat("p value 60% cutoff")
print(round(p.val60,4))
sink()



