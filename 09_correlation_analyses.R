# rm(list=ls())
setwd("~/Desktop/WCVP_special_issue/Eve_MyrtalesPAFTOL/myrtales")
#################################################################################################
#---------------------------
#---------------------------
# Correlation analyses
# Let's load everything first
library(phylolm)
library(ape)
library(phytools)
library(MuMIn)

# Tree
tree <- read.tree("tree/myrtales_pruned.tre")
tree$tip.label <- gsub("^[^_]*_","", tree$tip.label)
# plot(tree)
# axisPhylo()

# Traits
master_table <- read.csv("2022-07-06_master_table.csv") 
# Niche
summary_niches <- read.csv("niche_summaries.csv")
colnames(master_table)[which(colnames(master_table)=="Genus")] <- "genus"

# making mega table
combined_table <- merge(master_table, summary_niches, by="genus")
row.names(combined_table) <- combined_table$genus

# clade specifics:
myrtaceae_table <- subset(combined_table, combined_table$Family=="Myrtaceae")
onagraceae_table <- subset(combined_table, combined_table$Family=="Onagraceae")
combretaceae_table <- subset(combined_table, combined_table$Family=="Combretaceae")
melascapclade_table <- subset(combined_table, combined_table$Family%in%c("Melastomataceae","Alzateaceae","Crypteroniaceae","Penaeaceae"))
vochysiaceae_table <- subset(combined_table, combined_table$Family=="Vochysiaceae")
lythraceae_table <- subset(combined_table, combined_table$Family=="Lythraceae")

clade_specific_list <- list(myrtaceae_table, onagraceae_table, combretaceae_table, melascapclade_table, vochysiaceae_table, lythraceae_table)
names(clade_specific_list) <- c("Myrtaceae", "Onagraceae", "Combretaceae", "Melastomataceae_CAPclade","Vochysiaceae","Lythraceae")

##### copied and pasted from abstract:
# In this study we will address four hypotheses: 
# (1) key traits and environment are correlated and determine range size, geographic distribution, and speciation/extinction rates of extant taxa at a global scale. 


# a. niche breadth as a response var

# Global
model_1 <- phylolm(div_rate_eps0~fm_scoring_seed_number+CHELSA_bio10_12+GLOBAL_SLOPE_10MIN, data=combined_table, tree)

model_1 <- phylolm(vol~div_rate_eps0, data=combined_table, tree)
test <- dredge(model_1)


# b. habitat as a response var

x <- phylolm(div_rate_eps0.9~CHELSA_bio10_02+fm_scoring_fruit, data=combined_table, tree)
summary(x)

# c. div rates as a response var 
 

model_1 <- phylolm(div_rate_eps0.9~vol+fm_scoring_seed_number+most_common_life_form, data=combined_table, tree)

summary(model_1)
test <- dredge(model_1)
summary(test)

# (2) species with fleshy fruits and/or species with fruits with larger, fewer seeds are more species-rich in closed habitats. 


# (3) wind pollinated species and/or those with larger, more conspicuous flowers are more species-rich in open habitats. 


# (4) range expansion rates are positively correlated with speciation rates. 

# Clade specific

for(i in 1:length(clade_specific_list)) {
  one_subset <- clade_specific_list[[1]]
  one_clade <- names(clade_specific_list)[1]
  
  
  
}


#----------------------------------------------
# Correlations

traits_continuous <- c("fm_scoring_seed_number")
traits_discrete <- c("most_common_life_form","fm_scoring_fruit","fm_scoring_pollination1","fm_scoring_pollination2")
environment <- c("CHELSA_bio10_01","CHELSA_bio10_12","CHELSA_bio10_15","GLOBAL_SLOPE_10MIN","depthtobedrock2",
                 "meancarbon","meanpH","meanwatercap")
niche_breadth <- c("vol")
div_rates <- c("div_rate_eps0","div_rate_eps0.5","div_rate_eps0.9")

################################################
# (1) Model: div rates ~ niche breadth

combined_table_cor1 = combined_table
combined_table_cor1$niche_expansion <- combined_table_cor1$vol / combined_table_cor1$age # Following that paper, range expansion is range size divided by age of the group
combined_table_cor1 <- subset(combined_table_cor1, combined_table_cor1$div_rate_eps0!=0) # many genera have div rates of 0 because they are monotypic. let's take remove them for one first analysis
combined_table_cor1 <- subset(combined_table_cor1, !is.na(combined_table_cor1$niche_expansion)) 
combined_table_cor1 <- subset(combined_table_cor1, combined_table_cor1$niche_expansion!=0) 

combined_table_cor1$div_rate_eps0 <- log(combined_table_cor1$div_rate_eps0)
combined_table_cor1$niche_expansion <- log(combined_table_cor1$niche_expansion)

pruned_tree <- keep.tip(tree, combined_table_cor1$genus)
# Now let's run the phylolm, a phylogenetic linear model
model_1 <- phylolm(div_rate_eps0~niche_expansion, data=combined_table_cor1, pruned_tree, model="BM")
summary(model_1)


for(i in 1:length(clade_specific_list)) {
  one_subset <- clade_specific_list[[i]]
  one_clade <- names(clade_specific_list)[i]
  
  combined_table_cor1 = one_subset
  combined_table_cor1$niche_expansion <- combined_table_cor1$vol / combined_table_cor1$age # Following that paper, range expansion is range size divided by age of the group
  combined_table_cor1 <- subset(combined_table_cor1, combined_table_cor1$div_rate_eps0!=0) # many genera have div rates of 0 because they are monotypic. let's take remove them for one first analysis
  combined_table_cor1 <- subset(combined_table_cor1, !is.na(combined_table_cor1$niche_expansion)) 
  combined_table_cor1 <- subset(combined_table_cor1, combined_table_cor1$niche_expansion!=0) 
  
  combined_table_cor1$div_rate_eps0 <- log(combined_table_cor1$div_rate_eps0)
  combined_table_cor1$niche_expansion <- log(combined_table_cor1$niche_expansion)
  
  pruned_tree <- keep.tip(tree, combined_table_cor1$genus)
  # Now let's run the phylolm, a phylogenetic linear model
  model_1 <- phylolm(div_rate_eps0~niche_expansion, data=combined_table_cor1, pruned_tree, model="BM")
  summary(model_1)
  
  plot(combined_table_cor1$div_rate_eps0~combined_table_cor1$niche_expansion)
  abline(model_1)

}





################################################
# (2) Model: div rates ~ env factors
combined_table_cor2 = combined_table
combined_table_cor2 <- subset(combined_table_cor2, combined_table_cor2$div_rate_eps0!=0) # many genera have div rates of 0 because they are monotypic. let's take remove them for one first analysis

combined_table_cor2$div_rate_eps0 <- log(combined_table_cor2$div_rate_eps0)
combined_table_cor2$CHELSA_bio10_01 <- log(combined_table_cor2$CHELSA_bio10_01)
combined_table_cor2$CHELSA_bio10_12 <- log(combined_table_cor2$CHELSA_bio10_12)
combined_table_cor2$meanpH <- log(combined_table_cor2$meanpH)
combined_table_cor2$meanwatercap <- log(combined_table_cor2$meanwatercap)
combined_table_cor2$meancarbon <- log(combined_table_cor2$meancarbon)

pruned_tree <- keep.tip(tree, combined_table_cor2$genus)

# Now let's run the phylolm, a phylogenetic linear model
model_2.1 <- phylolm(div_rate_eps0~CHELSA_bio10_12, data=combined_table_cor2, pruned_tree, model="BM")
summary(model_2.1)
model_2.2 <- phylolm(div_rate_eps0~CHELSA_bio10_01, data=combined_table_cor2, pruned_tree, model="BM")
summary(model_2.2)
model_2.3 <- phylolm(div_rate_eps0~meanpH, data=combined_table_cor2, pruned_tree, model="BM")
summary(model_2.3)
model_2.4 <- phylolm(div_rate_eps0~meanwatercap, data=combined_table_cor2, pruned_tree, model="BM")
summary(model_2.4)
model_2.5 <- phylolm(div_rate_eps0~meancarbon, data=combined_table_cor2, pruned_tree, model="BM")
summary(model_2.5)

model_2.6 <- phylolm(div_rate_eps0~meancarbon+meanpH, data=combined_table_cor2, pruned_tree, model="BM")
summary(model_2.6)
model_2.7 <- phylolm(div_rate_eps0~CHELSA_bio10_01+CHELSA_bio10_12, data=combined_table_cor2, pruned_tree, model="BM")
summary(model_2.7)

model <- phylolm(div_rate_eps0~mean, data=combined_table_no0, lythraceae_tree)
summary(model)

plot(combined_table_no0$div_rate_eps0~combined_table_no0$mean)
abline(model, col="blue")

#----------------------------------------------
# (3) Model: div rates ~ traits
combined_table_cor3 = combined_table
combined_table_cor3 <- subset(combined_table_cor3, combined_table_cor3$div_rate_eps0!=0) # many genera have div rates of 0 because they are monotypic. let's take remove them for one first analysis
combined_table_cor3 <- subset(combined_table_cor3, !is.na(combined_table_cor3$fm_scoring_fruit)) 

combined_table_cor3$div_rate_eps0 <- log(combined_table_cor3$div_rate_eps0)
combined_table_cor3$CHELSA_bio10_12 <- log(combined_table_cor3$CHELSA_bio10_12)

pruned_tree <- keep.tip(tree, combined_table_cor3$genus)

#combined_table_cor3$fm_scoring_fruit_binary <- NA
#combined_table_cor3$fm_scoring_fruit_binary[which(combined_table_cor3$fm_scoring_fruit=="Fleshy")] <- 0
#combined_table_cor3$fm_scoring_fruit_binary[which(combined_table_cor3$fm_scoring_fruit=="Dry")] <- 1

# Now let's run the phylolm, a phylogenetic linear model
#model_3.1 <- phyloglm(fm_scoring_fruit_binary~div_rate_eps0, data=combined_table_cor3, pruned_tree)

rate <- combined_table_cor3$div_rate_eps0
env <- combined_table_cor3$CHELSA_bio10_17
trait1 <- combined_table_cor3$fm_scoring_fruit
trait2 <- combined_table_cor3$fm_scoring_pollination2


names(rate) <- names(env) <- names(trait) <- combined_table_cor3$genus

model_3.1 <- phylolm(div_rate_eps0~CHELSA_bio10_17+fm_scoring_fruit, data=combined_table_cor3, pruned_tree)

#phylANOVA(pruned_tree, trait, env)

summary(model_3.1)

#----------------------------------------------
# (4) Model: niche breadth ~ traits


#----------------------------------------------
# (5) Model: env_factors ~ traits

# precipitation and fruit type


#----------------------------------------------
# Correlation (2): div rates ~ niche breadth
# One thing that I just thought is that we could see also if the niche breath is correlated with diversification rates, so to see if the lineages that diversified more are those that were more flexible or those who were merely exploring a particular niche with a larger area.
model <- phylolm(div_rate_eps0~sd, data=combined_table, lythraceae_tree)
summary(model)

plot(combined_table$sd~combined_table$div_rate_eps0)
abline(model, col="blue")


# Correlation (4): fruit type and other things
model <- phylolm(div_rate_eps0~tv_scoring_fruit, data=combined_table, lythraceae_tree) # uncorrelated
summary(model)
model <- phylolm(tv_scoring_corolla_diam2~tv_scoring_fruit, data=combined_table, lythraceae_tree) # correlated
summary(model)
# boxplot(lythraceae_traits$tv_scoring_corolla_diam2 ~  lythraceae_traits$tv_scoring_fruit)
model <- phylolm(mean~tv_scoring_fruit, data=combined_table, lythraceae_tree) # uncorrelated
summary(model)
model <- phylolm(tv_scoring_seed_number2~tv_scoring_fruit, data=combined_table, lythraceae_tree) # uncorrelated
summary(model)
model <- phylolm(tv_scoring_seed_size2~tv_scoring_fruit, data=combined_table, lythraceae_tree) # uncorrelated
summary(model)

# Correlation (4): pollination syndrome and other things
model <- phylolm(div_rate_eps0~tv_scoring_pollination1, data=combined_table, lythraceae_tree) # uncorrelated
summary(model)
model <- phylolm(tv_scoring_corolla_diam2~tv_scoring_pollination1, data=combined_table, lythraceae_tree) # uncorrelated
summary(model)
model <- phylolm(mean~tv_scoring_fruit, data=tv_scoring_pollination1, lythraceae_tree) # uncorrelated
summary(model)
model <- phylolm(tv_scoring_seed_number2~tv_scoring_pollination1, data=combined_table, lythraceae_tree) # uncorrelated
summary(model)
model <- phylolm(tv_scoring_seed_size2~tv_scoring_pollination1, data=combined_table, lythraceae_tree) # uncorrelated
summary(model)

#---------------------------
# Table with:
# Models        r2  Slope ln likelihood AIC p-value
# Terrestrial 0.0268 + 154.1457 −300.2913 0.1728

# Plots with all:
