############################################
############################################
# rm(list=ls())
setwd("~/Desktop/myrtales")

library(ape)
library(phylolm)
############################################

# Load tree
tree <- read.tree("tree/myrtales_pruned.tre")

# Load master table
master_table <- readRDS("datasets/Myrtales_full_dataset.Rdata") 

############################################
############################################
master_table <- subset(master_table, master_table$div_rate_eps0.9!=0) # 0.9 extinction fraction

#
# The issues raised by the review are related to age and species-richness being confounding factors 
# in the correlation between niche breadth and diversification rates. As the review points out, 
# ascribing cause in these types of studies are dubious at best, but we agree that it is important 
# to consider other variables which may lead to spurious correlations. 

# We have some hesitation about including species-richness and age as a potential correlate with 
# diversification rates because these are an explicit part of the calculation in the method-of-moments 
# (Magallon and Sanderson, 2001; referenced in the manuscript), and we suspect that the positive 
# correlation between these could be a statistical artefact. 

# Nevertheless, to assess the potential for the confounding variables, we run regressions of diversification
# rates and niche breadth as predicted by species richness and clade age:

model1 <- phylolm(niche_through_time~n_species, data=master_table, phy=tree)
summary(model1)
model2 <- phylolm(niche_through_time~age, data=master_table, phy=tree)
summary(model2)
model3 <- phylolm(div_rate_eps0.9~n_species, data=master_table, phy=tree)
summary(model3)
model4 <- phylolm(div_rate_eps0.9~age, data=master_table, phy=tree)
summary(model4)

# Finally, to account for the possibility that niche-breadth expansion and diversification rates are 
# independently related to each other even when accounting for the confounding factor of species-richness 
# or clade age, we included species-richness and clade age as covariates in our regression. Running 
#this regression structure, we found only slight increases in R2. 

model5 <- phylolm(div_rate_eps0.9~niche_through_time+n_species, data=master_table, phy=tree)
summary(model5)
model6 <- phylolm(div_rate_eps0.9~niche_through_time+age, data=master_table, phy=tree)
summary(model6)
model7 <- phylolm(niche_through_time~div_rate_eps0.9+n_species, data=master_table, phy=tree)
summary(model7)
model8 <- phylolm(niche_through_time~div_rate_eps0.9+age, data=master_table, phy=tree)
summary(model8)


