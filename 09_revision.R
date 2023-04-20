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
master_table <- subset(master_table, master_table$div_rate_eps0.9!=0) # Excluding 0 div rates from monotypic genera (cant be logged, same procudure of all analyses)
master_table$div_rate_eps0.9 <- log(master_table$div_rate_eps0.9) 

#summary(phylolm(div_rate_eps0.9~n_species, data=master_table, phy=tree))

# Question 2 
# Let's run the same as in the paper:
model_div_niche_expansion <- phylolm(div_rate_eps0.9~age, data=master_table, phy=tree)
summary(model_div_niche_expansion)
#Call:
#  phylolm(formula = div_rate_eps0.9 ~ niche_through_time, data = master_table, 
#          phy = tree)
#
# AIC logLik 
# 687.0 -340.5 
#
# Raw residuals:
#  Min      1Q  Median      3Q     Max 
# -2.8698 -0.0368  0.5748  1.2439  4.0058 
#
# Mean tip height: 128.7
# Parameter estimate(s) using ML:
#  sigma2: 0.02494881 
#
# Coefficients:
#  Estimate    StdErr t.value   p.value    
# (Intercept)        -3.730089  0.589132 -6.3315 1.228e-09 ***
#  niche_through_time  0.357776  0.026378 13.5634 < 2.2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# R-squared: 0.4401	Adjusted R-squared: 0.4378 
####

model_div_niche_expansion_nspecies <- phylolm(div_rate_eps0.9~niche_through_time+n_species, data=master_table, phy=tree)
summary(model_div_niche_expansion_nspecies)


model_niche_expansion_div_nspecies <- phylolm(niche_through_time~div_rate_eps0.9+n_species, data=master_table, phy=tree)
summary(model_niche_expansion_div_nspecies)

# Call:
#  phylolm(formula = niche_through_time ~ div_rate_eps0.9 + n_species, 
#          data = master_table, phy = tree)
#
# AIC logLik 
# 975.1 -483.6 
#
# Raw residuals:
#  Min      1Q  Median      3Q     Max 
# -9.4634 -1.5237 -0.3272  0.5074  3.4076 
#
# Mean tip height: 128.7
# Parameter estimate(s) using ML:
#  sigma2: 0.0838665 
#
# Coefficients:
#  Estimate    StdErr t.value   p.value    
# (Intercept)     3.5978977 1.1601561  3.1012  0.002165 ** 
#  div_rate_eps0.9 1.1391505 0.0981404 11.6074 < 2.2e-16 ***
#  n_species       0.0014408 0.0006237  2.3100  0.021762 *  
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# R-squared: 0.4527	Adjusted R-squared: 0.448 
#

# We can show the differences between these in a scatterplot:
#plot(master_table$n_species, master_table$div_rate_eps0.9)
#plot(master_table$n_species, exp(master_table$niche_through_time))

# Create scatterplot of y vs x
#plot(master_table$niche_through_time, master_table$div_rate_eps0.9, 
#     xlab = "log(niche-breadth expansion)", ylab = "log(div. rates)",
#     pch=19, col="darkgray", cex=0.7)

# Add linear regression line
#abline(model_div_niche_expansion, col = "blue")
#abline(model_div_niche_expansion_nspecies, col = "blue", lty=2)

###################################################
# Question 3

model_div_niche_expansion_nspecies <- phylolm(div_rate_eps0.9~niche_through_time+age, data=master_table, phy=tree)
model_div_niche_expansion_nspecies <- phylolm(niche_through_time~div_rate_eps0.9+age, data=master_table, phy=tree)

summary(model_div_niche_expansion_nspecies)

model_niche_expansion_age <- phylolm(niche_through_time~age, data=master_table, phy=tree)
model_niche_diversification_age <- phylolm(div_rate_eps0.9~age, data=master_table, phy=tree)





abline(lm(y ~ x + z), col = "red")