# rm(list=ls())
# setwd("~/Desktop/WCVP_special_issue/Eve_MyrtalesPAFTOL/myrtales")
#################################################################################################
#--------------------------------------------
#--------------------------------------------
# (1) Dated tree
library(ape)
whole_myrtales_tree <- read.tree("tree/myrtales_pruned.tre")
lyt_tree <- keep.tip(whole_myrtales_tree, grep("Lythraceae", whole_myrtales_tree$tip.label))
write.tree(lyt_tree, file="pilot/pilot_tree.tre")

#################################################################################################
#--------------------------------------------
#--------------------------------------------
# (2) Distribution points
# Using the ones I cleaned for Myrtales

#################################################################################################
#--------------------------------------------
#--------------------------------------------
# (3) Trait data
# Table with most common life form:
dist_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")
all_vars <- subset(all_vars, all_vars$family %in% "Lythraceae") # subsetting it for the pilot
all_vars$species_name_full <- paste(all_vars$taxon_name, all_vars$taxon_authors, sep=" ")
# excluding genera only
all_vars <- subset(all_vars, grepl(" ", all_vars$taxon_name))

# Getting most common life form for each genus
all_genera <- unique(all_vars$genus)
life_forms <- data.frame(genera=all_genera, most_common_life_form=NA)
for(i in 1:length(all_genera)) {
  one_genus <- all_genera[i]
  one_subset <-  subset(all_vars, all_vars$genus==one_genus)
  one_subset <-  subset(one_subset, one_subset$lifeform_description!="")
  if(nrow(one_subset)>0) {
    life_forms[i,2] <- names(sort(table(one_subset$lifeform_description), decreasing=T))[1]
  }
  cat(i, "\r")
}
write.csv(life_forms, file="pilot/pilot_most_common_life_form.csv",row.names = F)
View(life_forms)
# proportion NA 0.06896552
# length(which(is.na(life_forms$most_common_life_form))) /length(all_genera)

#################################################################################################
#--------------------------------------------
#--------------------------------------------
# (4) WCVP and TWDG data
library(data.table)
library(maptools)
data("wrld_simpl")

# Load distribution data
source("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/GetRanges.R")
source("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/miscellaneous.R")

lythraceae_points <- fread("myrtales_gbif/Lythraceae_cleaned_points.csv")
ranges_lythraceae <- GetRanges(lythraceae_points, species="genus", lat="decimalLatitude", lon="decimalLongitude", threshold=0.75, buffer=25, res=10)
save(ranges_lythraceae, file="Pilot/ranges_lythraceae.Rsave")
range_sizes <- as.data.frame(unlist(lapply(ranges_lythraceae, "[[","range_size")))
colnames(range_sizes) <- c("range_size")
range_sizes$genus <- row.names(range_sizes)
range_sizes <- range_sizes[,c(2,1)]
range_sizes <- subset(range_sizes, range_sizes$genus!="")
write.csv(range_sizes, file="Pilot/lythraceae_range_sizes.csv", row.names=F)

# let's see how they look
plot(wrld_simpl)
plot(ranges_lythraceae$Decodon$range, add=T, col="blue", legend=F)
plot(ranges_lythraceae$Lafoensia$range, add=T, col="blue", legend=F)

# Now let's get some data from these ranges

all_points <- as.data.frame(matrix(nrow=0,ncol=3))
colnames(all_points) <- c("genus","lon","lat")
for(i in 1:length(ranges_lythraceae)){
  one_result <- ranges_lythraceae[[i]]
  points_tmp <- rasterToPoints(one_result$range)[,c(1,2)]
  points_tmp <- as.data.frame(cbind(rep(one_result$species_name, nrow(points_tmp)), points_tmp))
  colnames(points_tmp) <- colnames(all_points)
  all_points <- rbind(all_points, points_tmp)
}

layer <- raster("~/Desktop/rangers/example_datasets/bio12.bil")
all_biomes <- localityToBiome(all_points) 
all_biomes <- subset(all_biomes, !is.na(all_biomes$biome))

all_points$lon <- as.numeric(all_points$lon)
all_points$lat <- as.numeric(all_points$lat)

summary_biomes <- getBiomes(all_biomes, species="genus")
climate_dataset <- GetClimateSummStats(all_points, layer, species="genus")
colnames(climate_dataset)[1] <- "genus"
write.csv(climate_dataset, file="Pilot/lythraceae_niche.csv", row.names=F)

#################################################################################################
#--------------------------------------------
#--------------------------------------------
# (5) Div rates
library(geiger)
library(ape)
library(phangorn)

get.node.age <- function (phy) {
  root.node <- length(phy$tip.label)+1
  seq.nodes <- phy$edge
  dists <- phy$edge.length
  res <- numeric(max(phy$edge))
  for (i in seq_len(nrow(seq.nodes))) {
    res[seq.nodes[i, 2]] <- res[seq.nodes[i,1]] + dists[i]
  }
  ages <- abs(round(res,3)-round(max(res),3))
  return(ages)
} # fun?ao pra pegar os node ages

# Getting a table with genera age from dated phylogeny
myrtales_tree <- read.tree("tree/myrtales_pruned.tre")
all_node_ages <- get.node.age(myrtales_tree)
myrtales_genera_ages <- data.frame(genus=myrtales_tree$tip.label, age=NA)
for(i in 1:length(myrtales_tree$tip.label)){
  myrtales_genera_ages[i,2] <- all_node_ages[Ancestors(myrtales_tree, i, type="parent")]
}
myrtales_genera_ages$genus <- gsub("^[^_]*_","", myrtales_genera_ages$genus)

# Now getting estimates of species number per genus from WCVP
names_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

# Merge them in one big table
myrtales_families <- c("Lythraceae")

all_myrtales_species <- subset(names_sample, names_sample$family %in% myrtales_families)
all_myrtales_species$species_name_full <- paste(all_myrtales_species$taxon_name, all_myrtales_species$taxon_authors, sep=" ")
# excluding genera only
all_myrtales_species <- subset(all_myrtales_species, grepl(" ", all_myrtales_species$taxon_name))
all_myrtales_species <- subset(all_myrtales_species, !duplicated(all_myrtales_species$species_name_full))
all_myrtales_species <- subset(all_myrtales_species, all_myrtales_species$taxon_status=="Accepted") #only accepted

myrtales_species_total <- as.data.frame(table(all_myrtales_species$genus))
#myrtales_species_total <- myrtales_species_total[-1,]
colnames(myrtales_species_total) <- c("genus","n_species")

# Merge
myrtales_final <- merge(myrtales_genera_ages, myrtales_species_total, by="genus", all = T)
write.csv(myrtales_final, file="myrtales_div_table_full.csv", row.names = F)

# Preliminary
myrtales_final <- merge(myrtales_genera_ages, myrtales_species_total, by="genus")

myrtales_final$div_rate_eps0 <- NA
myrtales_final$div_rate_eps0.5 <- NA
myrtales_final$div_rate_eps0.9 <- NA

for(u in 1:nrow(myrtales_final)) {
  myrtales_final$div_rate_eps0[u] <- round(bd.ms(phy=NULL, myrtales_final$age[u], myrtales_final$n_species[u], crown=F, epsilon=0),3)
  myrtales_final$div_rate_eps0.5[u] <- round(bd.ms(phy=NULL, myrtales_final$age[u], myrtales_final$n_species[u], crown=F, epsilon=0.5),3)
  myrtales_final$div_rate_eps0.9[u] <- round(bd.ms(phy=NULL, myrtales_final$age[u], myrtales_final$n_species[u], crown=F, epsilon=0.9),3)
}
write.csv(myrtales_final, file="Pilot/lythraceae_div_table_preliminary.csv", row.names = F)

#################################################################################################
#---------------------------
#---------------------------
# (6) PGLS and other correlation analyses
# Let's load everything first:
library(phylolm)

# Tree
lythraceae_tree <- read.tree("Pilot/pilot_tree.tre")
lythraceae_tree$tip.label <- gsub("^[^_]*_","", lythraceae_tree$tip.label)
# plot(lythraceae_tree)
# axisPhylo()

# Traits
lythraceae_traits <- read.csv("Pilot/pilot_traits_scored.csv") 
lythraceae_traits <- lythraceae_traits[,c(1,2, grep("tv",colnames(lythraceae_traits)))]
colnames(lythraceae_traits)[2] <- "genus" # importante que sempre esteja como "genus"

# Range size
lythraceae_range <- read.csv("Pilot/lythraceae_range_sizes.csv")

# Niche
lythraceae_niche <- read.csv("Pilot/lythraceae_niche.csv")

# Div rates
lythraceae_divrates <- read.csv("Pilot/lythraceae_div_table_preliminary.csv")

# making mega table
combined_table <- merge(lythraceae_niche, lythraceae_divrates)
combined_table <- merge(combined_table, lythraceae_traits)
combined_table <- merge(combined_table, lythraceae_range)
row.names(combined_table) <- combined_table$genus

View(combined_table)

##### copied and pasted from abstract:
# In this study we will address four hypotheses: 
# (1) key traits and environment are correlated and determine range size, geographic distribution, and speciation/extinction rates of extant taxa at a global scale. 
# (2) species with fleshy fruits and/or species with fruits with larger, fewer seeds are more species-rich in closed habitats. 
# (3) wind pollinated species and/or those with larger, more conspicuous flowers are more species-rich in open habitats. 
# (4) range expansion rates are positively correlated with speciation rates. 

#----------------------------------------------
# Correlation (1): div rates ~ range expansion
combined_table$range_expansion <- combined_table$range_size / combined_table$age # Following that paper, range expansion is range size divided by age of the group
combined_table_no0 <- subset(combined_table, combined_table$div_rate_eps0!=0) # many genera have div rates of 0 because they are monotypic. let's take remove them for one first analysis
#combined_table_no0 <- combined_table

# Now let's run the phylolm, a phylogenetic linear model
model <- phylolm(div_rate_eps0~range_expansion, data=combined_table_no0, lythraceae_tree)
summary(model)
# Coefficients:
#  Estimate  StdErr t.value   p.value    
# (Intercept)     -41668   51570  -0.808 0.4336426    
# div_rate_eps0  1332920  282221   4.723 0.0003984 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-squared: 0.6318	Adjusted R-squared: 0.6035 

plot(combined_table_no0$range_expansion~combined_table_no0$div_rate_eps0)
abline(model, col="blue")

#----------------------------------------------
# Correlation (2): div rates ~ niche breadth
# One thing that I just thought is that we could see also if the niche breath is correlated with diversification rates, so to see if the lineages that diversified more are those that were more flexible or those who were merely exploring a particular niche with a larger area.
model <- phylolm(div_rate_eps0~sd, data=combined_table, lythraceae_tree)
summary(model)

plot(combined_table$sd~combined_table$div_rate_eps0)
abline(model, col="blue")

#----------------------------------------------
# Correlation (3): div rates ~ mean prec
model <- phylolm(div_rate_eps0~mean, data=combined_table_no0, lythraceae_tree)
summary(model)

plot(combined_table_no0$div_rate_eps0~combined_table_no0$mean)
abline(model, col="blue")

#----------------------------------------------
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
