# Getting data for pilot
# We will need:
#--------------------------------------------
#--------------------------------------------
# (1) Dated tree
library(ape)
whole_myrtales_tree <- read.tree("tree/myrtales_pruned.tre")
lyt_tree <- keep.tip(whole_myrtales_tree, grep("Lythraceae", whole_myrtales_tree$tip.label))
write.tree(lyt_tree, file="pilot/pilot_tree.tre")

#--------------------------------------------
#--------------------------------------------
# (2) Distribution points
source("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/WCVPtools_functions.R")
library(data.table)
library(maptools)
library(raster)

dist_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")
all_vars <- subset(all_vars, all_vars$family %in% "Lythraceae") # subsetting it for the pilot

# reference table for taxized names
reference_table <- list.files("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/taxized_reference_tables", full.names = T)
reference_table <- do.call(rbind, lapply(reference_table, read.csv))

gbif_data <- fread("pilot/0191898-210914110416597.csv") # load the table you downloaded from GBIF
reference_table <- subset(reference_table, reference_table$wcvp_name %in% unique(all_vars$taxon_name))

# Looking at the WCVP table and TDWG to clean GBIF points
cleaned_gbif <- FilterWCVP(points=gbif_data, all_vars, reference_table, path="/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/wgsrpd-master/level3/level3.shp") # This will filter the GBIF points acording to WCVP
write.csv(cleaned_gbif, file="pilot_cleaned_points.csv", row.names=F)
# "36126 points removed."


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

# proportion NA 0.06896552
# length(which(is.na(life_forms$most_common_life_form))) /length(all_genera)

#--------------------------------------------
#--------------------------------------------
# (4) WCVP and TWDG data

# Load distribution data
source("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/GetRanges.R")
lythraceae_points <- fread("myrtales_gbif/Lythraceae_cleaned_points.csv")
ranges_lythraceae <- GetRanges(lythraceae_points, species="genus", lat="decimalLatitude", lon="decimalLongitude", threshold=0.75, buffer=25, res=10)
save(ranges_lythraceae, file="Pilot/ranges_lythraceae.Rsave")
range_sizes <- as.data.frame(unlist(lapply(ranges_lythraceae, "[[","range_size")))
colnames(range_sizes) <- c("range_size")
range_sizes$genus <- row.names(range_sizes)
range_sizes <- range_sizes[,c(2,1)]
range_sizes <- subset(range_sizes, range_sizes$genus!="")
write.csv(range_sizes, file="Pilot/lythraceae_range_sizes.csv", row.names=F)

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

#---------------------------
#---------------------------
# (6) PGLS and other correlation analyses
library(phylolm)
tree <- read.tree("Pilot/pilot_tree.tre")
tree$tip.label <- gsub("^[^_]*_","", tree$tip.label)

#---------------------------
# load traits
lythraceae_traits <- read.csv("Pilot/pilot_traits_scored.csv") 

#---------------------------
# Correlation (1): div rates ~ range expansion--------
div_rates <- read.csv("Pilot/lythraceae_div_table_preliminary.csv")
range_sizes <- read.csv("Pilot/lythraceae_range_sizes.csv")

comb_table <- merge(range_sizes, div_rates)
comb_table$range_expansion <- comb_table$range_size / comb_table$age
comb_table_no0 <- subset(comb_table, comb_table$div_rate_eps0!=0)

row.names(comb_table_no0) <- comb_table_no0$genus
model <- phylolm(range_expansion~div_rate_eps0, data=comb_table_no0, tree)

# Coefficients:
# Estimate StdErr t.value p.value  
# (Intercept)      28633  42123  0.6797 0.50859  
# div_rate_eps0   673766 230520  2.9228 0.01188 *
# R-squared: 0.3966	Adjusted R-squared: 0.3501 

plot(comb_table_no0$range_expansion~comb_table_no0$div_rate_eps0)
abline(model, col="blue")
#---------------------------

#---------------------------
# Table with:
# Models        r2  Slope ln likelihood AIC p-value
# Terrestrial 0.0268 + 154.1457 −300.2913 0.1728

# Plots with all:
