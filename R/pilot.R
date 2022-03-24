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


# (5) Div rates


# (6) PGLS and other correlation analyses



