# cleaning Myrtales points
# rm(list=ls())

library(data.table)
library(maptools)
library(raster)
library(sp)
library(rgeos)
library(rworldmap)
data("wrld_simpl")

#-----------------------------
# If local
# setwd("~/Desktop/WCVP_special_issue/Eve_MyrtalesPAFTOL/myrtales")
source("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/WCVPtools_functions.R")
dist_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
#-----------------------------
# If labcomputer
# setwd("~/myrtales")
source("../WCVPtools/WCVPtools_functions.R")
dist_sample <- read.table("../life_history_houwie/wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("../life_history_houwie/wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
#-----------------------------
# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")
myrtales_families <- c("Alzateaceae",
                       "Combretaceae",
                       "Crypteroniaceae",
                       "Lythraceae",
                       "Melastomataceae",
                       "Myrtaceae",
                       "Onagraceae",
                       "Penaeaceae",
                       "Vochysiaceae")
all_vars <- subset(all_vars, all_vars$family %in% myrtales_families)

# reference table for taxized names
#-----------------------------
# If local
reference_table <- list.files("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/taxized_reference_tables", full.names = T)
reference_table <- do.call(rbind, lapply(reference_table, read.csv))
#-----------------------------
# If labcomputer
reference_table <- list.files("../WCVPtools/taxized_reference_tables", full.names = T)
reference_table <- do.call(rbind, lapply(reference_table, read.csv))

# Reading gbif file
gbif_data <- fread("myrtales_gbif/0188513-210914110416597.csv") # load the table you downloaded from GBIF

# Looking at the WCVP table and TDWG to clean GBIF points
#-----------------------------
# If local
path="/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/wgsrpd-master/level3/level3.shp"
#-----------------------------
# If labcomputer
path="../WCVPtools/wgsrpd-master/level3/level3.shp"
#-----------------------------

twgd_data <- suppressWarnings(maptools::readShapeSpatial(path))
issues_to_remove <- read.csv("gbif_issues_to_remove.csv")
cultivated <- read.csv("cultivated_species.csv")
for(family_index in 1:length(myrtales_families)){
  one_subset <- subset(gbif_data, gbif_data$family == myrtales_families[family_index]) #subsetting to make it more manageble
  cleaned_points <- subset(one_subset, one_subset$basisOfRecord == "PRESERVED_SPECIMEN")
  cleaned_points <- subset(cleaned_points, cleaned_points$scientificName!="")
  cleaned_points <- FilterWCVP_genus(cleaned_points, all_vars, twgd_data)
  for(issue_index in 1:nrow(issues_to_remove)) {
    cleaned_points <- subset(cleaned_points, !grepl(issues_to_remove$issues_to_remove[issue_index], cleaned_points$issue))
  }
  for(cultivated_index in 1:nrow(cultivated)) {
    cleaned_points <- subset(cleaned_points, !grepl(cultivated$cultivated_species[cultivated_index], cleaned_points$species))
  }
  subset_reference_table <- subset(reference_table, reference_table$gbif_name %in% unique(cleaned_points$scientificName))
  if(nrow(subset_reference_table)>0){
    cleaned_points <- FilterWCVP(cleaned_points, all_vars, subset_reference_table, twgd_data) # This will filter the GBIF points acording to WCVP for species
  }
  # Cleaning common problems:
  cleaned_points <- RemoveNoDecimal(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  cleaned_points <- RemoveCentroids(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  cleaned_points <- RemoveDuplicates(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  cleaned_points <- RemoveOutliers(cleaned_points, species="scientificName", lon="decimalLongitude", lat="decimalLatitude")
  cleaned_points <- RemoveZeros(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  cleaned_points <- RemoveSeaPoints(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  write.csv(cleaned_points, file=paste0("myrtales_gbif/", myrtales_families[family_index], "_cleaned_points.csv"), row.names=F)
}

#------------------------
all_cleaned_points_files <- list.files("myrtales_gbif", full.names = T)
all_cleaned_points_files <- all_cleaned_points_files[grep("cleaned", all_cleaned_points_files)]
labels <- gsub(paste0(c("myrtales_gbif/","_cleaned_points.csv"), collapse="|"),"", all_cleaned_points_files)
all_cleaned_points <- lapply(all_cleaned_points_files, read.csv)
names(all_cleaned_points) <- labels

# Plotting to inspect distributions
{; for(family_index in 1:length(all_cleaned_points)) {
  points_cleaned <- all_cleaned_points[[family_index]]
  genera <- unique(points_cleaned$genus)
  genera <- subset(genera, genera!="")
  pdf(paste0("plots/", names(all_cleaned_points)[family_index], "_points.pdf"))
  for(genus_index in 1:length(genera)){
    tmp_subset <- as.data.frame(points_cleaned[points_cleaned$genus==genera[genus_index],])
    coord <- tmp_subset[,c("decimalLongitude","decimalLatitude")]
    coordinates(coord) <- ~ decimalLongitude + decimalLatitude
    plot(wrld_simpl)
    plot(coord, col="red", add=T)
    title(genera[genus_index])
    print(genus_index)
  }
  dev.off()
}
  beepr::beep("fanfare"); } 


