# cleaning Myrtales points
# rm(list=ls())

library(data.table)
library(maptools)
library(raster)

# Local
# setwd("~/Desktop/WCVP_special_issue/Eve_MyrtalesPAFTOL/myrtales")
source("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/WCVPtools_functions.R")
dist_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

# Labcomputer
# setwd("~/myrtales")
source("../WCVPtools/WCVPtools_functions.R")
dist_sample <- read.table("../life_history_houwie/wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("../life_history_houwie/wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")
myrtales_families <- c("Alzateaceae","Combretaceae","Crypteroniaceae","Lythraceae",
                       "Melastomataceae","Myrtaceae","Onagraceae","Penaeaceae","Vochysiaceae")
all_vars <- subset(all_vars, all_vars$family %in% c(myrtales_families))

# reference table for taxized names
# Local
reference_table <- list.files("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/taxized_reference_tables", full.names = T)
reference_table <- do.call(rbind, lapply(reference_table, read.csv))
# Labcomputer
reference_table <- list.files("../WCVPtools/taxized_reference_tables", full.names = T)
reference_table <- do.call(rbind, lapply(reference_table, read.csv))

# Reading gbif file
gbif_data <- fread("myrtales_gbif/0188513-210914110416597.csv") # load the table you downloaded from GBIF

reference_table <- subset(reference_table, reference_table$wcvp_name %in% unique(all_vars$taxon_name))
# Looking at the WCVP table and TDWG to clean GBIF points
# Local
path="/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/wgsrpd-master/level3/level3.shp"
# Labcomputer
path="../WCVPtools/wgsrpd-master/level3/level3.shp"

cleaned_gbif <- FilterWCVP(gbif_data, all_vars, reference_table, path=path) # This will filter the GBIF points acording to WCVP
write.csv(cleaned_gbif, file="cleaned_gbif_wcvp.csv", row.names=F)
