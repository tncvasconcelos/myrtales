# rm(list=ls())
     
#################################################################################################
#---------------------------
library(ape)
library(phytools)

# Traits
setwd("~/Desktop/Pubs_inprep/WCVP_special_issue/Eve_MyrtalesPAFTOL/myrtales")

master_table <- read.csv("Supplementary Table 2/2022-25-04_Myrtales-master-tabl-Table 1.csv") 
master_table <- subset(master_table, master_table$Family_Genus!="")

# here we keep only columns that will be used in the analyses
cols_to_keep <- c("Family_Genus","Family","Genus",
                  "Most.Common.Life.Form","Dry.or.Fleshy.Fruit","Mean.Seed.Number.per.Fruit",
                  "Mean.Seed.Length","Mean.Corolla.Diameter")

master_table <- master_table[,cols_to_keep]
colnames(master_table)[which(colnames(master_table)=="Genus")] <- "genus"

# fixing typos
master_table$Family_Genus[grep("Myrtaceae_Astereomyrtus", master_table$Family_Genus)] <- "Myrtaceae_Asteromyrtus"
master_table$Family_Genus[grep("Melastomataceae_Potheranthera", master_table$Family_Genus)] <- "Melastomataceae_Poteranthera"
master_table$genus[grep("Potheranthera", master_table$genus)] <- "Poteranthera"

# Niche
summary_niches <- read.csv("datasets/myrtales_niche_summaries.csv")

# Habitat
habitat <- read.csv("datasets/myrtales_habitat.csv")
habitat <- habitat[,-1] #deleting family column
colnames(habitat)[which(colnames(habitat)=="genera")] <- "genus"

# Div rates
div_rates <- read.csv("datasets/myrtales_div_table_full.csv")
div_rates <- div_rates[,-4] #deleting family column

# MERGING TABLES
combined_table <- merge(master_table, summary_niches, by="genus")
combined_table <- merge(combined_table, habitat, by="genus")
combined_table <- merge(combined_table, div_rates, by="genus")

row.names(combined_table) <- combined_table$Family_Genus
#combined_table <- combined_table[,-2]

# # Adding mean seed length from min and max values
# combined_table$seed.length.mean <- NA
# # Getting a mean value for seed length
# for(i in 1:nrow(combined_table)) {
#   seed_min <- combined_table$seed.length.min[i]
#   seed_max <- combined_table$seed.length.max[i]
#   combined_table$seed.length.mean[i] <- mean(c(seed_min, seed_max))
# }

# Lumping life forms
combined_table$Most.Common.Life.Form[which(combined_table$Most.Common.Life.Form=="annual")] <- "herbaceous"
combined_table$Most.Common.Life.Form[which(combined_table$Most.Common.Life.Form=="woody perennial")] <- "woody"
combined_table$Most.Common.Life.Form[which(combined_table$Most.Common.Life.Form=="herbaceous perennial")] <- "herbaceous"
combined_table$Most.Common.Life.Form[which(combined_table$Most.Common.Life.Form=="epiphyte")] <- "herbaceous"

# Removing 0s
# combined_table <- subset(combined_table, combined_table$div_rate_eps0.9!=0)

# Adding niche through time
combined_table$niche_through_time <- combined_table$Vol / combined_table$age

# Transforming temperature in Kelvin before logging
combined_table$CHELSA_bio10_01 <- combined_table$CHELSA_bio10_01+273.15
combined_table$CHELSA_bio10_10 <- combined_table$CHELSA_bio10_10+273.15
combined_table$CHELSA_bio10_11 <- combined_table$CHELSA_bio10_11+273.15

# Logging all continuous vars
#combined_table$div_rate_eps0.9 <- log(combined_table$div_rate_eps0.9)
combined_table$Vol <- log(combined_table$Vol)
combined_table$niche_through_time <- log(combined_table$niche_through_time)

# Abiotic:
combined_table$CHELSA_bio10_01 <- log(combined_table$CHELSA_bio10_01)
combined_table$CHELSA_bio10_02 <- log(combined_table$CHELSA_bio10_02)
combined_table$CHELSA_bio10_10 <- log(combined_table$CHELSA_bio10_10)
combined_table$CHELSA_bio10_11 <- log(combined_table$CHELSA_bio10_11)
combined_table$CHELSA_bio10_12 <- log(combined_table$CHELSA_bio10_12)
combined_table$CHELSA_bio10_15 <- log(combined_table$CHELSA_bio10_15)
combined_table$CHELSA_bio10_16 <- log(combined_table$CHELSA_bio10_16)
combined_table$CHELSA_bio10_17 <- log(combined_table$CHELSA_bio10_17)
combined_table$GLOBAL_SLOPE_10MIN <- log(combined_table$GLOBAL_SLOPE_10MIN)
combined_table$depthtobedrock2 <- log(combined_table$depthtobedrock2)
combined_table$meancarbon <- log(combined_table$meancarbon)
combined_table$meanpH <- log(combined_table$meanpH)
combined_table$meanwatercap <- log(combined_table$meanwatercap)

# traits:

combined_table$Mean.Seed.Number.per.Fruit <- log(combined_table$Mean.Seed.Number.per.Fruit)
combined_table$Mean.Seed.Length <- log(combined_table$Mean.Seed.Length)
combined_table$Mean.Corolla.Diameter <- log(combined_table$Mean.Corolla.Diameter)

# SAVE
saveRDS(combined_table, file="datasets/Myrtales_full_dataset.Rdata")
# for sup mat:
write.csv(combined_table, file="Myrtales_full_dataset.csv")

# clade specifics:
myrtaceae_table <- subset(combined_table, combined_table$Family=="Myrtaceae")
onagraceae_table <- subset(combined_table, combined_table$Family=="Onagraceae")
combretaceae_table <- subset(combined_table, combined_table$Family=="Combretaceae")
melascapclade_table <- subset(combined_table, combined_table$Family%in%c("Melastomataceae","Alzateaceae","Crypteroniaceae","Penaeaceae"))
vochysiaceae_table <- subset(combined_table, combined_table$Family=="Vochysiaceae")
lythraceae_table <- subset(combined_table, combined_table$Family=="Lythraceae")

clade_specific_list <- list(myrtaceae_table, onagraceae_table, combretaceae_table, melascapclade_table, vochysiaceae_table, lythraceae_table)
names(clade_specific_list) <- c("Myrtaceae", "Onagraceae", "Combretaceae", "Melastomataceae_CAPclade","Vochysiaceae","Lythraceae")

saveRDS(clade_specific_list, file="datasets/Myrtales_by_clade_dataset.Rdata")

