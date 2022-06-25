# rm(list=ls())
setwd("~/Desktop/WCVP_special_issue/Eve_MyrtalesPAFTOL/myrtales")
#################################################################################################
#---------------------------
library(ape)
library(phytools)

# Traits
master_table <- read.csv("traits_myrtales_25.06.csv") 
cols_to_keep <- c("X","Family","Genus","most_common_life_form","seed.min","seed.max","fm_scoring_seed_number",
                  "seed.length.min","seed.length.max","seed.width.min","seed.length.max","seed.width.min",
                  "seed.width.max","fm_scoring_pollination1","fm_scoring_pollination2")

master_table <- master_table[,cols_to_keep]
colnames(master_table)[which(colnames(master_table)=="Genus")] <- "genus"

# Niche
summary_niches <- read.csv("myrtales_niche_summaries.csv")

# Habitat
habitat <- read.csv("myrtales_habitat.csv")
habitat <- habitat[,-1] #deleting family column
colnames(habitat)[which(colnames(habitat)=="genera")] <- "genus"

# Div rates
div_rates <- read.csv("myrtales_div_table_full.csv")
div_rates <- div_rates[,-4] #deleting family column

# MERGING TABLES
combined_table <- merge(master_table, summary_niches, by="genus")
combined_table <- merge(combined_table, habitat, by="genus")
combined_table <- merge(combined_table, div_rates, by="genus")

row.names(combined_table) <- combined_table$X
combined_table <- combined_table[,-2]

saveRDS(combined_table, file="Myrtales_full_dataset.Rdata")

# clade specifics:
myrtaceae_table <- subset(combined_table, combined_table$Family=="Myrtaceae")
onagraceae_table <- subset(combined_table, combined_table$Family=="Onagraceae")
combretaceae_table <- subset(combined_table, combined_table$Family=="Combretaceae")
melascapclade_table <- subset(combined_table, combined_table$Family%in%c("Melastomataceae","Alzateaceae","Crypteroniaceae","Penaeaceae"))
vochysiaceae_table <- subset(combined_table, combined_table$Family=="Vochysiaceae")
lythraceae_table <- subset(combined_table, combined_table$Family=="Lythraceae")

clade_specific_list <- list(myrtaceae_table, onagraceae_table, combretaceae_table, melascapclade_table, vochysiaceae_table, lythraceae_table)
names(clade_specific_list) <- c("Myrtaceae", "Onagraceae", "Combretaceae", "Melastomataceae_CAPclade","Vochysiaceae","Lythraceae")
