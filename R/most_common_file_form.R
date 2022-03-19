#-------------------------------------------------
#-------------------------------------------------
#-------------------------------------------------
# Set wd as the repo
setwd("~/Desktop/WCVP_special_issue/Eve_MyrtalesPAFTOL/myrtales/")

#----------------
# (1) Filter Myrtales from the WCVP
dist_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")

myrtales_families <- c("Alzateaceae","Combretaceae","Crypteroniaceae","Lythraceae",
                       "Melastomataceae","Myrtaceae","Onagraceae","Penaeaceae","Vochysiaceae")

all_vars_myrtales <- subset(all_vars, all_vars$family %in% myrtales_families)
all_vars_myrtales$species_name_full <- paste(all_vars_myrtales$taxon_name, all_vars_myrtales$taxon_authors, sep=" ")
# excluding genera only
all_vars_myrtales <- subset(all_vars_myrtales, grepl(" ", all_vars_myrtales$taxon_name))

# Getting most common life form for each genus
all_genera <- unique(all_vars_myrtales$genus)
life_forms <- data.frame(genera=all_genera, most_common_life_form=NA)
for(i in 1:length(all_genera)) {
   one_genus <- all_genera[i]
   one_subset <-  subset(all_vars_myrtales, all_vars_myrtales$genus==one_genus)
   one_subset <-  subset(one_subset, one_subset$lifeform_description!="")
   if(nrow(one_subset)>0) {
        life_forms[i,2] <- names(sort(table(one_subset$lifeform_description), decreasing=T))[1]
   }
   cat(i, "\r")
}
write.csv(life_forms, file="most_common_life_form.csv",row.names = F)

# proportion NA 37.1%
length(which(is.na(life_forms$most_common_life_form))) /length(all_genera)

