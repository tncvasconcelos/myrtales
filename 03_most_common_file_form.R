#-------------------------------------------------
# Set wd as the repo
setwd("~/Desktop/Pubs_inprep/WCVP_special_issue/Eve_MyrtalesPAFTOL/myrtales/")

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
scoring <- read.csv("datasets/lifeform_mapping.csv")
life_forms <- data.frame(family=NA,genera=all_genera, most_common_life_form=NA)
for(i in 1:length(all_genera)) {
   one_genus <- all_genera[i]
   one_subset <-  subset(all_vars_myrtales, all_vars_myrtales$genus==one_genus)
   one_family <- unique(one_subset$family)
   one_subset <-  subset(one_subset, one_subset$lifeform_description!="")
   if(nrow(one_subset)>0) {
      all_life_forms <- one_subset$lifeform_description
      lumped_life_forms <- c()
      for(life_form_index in 1:length(all_life_forms)){
         lumped_life_forms[life_form_index] <- scoring$humphreys_lifeform[scoring$lifeform_description==all_life_forms[life_form_index]]
      }
      life_forms[i,3] <- names(sort(table(lumped_life_forms), decreasing=T))[1]
   }
   life_forms[i,1] <- one_family
   cat(i, "\r")
}

# cross-checking against genera sampled in the tree
tree_genera <- read.csv("datasets/prevalent_life_form.csv")
tree_genera$keep <- "YES"
life_forms <- merge(life_forms, tree_genera, by.x="genera", by.y="Genus", all=T)
life_forms <- subset(life_forms, !is.na(life_forms$keep))
write.csv(life_forms, file="life_forms_myrtales.csv", row.names = F)




