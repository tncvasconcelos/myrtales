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

# Now getting estimates of species number per genus from WCVP
names_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

# Merge them in one big table
myrtales_families <- c("Alzateaceae","Combretaceae","Crypteroniaceae","Lythraceae",
                       "Melastomataceae","Myrtaceae","Onagraceae","Penaeaceae","Vochysiaceae")

all_myrtales_species <- subset(names_sample, names_sample$family %in% myrtales_families)
all_myrtales_species$species_name_full <- paste(all_myrtales_species$taxon_name, all_myrtales_species$taxon_authors, sep=" ")
# excluding genera only
all_myrtales_species <- subset(all_myrtales_species, grepl(" ", all_myrtales_species$taxon_name))
all_myrtales_species <- subset(all_myrtales_species, !duplicated(all_myrtales_species$species_name_full))
all_myrtales_species <- subset(all_myrtales_species, all_myrtales_species$taxon_status=="Accepted") #only accepted

myrtales_species_total <- as.data.frame(table(all_myrtales_species$genus))
myrtales_species_total <- myrtales_species_total[-1,]
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
write.csv(myrtales_final, file="myrtales_div_table_preliminary.csv", row.names = F)



