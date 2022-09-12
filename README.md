Codes used in: The confluence of traits and environmental factors driving diversification and niche expansion in the globally distributed order Myrtales (in prep.)

----
Description of folders: 
 
- **datasets/** 

Tables with traits to be analyzed. 

- **plots/** 

Folder to store plots. 

- **results/** 

Results organized in folders that correspond to the three hypotheses presented in the main manuscript (h1, h2, h3)

- **tree/** 

Folder to store trees.  

----
Scripts:

> 00_pruning_tree.R

Organizes tip labels in the Myrtales time-calibrated tree so that only one tip per genus is kept and tip labels match with those in trait dataset.

> 01_cleaning_gbif.R 

Filters occurence points downloaded from GBIF following a series of criteria and the distribution information available from the WCVP.

> 02_div_rates.R

Uses the stem age estimates and the number of species for each Myrtales genera to calculate net-diversification rates using the method-of-moments estimator from Magallon & Sanderson (2001)[(link)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.0014-3820.2001.tb00826.x).

> 03_most_common_life_form.R 

Calculates most common life form in each Myrtales genus based on the WCVP dataset.

> 04_Niche.R 

Script used to calculate niche volume for each genus (by Samuel Pironon)

> 05_get_habitat.R 

Uses shapefiles from WWF to calculate the proportion of dustribution points in each Myrtales genus that falls within closed or open canopy biomes.

> 06_organizing_tables_for_regression.R 

Organizes traits, diversification rates, niche, habitat, and life form data in a single table for subsequent analyses.

> 07.1_h1_individual_regressions.R 

Runs individual phylogenetic regression analyses (hypothesis 1)

> 07.2_h2_dredging_models.R 



> 07.3_phylanova_habitat.R 

Script used to run phyloanovas and corhmm (hypothesis 3).

> 07.4_testing_colinearity.R 

Short script to test collinearity between environmental variables used in the multiple regressions analysis (hypothesis 2).

> 08_tree_plots.R 

Script to make tree plots for Figure 1 in the main paper (Myrtales tree and trait data).


----
Distribution data used in niche analyses comes from:

GBIF.org (19 March 2022) GBIF Occurrence Download https://doi.org/10.15468/dl.ba76td
