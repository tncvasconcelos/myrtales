# Regression analyses H1: individual regressions
# rm(list=ls())
setwd("~/Desktop/Pubs_inprep/WCVP_special_issue/Eve_MyrtalesPAFTOL/myrtales")
#################################################################################################
library(phylolm)
library(ape)

# Tree
tree <- read.tree("tree/myrtales_pruned.tre")

# Master table
master_table <- readRDS("datasets/Myrtales_full_dataset.Rdata") 
master_table_by_clade <- readRDS("datasets/Myrtales_by_clade_dataset.Rdata") 

#boxplot(master_table$div_rate_eps0.9~master_table$most_common_life_form)

############################################
############################################
# Vars to analyse
div_model_var_to_keep <- c("most_common_life_form","fm_scoring_fruit","seed.length.mean",
                           "niche_through_time","fm_scoring_seed_number","fm_scoring_corolla_diam",
                           "Vol", "CHELSA_bio10_01","CHELSA_bio10_02","CHELSA_bio10_10",
                           "CHELSA_bio10_11","CHELSA_bio10_12","CHELSA_bio10_15","CHELSA_bio10_16",
                           "CHELSA_bio10_17","GLOBAL_SLOPE_10MIN","depthtobedrock2", "meancarbon",
                           "meanpH","meanwatercap")

# Running individual models
subset_master_table <- master_table
results_div <- matrix(nrow=0, ncol=7)
for(var_index_div in 1:length(div_model_var_to_keep)) {
  one_var_subset <- subset_master_table[,c("genus","Family","div_rate_eps0.9",div_model_var_to_keep[var_index_div])]
  colnames(one_var_subset)[4] <- "var_to_test"
  fail <- F
  tryCatch(model_div_full <- phylolm(div_rate_eps0.9~ var_to_test, data=one_var_subset, phy=tree), error = function(e) { fail <<- TRUE})
  if(fail) {
    var <- div_model_var_to_keep[var_index_div]
    r.sq <- p <- slope <- n_points <- "fail"
    n_points <- length(model_div_full$residuals)
  } else {
    
    #div <- subset_master_table$div_rate_eps0.9
    #lifeform <- subset_master_table$most_common_life_form
    #names(div) <- names(lifeform) <- rownames(subset_master_table)
    #pruned_tree <- keep.tip(tree, rownames(subset_master_table))
    #x <- phylANOVA(pruned_tree, lifeform, div)
    #boxplot(subset_master_table$div_rate_eps0.9~subset_master_table$most_common_life_form)
    
    var <- div_model_var_to_keep[var_index_div]
    r.sq <- round(model_div_full$r.squared,3)
    p <- round(unname(summary(model_div_full)$coefficients[,4][2]),3)
    slope <- unname(summary(model_div_full)$coefficients[,1][2])
    slope <- ifelse(slope>0,"+","-")
    n_points <- length(model_div_full$residuals) 
    
    if(class(one_var_subset$var_to_test)=="numeric") {
      #-------------------------
      # Plots:
      pdf(paste0("results/h1/Myrtales_div_",var,".pdf"))
      par(mar=c(3,2,1,0.5))
      par(lwd=.3)
      plot(one_var_subset$div_rate_eps0.9~one_var_subset$var_to_test, bty="n", pch=21, bg="grey", xaxt="n", yaxt="n",xlab="", ylab="",  cex=1, lwd=0.5)
      # adding axes
      par(lwd=.6)
      axis(side=1, tick = TRUE, line = 0, lwd = .4, cex.axis=.6, tcl=NA, mgp=c(0,0.1,0),las=1, cex.lab=.4)
      axis(side=2, tick = TRUE, line = 0, lwd = .4, cex.axis=.6, tcl=NA, mgp=c(2.5,0.2,0),las=1)
      
      mtext(text=paste0(var), side=1, line=1, cex=.5)
      par(las=0)
      mtext(text="(log) net-diversification rate", side=2, line=1.5, cex=.5)
      par(lwd=.5)
      abline(model_div_full)
      dev.off()
    }
  }
  results_div <- rbind(results_div,c("Myrtales", n_points, "div.rates", var, r.sq, p, slope))
}

colnames(results_div) <- c("group","n_points","d_var","predictor","r.sqr","p","slope")
write.csv(results_div, file="results/h1/Myrtales_results_div.csv", row.names=F)



#################
# Clade specific

full_results_div <- list()
for(i in 1:length(master_table_by_clade)) {
  subset_master_table <- master_table_by_clade[[i]]
  one_clade <- names(master_table_by_clade)[i]
  
  results_div <- matrix(nrow=0, ncol=7)
  for(var_index_div in 1:length(div_model_var_to_keep)) {
    one_var_subset <- subset_master_table[,c("genus","Family","div_rate_eps0.9",div_model_var_to_keep[var_index_div])]
    colnames(one_var_subset)[4] <- "var_to_test"
    fail <- F
    tryCatch(model_div_full <- phylolm(div_rate_eps0.9~ var_to_test, data=one_var_subset, phy=tree), error = function(e) { fail <<- TRUE})
    if(fail) {
      var <- div_model_var_to_keep[var_index_div]
      r.sq <- p <- slope <- n_points <- "fail"
    } else {
      var <- div_model_var_to_keep[var_index_div]
      r.sq <- round(model_div_full$r.squared,3)
      p <- round(unname(summary(model_div_full)$coefficients[,4][2]),3)
      slope <- unname(summary(model_div_full)$coefficients[,1][2])
      slope <- ifelse(slope>0,"+","-")
      n_points <- length(model_div_full$residuals) 
      
      if(class(one_var_subset$var_to_test)=="numeric") {
        #-------------------------
        # Plots:
        pdf(paste0("results/h1/",one_clade,"_div_",var,".pdf"))
        par(mar=c(3,2,1,0.5))
        par(lwd=.3)
        plot(one_var_subset$div_rate_eps0.9~one_var_subset$var_to_test, bty="n", pch=21, bg="grey", xaxt="n", yaxt="n",xlab="", ylab="",  cex=1, lwd=0.5)
        # adding axes
        par(lwd=.6)
        axis(side=1, tick = TRUE, line = 0, lwd = .4, cex.axis=.6, tcl=NA, mgp=c(0,0.1,0),las=1, cex.lab=.4)
        axis(side=2, tick = TRUE, line = 0, lwd = .4, cex.axis=.6, tcl=NA, mgp=c(2.5,0.2,0),las=1)
        
        mtext(text=paste0(var), side=1, line=1, cex=.5)
        par(las=0)
        mtext(text="(log) net-diversification rate", side=2, line=1.5, cex=.5)
        par(lwd=.5)
        abline(model_div_full)
        dev.off()
      }
    }
    results_div <- rbind(results_div,c(one_clade, n_points, "div.rates", var, r.sq, p, slope))
  }
  
  tmp_full_results_div <- as.data.frame(results_div)
  colnames(tmp_full_results_div) <- c("family","n_points","ind_var","predictor","r.sqr","p","slope")
  full_results_div[[i]] <- tmp_full_results_div
}

full_results_div <- do.call(rbind, full_results_div)
write.csv(full_results_div, file="results/h1/clade_specifc_results_div.csv", row.names=F)
