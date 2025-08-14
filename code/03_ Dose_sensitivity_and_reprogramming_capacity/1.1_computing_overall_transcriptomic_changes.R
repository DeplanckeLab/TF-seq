# title: "Quantification of overall transcriptomic changes upon TF overexpression"
# input dta: "D0 regressed, G1-phase TF atlas"
# author: "Wangjie Liu"



setwd("./")


library(Seurat)
library(tidyverse)
library(dplyr)

# functions
calculate_mean_D0s <- function(pca, seu, Phase){
  mean_D0s <- colMeans(pca[[Phase]][colnames(seu[[Phase]])[seu[[Phase]]$TF %in% c("D0","D0_confluent")],]) 
  return(mean_D0s)
} 

add_meanD0_to_pca_matrix <- function(pca, meanD0, Phase){
  matrix_new <- rbind(pca[[Phase]], meanD0[[Phase]])
  rownames(matrix_new)[nrow(matrix_new)] <- "mean_D0s"
  return(matrix_new)
} 



# load TF atlas
seu <- readRDS("results/C3H10_10X_all_exps_D0regressed_integrated_dosealigned.rds")

# split cells by phase
Phase_to_use <- c("G1","G2M","S")
seu.list <- lapply(Phase_to_use, function(x){
  seu.sub <- subset(seu, Phase_corrected == x)
  return(seu.sub)
})
names(seu.list) <- Phase_to_use

# compute correlation
pca.list <- lapply(seu.list, function(x) return(x@reductions$pca@cell.embeddings[, 1:200]))
mean_D0s <- lapply(Phase_to_use, function(x) calculate_mean_D0s(pca = pca.list, seu = seu.list, Phase = x ))
names(mean_D0s) <- Phase_to_use

matri_for_corr <- lapply(Phase_to_use, function(x) add_meanD0_to_pca_matrix(pca = pca.list, meanD0 = mean_D0s, Phase = x))
names(matri_for_corr) <- Phase_to_use
cor_res <- lapply(matri_for_corr, function(x) return(cor(t(x))[, "mean_D0s"]))
names(cor_res) <- Phase_to_use
cor_res$G1[length(cor_res$G1)]

# output the correlation, normalize to control, and transform to overall transcriptomic changes
df <- data.frame(TF = seu.list$G1$TF)
df$TF <- as.factor(df$TF)
df$correlation <- cor_res$G1[rownames(df)]
mean_D0_acti <- mean(-df$correlation[df$TF %in% c("D0","D0_confluent")])
df$Overall_transcriptomic_change <- -df$correlation
df$Overall_transcriptomic_change <- df$Overall_transcriptomic_change - mean_D0_acti

saveRDS(df, file = "results/df_overall_transcriptomic_changes.rds")
