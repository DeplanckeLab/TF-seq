# title: "Differential expression (DE) analysis for the inflammatory cluster and TF cells in G1-phase functional TF atlas"
# author: "Wangjie Liu"
# date: "2024/10/23"


setwd("./")


library(tidyverse)
library(Seurat)
library(snowfall)


# load gene annotations
data.annot <- read.table("data/GRCm38.96_Vector_data.annot.txt", sep = "\t")
rownames(data.annot) <- data.annot$ens_id
colnames(data.annot)

# load G1-phase functional atlas
seu <- readRDS("results/C3H10_10X_all_exps_integrated_functional_TF_atlas_G1_Phase.rds")


## -------------------------------------------DE for clusters, e.g., cluster 8 with inflammatory geen expressions
Idents(seu) <- seu$integrated_snn_res.0.2
table(seu$integrated_snn_res.0.2[seu$TF %in% c("D0","D0_confluent")])
clusters.ctr <- table(seu$integrated_snn_res.0.2[seu$TF %in% c("D0","D0_confluent")])/sum(table(seu$integrated_snn_res.0.2[seu$TF %in% c("D0","D0_confluent")]))
clusters.ctr <- names(clusters.ctr)[clusters.ctr > 0.1]

seu_diffexp <- subset(seu, integrated_snn_res.0.2 %in% c(8, clusters.ctr))
table(seu_diffexp$integrated_snn_res.0.2)
rm(seu)
seu_diffexp$group <- seu_diffexp$integrated_snn_res.0.2 %>% as.vector()
seu_diffexp$group[seu_diffexp$group %in% clusters.ctr] <- "Ctr"
seu_diffexp$group[seu_diffexp$group == 8] <- "Cluster_8"
print(table(seu_diffexp$group))
seu_diffexp$group <- factor(seu_diffexp$group, levels = c("Ctr", "Cluster_8")) # set the order of Ctr and TF in order to use Ctr as contrast

counts_diffexp <- seu_diffexp@assays$RNA@counts
metadata_diffexp <- seu_diffexp@meta.data %>%
  dplyr::select(c(group, batch_overall))
# dplyr::select(c(TF, batch, log_vector))

metadata_diffexp$batch_overall <- as.factor(metadata_diffexp$batch_overall)
dge <- edgeR::DGEList(counts = counts_diffexp)
dge <- edgeR::calcNormFactors(dge) # calculate scaling factors to convert raw library size to effective library size
metadata_diffexp$cdr <- scale(Matrix::colMeans(counts_diffexp > 0)) # center or scale the columns of matrix

# cdr -> scaled gene expression mean
formula <- "~cdr + group" # tried with vector, but non-monotonicity messes the linear regression up
model_batch <- T
# add extra covariates if necessary
if(model_batch && length(levels(metadata_diffexp$batch_overall)) > 1) {
  formula <- paste0(formula, " + batch_overall")
}

# create design
design <- model.matrix(as.formula(formula), metadata_diffexp)
colnames(design)
## notice that the first level of TF group like Ctr is taking as base automatically.
#colnames(design)[startsWith(colnames(design), x)] <- x
dge <- edgeR::estimateDisp(dge, design = design)# maximize negative binomial likelihood
# to give the estimate the common, trended and tagwise dispersions across all tags
print(colnames(design))
Coef <- colnames(design)[c(3)] # group Cluster_8

fit <- edgeR::glmQLFit(dge, design = design)
qlf <- edgeR::glmQLFTest(fit, coef = Coef)
tt <- edgeR::topTags(qlf, n = Inf)
score <- as.data.frame(tt$table)
score$ens_id <- rownames(score)
score$gene_short_name <- data.annot[rownames(score), "gene_short_name"]
score$group <- "Cluster_8_versus_clusters_ctr_at_R0.2"
print(colnames(score))


saveRDS(score, 
        file = "results/DEGs_glmQLFtest_by_Cluster_8_VS_Clusters_Ctr_at_R0.2_G1_Phase_batchcovariate.rds")



##--------------------------------------------DE for TF cells
Idents(seu) <- seu$TF
TFois <- table(seu$TF) %>% names() 
TFois <- TFois[!TFois %in% c("D0","D0_confluent")]

sfInit(parallel = TRUE, cpus = 16)
sfLibrary(Seurat)
sfLibrary(tidyverse)
sfLibrary(ggplot2)

sfExport("data.annot", local=FALSE)
sfExport("seu", local=FALSE)
sfExport("TFois", local=FALSE)

glmQLFfit_test <- function(x){
  seu_diffexp <- subset(seu, TF %in% c("D0", "D0_confluent", x))
  rm(seu)
  cat(paste0("D0 + D0_conlfuent + ",x," have ", ncol(seu_diffexp), " cells"))
  
  seu_diffexp$TF[seu_diffexp$TF %in% c("D0","D0_confluent")] <- "Ctr"
  print(table(seu_diffexp$TF))
  seu_diffexp$TF <- factor(seu_diffexp$TF, levels = c("Ctr", x)) # set the order of Ctr and TF in order to use Ctr as contrast
  
  counts_diffexp <- seu_diffexp@assays$RNA@counts
  metadata_diffexp <- seu_diffexp@meta.data %>%
    dplyr::select(c(TF, batch_overall))
  
  metadata_diffexp$batch_overall <- as.factor(metadata_diffexp$batch_overall)
  dge <- edgeR::DGEList(counts = counts_diffexp)
  dge <- edgeR::calcNormFactors(dge) # calculate scaling factors to convert raw library size to effective library size
  metadata_diffexp$cdr <- scale(Matrix::colMeans(counts_diffexp > 0)) # center or scale the columns of matrix
  
  # cdr -> scaled gene expression mean
  formula <- "~cdr + TF" # tried with vector, but non-monotonicity messes the linear regression up
  model_batch <- T
  # add extra covariates if necessary
  if(model_batch && length(levels(metadata_diffexp$batch_overall)) > 1) {
    formula <- paste0(formula, " + batch_overall")
  }
  
  # create design
  design <- model.matrix(as.formula(formula), metadata_diffexp)
  colnames(design)
  ## notice that the first level of TF group like Ctr is taking as base automatically.
  #colnames(design)[startsWith(colnames(design), x)] <- x
  dge <- edgeR::estimateDisp(dge, design = design)# maximize negative binomial likelihood
  # to give the estimate the common, trended and tagwise dispersions across all tags
  print(colnames(design))
  Coef <- colnames(design)[c(3)]
  
  fit <- edgeR::glmQLFit(dge, design = design)
  qlf <- edgeR::glmQLFTest(fit, coef = Coef)
  tt <- edgeR::topTags(qlf, n = Inf)
  score <- as.data.frame(tt$table)
  score$ens_id <- rownames(score)
  score$gene_short_name <- data.annot[rownames(score), "gene_short_name"]
  score$TF <- x
  print(colnames(score))
  return(score)
}


df_DEGs.TFois <- sfLapply(TFois, glmQLFfit_test)
df_DEGs.TFois <- data.table::rbindlist(df_DEGs.TFois)
write.csv(df_DEGs.TFois, file = "results/DEGs_glmQLFtest_by_TF_batchoverall_allD0_functionalCells.csv", quote = F)

sfStop(nostop = F)