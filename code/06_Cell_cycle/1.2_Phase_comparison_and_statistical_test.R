# title: "Cell cycle analysis"
# input: TF atlas with phase corrected
# author: "Wangjie Liu"
# date: "2024/10/23"

setwd("./")

library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggdensity)
library(viridis)
library(ks)
library(Hmisc)
library(dplyr)

# load functions
source("code/05_Cell_cycle/1.1_cell_cycle_utils.R")

# load data, annotation, and add metadata
data.annot <- fread("metadata/GRCm38_gene_annot.tsv", data.table = F) 
rownames(data.annot) <- data.annot$ensembl_id

seu <- readRDS("results/C3H10_10X_all_exps_merged_genefiltered_integrated")
# exclude TF combinations adn reference cells
seu <- subset(seu, TF %nin% c("Pparg-Runx2",
                                "Mycn-Runx2",
                                "Mycn-Myog",
                                "Mycn-Pparg",
                                "Cebpa-Mycn",
                                "Cebpa-Myog",
                                "Cebpa-Pparg",
                                "Adipo_ref",
                                "Myo_ref")) 

seu$dose <- seu$Log_Vector_UMI
cc.meta <- seu@meta.data[, c("TF","dose","Phase","Phase_corrected","S.Score","G2M.Score")]



Phase_to_use <- "Phase_corrected" # Figure 6b
# Phase_to_use <- "Phase" # for Supplementary Figure 7

##---------------------------------1. Barplot. proportion of cycling cells
P_ALLCells <- barplot_TopCyclingTF(cc.meta, Phase_to_use) 

#  plot significance in labels of y.axis
FDR.threshold <- 0.25
Only_Significant <- T
P_ALLCells_new <- barplot_TopCyclingTF_step2(P_ALLCells$df_PhaseperTF, 
                                             q_value_threshold = FDR.threshold,
                                             plot_only_significant = Only_Significant)

ggsave(P_ALLCells_new$p_h, height = 4.79, width = 15, 
       filename = paste0("figures/Fig6n-", Phase_to_use, "_FisherTest_FDR",FDR.threshold,"_Sig-", Only_Significant,"_total_nCell.pdf"))



##---------------------------------2. 1D, 2D density plot and Scatter plot of S.score VS G2M.score

### kde.test
TFois <- table(cc.meta$TF)[table(cc.meta$TF) > 50] %>% names()
pvalues.kde <- lapply(TFois[TFois!= "D0_confluent"], function(x) 
  Pvalue.KDE_TEST(TFoi = x, data = cc.meta, D0_to_use = "D0_confluent"))
names(pvalues.kde) <- TFois[TFois!= "D0_confluent"]
pvalues.kde <- unlist(pvalues.kde)
pvalues.kde.adj <- p.adjust(pvalues.kde, method = "fdr")
pvalues.kde.adj
df.kdeTest <- as.data.frame(pvalues.kde.adj)
df.kdeTest$TF <- rownames(df.kdeTest)
df.kdeTest$pvalues.kde <- pvalues.kde
view(df.kdeTest)
df.kdeTest[df.kdeTest$pvalues.kde.adj < 0.25,]$TF


### plot 2-D contour density for TFs of interest
Plots.contour <- lapply(TFois[TFois!= "D0_confluent"], function(x) 
  Plot_2D_contour(TFoi = x, data = cc.meta, D0_to_use =  "D0_confluent", OutputDir = "figures/"))



### scatter plots of cellcylescore for TFs of interest
lapply(TFois, function(x) 
  Plot_Scatter_CellCycleScore(TFoi = x, data = cc.meta, OutputDir = "figures/", plot_vector = T, plot_phase = F ))

# for D0 and D0_confluent (color by phase)
lapply(c("D0","D0_confluent"), function(x) 
  Plot_Scatter_CellCycleScore(TFoi = x, data = cc.meta, OutputDir = "figures/", plot_vector = F, plot_phase = T, Phase.use = "Phase"))

lapply(c("Yap1","Atf3"), function(x) 
  Plot_Scatter_CellCycleScore(TFoi = x, data = cc.meta, OutputDir = "figures/", plot_vector = F, plot_phase = T, Phase.use = "Phase_corrected"))



### 1-D density plots of cellcylescore for TFs of interest
lapply(TFois[TFois!= "D0_confluent"], function(x) 
  Plot_1D_density(TFoi = x, data = cc.meta, OutputDir = "figures/"))


### wilcox test on cellcyclescore
df.pvalues_wilcox <- lapply(TFois[TFois!= "D0_confluent"], function(x)
  Pvalue.WILCOX_TEST(TFoi = x, data = cc.meta, D0_to_use = "D0_confluent"))

df.pvalues_wilcox <- data.table::rbindlist(df.pvalues_wilcox)
df.pvalues_wilcox$padj.wilcox_S <- p.adjust(df.pvalues_wilcox$pvalue.wilcox_S, method = "fdr")
df.pvalues_wilcox$padj.wilcox_G2M <- p.adjust(df.pvalues_wilcox$pvalue.wilcox_G2M, method = "fdr")
rownames(df.pvalues_wilcox) <- df.pvalues_wilcox$TF
df.pvalues_wilcox$TF <- NULL

### combine all info together 
df.final <- cbind(df.kdeTest, df.pvalues_wilcox)
df.final <- df.final[, c(2, 1, 3:7)]
colnames(df.final) <- c("TF",
                        "padj.kde",
                        "p_value.kde",
                        "p_value.wilcox.S",
                        "p_value.wilcox.G2M",
                        "padj.wilcox_S",
                        "padj.wilcox_G2M")

### add information of nCell, proportion of each phase
df.final <- df.final %>% 
  group_by(TF) %>% 
  mutate(Total_number_of_cells = table(seu$TF)[TF]) %>%
  mutate(Number_of_G1.phase = length(seu$Phase[seu$TF == TF & seu$Phase == "G1"])) %>%
  mutate(Number_of_S.phase = length(seu$Phase[seu$TF == TF & seu$Phase == "S"])) %>%
  mutate(Number_of_G2M.phase = length(seu$Phase[seu$TF == TF & seu$Phase == "G2M"])) %>%
  mutate(Number_of_G1.phase_corrected = length(seu$Phase[seu$TF == TF & seu$Phase_corrected == "G1"])) %>%
  mutate(Number_of_S.phase_corrected = length(seu$Phase[seu$TF == TF & seu$Phase_corrected == "S"])) %>%
  mutate(Number_of_G2M.phase_corrected = length(seu$Phase[seu$TF == TF & seu$Phase_corrected == "G2M"]))

write.csv(df.final, file = "results/Table_1D_2D_statistical_test_for_cellcycle.csv",
          quote = F)


##---------------------------------3. Interplay between adipogenesis and cell cycle (Figure 6e)
adipo.markers <-  c("Fabp4", "Lpl", "Pparg", "Lipe", "Adipoq", "Cd36",
                    "Plin4", "Plin2", "Plin1", "Cebpa", "Cebpb",
                    "Cidec", "Cidea")
adipo.markers <-  data.annot[data.annot$gene_short_name %in% adipo.markers ,"ens_id"]

DefaultAssay(seu) <- "RNA"
seu <- AddModuleScore(seu, features = list(adipo.markers), name = "adiposcore")

# DE between all cells of Mycn, Cebpa, and Pparg
seu.adipo <- subset(seu, TF %in% c("Mycn","Cebpa","Pparg"))
DefaultAssay(seu.adipo) <- "RNA"
Idents(seu.adipo) <- "TF"
DEG.adipo <- FindAllMarkers(seu.adipo, logfc.threshold = 0.25, min.pct = 0.25, only.pos = T)
DEG.adipo$gene_name <- data.annot[DEG.adipo$gene, "gene_short_name"]
DEG.adipo$diff.pct <- abs(DEG.adipo$pct.1 - DEG.adipo$pct.2)
DEG.adipo <- DEG.adipo[DEG.adipo$p_val_adj < 0.05,]

# plot normalized expression of p21 (all Mycn cells versus all Cebpa or Pparg cells)
seu.adipo_D0 <- subset(seu, TF %in% c("D0_confluent","Cebpa","Pparg","Mycn"))
data.annot[data.annot$gene_short_name == "Cdkn1a",] # p21
table(seu.adipo_D0$TF)
seu.adipo_D0$TF <- factor(seu.adipo_D0$TF, levels = c("D0","Cebpa","Pparg","Mycn"))

seu.adipo_D0$p21 <- seu.adipo_D0@assays$RNA@data["ENSMUSG00000023067",]
pairwise.t.test(seu.adipo_D0$p21, g = seu.adipo_D0$TF, p.adjust.method = "fdr")

# DE between S-G2M cells of Mycn, Cebpa, and Pparg
cells.oi <- colnames(seu.adipo)[seu.adipo$Phase_corrected != "G1"]
seu.adipo.CC  <- subset(seu.adipo, cells = cells.oi)
Idents(seu.adipo.CC) <- "TF"
DEG.adipo.CC <- FindAllMarkers(seu.adipo, only.pos = T)
DEG.adipo.CC$gene_name <- data.annot[rownames(DEG.adipo.CC), "gene_short_name"]


tfoi <- "Mycn"
tfoi <- "Cebpa"
tfoi <- "Pparg"

gene.oi <- "Cdkn1a" # Mycn - CC -adipo
Phase_to_use <- "Phase_corrected"
phase_combine <- T
scale_gene <- F
check_gene <- T # check gene of interest (T) or adiposcore (F)

if(phase_combine == T){
  seu.oi <- subset(seu, TF == c("D0_confluent",tfoi) ) # & Phase != "G1"
  seu.oi$TF_dosage <- seu.oi$dose
  
  # DefaultAssay(seu.oi) <- "RNA"
  # seu.oi <- AddModuleScore(seu.oi, features = list(adipo.markers), name = "adiposcore")
  
  plotdata.gene <- FetchData(seu.oi, vars = c("TF_dosage",
                                              Phase_to_use,
                                              "S.Score",
                                              "G2M.Score", 
                                              data.annot[data.annot$gene_short_name == gene.oi, "ens_id"], 
                                              "adiposcore1"))
  
  colnames(plotdata.gene)[2] <- "Phase"
  colnames(plotdata.gene)[ncol(plotdata.gene)-1] <- "gene.oi"
  colnames(plotdata.gene)[ncol(plotdata.gene)] <- "adiposcore"
  
  # Bins: take max and min and then divide each dimension in 5 bins
  # Missing values: if a bin doesn't contain at least 3 cells => NaN
  
  
  n_bins <- 5
  # Bin s_score and g2m_score, dosage into discrete intervals
  plotdata.gene$dosage_bin <- cut(plotdata.gene$TF_dosage, 
                                  breaks = c(-Inf, log1p(1), seq(log1p(1), max(plotdata.gene$TF_dosage), length.out = n_bins)[-1]), 
                                  include.lowest = T,
                                  right = F)
  
  
  # Calculate mean expression level in cells of individual TF, for each bin and phase (G1 or CC) 
  plotdata.gene$Phase[plotdata.gene$Phase !="G1"] <- "CC"
  plotdata.gene$Phase <- factor(plotdata.gene$Phase, levels = c("G1","CC"))
  
  if (check_gene == T){
    df_binned <- plotdata.gene %>%
      group_by(Phase, dosage_bin) %>%
      filter(n() >= 3) %>%
      summarise(mean_expression = mean(gene.oi, na.rm = TRUE))
    df_binned$mean_expression_scale <- df_binned$mean_expression / max(df_binned$mean_expression)
  } else {
    df_binned <- plotdata.gene %>%
      group_by(Phase, dosage_bin) %>%
      filter(n() >= 3) %>%
      summarise(mean_expression = mean(adiposcore, na.rm = TRUE))
    df_binned$mean_expression_scale <- df_binned$mean_expression / max(df_binned$mean_expression)
  }
  
  
  if (scale_gene == F){
    p1 <- ggplot(df_binned, aes(x = dosage_bin, y = Phase,fill = mean_expression)) +
      geom_tile() +
      scale_fill_viridis_c() +
      cowplot::theme_cowplot()+
      labs(subtitle = paste0(tfoi, " cells"))
    print(p1)
  } else {
    p1 <- ggplot(df_binned, aes(x = dosage_bin, y = Phase, fill = mean_expression_scale)) +
      geom_tile() +
      scale_fill_viridis_c() +
      cowplot::theme_cowplot()+
      labs(subtitle = paste0(tfoi, " cells"))
    print(p1)
  }
}



if(check_gene == T){
  ggsave(p1, width = 6, height = 2,
         filename = paste0("output_version2_2023_July7/9-cell-cycle/Heatmap_lineage_genes_scores/",gene.oi, "_in_", tfoi,"-",Phase_to_use,".pdf"))
} else {
  ggsave(p1, width = 6, height = 2,
         filename = paste0("output_version2_2023_July7/9-cell-cycle/Heatmap_lineage_genes_scores/Adiposcore_in_", tfoi,"-",Phase_to_use,".pdf"))
}






##-----------------------------------------------------Figure 6e
# combine TFs of interest together
Phase_to_use <- "Phase_corrected"
phase_combine <- T
scale_gene <- F
check_gene <- T

tfois <- c("Mycn","Cebpa","Pparg")
gene.oi <- "Cdkn1a"

if(phase_combine == T){
  seu.oi <- subset(seu, TF == c("D0_confluent",tfois) ) 
  seu.oi$TF_dosage <- seu.oi$dose
  
  plotdata.gene <- FetchData(seu.oi, vars = c("TF",
                                              "TF_dosage",
                                              Phase_to_use,
                                              "S.Score",
                                              "G2M.Score", 
                                              data.annot[data.annot$gene_short_name == gene.oi, "ens_id"], 
                                              "adiposcore1"))
  
  colnames(plotdata.gene)[1] <- "TF" 
  colnames(plotdata.gene)[3] <- "Phase"
  colnames(plotdata.gene)[ncol(plotdata.gene)-1] <- gene.oi
  colnames(plotdata.gene)[ncol(plotdata.gene)] <- "adiposcore"
  
  # Bins: take max and min and then divide each dimension in 5 bins
  # Missing values: if a bin doesn't contain at least 3 cells => NaN
  
  plotdata.gene$Phase[plotdata.gene$Phase !="G1"] <- "CC"
  plotdata.gene$Phase <- as.vector(plotdata.gene$Phase)
  plotdata.gene$Phase <- factor(plotdata.gene$Phase, levels = c("G1","CC"))
  
  plotdata.gene.ctr <- plotdata.gene[plotdata.gene$TF == "D0_confluent",]
  plotdata.gene.ctr$dosage_bin <- as.integer(0)
  plotdata.gene.adipo <- plotdata.gene[plotdata.gene$TF != "D0_confluent",]
  
  n_bins <- 5
  # Bin s_score, g2m_score, dose into discrete intervals for each TF
  
  
  plotdata.gene.adipo <- plotdata.gene.adipo %>% 
    group_by(TF) %>% 
    mutate(dosage_bin = cut(TF_dosage, 
                            breaks = c( seq(log1p(1), max(TF_dosage), length.out = n_bins)), 
                            include.lowest = T,
                            right = F, 
                            labels = F))
  
  
  # Calculate mean expression level for each bin and phase (G1 or CC) for each TF-Ctr pair
  df <- lapply(tfois, function(x){
    df.tfoi_ctr <- rbind(plotdata.gene.adipo[plotdata.gene.adipo$TF == x,], plotdata.gene.ctr)
    df.tfoi_ctr$meta <- paste0(x, "-", df.tfoi_ctr$Phase)
    return(df.tfoi_ctr)
  })
  df <- data.table::rbindlist(df)
  
  
  if (check_gene == T){
    df_binned <- df %>%
      group_by(meta, dosage_bin) %>%
      filter(n() >= 3) %>%
      summarise(mean_expression = mean(gene.oi, na.rm = TRUE))
    df_binned$mean_expression_scale <- df_binned$mean_expression / max(df_binned$mean_expression)
  } else {
    df_binned <- df %>%
      group_by(meta, dosage_bin) %>%
      filter(n() >= 3) %>%
      summarise(mean_expression = mean(adiposcore, na.rm = TRUE))
    df_binned$mean_expression_scale <- df_binned$mean_expression / max(df_binned$mean_expression)
  }
  
  
  if (scale_gene == F){
    p1 <- ggplot(df_binned, aes(x = dosage_bin, y = meta, fill = mean_expression)) +
      geom_tile() +
      scale_fill_viridis_c() +
      cowplot::theme_cowplot()+
      labs(subtitle = gene.oi)
    print(p1)
  } else {
    p1 <- ggplot(df_binned, aes(x = dosage_bin, y = meta, fill = mean_expression_scale)) +
      geom_tile() +
      scale_fill_viridis_c() +
      cowplot::theme_cowplot()+
      labs(subtitle = gene.oi)
    print(p1)
  }
}



##-----------------------------------------------------Figure 6d & ED_Fig8e
# Fraction of cycling cells across dose bins
Phase_to_use <- "Phase_corrected"

plot_fraction_per_dose_bin <- function(tfoi, n_bins = 5, phase){
  seu.oi <- subset(seu, TF %in% c(tfoi, "D0_confluent") )
  seu.oi$dose[seu.oi$TF == "D0_confluent" ] <- 0
  
  plotdata.density <- FetchData(seu.oi, vars = c("TF",phase,"dose"))
  
  plotdata.density$dosage_bin <- cut(plotdata.density$dose, 
                                     breaks = c(-Inf, log1p(1), seq(log1p(1), max(plotdata.density$dose), length.out = n_bins)[-1]), 
                                     include.lowest = T,
                                     right = F)
  
  plotdata.density$is.CC_phase <- plotdata.density[[phase]]
  
  
  plotdata.density <- plotdata.density %>% 
    group_by(dosage_bin) %>%
    mutate(nCell_total = n()) %>%
    ungroup(dosage_bin) %>%
    group_by(dosage_bin, is.CC_phase) %>%
    mutate(nCell = n())
  
  plotdata.density$Fraction_cc <- plotdata.density$nCell/plotdata.density$nCell_total*100
  
  plotdata.density$is.CC_phase <- factor(plotdata.density$is.CC_phase, levels = c("G1","G2M","S"))
  
  
  plotdata.density <- plotdata.density[, c("dosage_bin","Fraction_cc","is.CC_phase","nCell_total")]
  plotdata.density <- plotdata.density[!duplicated(plotdata.density),]
  
  # if there are less than 9 total cells in a bin, set that dose bin to NA/empty
  df.nCell.bin <- plotdata.density[,c("dosage_bin","Fraction_cc","nCell_total")]
  df.nCell.bin <- df.nCell.bin[!duplicated(df.nCell.bin),]
  bins_fewCells <- df.nCell.bin$dosage_bin[df.nCell.bin$nCell_total <9]
  
  plotdata.density$Fraction_cc[plotdata.density$dosage_bin %in% bins_fewCells] <- NA
  print(plotdata.density)
  p.fraction <- ggplot(plotdata.density, mapping = aes(x = dosage_bin, y = Fraction_cc, fill = is.CC_phase))+
    geom_bar( position = "stack", stat = "identity")+
    scale_fill_manual(values = c("gray90", "#a6cee3","#1f78b4"))+
    cowplot::theme_cowplot()+
    labs(title = tfoi)
  p.fraction
  return(p.fraction)
}
plot_fraction_per_dose_bin("Myod1", n_bins = 5, phase = Phase_to_use)

# all TFs :
# 1. that have pro-cell cycle effect 
# 2. and tested significant on all cells
# 3. and have more than 
tfois <- c("E2f2","T","Mycn","Runx2","Rela","Klf4","Zfp78","Isx","Ahrr","Tfec","Fosl1","Myod1","Zkscan14","Dlx3","Pax9")
P_list.fraction <- lapply(tfois, function(x) plot_fraction_per_dose_bin(x, n_bins = 5, phase = Phase_to_use))
names(P_list.fraction) <- tfois
P_list.fraction

pdf("output_version2_2023_July7/9-cell-cycle/Fraction_of_SG2M_phaseAdjusted_5_dose_bins.pdf", width = 10) 
gridExtra::marrangeGrob(P_list.fraction, ncol = 2, nrow = 2)
dev.off()



