# title: "Downstream QC and plotting on the single-cell fluorescence quantification of RNAscope images - KLF4, ESR2"
# input data: "RNAscope_measurements_Fig5_EDFig6.tsv"
# author: "Wangjie Liu"


## Non-confocal image data of KLF4, ESR2 (rep2 - 20240918)
## best Z

library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(gridExtra)

setwd("./")

measurements <- read_delim("data/RNAscope_measurements_Fig5_EDFig6.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(measurements) <- measurements$`Object ID`

ggplot(measurements)+geom_density(mapping = aes(x =`Nucleus: Area µm^2`))+geom_vline(xintercept = 100)
# need to decide empirically based on the data
ggplot(measurements)+geom_density(mapping = aes(x =`Channel_1: Nucleus: Mean`)) + geom_vline(xintercept = 800)# DAPI signal
ggplot(measurements)+geom_density(mapping = aes(x =`Channel_1: Nucleus: Median`)) # DAPI signal
ggplot(measurements)+geom_density(mapping = aes(x =`Channel_1: Cytoplasm: Mean`)) # DAPI signal
ggplot(measurements)+geom_density(mapping = aes(x = log10(`Channel_1: Cell: IntDen Corrected`)))  #  DAPI signal
ggplot(measurements)+geom_density(mapping = aes(x = `Nucleus/Cell area ratio`)) # high Nucleus/Cell area ratio as an indicative of cell clumps or multinucleated cells

# filter out outlier with very small cell size  
measurements <- measurements[measurements$`Nucleus: Area µm^2` >= 100 ,] # size

# filter out outlier with weak DAPI signals (wrong cell segmentation) 
## for best Z
measurements <- measurements[measurements$`Channel_1: Nucleus: Mean` >= 800,] # DAPI signals

# exclude wrong segmentation or cells without intact nuclei indicated by low sum intensity of DAPI
## for best Z
measurements <- measurements[log10(measurements$`Channel_1: Cell: IntDen Corrected`) >= 6.5, ] # sum of background corrected DAPI signals

# filter out cells that in cell clumps indicated by high Nucleus/Cell area ratio & filter out artifacts that come from dirt indicated by high median DAPI signals in cytoplasm
measurements <- measurements[measurements$`Nucleus/Cell area ratio` <= 0.75 & measurements$`Channel_1: Cytoplasm: Mean` <= 1000, ]

measurements$analysis <- "Fluorescence"
measurements$analysis[grep("DPC", measurements$Image)] <- "DPC" 
table(measurements$analysis)
# DPC Fluorescence 
# 6417         6841 
measurements2 <- measurements[measurements$analysis != "DPC",] 
measurements2 <- measurements2[, c(1, 2, grep("Cell: Mean",colnames(measurements2)), grep("Cell: IntDen Corrected",colnames(measurements2)),  grep("Cell: ROI Area",colnames(measurements2)))]
measurements2$Well <- substring(measurements2$Image, 33, 39)
table(measurements2$Well)

# R02-C02 R02-C03 R02-C04 R02-C05 R03-C02 R03-C03 R03-C04 R03-C05 R04-C02 R04-C03 R04-C04 R05-C02 R05-C03 R05-C04 R06-C02 R06-C03 R06-C04 
# 312     419     410     356     454     431     419     385     323     489     512     418     356     481     346     415     315 


measurements2$condition <- "Pos_neg_controls"
measurements2$condition[measurements2$Well %in% c(paste0("R02-C0",2:4))] <- "Ch2=Net1; Ch3=WPRE; Ch4=Esr2-ORF"
measurements2$condition[measurements2$Well %in% c(paste0("R03-C0",2:4))] <- "Ch2=Net1; Ch3=WPRE; Ch4=Cdsn"
measurements2$condition[measurements2$Well %in% c(paste0("R04-C0",2:4))] <- "Ch2=Gng12; Ch3=WPRE; Ch4=Aspn"
measurements2$condition[measurements2$Well %in% c(paste0("R05-C0",2:4))] <- "Ch2=Postn; Ch3=WPRE; Ch4=Glul"
measurements2$condition[measurements2$Well %in% c(paste0("R06-C0",2:4))] <- "Ch2=Postn; Ch3=WPRE; Ch4=Glul"
table(measurements2$condition)


measurements2$cell_group <- "WT-C3H10 cells"
measurements2$cell_group[measurements2$Well %in% c(paste0("R0",2:4,"-C02"))] <- "Esr2 cells"
measurements2$cell_group[measurements2$Well %in% c(paste0("R0",5:6,"-C02"))] <- "Klf4 cells"
measurements2$cell_group[measurements2$Well %in% c(paste0("R0",2:6,"-C03"))] <- "mCherry-noDox cells"
measurements2$cell_group[measurements2$Well %in% c(paste0("R0",2:6,"-C04"))] <- "WT-C3H10 cells"
table(measurements2$cell_group)



condition.ois <- unique(measurements2$condition)
condition.ois <- condition.ois[condition.ois != "Pos_neg_controls"]
lapply(condition.ois, function(x){
  condition.oi <- x
  print(condition.oi)
  plot_data <- measurements2[measurements2$condition == condition.oi,]
  intensity.oi <- "Mean BgCorrected"
  colnames(plot_data) <- sub(paste0(": Cell: ", intensity.oi),"",colnames(plot_data))
  colnames(plot_data)
  
  # Channel_0 = DPC
  # Channel_1 = DAPI
  # Channel_2 = Vivid520
  # Channel_3 = Vivid570
  # Channel_4 = Vivid650
  # Channel_5 = BF
  
  plot_data <- as.data.frame(plot_data)
  
  
  # for each well and each channel, remove 1% of the upper intensity (extreme values/outliers)
  cells.out <- c()
  for (Ch.oi in paste0("Channel_", 2:4)){
    for (Well.oi in unique(plot_data$Well)){
      df <- plot_data[plot_data$Well == Well.oi, ]
      cells.out <- c(cells.out, df$`Object ID`[df[, Ch.oi] >= quantile(df[, Ch.oi], 0.99)], df$`Object ID`[df[, Ch.oi] <= quantile(df[, Ch.oi], 0.01)]) 
    }
  }
  
  plot_data <- plot_data[!plot_data$`Object ID` %in% cells.out, ]
  
  ## contamination correction (before log transformation)
  
  model_Ch2 <- lm(data = plot_data[plot_data$cell_group %in% c("mCherry-noDox cells"),], formula = Channel_2 ~ Channel_3)
  model_Ch2
  model_Ch4 <- lm(data = plot_data[plot_data$cell_group %in% c("mCherry-noDox cells"),], formula = Channel_4 ~ Channel_3)
  model_Ch4
  if(model_Ch2$coefficients[2] > 0.01){
    plot_data$Channel_2 <- plot_data$Channel_2 - predict(model_Ch2, newdata = plot_data)
    p1 <- ggplot() +
      geom_point(plot_data, mapping = aes(x=(`Channel_3`), y=(`Channel_2`), col=cell_group), alpha = 0.5,) +
      facet_grid(~cell_group)+
      geom_smooth(plot_data[plot_data$cell_group != "WT-C3H10 cells",], mapping = aes(x=(`Channel_3`), y=(`Channel_2`)))+
      facet_grid(~cell_group)+
      labs(title = plot_data$condition)
  }
  if(model_Ch4$coefficients[2] > 0.01){
    plot_data$Channel_4 <- plot_data$Channel_4 - predict(model_Ch4, newdata = plot_data)
    p2 <-  ggplot() +
      geom_point(plot_data, mapping = aes(x=(`Channel_3`), y=(`Channel_4`), color=cell_group), alpha = 0.5) +
      facet_grid(~cell_group)+
      geom_smooth(plot_data[plot_data$cell_group != "WT-C3H10 cells",], mapping = aes(x=(`Channel_3`), y=(`Channel_4`)))+
      facet_grid(~cell_group)
  }
  grid.arrange(p1,p2)
  
  
  ## log transformation after contamination correction
  plot_data.new <- plot_data
  base.oi <- 10
  plot_data.new[plot_data.new < base.oi] <- base.oi 
  for (i in paste0("Channel_", 2:4)){
    plot_data.new[, i] <- log(plot_data.new[, i], base = base.oi) 
  }
  
  
  # if plot WT and TF cells together
  plot_data.sub <- plot_data.new[plot_data.new$cell_group != "mCherry-noDox cells",]
  plot_data.sub$cell_group |> table() |> print()
  cell_group <- names(table(plot_data.sub$cell_group))
  cell_count <- table(plot_data.sub$cell_group)
  p1 <- ggplot() +
    geom_point(plot_data.sub, mapping = aes(x=(`Channel_3`), y=(`Channel_2`), col=cell_group), alpha = 0.5) +
    geom_smooth(plot_data.sub, mapping = aes(x=(`Channel_3`), y=(`Channel_2`)))+
    labs(title = plot_data.sub$condition, subtitle = paste0(cell_group[1], "=", cell_count[1], "cells;", cell_group[2], "=", cell_count[2]))+
    theme_cowplot()+
    scale_x_continuous(breaks = seq(1, 4, by = 0.5), labels = c(1, 1.5, 2, 2.5,  3, 3.5,  4))+
    scale_y_continuous(breaks = seq(1, 4, by = 0.5), labels = c(1, 1.5, 2, 2.5,  3, 3.5,  4))
  p2 <-  ggplot() +
    geom_point(plot_data.sub, mapping = aes(x=(`Channel_3`), y=(`Channel_4`), col=cell_group), alpha = 0.5) +
    geom_smooth(plot_data.sub, mapping = aes(x=(`Channel_3`), y=(`Channel_4`)))+
    theme_cowplot()+
    scale_x_continuous(breaks = seq(1, 4, by = 0.5), labels = c(1, 1.5, 2, 2.5,  3, 3.5,  4))+
    scale_y_continuous(breaks = seq(1, 4, by = 0.5), labels = c(1, 1.5, 2, 2.5,  3, 3.5,  4))
  grid.arrange(p1,p2)
  
  
  p <- grid.arrange(p1, p2)
  p
  
  Output_name <- gsub(";","_",condition.oi)
  Output_name <- gsub("=","-", Output_name)
  ggsave(p, filename = paste0("figures/",Output_name, ".pdf"), height = 8, width = 6)
  
  
})
