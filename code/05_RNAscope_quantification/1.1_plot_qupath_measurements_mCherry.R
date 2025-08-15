# title: "Downstream QC and plotting on the single-cell fluorescence quantification of RNAscope images - mCherry"
# input data: "RNAscope_measurements_Fig1_EDFig1.tsv"
# author: "Wangjie Liu"


# Non-confocoal image data of mCherry (rep - 20241021)
# Best Z

library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(gridExtra)

setwd("./")

measurements <- read_delim("data/RNAscope_measurements_Fig1_EDFig1.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>% as.data.frame()
rownames(measurements) <- measurements$`Object ID`

measurements$analysis <- "Fluorescence"
measurements$analysis[grep("DPC", measurements$Image)] <- "DPC"
table(measurements$analysis)
# DPC Fluorescence 
# 2569         2412 

df <- measurements[measurements$analysis == "Fluorescence",] 

ggplot(df)+geom_density(mapping = aes(x =`Nucleus: Area µm^2`))+geom_vline(xintercept = 600)
# need to decide empirically based on the data
ggplot(df)+geom_density(mapping = aes(x =`Channel_1: Nucleus: Mean`)) + geom_vline(xintercept = 2000)# DAPI signal 
ggplot(df)+geom_density(mapping = aes(x =`Channel_1: Nucleus: Median`)) # DAPI signal
## this batch have some cells  have scattered and weak DAPI signal, exclude by adding cutoff for `Channel_1: Nucleus: Max`
ggplot(df)+geom_density(mapping = aes(x =`Channel_1: Nucleus: Max`)) + geom_vline(xintercept = 12000) # DAPI signal
ggplot(df)+geom_density(mapping = aes(x =`Channel_1: Cytoplasm: Mean`)) # DAPI signal
ggplot(df)+geom_density(mapping = aes(x = log10(`Channel_1: Cell: IntDen Corrected`)))  #  DAPI signal
ggplot(df)+geom_density(mapping = aes(x = `Nucleus/Cell area ratio`)) # high Nucleus/Cell area ratio is a indicative of cell clumps or multinucleated cells




# filter out outlier with very small cell size  
measurements <- measurements[measurements$`Nucleus: Area µm^2` >= 600 ,] # size

# filter out outlier with weak DAPI signals (wrong cell segmentation) 
measurements <- measurements[measurements$`Channel_1: Nucleus: Mean` >= 2000 & measurements$`Channel_1: Nucleus: Max` >=12000,] # DAPI signals

# filter out outlier with very strong DAPI signals (wrong cell segmentation) 
measurements <- measurements[measurements$`Channel_1: Nucleus: Mean` <= 15000 & measurements$`Channel_1: Nucleus: Max` <=45000,] # DAPI signals

# exclude cells without intact nuclei indicated by low sum intensity of DAPI
measurements <- measurements[measurements$`Channel_1: Cell: IntDen Corrected` >= 10^6.5, ] # sum of background corrected DAPI signals

measurements <- measurements[measurements$`Nucleus/Cell area ratio` >= 0.3 & measurements$`Nucleus/Cell area ratio` <= 0.75 & measurements$`Channel_1: Cytoplasm: Mean` <= 4000, ]

table(measurements$analysis)


# subset to Cell Mean Intensity; Cell Mean intensity with background corrected; Cell integrated intensity with background corrected
measurements2 <- measurements[measurements$analysis != "DPC",] 
measurements2 <- measurements2[, c(1, 2, grep("Cell: Mean",colnames(measurements2)), grep("Cell: IntDen Corrected",colnames(measurements2)),  grep("Cell: ROI Area",colnames(measurements2)))]
measurements2$Well <- substring(measurements2$Image, 29, 35)
table(measurements2$Well)
# R02-C02 R02-C03 R02-C04 R02-C05 R03-C02 R03-C03 R03-C04 R03-C05 R04-C02 R04-C03 R04-C04 R04-C05 
# 275     211     122     168     223     270     111     184     207     379     137     125 


measurements2$condition <- "Pos_neg_controls"
measurements2$condition[measurements2$Well %in% c(paste0("R02-C0",2:4))] <- "Ch2=mCherry; Ch3=WPRE"
measurements2$condition[measurements2$Well %in% c(paste0("R03-C0",2:4))] <- "Ch2=mCherry; Ch3=WPRE"
measurements2$condition[measurements2$Well %in% c(paste0("R04-C0",2:4))] <- "Ch2=mCherry; Ch3=WPRE"

table(measurements2$condition)
# Ch2=mCherry; Ch3=WPRE      Pos_neg_controls 
# 1935                   477 

# replicates:
measurements2$replicate <- "Others"
measurements2$replicate[measurements2$Well %in% c(paste0("R02-C0",2:4))] <- "Rep1"
measurements2$replicate[measurements2$Well %in% c(paste0("R03-C0",2:4))] <- "Rep2"
measurements2$replicate[measurements2$Well %in% c(paste0("R04-C0",2:4))] <- "Rep3"
table(measurements2$replicate)
# Others   Rep1   Rep2   Rep3 
# 477    608    604    723 

measurements2$cell_group <- "WT-C3H10 cells"
measurements2$cell_group[measurements2$Well %in% c(paste0("R0",2:4,"-C02"))] <- "mCherry cells"
measurements2$cell_group[measurements2$Well %in% c(paste0("R0",2:4,"-C03"))] <- "mCherry-noDox cells"
measurements2$cell_group[measurements2$Well %in% c(paste0("R0",2:4,"-C04"))] <- "WT-C3H10 cells"
table(measurements2$cell_group)
# mCherry-noDox cells       mCherry cells      WT-C3H10 cells 
# 860                 705                 847 

measurements2$cell_group <- factor(measurements2$cell_group, levels = c("mCherry cells", "mCherry-noDox cells", "WT-C3H10 cells"))



condition.oi <- "Ch2=mCherry; Ch3=WPRE"
plot_data <- measurements2[measurements2$condition == condition.oi,]
intensity.oi <- "Mean BgCorrected"
colnames(plot_data) <- sub(paste0(": Cell: ", intensity.oi),"",colnames(plot_data))
colnames(plot_data)

# Channel_0 = DPC
# Channel_1 = DAPI
# Channel_2 = mCherry
# Channel_3 = Vivid650
# Channel_4 = BF


plot_data <- as.data.frame(plot_data)



  
#  plot_data <- plot_data[plot_data$replicate == "Rep1",]
  # for the first look
  ggplot() +
    geom_point(plot_data, mapping = aes(x=(`Channel_3`), y=(`Channel_2`), col=replicate), alpha = 0.5,) +
    facet_grid(~cell_group)+
    geom_smooth(plot_data[plot_data$cell_group != "WT-C3H10 cells",], mapping = aes(x=(`Channel_3`), y=(`Channel_2`)))+
    facet_grid(~cell_group)+
    labs(title = plot_data$condition)
  
  
  # for each well and each channel, remove 1% of the upper/lower intensity (extreme values/outliers)
  cells.out <- c()
  for (Ch.oi in paste0("Channel_", 2:3)){
    for (Well.oi in unique(plot_data$Well)){
      df <- plot_data[plot_data$Well == Well.oi, ]
      cells.out <- c(cells.out, df$`Object ID`[df[, Ch.oi] >= quantile(df[, Ch.oi], 0.99)], df$`Object ID`[df[, Ch.oi] <= quantile(df[, Ch.oi], 0.01)]) 
    }
  }
  
  plot_data <- plot_data[!plot_data$`Object ID` %in% cells.out, ]
  
  
  ggplot() +
    geom_point(plot_data, mapping = aes(x=(`Channel_3`), y=(`Channel_2`), col=replicate), alpha = 0.5,) +
    facet_grid(~cell_group)+
    geom_smooth(plot_data[plot_data$cell_group != "WT-C3H10 cells",], mapping = aes(x=(`Channel_3`), y=(`Channel_2`)))+
    facet_grid(~cell_group)+
    labs(title = plot_data$condition)
  
  ## since there is intercorrelation between mCherry and WPRE even in mCherry-noDox cells because of mCherry leaky expression,
  ## mCherry-noDox cannot be used for contamination correction, instead, use WT-C3H10 cells to correct spillover
  
  model_Ch2 <- lm(data = plot_data[plot_data$cell_group %in% c("WT-C3H10 cells"),], formula = Channel_2 ~ Channel_3)
  model_Ch2
  if(model_Ch2$coefficients[2] > 0.01){
    plot_data$Channel_2 <- plot_data$Channel_2 - predict(model_Ch2, newdata = plot_data)
    p1 <- ggplot() +
      geom_point(plot_data, mapping = aes(x=(`Channel_3`), y=(`Channel_2`), col=replicate), alpha = 0.5,) +
      facet_grid(~cell_group)+
      geom_smooth(plot_data[plot_data$cell_group != "WT-C3H10 cells",], mapping = aes(x=(`Channel_3`), y=(`Channel_2`)))+
      facet_grid(~cell_group)+
      labs(title = plot_data$condition)
  }
  p1
  

  
  ## log transformation after contamination correction
  plot_data.new <- plot_data
  base.oi <- 10
  plot_data.new[plot_data.new <= base.oi] <- base.oi 
  for (i in paste0("Channel_", 2:3)){
    plot_data.new[, i] <- log(plot_data.new[, i], base = base.oi) 
  }
  
  
  ggplot() +
    geom_point(plot_data.new, mapping = aes(x=(`Channel_3`), y=(`Channel_2`), col=replicate), alpha = 0.5,) +
    facet_grid(~cell_group)+
    geom_smooth(plot_data.new[plot_data.new$cell_group != "WT-C3H10 cells",], mapping = aes(x=(`Channel_3`), y=(`Channel_2`)))+
    facet_grid(~cell_group)+
    labs(title = plot_data.new$condition)
  
  
  # if plot WT and TF cells together
  plot_data.sub <- plot_data.new[plot_data.new$cell_group != "mCherry-noDox cells",]
  plot_data.sub <- plot_data.sub[plot_data.sub$Channel_3 >= 1.75, ] # exclude remaining outliers after checking the images
  plot_data.sub$cell_group |> table() |> print()
  cell_group <- names(table(plot_data.sub$cell_group))
  cell_count <- table(plot_data.sub$cell_group)
  p1 <- ggplot() +
    geom_point(plot_data.sub, mapping = aes(x=(`Channel_3`), y=(`Channel_2`), col=cell_group), alpha = 0.8) +
    scale_color_manual(values = c("#F3766D","#D3D3D3"))+
    geom_smooth(plot_data.sub, mapping = aes(x=(`Channel_3`), y=(`Channel_2`)), span = 0.8,)+
    labs(title = plot_data.sub$condition, subtitle = paste0(cell_group[1], "=", cell_count[1], "cells;", cell_group[3], "=", cell_count[3]))+
    theme_cowplot()+
    scale_x_continuous(breaks = seq(1, 4, by = 0.5), labels = c(1, 1.5, 2, 2.5,  3, 3.5,  4))+
    scale_y_continuous(breaks = seq(1, 4, by = 0.5), labels = c(1, 1.5, 2, 2.5,  3, 3.5,  4))
  p1

  Output_name <- gsub(";","_",condition.oi)
  Output_name <- gsub("=","-", Output_name)
  ggsave(p1, filename = paste0("figures/",Output_name, ".pdf"), height = 4, width = 6)