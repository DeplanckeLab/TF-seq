# title: "Modeling overall transcriptomic changes in function of TF dose"
# author: "Wangjie Liu"
# date: "2024/10/23"

setwd("./")

library(Seurat)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(dplyr)

# function
mod.nls.nostart <- function(data, TFoi){
  mod.nls <- nls(Overall_transcriptomic_change ~ SSlogis(Dose,  Asym, xmid, scal), 
                 data = subset(df, TF == TFoi), trace = T)
  print(coef(mod.nls))
  print(summary(mod.nls))
  return(mod.nls)
}

mod.nls.start.range <- function(data, TFoi, init, upper.value, lower.value){
  mod.nls <- nls(Overall_transcriptomic_change ~ SSlogis(Dose,  Asym, xmid, scal), 
                 data = subset(df, TF == TFoi),
                 start = init,
                 upper = upper.value,
                 lower = lower.value,
                 algorithm = "port", trace = T)
  print(coef(mod.nls))
  print(summary(mod.nls))
  return(mod.nls)
}


plot.TF.activities <- function(tfoi, data, col.TFs = NULL, model = NULL, alpha_given = 0.5){
  df.sub <- data[data$TF %in% tfoi,]
  if (length(tfoi) == 1){
    p <- ggplot(df.sub, mapping = aes(x = Dose, y = Overall_transcriptomic_change, col = TF)) + 
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(size =1, alpha = alpha_given)+
      theme_cowplot()
    if (is.null(model) == F){
      p <- p+geom_line(data = df.sub, mapping = aes(x = Dose, y = predict(model, data.frame(x = Dose)), col = TF), size =1)
    }
  } else {
    p <- ggplot(df.sub, mapping = aes(x = Dose, y = Overall_transcriptomic_change, col = TF)) + 
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(size =1, alpha = alpha_given)+
      theme_cowplot()+ggtitle(label = "TF regulatory activities")+
      scale_color_manual(values = col.TFs) # +
    if (is.null(model) == F){
      p <- p+geom_line(data = df.sub[df.sub$TF == tfoi[1],], mapping = aes(x = Dose, y = predict(model, data.frame(x = Dose))), size =1 )
    }
  }
  return(p)
}




# load data
seu <- readRDS("results/C3H10_10X_all_exps_D0regressed_integrated_dosealigned.rds") 
df <- readRDS("output_version3_2024_June/7-TF-potency-capacity/df_overall_transcriptomic_changes.rds")

# add metadata
seu <- subset(seu, subset = Phase_corrected == "G1")
seu$Overall_transcriptomic_changes <- df[colnames(seu), "Overall_transcriptomic_changes"]

df$Dose <- seu$Dose_aligned[rownames(df)]
df$TF[df$TF %in% c("D0","D0_confluent")] <- "D0"
df$Dose_unaligned[df$TF %in% c("D0")] <- 0
df$cell_barcode <- rownames(df)
df$TF <- as.character(df$TF)
df$batch <- seu$batch_overall[rownames(df)]

# exclude controls, references, and TF combinations
df <- subset(df, !TF %in% c("Mycn-Runx2","Mycn-Myog","Pparg-Runx2",
                            "Cebpa-Pparg","Cebpa-Myog","Cebpa-Mycn",
                            "Mycn-Pparg",
                            "Myo_ref",
                            "Adipo_ref",
                            "D0"
))
all_TFs <- names(table(df$TF))

##-------------------------------------quality control before modeling
df <- df %>% group_by(TF) %>% mutate(nCell = length(TF))
# 1) filter by number of G1 cells
df <- df %>% filter(nCell >= 30)
TFs_out.nCell <- all_TFs[!all_TFs %in% unique(df$TF)] 

# 2) control dose range 
df %>% ungroup() %>% summarise(min = min(Dose), max = max(Dose))

max.dose <- max(df$Dose)
df <- df %>% ungroup() %>% 
  mutate(dose_bin = cut(Dose, 
                        breaks = c(-Inf, log1p(1), seq(log1p(1), max(Dose), length.out = round(max.dose)/1)[-1]),
                        include.lowest = T,
                        right = F))
df <- df %>% group_by(TF) %>% filter(min(Dose) < 1.68) 

df.nCell.bin <- df %>%
  group_by(TF, dose_bin) %>%
  dplyr::count(TF)
length(table(df.nCell.bin$TF)) 


df.lowdose <- df.nCell.bin[df.nCell.bin$dose_bin %in% c("[0.693,1.68)"),] 
TFs_out.nCell_lowdose <- c((df.lowdose[df.lowdose$n < 3, ]$TF)) # exclude TFs that have less than 3 cells at the first bin (dose < 1.68)
TFs_out.nCell_lowdose 
df <- df[!df$TF %in% TFs_out.nCell_lowdose,]

df <- df %>% group_by(TF) %>% filter(max(Dose) > 3.5)  # exclude TFs that do not have an overall high dose


##------------------------------------- run logistic model
df <- df[order(df$Dose),]
TFois <- unique(df$TF)
model_list <- list()
plot_list <- list()
NotApplic <- c()

for (i in seq_along(TFois)){ 
  print(TFois[i])
  df.tfoi <- df[df$TF == TFois[i],]
  Init <- list(Asym = median(df.tfoi$Overall_transcriptomic_change), xmid = median(df.tfoi$Dose), scal = 0.1)
  Upper.value <- list(Asym = max(df.tfoi$Overall_transcriptomic_change), xmid = max(df.tfoi$Dose), scal = 10)
  Lower.value <- list(Asym = min(df.tfoi$Overall_transcriptomic_change), xmid = 0, scal = 0.001)
  model_i <- try(mod.nls.start.range(df.tfoi, TFois[i], Init, Upper.value, Lower.value))
  # check if the model was successful
  if(inherits(model_i, "try-error")){
    model_list[[i]] <- FALSE
    plot_list[[i]] <- NULL
    NotApplic <- c(NotApplic, TFois[i])
  } else {
    model_list[[i]] <- model_i
    
    # generate predicted values
    pred_data <- data.frame(Dose = df[df$TF == TFois[i],]$Dose)
    pred_data$Overall_transcriptomic_change <- predict(model_i, newdata = pred_data)
    
    plot_i <- plot.TF.activities(TFois[i], df, model = model_i)
    print(plot_i)
    plot_list[[i]] <- plot_i
  }
}
names(model_list) <- TFois
plot_list <- plot_list[!sapply(plot_list, is.null)] 
TFs.notapplic  <- model_list[sapply(model_list, isFALSE)] %>% names()


##------------------------------------- run logistic model for TFs.notapplicable with new init scal -> 1
model_list_2 <- list()
plot_list_2 <- list()
NotApplic_2 <- c()
for (i in seq_along(TFs.notapplic)){
  print(TFs.notapplic[i])
  df.tfoi <- df[df$TF == TFs.notapplic[i],]
  Init <- list(Asym = median(df.tfoi$Overall_transcriptomic_change), xmid = median(df.tfoi$Dose), scal = 1)
  Upper.value <- list(Asym = max(df.tfoi$Overall_transcriptomic_change), xmid = max(df.tfoi$Dose), scal = 10)
  Lower.value <- list(Asym = min(df.tfoi$Overall_transcriptomic_change), xmid = 0, scal = 0.001)
  model_i <- try(mod.nls.start.range(df.tfoi, TFs.notapplic[i], Init, Upper.value, Lower.value))
  # check if the model was successful
  if(inherits(model_i, "try-error")){
    model_list_2[[i]] <- FALSE
    plot_list_2[[i]] <- NULL
    NotApplic_2 <- c(NotApplic_2, TFs.notapplic[i])
  } else {
    model_list_2[[i]] <- model_i
    
    # generate predicted values
    pred_data <- data.frame(Dose = df[df$TF == TFs.notapplic[i],]$Dose)
    pred_data$Overall_transcriptomic_change <- predict(model_i, newdata = pred_data)
    
    plot_i <- plot.TF.activities(TFs.notapplic[i], df, model = model_i)
    print(plot_i)
    plot_list_2[[i]] <- plot_i
  }
}
names(model_list_2) <- TFs.notapplic
plot_list_2 <- plot_list_2[!sapply(plot_list_2, is.null)] 
TFs.notapplic_2 <- model_list_2[sapply(model_list_2, isFALSE)] %>% names()


##------------------------------------- Combine all runs and extract model parameters
model_list <- model_list[!sapply(model_list, isFALSE)]
model_list_2 <- model_list_2[!sapply(model_list_2, isFALSE)]
model_list.com <- c(model_list, model_list_2)
names(model_list.com) <- c(names(model_list), names(model_list_2))
saveRDS(model_list.com, file = "results/nls_models.rds") 

nls.Par <- lapply(names(model_list.com), function(x){
  nls.Par.tfoi <- data.frame(round(coef(model_list.com[[x]]), 3))
  colnames(nls.Par.tfoi) <- x
  nls.Par.tfoi <- t(nls.Par.tfoi) %>% as.data.frame()
  nls.Par.tfoi$TF <- x
  nls.Par.tfoi$Dose.max <- max(df$Dose[df$TF == x]) %>% round(digits = 3)
  nls.Par.tfoi$Dose.min <- min(df$Dose[df$TF == x]) %>% round(digits = 3)
  nls.Par.tfoi$startpoint <- fitted(model_list.com[[x]])[1] %>% round(digits = 2)
  nls.Par.tfoi$endpoint <- fitted(model_list.com[[x]]) %>% tail(1) %>% round(digits = 2)
  return(nls.Par.tfoi)
})
names(nls.Par) <- names(model_list.com)
nls.Par <- data.table::rbindlist(nls.Par) %>% as.data.frame()
rownames(nls.Par) <- nls.Par$TF


##------------------------------------- categorizing TFs based on model fit

nls.Par$category <- "-"

# 1. define high-capacity and low-capacity groups
nls.Par$category[(nls.Par$endpoint >= 0.23)] <- "high-capacity"
nls.Par$category[(nls.Par$endpoint < 0.23)] <- "low-capacity"
table(nls.Par$category)

mean.Dose <- round(mean(df$Dose),0) # 3
predict.value_meanDose <- lapply(names(model_list.com),function(x){
  predict.value <- predict(model_list.com[[x]], newdata = data.frame(Dose = mean.Dose)) # Dose ==3
  predict.value <- data.frame(TF = x, predicted_value_meanDose = as.numeric(predict.value))
  return(predict.value)
})
predict.value_meanDose <- data.table::rbindlist(predict.value_meanDose)
nls.Par$predict.value_meanDose <- predict.value_meanDose$predicted_value_meanDose[predict.value_meanDose$TF %in% rownames(nls.Par)]
nls.Par$predict.value_meanDose <- round(nls.Par$predict.value_meanDose, 2)
summary(nls.Par[nls.Par$category %in% c("high-capacity"),]$predict.value_meanDose)
nls.Par$category[nls.Par$category == "high-capacity" &
                   nls.Par$predict.value_meanDose >= 0.23] <- "high-capacity & high-sensitivity"
nls.Par$category[nls.Par$category == "high-capacity" &
                   nls.Par$predict.value_meanDose < 0.23 ] <- "high-capacity & low-sensitivity"

# save the results
write.csv(nls.Par, file = "data/TF_categories_on_potency_capacity_dosealigned.csv", row.names = T) # Supplementary_table_4


##------------------------------------- plotting

# Figure 3a
p <- FeaturePlot(seu, features = c("Overall_transcriptomic_changes"), order = T)+scale_color_viridis_c(option = "B")#+
#theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.text = element_blank())
p


# Figure 3g
p <- ggplot()+
  geom_point(df[df$TF == "Hoxa6",], mapping = aes(x = Dose, y = Overall_transcriptomic_change), color = "#9EC8EA",alpha = 0.5,size = 0.5)+
  geom_point(df[df$TF == "Pou5f1",], mapping = aes(x = Dose, y = Overall_transcriptomic_change), color = "#CA99C5",alpha = 0.5, size = 0.5)+
  geom_point(df[df$TF == "Vdr",],mapping = aes(x = Dose, y = Overall_transcriptomic_change), color = "#98CB67",alpha = 0.5, size = 0.5)+
  geom_line(data = df[df$TF == "Hoxa6",], mapping = aes(x = Dose, y = predict(model_list.com$Hoxa6, data.frame(x = Dose))), size =1, color = "#9EC8EA")+
  geom_line(data = df[df$TF == "Pou5f1",], mapping = aes(x = Dose, y = predict(model_list.com$Pou5f1, data.frame(x = Dose))), size =1, color = "#CA99C5")+
  geom_line(data = df[df$TF == "Vdr",], mapping = aes(x = Dose, y = predict(model_list.com$Vdr, data.frame(x = Dose))), size =1, color = "#98CB67")+
  theme_minimal_grid(line_size = 0.5, color = "gray90")+
  facet_wrap("TF")
p


# Extended_data_figure_9e
P <- ggplot()+
  geom_point(df[df$TF == "Cebpa",], mapping = aes(x = Dose, y = Overall_transcriptomic_change), color = "#377EB8",alpha = 0.5, size = 0.5)+
  geom_point(df[df$TF == "Mycn",], mapping = aes(x = Dose, y = Overall_transcriptomic_change), color = "#FB8072",alpha = 0.5,size = 0.5)+
  geom_point(df[df$TF == "Myog",], mapping = aes(x = Dose, y = Overall_transcriptomic_change), color = "#CAB2D6",alpha = 0.5, size = 0.5)+
  geom_point(df[df$TF == "Runx2",],mapping = aes(x = Dose, y = Overall_transcriptomic_change), color = "#996600",alpha = 0.5, size = 0.5)+
  geom_point(df[df$TF == "Pparg",],mapping = aes(x = Dose, y = Overall_transcriptomic_change), color = "#A6CEE3",alpha = 0.5, size = 0.5)+
  geom_line(data = df[df$TF == "Cebpa",], mapping = aes(x = Dose, y = predict(model_list.com$Cebpa, data.frame(x = Dose))), size =1, color = "#377EB8")+
  geom_line(data = df[df$TF == "Mycn",], mapping = aes(x = Dose, y = predict(model_list.com$Mycn, data.frame(x = Dose))), size =1, color = "#FB8072")+
  geom_line(data = df[df$TF == "Myog",], mapping = aes(x = Dose, y = predict(model_list.com$Myog, data.frame(x = Dose))), size =1, color = "#CAB2D6")+
  geom_line(data = df[df$TF == "Runx2",], mapping = aes(x = Dose, y = predict(model_list.com$Runx2, data.frame(x = Dose))), size =1, color = "#996600")+
  geom_line(data = df[df$TF == "Pparg",], mapping = aes(x = Dose, y = predict(model_list.com$Pparg, data.frame(x = Dose))), size =1, color = "#A6CEE3")+
  theme_minimal_grid()
P


