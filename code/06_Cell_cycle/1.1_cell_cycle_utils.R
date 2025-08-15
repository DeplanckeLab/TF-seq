# functions for cell cycle analysis

library(ks)

mashaGgplot2Theme <- list(
  theme_classic(base_size = 14) + 
    theme(text = element_text(size = 14)) +
    theme(axis.line.x = element_line(colour = 'black', size = 0.5,
                                     linetype = 'solid'),
          axis.line.y = element_line(colour = 'black', size=0.5,
                                     linetype ='solid'),
          panel.grid.minor = element_line(colour = "white", size = 0.5,
                                          linetype = 2))
)

barplot_TopCyclingTF <- function(data, Phase.select){ # data is ALL_func from input
  # quality control: filter out TFs with less than 50 cells 
  TF.filter <- table(data$TF)[table(data$TF) > 50] %>% names()
  data <- data[data$TF %in% TF.filter,]
  df_PhaseperTF <- reshape2::melt(table(data$TF, data[[Phase.select]]))
  df_PhaseperTF$Fraction <- apply(df_PhaseperTF, 1, 
                                  function(x) return(as.numeric(x["value"])/sum(data$TF == x["Var1"])))
  colnames(df_PhaseperTF) <- c("TF", "Phase", "Counts", "Fraction")
  #detach(package:plyr)
  df_PhaseperTF <- as.data.frame(df_PhaseperTF)
  TFoi <- unique(df_PhaseperTF$TF) %>% as.character()
  
  # add number of S+G2M cells as label for each TF
  df_PhaseperTF.sub <- lapply(TFoi, function(x){
    df_PhaseperTF.sub <- subset(df_PhaseperTF, TF == x)
    df_PhaseperTF.sub$Total_nCell <- sum(df_PhaseperTF.sub$Counts)
    S.counts <- df_PhaseperTF.sub$Counts[df_PhaseperTF.sub$Phase == "S"]
    G2M.counts <- df_PhaseperTF.sub$Counts[df_PhaseperTF.sub$Phase == "G2M"]
    df_PhaseperTF.sub$label <- rep((S.counts+G2M.counts), 3)
    df_PhaseperTF.sub$label <- round(df_PhaseperTF.sub$label, 1)
    
    # if label TFs with total number of cells
    df_PhaseperTF.sub$label <- paste0(x, " (", df_PhaseperTF.sub$Total_nCell, ")")
    return(df_PhaseperTF.sub)
  })
  df_PhaseperTF <- data.table::rbindlist(df_PhaseperTF.sub)
  df_PhaseperTF$TF <- as.factor(df_PhaseperTF$TF)
  #Reorder from highest sum of cycling cells 
  cyclingFreq <- subset(df_PhaseperTF, Phase %in% c("S", "G2M"))
  ord <- as.data.frame(cyclingFreq %>% group_by(label) %>% summarise(SumCyclingFreq = sum(Fraction)))
  ord <- ord[order(ord$SumCyclingFreq, decreasing = F),]
  
  df_PhaseperTF$label <- as.character(df_PhaseperTF$label)
  df_PhaseperTF$label <- factor(df_PhaseperTF$label, levels = ord$label) 
  print(df_PhaseperTF)
  p_v <- ggplot(df_PhaseperTF, aes(fill = Phase, x = Fraction, y = label, alpha = Phase)) + 
    geom_bar(position="stack", stat="identity") + 
    scale_fill_manual(values = c("gray", "#a6cee3", "#1f78b4")) + 
    scale_alpha_manual(values = c(0.2,1,1)) +
    mashaGgplot2Theme + 
    theme(
      axis.text.y = element_text(size = 9),
      panel.grid.major.x = element_line(size = 0.1, color = "grey"),
      legend.position = "top"
    ) 
  
  
  df_PhaseperTF$label <- as.character(df_PhaseperTF$label)
  df_PhaseperTF$label <- factor(df_PhaseperTF$label, levels = rev(ord$label))
  p_h <- ggplot(df_PhaseperTF, aes(fill = Phase, x = label, y = Fraction, alpha = Phase)) + 
    geom_bar(position="stack", stat="identity") + 
    scale_fill_manual(values = c("gray", "#a6cee3", "#1f78b4")) + 
    scale_alpha_manual(values = c(0.2,1,1)) +
    mashaGgplot2Theme + 
    theme(
      axis.text.x = element_text(size = 9, hjust = 1, angle = 90),
      panel.grid.major.y = element_line(size = 0.5, color = "grey"),
      legend.position = "top"
    ) #+
  #scale_y_discrete(guide = guide_axis(n.dodge = 4))
  P <-  list(vertical = p_v, 
             horizontal = p_h, df_PhaseperTF)
  names(P) <- c("p_v","p_h","df_PhaseperTF")
  
  return(P)
}


## Fisher test or Chi-square test on number of cycling cells (TF vs D0)
# label whether significant or not (TF vs D0)
barplot_TopCyclingTF_step2 <- function(df, q_value_threshold, plot_only_significant){ 
  TFoi <- as.character(unique(df$TF))
  TFoi <- TFoi[TFoi != "D0_confluent"]
  
  P_values <- lapply(TFoi, function(x){
    df.TFoi <- subset(df, TF %in% c(x,"D0_confluent"))
    df.TFoi <- df.TFoi[, c(1:3)]
    df.TFoi$Counts <- as.numeric(df.TFoi$Counts)
    
    # for fisher.test
    df.test <- data.frame(TF = c("D0_confluent", x),
                          G1 = c(df.TFoi$Counts[df.TFoi$TF == "D0_confluent" & df.TFoi$Phase == "G1"], df.TFoi$Counts[df.TFoi$TF == x & df.TFoi$Phase == "G1"]),
                          G2MS =  c(sum(df.TFoi$Counts[df.TFoi$TF == "D0_confluent" & df.TFoi$Phase %in% c("G2M", "S")]), sum(df.TFoi$Counts[df.TFoi$TF == x & df.TFoi$Phase %in% c("G2M","S")]))
    )
    
    
    rownames(df.test) <- df.test$TF
    df.test$TF <-NULL
    df.test <- as.matrix(df.test)
    print(fisher.test(df.test))
    
    # combine S and G2M as one cell cycle group, do fisher.test (2x2contingency table)
    p.value <- data.frame(TF = x, p.value = fisher.test(df.test)$p.value, odds_ratio = fisher.test(df.test)$estimate["odds ratio"], Side = fisher.test(df.test)$alternative)
    rownames(p.value) <- x
    return(p.value)
  })
  names(P_values) <- TFoi
  P_values <- data.table::rbindlist(P_values)
  P_values$q.value <- p.adjust(P_values$p.value, method = "fdr") 
  P_values$significant <- F
  P_values$significant[P_values$q.value < q_value_threshold] <- T
  P_values$significant <- as.factor(P_values$significant)
  # add significance information in data frame for plotting
  df$is.significant <- F
  df.new <- lapply(TFoi, function(x){
    df.TFoi <- subset(df, TF == x)
    df.TFoi$is.significant <- rep(P_values$significant[P_values$TF == x], 3)
    return(df.TFoi)
  })
  names(df.new) <- TFoi
  df.new <- data.table::rbindlist(df.new) %>% as.data.frame()
  df.new <- rbind(df.new, df[df$TF == "D0_confluent",])
  df.new$TF <- as.factor(df.new$TF)
  df.new$is.significant <- as.factor(df.new$is.significant)
  
  #Reorder from highest sum of cycling cells 
  cyclingFreq <- subset(df.new, Phase %in% c("S", "G2M"))
  ord <- as.data.frame(cyclingFreq %>% group_by(label) %>% summarise(SumCyclingFreq = sum(Fraction)))
  ord <- ord[order(ord$SumCyclingFreq, decreasing = T),]
  
  df.new$label <- as.character(df.new$label)
  df.new$label <- factor(df.new$label, levels = ord$label) 
  df.new$is.significant <- factor(df.new$is.significant)
  cols <- ifelse(df.new$is.significant == T, "darkred","black")
  names(cols) <- df.new$label
  
  
  rownames(ord) <- ord$label
  ord$col <- cols[rownames(ord)] %>% as.factor()
  cols <- ord$col
  
  if (plot_only_significant == T){
    df.new <- df.new[df.new$is.significant == T | df.new$TF == "D0_confluent" , ]
    p_h <- ggplot(df.new, aes(x = label, y = Fraction)) + 
      geom_bar(mapping = aes(fill = Phase, alpha = Phase), position="stack", stat="identity") + 
      scale_fill_manual(values = c("lightgrey", "#a6cee3", "#1f78b4")) + 
      scale_alpha_manual(values = c(0.2,1,1)) +
      mashaGgplot2Theme + 
      theme(
        axis.text.x = element_text(size = 9, angle = 45, colour = "black", vjust = 0.9, hjust = 0.9),
        panel.grid.major.y = element_line(size = 0.1, color = "grey"),
        legend.position = "left")
  } else {
    p_h <- ggplot(df.new, aes(x = label, y = Fraction)) + 
      geom_bar(mapping = aes(fill = Phase, alpha = Phase), position="stack", stat="identity") + 
      scale_fill_manual(values = c("lightgrey", "#a6cee3", "#1f78b4")) + 
      scale_alpha_manual(values = c(0.2,1,1)) +
      mashaGgplot2Theme + 
      theme(
        axis.text.x = element_text(size = 9, colour = cols, angle = 45, vjust = 0.9,  hjust = 0.9 ),
        panel.grid.major.y = element_line(size = 0.1, color = "grey"),
        legend.position = "left") 
    
  }
  
  P <- list(p_h, P_values, df.new)
  names(P) <- c("p_h","P_values","df_plot")
  return(P)
}



# # plot the G2M score vs S score of all cells (D0, TFs of interests)
# #plot G2M score or S.Score vs Vector expression
# # plot phase density of TF cells and D0 cells for each cell cycle score 

Plot_1D_density <- function(TFoi, data, OutputDir, D0_to_use = "D0_confluent"){
  plotdata <- data[data$TF %in% c(D0_to_use,TFoi),]
  plotdata$TF <- factor(plotdata$TF, levels = c(TFoi, D0_to_use))
  #S score for all phase cells
  p1 <- ggplot(plotdata)+
    geom_density(aes(x = S.Score, col = TF), position = "identity")+theme_bw()+
    labs(x = "S.Score_all", y = "density")
  #G2M score for all phase cells
  p2 <- ggplot(plotdata )+
    geom_density(aes(x = G2M.Score, col = TF), position = "identity")+theme_bw()+
    labs(x = "G2M.Score_all", y = "density")
  p_1D_density <- gridExtra::grid.arrange(p1,p2, ncol = 2, nrow = 1)
  ggsave(p_1D_density,
         width = 8, height = 4,
         units = "in",
         filename = paste0(OutputDir,"Density_PhaseScore_all_cells_",D0_to_use,"-",TFoi,".pdf"))
}


## perform statistical test on cell cycle scores of two groups
Pvalue.WILCOX_TEST <- function(TFoi, data, D0_to_use = "D0_confluent"){
  X1.S <- data[data$TF == TFoi, c("S.Score")] 
  X2.S <- data[data$TF == D0_to_use, c("S.Score")] 
  pvalue.S <- wilcox.test(X1.S, X2.S)$p.value
  X1.G2M <- data[data$TF == TFoi, c("G2M.Score")] 
  X2.G2M <- data[data$TF == D0_to_use, c("G2M.Score")] 
  pvalue.G2M <- wilcox.test(X1.G2M, X2.G2M)$p.value
  df.pvalues <- data.frame(TF = TFoi, pvalue.wilcox_S = pvalue.S, pvalue.wilcox_G2M = pvalue.G2M)
  return(df.pvalues)
}

## perform 2-D kde test to compare kernal density of two groups
Pvalue.KDE_TEST <- function(TFoi, data, D0_to_use = "D0_confluent"){
  X1 <- data[data$TF == TFoi, c("S.Score","G2M.Score")] %>% as.matrix()
  X2 <- data[data$TF == D0_to_use, c("S.Score","G2M.Score")] %>% as.matrix()
  pvalue.kde <- kde.test(x1 = X1, x2 = X2)$pvalue
  # The null hypothesis is H0 : f1 ≡ f2 where f1, f2 are the respective density functions
  # This test statistic has a null distribution which is asymptotically normal
  return(pvalue.kde)
}


# plot the 2-D density contour plot 
Plot_2D_contour <- function(TFoi, data, D0_to_use = "D0_confluent", OutputDir){
  data.TFoi <- data[data$TF == TFoi, ] 
  data.ctr <- data[data$TF == D0_to_use,] 
  
  p_ctr_density <- ggplot(data.ctr, mapping = aes(x = S.Score, y = G2M.Score))+
    geom_point(alpha = 0.4, color = "lightblue", size = 2)+
    geom_hdr_lines(aes(colour = after_stat(probs)), method = "kde") +
    scale_fill_viridis_c()+
    geom_abline(intercept = 0, slope = 1, "darkgray")+
    geom_vline(xintercept = 0, color = "darkgray" )+
    geom_hline(yintercept = 0, color = "darkgray")+
    geom_vline(xintercept = 0.1 , color = "lightgray")+
    geom_hline(yintercept = 0.1, color = "lightgray")+
    labs(subtitle = D0_to_use)+
    cowplot::theme_cowplot()
  # p_ctr_density
  
  p_tf_density <- ggplot(data.TFoi,  aes(x = S.Score, y = G2M.Score))+
    geom_point(alpha = 0.4, color = "lightblue", size = 2)+
    geom_hdr_lines(aes(colour = after_stat(probs)), method = "kde", show.legend = T) +
    scale_fill_viridis_c()+
    geom_abline(intercept = 0, slope = 1, color = "darkgray")+
    geom_vline(xintercept = 0, color = "darkgray" )+
    geom_hline(yintercept = 0, color = "darkgray")+
    geom_vline(xintercept = 0.1 , color = "lightgray")+
    geom_hline(yintercept = 0.1, color = "lightgray")+
    labs(subtitle = TFoi)+
    cowplot::theme_cowplot()
  # p_tf_density
  
  P.com <- p_ctr_density+p_tf_density
  ggsave(P.com,
         width = 14, height = 6,
         filename = paste0(OutputDir,"Density_contour_D0con-",TFoi,".pdf"))
  return(P.com)
}


# # scatter plots of cell cycle scores (colored by vector expression or Phase group) for TFs of interest
Plot_Scatter_CellCycleScore <- function(TFoi, data, plot_vector = T, plot_phase = F, Phase.use = "Phase", OutputDir){
  plotdata <- data[data$TF == TFoi,]
  plotdata <- plotdata %>% arrange(dose)
  if(plot_vector == T & TFoi %nin% c("D0","D0_confluent")){
    p_vector <- ggplot(data = plotdata) +
      geom_point(mapping = aes(x = S.Score, y = G2M.Score, color = dose), position = 'jitter', size = 1) +
      scale_color_viridis_c()+
      geom_abline(intercept = 0, slope = 1, color = "darkgray")+
      geom_vline(xintercept = 0, color = "darkgray" )+
      geom_hline(yintercept = 0, color = "darkgray")+
      geom_vline(xintercept = 0.1 , color = "lightgray")+
      geom_hline(yintercept = 0.1, color = "lightgray")+
      labs(subtitle = TFoi)+
      cowplot::theme_cowplot()
    p_vector
    ggsave(p_vector, 
           width = 7, height = 6,
           filename = paste0(OutputDir,"Scatterplot_",TFoi,"-vector.pdf"))
  }
  
  if(plot_phase == T & Phase.use == "Phase"){
    p_phase <- ggplot(data = plotdata) +
      geom_point(mapping = aes(x = S.Score, y = G2M.Score, col = Phase), position = 'jitter', size = 1) +
      scale_color_manual(values = c("lightgrey", "#a6cee3", "#1f78b4"))+
      geom_abline(intercept = 0, slope = 1, color = "darkgray")+
      geom_vline(xintercept = 0, color = "darkgray" )+
      geom_hline(yintercept = 0, color = "darkgray")+
      geom_vline(xintercept = 0.1 , color = "lightgray")+
      geom_hline(yintercept = 0.1, color = "lightgray")+
      labs(subtitle = TFoi)+
      cowplot::theme_cowplot()
    p_phase
    ggsave(p_phase, 
           width = 7, height = 6,
           filename = paste0(OutputDir, "Scatterplot_",TFoi,"-",Phase.use,".pdf"))
  } else if(plot_phase == T & Phase.use == "Phase_corrected"){
    p_phase <- ggplot(data = plotdata) +
      geom_point(mapping = aes(x = S.Score, y = G2M.Score, col = Phase_corrected), position = 'jitter', size = 1) +
      scale_color_manual(values = c("lightgrey", "#a6cee3", "#1f78b4"))+
      geom_abline(intercept = 0, slope = 1, color = "darkgray")+
      geom_vline(xintercept = 0, color = "darkgray" )+
      geom_hline(yintercept = 0, color = "darkgray")+
      geom_vline(xintercept = 0.1 , color = "lightgray")+
      geom_hline(yintercept = 0.1, color = "lightgray")+
      labs(subtitle = TFoi)+
      cowplot::theme_cowplot()
    p_phase
    ggsave(p_phase, 
           width = 7, height = 6,
           filename = paste0(OutputDir, "Scatterplot_",TFoi,"-",Phase.use,".pdf"))
  }
}





