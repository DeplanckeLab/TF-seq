library(Seurat)
library(tidyverse)
library(edgeR)
library(ggplot2)
library(viridis)

setwd("~/NAS2/TF-seq/Wangjie/TF_resource_paper/")
source("code/12-combinations/functions-diffexp.R")

combinations
# combination_ix <- 5 # Mycn-Myog
combination_ix <- 3 # Cebpa-Mycn
# combination_ix <- 2 #Cebpa-Pparg
# combination_ix <- 6
combination <- combinations[[combination_ix]]
combination
for (combination_ix in c(2, 3, 5, 6)) {
  combination <- combinations[[combination_ix]]
  linkage <- linkages[[combination_ix]]
  seu_diffexp <- seu_diffexps[[combination_ix]]
  enrichment <- enrichments[[combination_ix]]
  condition_oi <- paste0(combination[[1]], "-", combination[[2]])
  
  cells_oi <- linkage %>% 
    filter(condition_from == condition_oi) %>% 
    group_by(cell_from) %>% 
    summarize(percentage = mean(condition_to == condition_oi)) %>% 
    filter(percentage >= 0.5) %>% 
    pull(cell_from)
  length(cells_oi)
  
  cells_ref <- map(combination, function(tf) {
    dosage_oi <- log1p(seu_diffexp@assays$oeumi@counts[tf,cells_oi])
    cells_ref_all <- linkage %>% filter(condition_from == tf) %>% pull(cell_from) %>% unique()
    dosage_ref_all <- log1p(seu_diffexp@assays$oeumi@counts[tf,cells_ref_all])
    library(FNN)
    knn <- knnx.index(dosage_ref_all, dosage_oi, k = 10)
    cells_ref <- unique(names(dosage_ref_all)[knn])
  }) %>% unlist() %>% unique()
  seu_diffexp$oi <- seu_diffexp$cell %in% cells_oi
  seu_diffexp$ref <- seu_diffexp$cell %in% cells_ref
  DimPlot(seu_diffexp, group.by = c("condition", "ref", "oi"))
  
  markers_ <- map(combination, function(tf) {
    dosage_oi <- log1p(seu_diffexp@assays$oeumi@counts[tf,cells_oi])
    cells_ref_all <- linkage %>% filter(condition_from == tf) %>% pull(cell_from) %>% unique()
    dosage_ref_all <- log1p(seu_diffexp@assays$oeumi@counts[tf,cells_ref_all])
    library(FNN)
    knn <- knnx.index(dosage_ref_all, dosage_oi, k = 10)
    cells_ref <- unique(names(dosage_ref_all)[knn])
    
    markers <- FindMarkers(seu_diffexp, ident.1 = cells_oi, ident.2 = cells_ref, logfc.threshold = 0.0)
    markers %>% mutate(tf = tf) %>% tibble::rownames_to_column("gene")
  })
  
  markers <- full_join(markers_[[1]], markers_[[2]], by = "gene", suffix = c("_tf1", "_tf2"))
  markers <- markers %>% filter(!is.na(avg_log2FC_tf1) & !is.na(avg_log2FC_tf2))
  
  markers$significant_1 <- (markers$p_val_adj_tf1 < 5e-2) & (abs(markers$avg_log2FC_tf1) > log(1.5))
  markers$significant_2 <- (markers$p_val_adj_tf2 < 5e-2) & (abs(markers$avg_log2FC_tf2) > log(1.5))
  gsub_combination <- function(text, combination) {
    text = gsub("TF1", toupper(combination[[1]]), text)
    text = gsub("TF2", toupper(combination[[2]]), text)
  }
  markers$group_ix <- ifelse(
    markers$significant_1 & markers$significant_2,
    4,
    ifelse(
      markers$significant_1,
      2,
      ifelse(
        markers$significant_2,
        3,
        1
      )
    )
  )
  groups <- c(
    "not differential",
    "TF1+TF2 vs TF1" %>% gsub_combination(combination),
    "TF1+TF2 vs TF2"  %>% gsub_combination(combination),
    "TF1+TF2 vs TF1 and TF1+TF2 vs TF2" %>% gsub_combination(combination)
  )
  markers$group <- factor(groups[markers$group_ix], levels = groups)
  
  lim <- 3
  plotdata <- markers %>% mutate(x = pmax(pmin(avg_log2FC_tf1, lim), -lim), y = pmax(pmin(avg_log2FC_tf2, lim), -lim))
  plotdata$symbol <- data.annot %>% dplyr::slice(match(plotdata$gene, gene)) %>% pull(gene_short_name)
  plotdata_oi <- bind_rows(
    plotdata %>% filter(x > 0, y > 0) %>% filter(significant_1 & significant_2) %>% filter(!stringr::str_starts(symbol, "Gm")) %>% filter(!stringr::str_starts(symbol, "Fam")) %>% arrange(desc(x*y)) %>% head(5),
    plotdata %>% filter(x < 0, y < 0) %>% filter(significant_1 & significant_2) %>% filter(!stringr::str_starts(symbol, "Gm")) %>% filter(!stringr::str_starts(symbol, "Fam")) %>% arrange(desc(x*y)) %>% head(5)
  )
  
  plot <- ggplot(plotdata) + geom_point(aes(x, y, color = group)) + 
    scale_color_manual(values = setNames(c(fourway_colors[1], fourway_colors[3], fourway_colors[2], "orange"), groups), name = NULL) +
    ggrepel::geom_label_repel(aes(x, y, label = symbol), data = plotdata_oi, fill = "#FFFFFF99") +
    scale_x_continuous(name = "log2 fold-change TF1+TF2 vs TF1" %>% gsub_combination(combination), limit = c(-lim, lim)) +
    scale_y_continuous(name = "log2 fold-change TF1+TF2 vs TF2" %>% gsub_combination(combination), limit = c(-lim, lim)) +
    theme_classic() +
    theme(legend.position="None") +
    # guides(color=guide_legend(ncol=3, nrow = 2)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_equal() +
    ggtitle(paste0(toupper(combination[[1]]), "+", toupper(combination[[2]])))
  plot
  ggsave(file.path("output", "12-combinations", "unique", paste0(condition_oi, ".pdf")), width = 5, height = 5, plot = plot)
  
  
  # enrichment
  markers_oi <- markers %>% 
    filter(
      (markers$p_val_adj_tf1 < 5e-2) & 
      (abs(markers$avg_log2FC_tf1) > log(1.5)) & 
     (markers$p_val_adj_tf1 < 5e-2) & 
     (abs(markers$avg_log2FC_tf1) > log(1.5))
    )
  
  features_oi <- markers_oi %>% filter(avg_log2FC_tf2 > 0) %>% arrange(desc(avg_log2FC_tf2*avg_log2FC_tf1)) %>% pull(gene)
  length(features_oi)
  FeaturePlot(
    seu_diffexp,
    features = features_oi %>% head(10)
  )
  
  symbols_oi <- data.annot %>% dplyr::slice(match(features_oi, gene)) %>% pull(gene_short_name)
  symbols_oi %>% head(10)
  
  library(enrichR)
  websiteLive <- getOption("enrichR.live")
  setEnrichrSite("Enrichr") # Human genes   
  
  dbs <- listEnrichrDbs()
  dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "PanglaoDB_Augmented_2021")
  enriched <- enrichr(symbols_oi, dbs)
  
  enriched$GO_Biological_Process_2023 %>% arrange(-Combined.Score) %>% head(30)
  enriched$PanglaoDB_Augmented_2021 %>% head(30)
  
  adipo_markers <-  c("Fabp4", "Lpl", "Pparg", "Lipe", "Adipoq", "Cd36",
                "Plin4", "Plin2", "Plin1", "Cebpa", "Cebpb",
                "Cidec")
  
  a <- sum(symbols_oi %in% adipo_markers)
  b <- length(symbols_oi) - a
  c <- length(markers) - a
  d <- nrow(plotdata)
  
  pval <- fisher.test(matrix(c(d, b, c, a), nrow = 2))
  
  
  plot <- plotEnrich(enriched$GO_Biological_Process_2023, showTerms = 5, numChar = 40, y = "Count", orderBy = "P.value", title = paste0(condition_oi, " GO BP enrichment"))
  ggsave(file.path("output", "12-combinations", "unique", paste0(condition_oi, "_enrichment_bp", ".pdf")), height = 2, width = 4)
  
  plot <- plotEnrich(enriched$PanglaoDB_Augmented_2021, showTerms = 5, numChar = 40, y = "Count", orderBy = "Ratio", title = paste0(condition_oi, " PanglaoDB enrichment"))
  plot
  ggsave(file.path("output", "12-combinations", "unique", paste0(condition_oi, "_enrichment_pan", ".pdf")), height = 2, width = 4)
}
