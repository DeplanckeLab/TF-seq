# title: "Dissecting identified reprogramming heterogeneity per TF of interest"
# input data: "TF atlas"
# author: "Wangjie Liu"

setwd("./")

suppressPackageStartupMessages(library(Seurat)) 
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2)) 
suppressPackageStartupMessages(library(cowplot)) 
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gridExtra)) 
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(monocle3))  # 1.0.0
suppressPackageStartupMessages(library(org.Mm.eg.db)) 
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(pheatmap))



#functions:

Run_UMAP_on_subset <- function(tfoi, data.seu, N.pcs = 50,  assay.oi = "corrected_integrated", Res = 0.5, fract_TF_in_clusters.cutoff = 0.6){
  batch.ois <- table(data.seu$batch_overall[data.seu$TF == tfoi])[table(data.seu$batch_overall[data.seu$TF == tfoi]) >= 30] |> names()
  cells.oi <- colnames(data.seu)[data.seu$TF %in% c(tfoi, "D0","D0_confluent") & data.seu$batch_overall %in% batch.ois]
  seu.oi <- subset(data.seu, cells = cells.oi)
  set.seed(10)
  # re-run UMAP and plot only the subset of cells
  DefaultAssay(seu.oi) <- assay.oi
  seu.oi <- RunPCA(seu.oi, npcs = N.pcs)
  seu.oi <- RunUMAP(seu.oi, dims = 1:N.pcs)
  
  seu.oi <- FindNeighbors(seu.oi, dims = 1:N.pcs)
  seu.oi <- FindClusters(seu.oi, resolution = Res)
  Idents(seu.oi) <- seu.oi$seurat_clusters
  clusters_TF.oi <- names(table(seu.oi$seurat_clusters)[table(seu.oi$seurat_clusters[seu.oi$TF == TFoi])/table(seu.oi$seurat_clusters) > fract_TF_in_clusters.cutoff])
  
  seu.oi$meta.new <- seu.oi$seurat_clusters |> as.character()
  seu.oi$meta.new[!seu.oi$seurat_clusters %in% clusters_TF.oi] <- "Controls" 
  return(seu.oi)
}


# monocle3 single-cell pipeline
infer_pseudotime <- function(seu.oi, data.gene_annot, Return_trajectory = FALSE){
  # transfer data and embeddings from seurat to monocle3 cds
  cell.meta <- seu.oi@meta.data[, c("TF", "Dose", "meta.new")]
  gene.meta <- data.frame(ens_id = rownames(seu.oi@assays$RNA@counts), gene_short_name = data.gene_annot[rownames(seu.oi@assays$RNA@counts),"gene_symbol"], gene_type = data.gene_annot[rownames(seu.oi@assays$RNA@counts), "biotype"])
  rownames(gene.meta) <- rownames(seu.oi@assays$RNA@counts)
  cds <- new_cell_data_set(expression_data = seu.oi@assays$RNA@counts, 
                           cell_metadata = cell.meta,
                           gene_metadata = gene.meta)
  
  # Set dimensionality reduction information from Seurat to Monocle3
  # reducedDims(cds)$PCA <- seu.oi@reductions$pca@cell.embeddings
  reducedDims(cds)$UMAP <- seu.oi@reductions$umap@cell.embeddings
  
  # Adding clustering information from Seurat
  cds@clusters$UMAP$clusters <- seu.oi@meta.data$meta.new
  
  # Assign partitions based on clusters if they exist
  cds@clusters$UMAP$partitions <- cds@clusters$UMAP$clusters
  
  # assign Named partitions
  cell_names <- colnames(cds)
  cds@clusters$UMAP$partitions <- setNames(rep(1, length(cell_names)), cell_names)
  
  # Learn the trajectory graph using the UMAP embedding
  cds <- learn_graph(cds, use_partition = TRUE)
  cds <- order_cells(cds)
  
  # Extract pseudotime and add it to seurat metadata
  pseudotime_values <- as.numeric(pseudotime(cds))
  names(pseudotime_values) <- colnames(cds) 
  seu.oi$Pseudotime <- pseudotime_values[rownames(seu.oi@meta.data)]
  if(Return_trajectory == T){
    return(list(seu.oi, cds))
  } else {
    return(seu.oi)
  }
}




# DE analysis
Calculate_DE <- function(seu.oi, data.gene_annot){
  Idents(seu.sub) <- seu.sub$meta.new
  DefaultAssay(seu.oi) <- "RNA"
  DEGs <- FindAllMarkers(seu.oi, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, slot = "data", test.use = "MAST", latent.vars = "batch" )
  DEGs$pct.diff <- abs(DEGs$pct.1 - DEGs$pct.2)
  DEGs$gene_name <- data.gene_annot[DEGs$gene, "gene_symbol"]
  DEGs <- DEGs[DEGs$p_val_adj < 0.05 ,] # 0.1
  dim(DEGs) |> print()
  return(DEGs)
}

# exclude overlapped GO terms
exclude_overlapped_GO <- function(percent.overlap = 0.5, df.GOenrich.res){
  # order GO terms by p.adjust
  # so that the intersection check is performed and 
  # GO terms that are both overlapped with others and have lower p values than overlapped terms 
  # will be excluded
  df.GOenrich.res <- df.GOenrich.res[order(df.GOenrich.res$p.adjust, decreasing = T),] 
  go.terms <- df.GOenrich.res$geneID
  go.terms <- str_split(go.terms, "/")
  names(go.terms) <- df.GOenrich.res$ID
  
  terms.overlap <- c()
  nTerms <- length(go.terms)
  if(nTerms > 1){
    for (i in 1:(nTerms-1)){
      Fractions <- lapply(go.terms[i:nTerms], function(x){
        overlap.fract <- length(intersect(x, go.terms[[i]]))/length(x) 
        return(overlap.fract)
      })
      names(Fractions) <- names(go.terms)[i:nTerms]
      Fractions <- unlist(Fractions)
      # if the term has more than 50% overlap with other terms
      if (table(Fractions > percent.overlap)["TRUE"] == 1){ #only itself
        terms.overlap[i] <- NA
      } else{ # overlap with other terms
        terms.overlap[i] <- names(go.terms)[i]
      }
    }
    # exclude overlapped GO terms
    df.new <- df.GOenrich.res[!df.GOenrich.res$ID %in% terms.overlap,]
  } else if (nTerms == 1){
    df.new <- df.GOenrich.res
  } else if (nTerms <1){ # when 0 terms is significantly enriched
    df.new <- NA
  }
  return(df.new)
}

# GO enrichment analysis and plot non-overlapped top GO terms 
Calculate_topGO_plot <- function(data, num_topGO = 10, deg.sig, cols){
  clusters.TFoi <- unique(as.character(data$meta.new))
  clusters.TFoi <- clusters.TFoi[clusters.TFoi != "Controls"]
  cat("Clusters of interest are:", clusters.TFoi, "\n")
  erich.summary.clusters <- lapply(clusters.TFoi, function(x){
    erich.go.BP <- enrichGO(gene = deg.sig[deg.sig$cluster == x, ]$gene_name,
                            OrgDb = org.Mm.eg.db,
                            keyType = "SYMBOL",
                            ont = "BP",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
    erich.summary <- as.data.frame(erich.go.BP@result)
    erich.summary$cluster <- x
    erich.summary <- erich.summary %>% filter(p.adjust < 0.06) 
    cat("check the summary of GO_BP enrichment for ", x, ": \n")
    print(erich.summary)
    return(erich.summary)
  })
  
  # exclude overlapped GO withih each cluster
  df.non_overlap <- lapply(erich.summary.clusters, function(x) exclude_overlapped_GO(percent.overlap = 0.5,
                                                                                     df.GOenrich.res = x))
  df.non_overlap <- rbindlist(df.non_overlap)
  
  cat("!!! check if enriched GO terms are unique (FALSE) across clusters:", "\n")
  print(table(duplicated(df.non_overlap$ID)))
  
  # exclude non-unique GO terms shared among the clusters
  cat("These descriptions are duplicated among clusters and thus exlcuded: ", df.non_overlap$Description[duplicated(df.non_overlap$Description)], "\n")
  Description.dup <- df.non_overlap$Description[duplicated(df.non_overlap$Description)]
  df.non_overlap <- df.non_overlap[!df.non_overlap$Description %in% Description.dup, ]
  
  results <- list()
  results[[1]] <- df.non_overlap
  
  df.non_overlap <- df.non_overlap %>% 
    group_by(cluster) %>% 
    arrange(p.adjust, .by_group = T) %>% 
    top_n(num_topGO, wt = -log10(p.adjust))
  
  data$meta.new <- fct_rev(data$meta.new)
  df.non_overlap$cluster <- factor(df.non_overlap$cluster, levels = levels(data$meta.new)[levels(data$meta.new) != "Controls"])
  df.non_overlap <- df.non_overlap %>% 
    group_by(cluster) %>% 
    arrange(desc(p.adjust), .by_group = T) 
  
  
  cat("These descriptions are duplicated among clusters and thus exlcuded: ", df.non_overlap$Description[duplicated(df.non_overlap$Description)], "\n")
  
  Description.order <- df.non_overlap$Description[!duplicated(df.non_overlap$Description)]
  
  df.non_overlap$Description <- factor(df.non_overlap$Description, levels = Description.order)
  
  p_GOenrich <- ggplot(df.non_overlap)+
    geom_bar(mapping = aes(x = -log10(p.adjust), y = Description, fill = cluster), 
             stat = "identity",  color = "white",
             width=0.5)+
    scale_fill_manual(values = cols)+
    geom_vline(xintercept = -log10(0.05), lty =2, color = "darkgray")+
    theme_cowplot()+
    scale_x_continuous(position = "top")+
    xlab(expression("-log10(adjusted "*italic("P")*" value)"))+
    ylab("Top 10 unique GO terms")
  print(p_GOenrich)
  results[[2]] <- p_GOenrich
  names(results) <- c("df_uniqueGO_all","plot")
  return(results)
} 

# set color scheme
cols <- c("#2E2591", "#338935", "#5DA899", "#7E2954", "#9F4A96", "#94CBE8", "#E69F00", "#DCC07D", "#C26A76")
names(cols) <- paste0("exp",1:9)


## ----------------------------------------------- load data, metadata, and annotation
# 1) annotation
data.annot <- fread("metadata/GRCm38_gene_annot.tsv", data.table = F) 
rownames(data.annot) <- data.annot$ensembl_id
# 2) aligned dose
Dose_aligned <- readRDS("data/Dose_aligned.rds")
# 3) functional TF atlas
seu_funct <- readRDS("results/C3H10_10X_all_exps_integrated_functional_TF_atlas_G1_Phase_corrected.rds")
dim(seu_funct)
## apply aligned dose
seu_funct$Dose_aligned <- Dose_aligned[colnames(seu_funct)]
seu_funct$Dose_unaligned <- seu_funct$Dose
seu_funct$Dose <- seu_funct$Dose_aligned
## rearrange batch
seu_funct$batch[seu_funct$batch == "exp05"] <- "exp1"
seu_funct$batch[seu_funct$batch == "exp06"] <- "exp2"
seu_funct$batch[seu_funct$batch == "exp07"] <- "exp3"
seu_funct$batch[seu_funct$batch == "exp08"] <- "exp4"
seu_funct$batch[seu_funct$batch == "exp09"] <- "exp5"
seu_funct$batch[seu_funct$batch == "exp10"] <- "exp6"
seu_funct$batch[seu_funct$batch == "exp11"] <- "exp7"
seu_funct$batch[seu_funct$batch == "exp12-13"] <- "exp8"
seu_funct$batch[seu_funct$batch == "exp14"] <- "exp9"
names(cols.batch) <- c(paste0("exp", 1:9))
## set clustering resolution
Idents(seu_funct) <- seu_funct$integrated_snn_res.1.2
# 4) TF atlas
seu <- readRDS("TF_resource_paper/output_version3_2024_June/2-D0-regression/C3H10_10X_all_exps_D0regressed10pc_50pc_integrated_dosealigned.rds")
## subset to G1 (Phase corrected) cells 
seu <- subset(seu, Phase_corrected == "G1")
DefaultAssay(seu) <- "corrected_integrated"
## apply aligned dose
seu$Dose_aligned <- Dose_aligned[colnames(seu)]
seu$Dose_unaligned <- seu$Dose
seu$Dose <- seu$Dose_aligned
## rearrange batch
seu$batch[seu$batch == "exp05"] <- "exp1"
seu$batch[seu$batch == "exp06"] <- "exp2"
seu$batch[seu$batch == "exp07"] <- "exp3"
seu$batch[seu$batch == "exp08"] <- "exp4"
seu$batch[seu$batch == "exp09"] <- "exp5"
seu$batch[seu$batch == "exp10"] <- "exp6"
seu$batch[seu$batch == "exp11"] <- "exp7"
seu$batch[seu$batch == "exp12_13"] <- "exp8"
seu$batch[seu$batch == "exp14"] <- "exp9"
names(cols.batch) <- c(paste0("exp", 1:9))





## ----------------------------------------------- Figure 4a, b - reproducibility script
DefaultAssay(seu) <- "RNA"

## gene list for adiposcore
Adipo <- c("Fabp4", "Lpl", "Pparg", "Lipe", "Adipoq", "Cd36",
           "Plin4", "Plin2", "Plin1", "Cebpa", "Cebpb",
           "Cidec", "Cidea")
Adipo <- data.annot[data.annot$gene_symbol %in% Adipo, "ensembl_id"]


TFoi <- "Cebpa"
seu.oi <- subset(seu, TF %in% c("D0_confluent","Cebpa"))


# calculate adiposcore
seu.oi <- AddModuleScore(seu.oi, features = list(Adipo), name = "adiposcore_")
plotdata_adiposcore <- subset(seu.oi, TF %in% c("D0_confluent","Cebpa"))@meta.data[,c("adiposcore_1","Dose")]
colnames(plotdata_adiposcore) <- c("Adiposcore","Cebpa_dose")
p_adiposcore <- ggplot(plotdata_adiposcore, aes(x = Cebpa_dose, y = Adiposcore)) +
  geom_point(color = "#CCCC66", alpha = 0.8)+
  geom_smooth(show.legend = T, linewidth = 1.5,alpha = 0.15, color = "#999900")+ # se = F,
  theme_cowplot()
p_adiposcore
ggsave(p_adiposcore, filename =  "figures/Fig_4a.pdf", 
       height = 5, width = 6)


cds <- GetAssayData(seu.oi, assay = "RNA", slot = "data")
gene.ois <- c("Pparg","Cebpd","Fabp5","C3")
plotdata_ois <- lapply(gene.ois, function(x){
  plotdata <- data.frame(dose = seu.oi$Dose, expression = cds[data.annot[data.annot$gene_symbol == x, "ensembl_id"], colnames(seu.oi)], gene = x)
})
plotdata_ois <- data.table::rbindlist(plotdata_ois)

p <- ggplot(plotdata_ois, aes(x = dose, y = expression, col = gene)) +
  scale_color_brewer(palette = "Dark2")+
  geom_smooth(show.legend = T, linewidth = 1.5,alpha = 0.15)+ # se = F,
  theme_cowplot()
p
ggsave(p, filename =  "figures/Fig_4b.pdf", 
       height = 5, width = 7)



## ----------------------------------------------- Extended_data_figure 5d - reproducibility script
TFois <- c("Klf4","Runx2","Etv1","Egr1","Grhl2","Esr2")
p_TFois <- lapply(TFois, function(x){
  p <- DimPlot(seu_funct, cells.highlight = colnames(seu_funct)[seu_funct$TF == x], cols = "gray85", cols.highlight = "red", sizes.highlight = 1, pt.size = 0.1)+
    theme(axis.line = element_blank(),
          legend.position = "none", 
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          axis.text=element_blank(), 
          axis.ticks=element_blank(), 
          plot.background=element_blank(), 
          panel.background=element_blank())+
    xlab(label = "") + ylab(label = "")+
    labs(title = x)
  return(p)
})
p_TFois

pdf("figures/Fig_S4.1d-1.pdf", width = 10, height = 12)
marrangeGrob(p_TFois, ncol = 2, nrow = 3)
dev.off()



TFois <- c("Meis2","Myog")
p_TFois <- lapply(TFois, function(x){
  p <- DimPlot(seu_funct, cells.highlight = colnames(seu_funct)[seu_funct$TF == x],cols = "gray85", cols.highlight = "darkgreen", sizes.highlight = 1, pt.size = 0.1)+
    theme(axis.line = element_blank(),
          legend.position = "none", 
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          axis.text=element_blank(), 
          axis.ticks=element_blank(), 
          plot.background=element_blank(), 
          panel.background=element_blank())+
    xlab(label = "") + ylab(label = "")+
    labs(title = x)
  return(p)
})
p_TFois

pdf("figures/Fig_S4.1d-2.pdf", width = 10, height = 4)
marrangeGrob(p_TFois, ncol = 2, nrow = 1)
dev.off()


## -----------------------------------------------Dissecting reprogramming heterogeneity (Klf4)
TFoi <- "Klf4"
seu.sub <- Run_UMAP_on_subset(tfoi = TFoi, 
                              data.seu = seu, 
                              N.pcs = 50, 
                              assay.oi = "corrected_integrated",
                              Res = 0.5,
                              fract_TF_in_clusters.cutoff = 0.6)
DimPlot(seu.sub, label = T)
seu.sub$meta.new <- factor(seu.sub$meta.new, levels = c("Controls", 2, 1, 4))

# check number of cells per cluster
table(seu.sub$meta.new)
# Controls        2        1        4 
#      670      229      266       68

# check number of TF cells per cluster
table(seu.sub$meta.new[seu.sub$TF == TFoi])
# Controls        2        1        4 
#       43      206      266       68 

p1 <- DimPlot(seu.sub, group.by = "batch", cols = cols.batch)+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by batch"))

ggsave(p1, filename = "figures/Fig_5-klf4-batch.pdf", width = 5, height = 4)


cols.tf <- c("gray85","gray70","#b3cde3")
names(cols.tf) <- c("D0","D0_confluent",TFoi)
p2 <- DimPlot(seu.sub, group.by = "TF", cols = cols.tf)+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by TF/control"))

ggsave(p2, filename = "figures/Fig_5-klf4-TF.pdf", width = 5.5, height = 4)

p3 <- FeaturePlot(seu.sub, features = "Dose") + 
  scale_color_viridis_c(breaks = c(0, 6), labels = c(0, 6))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by Dose"))

ggsave(p3, filename = "figures/Fig_5-klf4-dose.pdf", width = 4.5, height = 4)


cols.cluster <- c("lightgrey","#E69F00", "#56B4E9", "#009E73" )
names(cols.cluster) <- c("Controls", 2, 1, 4)

p4 <- DimPlot(seu.sub, group.by = "meta.new", cols = cols.cluster)+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by clusters"))
ggsave(p4, filename = "figures/Fig_5-klf4-cluster.pdf", width = 5, height = 4)


p5 <- VlnPlot(subset(seu.sub, TF == TFoi), group.by = "meta.new", features = "Dose", cols = cols.cluster, pt.size = 0)+
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        legend.position = "none",
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = paste0(TFoi, " dose"))+
  labs(title = NULL)
ggsave(p5, filename = "figures/Fig_5-klf4.pdf", width = 3.5, height = 5)

pairwise.t.test(seu.sub$Dose[seu.sub$TF == TFoi], g = seu.sub$meta.new[seu.sub$TF == TFoi], p.adjust.method = "fdr")

# infer pseudotime
seu.sub <- infer_pseudotime(seu.oi = seu.sub, data.gene_annot = data.annot)

# DE_heatmap
DEGs <- Calculate_DE(seu.sub, data.annot) # significant DEGs

DEGs.oi <- c(top_n(n = 10, x = DEGs[DEGs$cluster == "Controls", ], wt = avg_log2FC*pct.diff)$gene,
             top_n(n = 20, x = DEGs[DEGs$cluster == "2", ], wt = avg_log2FC*pct.diff)$gene,
             top_n(n = 20, x =  DEGs[DEGs$cluster == "1", ], wt = avg_log2FC*pct.diff)$gene,
             top_n(n = 20, x = DEGs[DEGs$cluster == "4", ], wt = avg_log2FC*pct.diff)$gene)

genes.heatmap <- unique(DEGs.oi)

cds <- GetAssayData(seu.sub, assay = "RNA", slot = "data")
heatmap_data <- cds[genes.heatmap,]
rownames(heatmap_data) <- data.annot[rownames(heatmap_data), "gene_symbol"]
Order1 <- seu.sub$Pseudotime[seu.sub$meta.new == "Controls"]
Order1 <- Order1[order(Order1, decreasing = F)] %>% names()
Order2 <- seu.sub$Pseudotime[seu.sub$meta.new == 2]
Order2 <- Order2[order(Order2, decreasing = F)] %>% names()
Order3 <- seu.sub$Pseudotime[seu.sub$meta.new == 1]
Order3 <- Order3[order(Order3, decreasing = F)] %>% names()
Order4 <- seu.sub$Pseudotime[seu.sub$meta.new == 4]
Order4 <- Order4[order(Order4, decreasing = F)] %>% names()
heatmap_data <- heatmap_data[, c(Order1, Order2, Order3, Order4)]
heatmap_data[1:3,1:3]

my_sample_col <- seu.sub@meta.data[,"meta.new", drop = F]
my_sample_col$Dose <- seu.sub$Dose
my_sample_col$TF <- seu.sub$TF

my_colour = list(
  meta.new = c(Controls = "lightgray", '2' = "#E69F00", '1' = "#56B4E9", '4' = "#009E73"),
  Dose = rocket(200, direction = -1),
  TF = c(D0 = "gray85", D0_confluent = "gray70", Klf4 = "#b3cde3")
)

heatmap_data <- as.matrix(heatmap_data)
new_name <- lapply(rownames(heatmap_data), function(x) bquote(italic(.(x))))
heatmap <- pheatmap(heatmap_data,
                    color = colorRampPalette(c("blue","white","red"))(200),
                    annotation_col = my_sample_col,
                    annotation_colors = my_colour,
                    border_color = NA, 
                    show_colnames = F, show_rownames = T,
                    labels_row = as.expression(new_name),
                    cluster_cols = F, cluster_rows = F, 
                    scale = "row", 
                    clustering_method = "ward.D2",
                    annotation_legend = T)
ggsave(heatmap, filename = "figures/Fig_5-klf4-DE_heatmap.pdf", width = 8, height = 12)


# GO enrichment
results.GO <- Calculate_topGO_plot(data = seu.sub,
                                   num_topGO = 10,
                                   deg.sig = DEGs, 
                                   cols = cols.cluster)

p_top10_unique_GO <- results.GO$plot
ggsave(p_top10_unique_GO, file = "figures/Fig_5-klf4-GOenrich.pdf", width = 12, height = 6)
saveRDS(results.GO$df_uniqueGO_all, paste0("figures/", TFoi, "_uniqueGO_all.rds"))


# gene expression (RNAscope probes, Figure 5i)
DefaultAssay(seu.sub) <- "RNA"
p1 <- FeatureScatter(seu.sub, feature1 = "Dose", feature2 = data.annot[data.annot$gene_symbol == "Glul", "ensembl_id"], cols = cols.cluster, slot = "data", group.by = "meta.new") +
  geom_smooth()+
  labs(title = "")+
  xlab("Klf4 dose")+
  ylab("Glul expression level")

p2 <- FeatureScatter(seu.sub, feature1 = "Dose", feature2 = data.annot[data.annot$gene_symbol == "Postn", "ensembl_id"], cols = cols.cluster, slot = "data", group.by = "meta.new") +
  geom_smooth()+
  labs(title = "")+
  xlab("Klf4 dose")+
  ylab("Postn expression level")

p <- p1+p2
ggsave(p, filename = "figures/Fig_S4_Glul_Postn_expression.pdf", width = 12, height = 4)


## -----------------------------------------------Dissecting reprogramming heterogeneity (Esr2) 
TFoi <- "Esr2"
seu.sub <- Run_UMAP_on_subset(tfoi = TFoi, 
                              data.seu = seu, 
                              N.pcs = 20, 
                              assay.oi = "corrected_integrated",
                              Res = 0.7,
                              fract_TF_in_clusters.cutoff = 0.6)
DimPlot(seu.sub, label = T, label.size = 8)+DimPlot(seu.sub, group.by = "TF")

seu.sub$meta.new <- factor(seu.sub$meta.new, levels = c("Controls", 3, 1, 4))

# check number of cells per cluster
table(seu.sub$meta.new)
# Controls        3        1        4 
#      664      162      181      105 

# check number of TF cells per cluster
table(seu.sub$meta.new[seu.sub$TF == TFoi])
# Controls        3        1        4 
#       22      131      181      105 

p1 <- DimPlot(seu.sub, group.by = "batch", cols = cols.batch)+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by rearranged batch"))

ggsave(p1, filename = "figures/ED_Fig6-Esr2-batch.pdf", width = 5, height = 4)


cols.tf <- c("gray85","gray70","#cab2d6")
names(cols.tf) <- c("D0","D0_confluent",TFoi)
p2 <- DimPlot(seu.sub, group.by = "TF", cols = cols.tf)+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by TF/control"))

ggsave(p2, filename = "figures/ED_Fig6-Esr2-TF.pdf", width = 5.5, height = 4)

p3 <- FeaturePlot(seu.sub, features = "Dose") + 
  scale_color_viridis_c(breaks = c(0, 6), labels = c(0, 6))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by Dose"))

ggsave(p3, filename = "figures/ED_Fig6-Esr2-dose.pdf", width = 4.5, height = 4)


cols.cluster <- c("lightgrey","#E69F00", "#56B4E9", "#009E73" )
names(cols.cluster) <- c("Controls", 3, 1, 4)

p4 <- DimPlot(seu.sub, group.by = "meta.new", cols = cols.cluster)+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by clusters"))
ggsave(p4, filename = "figures/ED_Fig6-Esr2-cluster.pdf", width = 5, height = 4)


p5 <- VlnPlot(subset(seu.sub, TF == TFoi), group.by = "meta.new", features = "Dose", cols = cols.cluster, pt.size = 0)+
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        legend.position = "none",
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = paste0(TFoi, " dose"))+
  labs(title = NULL)
ggsave(p5, filename = "figures/ED_Fig6-Esr2-dose-vlnplot.pdf", width = 3.5, height = 5)

pairwise.t.test(seu.sub$Dose[seu.sub$TF == TFoi], g = seu.sub$meta.new[seu.sub$TF == TFoi], p.adjust.method = "fdr")
# 	Pairwise comparisons using t tests with pooled SD
# 
# data:  seu.sub$Dose[seu.sub$TF == TFoi] and seu.sub$meta.new[seu.sub$TF == TFoi]
# 
#   Controls 3      1
# 3 0.091    -      -
# 1 <2e-16   <2e-16 -
# 4 <2e-16   <2e-16 <2e-16
# 
# P value adjustment method: fdr

# gene expression for RNAscope probes (Esr2)
DefaultAssay(seu.sub) <- "RNA"
p1 <- FeatureScatter(seu.sub, feature1 = "Dose", feature2 = data.annot[data.annot$gene_symbol == "Gng12", "ensembl_id"], cols = cols.cluster, slot = "data", group.by = "meta.new") +
  geom_smooth()+
  labs(title = "")+
  xlab("Esr2 dose")+
  ylab("Gng12 expression level")

p2 <- FeatureScatter(seu.sub, feature1 = "Dose", feature2 = data.annot[data.annot$gene_symbol == "Aspn", "ensembl_id"], cols = cols.cluster, slot = "data", group.by = "meta.new") +
  geom_smooth()+
  labs(title = "")+
  xlab("Esr2 dose")+
  ylab("Aspn expression level")

p <- p1+p2
ggsave(p, filename = "figures/ED_Fig6-Esr2_Glul_Postn_expression.pdf", width = 12, height = 4)




## -----------------------------------------------Dissecting reprogramming heterogeneity (Egr1)
TFoi <- "Egr1"
seu.sub <- Run_UMAP_on_subset(tfoi = TFoi, 
                              data.seu = seu, 
                              N.pcs = 30, 
                              assay.oi = "corrected_integrated",
                              Res = 0.3,
                              fract_TF_in_clusters.cutoff = 0.6)
seu.sub$meta.new <- factor(seu.sub$meta.new, levels = c("Controls", 1, 2))

# check number of cells per cluster
table(seu.sub$meta.new)
# Controls        1        2 
#      785      220       86 

# check number of TF cells per cluster
table(seu.sub$meta.new[seu.sub$TF == TFoi])
# Controls        1        2 
#      133      199       86 

p1 <- DimPlot(seu.sub, group.by = "batch", cols = cols.batch)+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by rearranged batch"))
ggsave(p1, filename = paste0("figures/ED_Fig6-",TFoi,"-batch.pdf"), width = 5, height = 4)

cols.tf <- c("gray85","gray70","#FB9A99")
names(cols.tf) <- c("D0","D0_confluent",TFoi)
p2 <- DimPlot(seu.sub, group.by = "TF", cols = cols.tf)+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by TF/control"))
ggsave(p2, filename = paste0("figures/ED_Fig6-",TFoi,"-TF.pdf"), width = 5.5, height = 4)

p3 <- FeaturePlot(seu.sub, features = "Dose") + 
  scale_color_viridis_c(breaks = c(0, 7), labels = c(0, 7))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by Dose"))
ggsave(p3, filename = paste0("figures/ED_Fig6-",TFoi,"-dose.pdf"), width = 4.5, height = 4)


cols.cluster <- c("lightgrey","#009E73", "#D55E00" )
names(cols.cluster) <- c("Controls", 1, 2)
p4 <- DimPlot(seu.sub, group.by = "meta.new", cols = cols.cluster)+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by clusters"))
ggsave(p4, filename = paste0("figures/ED_Fig6-",TFoi,"-cluster.pdf"), width = 5, height = 4)

p5 <- VlnPlot(subset(seu.sub, TF == TFoi), group.by = "meta.new", features = "Dose", cols = cols.cluster, pt.size = 0)+
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        legend.position = "none",
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = paste0(TFoi, " dose"))+
  labs(title = NULL)
ggsave(p5, filename = paste0("figures/ED_Fig6-",TFoi,"-dose-vlnplot.pdf"), width = 3.5, height = 5)
pairwise.t.test(seu.sub$Dose[seu.sub$TF == TFoi], g = seu.sub$meta.new[seu.sub$TF == TFoi], p.adjust.method = "fdr")
# 	Pairwise comparisons using t tests with pooled SD
# 
# data:  seu.sub$Dose[seu.sub$TF == TFoi] and seu.sub$meta.new[seu.sub$TF == TFoi]
# 
#   Controls 1
# 1 <2e-16   -
# 2 <2e-16   <2e-16
# 
# P value adjustment method: fdr



## -----------------------------------------------Dissecting reprogramming heterogeneity (Grhl2) 
TFoi <- "Grhl2"
seu.sub <- Run_UMAP_on_subset(tfoi = TFoi, 
                              data.seu = seu, 
                              N.pcs = 50, 
                              assay.oi = "corrected_integrated",
                              Res = 0.8,
                              fract_TF_in_clusters.cutoff = 0.6)
DimPlot(seu.sub, label = T, label.size = 8)+DimPlot(seu.sub, group.by = "TF")+DimPlot(seu.sub, group.by = "batch",cols = cols.batch)+DimPlot(seu.sub, group.by = "Phase")
seu.sub$meta.new <- factor(seu.sub$meta.new, levels = c("Controls", 2, 3))
# check number of cells per cluster
table(seu.sub$meta.new)
# Controls        2        3 
#      805      137       87 
# check number of TF cells per cluster
table(seu.sub$meta.new[seu.sub$TF == TFoi])
# Controls        2        3 
#      134      135       87 

p1 <- DimPlot(seu.sub, group.by = "batch", cols = cols.batch)+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by rearranged batch"))
ggsave(p1, filename = paste0("figures/ED_Fig6-",TFoi,"-batch.pdf"), width = 5, height = 4)

cols.tf <- c("gray85","gray70","#E6AB02")
names(cols.tf) <- c("D0","D0_confluent",TFoi)
p2 <- DimPlot(seu.sub, group.by = "TF", cols = cols.tf)+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by TF/control"))
ggsave(p2, filename = paste0("figures/ED_Fig6-",TFoi,"-TF.pdf"), width = 5.5, height = 4)

p3 <- FeaturePlot(seu.sub, features = "Dose") + 
  scale_color_viridis_c(breaks = c(0, 6), labels = c(0, 6))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by Dose"))
ggsave(p3, filename = paste0("figures/ED_Fig6-",TFoi,"-dose.pdf"), width = 4.5, height = 4)


cols.cluster <- c("lightgrey","#D55E00", "#33A02C" )
names(cols.cluster) <- c("Controls", 2, 3)
p4 <- DimPlot(seu.sub, group.by = "meta.new", cols = cols.cluster)+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by clusters"))
ggsave(p4, filename = paste0("figures/ED_Fig6-",TFoi,"-cluster.pdf"), width = 5, height = 4)

p5 <- VlnPlot(subset(seu.sub, TF == TFoi), group.by = "meta.new", features = "Dose", cols = cols.cluster, pt.size = 0)+
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        legend.position = "none",
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = paste0(TFoi, " dose"))+
  labs(title = NULL)
ggsave(p5, filename = paste0("figures/ED_Fig6-",TFoi,"-dose-vlnplot.pdf"), width = 3.5, height = 5)

pairwise.t.test(seu.sub$Dose[seu.sub$TF == TFoi], g = seu.sub$meta.new[seu.sub$TF == TFoi], p.adjust.method = "fdr")
# 	Pairwise comparisons using t tests with pooled SD
# 
# data:  seu.sub$Dose[seu.sub$TF == TFoi] and seu.sub$meta.new[seu.sub$TF == TFoi]
# 
#   Controls 2
# 2 <2e-16   -
# 3 <2e-16   <2e-16
# 
# P value adjustment method: fdr


## -----------------------------------------------Dissecting reprogramming heterogeneity (Meis2)
TFoi <- "Meis2"
seu.sub <- Run_UMAP_on_subset(tfoi = TFoi,
                              data.seu = seu,
                              N.pcs = 30,
                              assay.oi = "corrected_integrated",
                              Res = 1,
                              fract_TF_in_clusters.cutoff = 0.6)
seu.sub$meta.new <- factor(seu.sub$meta.new, levels =  c("Controls", 8, 5, 4, 1, 7))
# check number of cells per cluster
table(seu.sub$meta.new)
# Controls        8        5        4        1        7 
#      860       38       86      124      309       39 

# check number of TF cells per cluster
table(seu.sub$meta.new[seu.sub$TF == TFoi])
# Controls        8        5        4        1        7 
#      212       29       70      124      309       39  

p1 <- DimPlot(seu.sub, group.by = "batch", cols = cols.batch)+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by rearranged batch"))
ggsave(p1, filename = paste0("figures/ED_Fig_7d-",TFoi,"-batch.pdf"), width = 5, height = 4)

cols.tf <- c("gray85","gray70","#FB9A99")
names(cols.tf) <- c("D0","D0_confluent",TFoi)
p2 <- DimPlot(seu.sub, group.by = "TF", cols = cols.tf)+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by TF/control"))
ggsave(p2, filename = paste0("figures/ED_Fig_7d-",TFoi,"-TF.pdf"), width = 5.5, height = 4)

p3 <- FeaturePlot(seu.sub, features = "Dose") + 
  scale_color_viridis_c(breaks = c(0, 8), labels = c(0, 8))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by Dose"))
ggsave(p3, filename = paste0("figures/Fig_5k-",TFoi,"-dose.pdf"), width = 4.5, height = 4)


cols.cluster <- c("lightgrey","#F0E442","#009E73", "#56B4E9","#D55E00", "#6A3D9A")
names(cols.cluster) <- c("Controls", 8, 5, 4, 1, 7)
p4 <- DimPlot(seu.sub, group.by = "meta.new", cols = cols.cluster)+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by clusters"))
ggsave(p4, filename = paste0("figures/Fig_5k-",TFoi,"-cluster.pdf"), width = 5, height = 4)


p5 <- VlnPlot(subset(seu.sub, TF == TFoi), group.by = "meta.new", features = "Dose", cols = cols.cluster, pt.size = 0)+
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        legend.position = "none",
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = paste0(TFoi, " dose"))+
  labs(title = NULL)
ggsave(p5, filename = paste0("figures/Fig_5l-",TFoi,"-dose-vlnplot.pdf"), width = 3.5, height = 5)
pairwise.t.test(seu.sub$Dose[seu.sub$TF == TFoi], g = seu.sub$meta.new[seu.sub$TF == TFoi], p.adjust.method = "fdr")
# 	Pairwise comparisons using t tests with pooled SD 
# 
# data:  seu.sub$Dose[seu.sub$TF == TFoi] and seu.sub$meta.new[seu.sub$TF == TFoi] 
# 
#   Controls 8       5       4       1      
# 8 3.3e-08  -       -       -       -      
# 5 7.9e-07  0.066   -       -       -      
# 4 1.1e-10  0.086   0.687   -       -      
# 1 < 2e-16  0.173   2.3e-07 4.9e-09 -      
# 7 < 2e-16  4.9e-09 < 2e-16 < 2e-16 8.6e-12
# 
# P value adjustment method: fdr 

DefaultAssay(seu.sub) <- "RNA"
for (i in unique(seu.sub$meta.new)){
  genes.oi <- DEGs[DEGs$cluster == i,]$gene
  seu.sub <- AddModuleScore(seu.sub, features = list(genes.oi), name = paste0("gEX_score_",i,"_"))
}

# plot modulescore on UMAP
p0 <- FeaturePlot(seu.sub, features = "gEX_score_Controls_1", order = T, pt.size = 1)+
  scale_color_viridis_c(option = "B", breaks = c(min(seu.sub$gEX_score_8_1), max(seu.sub$gEX_score_8_1)), labels = c("min", "max"))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab("")+ylab("")+
  labs(title = "", subtitle = "Control module")

p1 <- FeaturePlot(seu.sub, features = "gEX_score_8_1", order = T, pt.size = 1)+
  scale_color_viridis_c(option = "B", breaks = c(min(seu.sub$gEX_score_8_1), max(seu.sub$gEX_score_8_1)), labels = c("min", "max"))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab("")+ylab("")+
  labs(title = "", subtitle = "Module 1")

p2 <- FeaturePlot(seu.sub, features = "gEX_score_5_1", order = T, pt.size = 1)+
  scale_color_viridis_c(option = "B", breaks = c(min(seu.sub$gEX_score_5_1), max(seu.sub$gEX_score_5_1)), labels = c("min", "max"))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab("")+ylab("")+
  labs(title = "", subtitle = "Module 2")

p3 <- FeaturePlot(seu.sub, features = "gEX_score_4_1", order = T, pt.size = 1)+
  scale_color_viridis_c(option = "B", breaks = c(min(seu.sub$gEX_score_4_1), max(seu.sub$gEX_score_4_1)), labels = c("min", "max"))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab("")+ylab("")+
  labs(title = "", subtitle = "Module 3")

p4 <- FeaturePlot(seu.sub, features = "gEX_score_1_1", order = T, pt.size = 1)+
  scale_color_viridis_c(option = "B", breaks = c(min(seu.sub$gEX_score_1_1), max(seu.sub$gEX_score_1_1)), labels = c("min", "max"))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab("")+ylab("")+
  labs(title = "", subtitle = "Module 4")


p5 <- FeaturePlot(seu.sub, features = "gEX_score_7_1", order = T, pt.size = 1)+
  scale_color_viridis_c(option = "B", breaks = c(min(seu.sub$gEX_score_7_1), max(seu.sub$gEX_score_7_1)), labels = c("min", "max"))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab("")+ylab("")+
  labs(title = "", subtitle = "Module 5")

pdf("figures/ED_Fig5e-Meis2_modulescore_UMAP.pdf", width = 14, height = 8)
marrangeGrob(list(p0,p1,p2,p3,p4,p5), ncol = 3, nrow = 2)
dev.off()



## Meis2 - exp4
TFoi <- "Meis2"
seu.sub.batch <- subset(seu, batch == "exp4" & TF == TFoi)
DefaultAssay(seu.sub.batch) <- "RNA"
N.pcs <- 10
seu.sub.batch <- seu.sub.batch |> NormalizeData() |> FindVariableFeatures() |> ScaleData() |> RunPCA(npcs = N.pcs) |> RunUMAP(dims = 1:N.pcs) |> FindNeighbors(dims = 1:N.pcs ) |> FindClusters(resolution = 0.8)
DimPlot(seu.sub.batch, label = T, label.size = 8)+DimPlot(seu.sub.batch, group.by = "TF")+DimPlot(seu.sub.batch, group.by = "Phase")+FeaturePlot(seu.sub.batch, features = "Dose", order = T)+scale_color_viridis_c()
VlnPlot(seu.sub.batch, features = "Dose")
pairwise.t.test(seu.sub.batch$Dose, g = seu.sub.batch$seurat_clusters, p.adjust.method = "fdr")
#	Pairwise comparisons using t tests with pooled SD 
# data:  seu.sub.batch$Dose and seu.sub.batch$seurat_clusters 
# 
#   0       1      
# 1 1.2e-05 -      
# 2 0.00011 0.98539
# 
# P value adjustment method: fdr 

DefaultAssay(seu.sub.batch) <- "RNA"
for (i in unique((DEGs$cluster))){
  print(i)
  genes.oi <- DEGs[DEGs$cluster == i,]$gene
  seu.sub.batch <- AddModuleScore(seu.sub.batch, features = list(genes.oi), name = paste0("gEX_score_",i,"_"))
}

p1 <- FeaturePlot(seu.sub.batch, features = "Dose", pt.size = 2) + 
  scale_color_viridis_c(breaks = c(1.1, 7), labels = c(1.1, 7))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = NULL,  subtitle = "Dose")

cols.cluster <- c("#009E73", "#D55E00", "#56B4E9","lightgrey","#F0E442","#6A3D9A")
names(cols.cluster) <- c(0,1,2)
p2 <- DimPlot(seu.sub.batch, group.by = "seurat_clusters", cols = c("#009E73","#D55E00","#56B4E9"), pt.size = 2)+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(subtitle = paste0("n = ", ncol(seu.sub.batch)), title = NULL)
pdf("figures/ED_Fig7g-Meis2_exp4-Dose-and-cluster.pdf", width = 6, height = 3)
marrangeGrob(list(p1,p2), ncol = 2, nrow = 1)
dev.off()

p <- VlnPlot(seu.sub.batch, features = "Dose", cols = cols.cluster, pt.size = 0)+
  theme(legend.position = "none")+
  xlab("")+
  ylab("Meis2 dose")
ggsave(p, filename = "figures/ED_Fig7h-Meis2_exp4-vlnplot-dose.pdf", width = 4, height = 3)

# plot modulescore on UMAP
p0 <- FeaturePlot(seu.sub.batch, features = "gEX_score_Controls_1", order = T, pt.size = 1)+
  scale_color_viridis_c(option = "B", breaks = c(min(seu.sub.batch$gEX_score_8_1), max(seu.sub.batch$gEX_score_8_1)), labels = c("min", "max"))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab("")+ylab("")+
  labs(title = "", subtitle = "Control module")

p1 <- FeaturePlot(seu.sub.batch, features = "gEX_score_8_1", order = T, pt.size = 1)+
  scale_color_viridis_c(option = "B", breaks = c(min(seu.sub.batch$gEX_score_8_1), max(seu.sub.batch$gEX_score_8_1)), labels = c("min", "max"))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab("")+ylab("")+
  labs(title = "", subtitle = "Module 1")

p2 <- FeaturePlot(seu.sub.batch, features = "gEX_score_5_1", order = T, pt.size = 1)+
  scale_color_viridis_c(option = "B", breaks = c(min(seu.sub.batch$gEX_score_5_1), max(seu.sub.batch$gEX_score_5_1)), labels = c("min", "max"))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab("")+ylab("")+
  labs(title = "", subtitle = "Module 2")

p3 <- FeaturePlot(seu.sub.batch, features = "gEX_score_4_1", order = T, pt.size = 1)+
  scale_color_viridis_c(option = "B", breaks = c(min(seu.sub.batch$gEX_score_4_1), max(seu.sub.batch$gEX_score_4_1)), labels = c("min", "max"))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab("")+ylab("")+
  labs(title = "", subtitle = "Module 3")

p4 <- FeaturePlot(seu.sub.batch, features = "gEX_score_1_1", order = T, pt.size = 1)+
  scale_color_viridis_c(option = "B", breaks = c(min(seu.sub.batch$gEX_score_1_1), max(seu.sub.batch$gEX_score_1_1)), labels = c("min", "max"))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab("")+ylab("")+
  labs(title = "", subtitle = "Module 4")


p5 <- FeaturePlot(seu.sub.batch, features = "gEX_score_7_1", order = T, pt.size = 1)+
  scale_color_viridis_c(option = "B", breaks = c(min(seu.sub.batch$gEX_score_7_1), max(seu.sub.batch$gEX_score_7_1)), labels = c("min", "max"))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab("")+ylab("")+
  labs(title = "", subtitle = "Module 5")

pdf("figures/ED_Fig7i-Meis2_all-modulescore_UMAP-exp4.pdf", width = 14, height = 8)
marrangeGrob(list(p0,p1,p2,p3,p4,p5), ncol = 3, nrow = 2)
dev.off()


## Mesie- exp 9
set.seed(10)
TFoi <- "Meis2"
seu.sub.batch <- subset(seu, batch == "exp9" & TF == TFoi)
DefaultAssay(seu.sub.batch) <- "RNA"
N.pcs <- 15
seu.sub.batch <- seu.sub.batch |> NormalizeData() |> FindVariableFeatures() |> ScaleData() |> RunPCA(npcs = N.pcs) |> RunUMAP(dims = 1:N.pcs) |> FindNeighbors(dims = 1:N.pcs ) |> FindClusters(resolution = 0.8)

pairwise.t.test(seu.sub.batch$Dose, g = seu.sub.batch$seurat_clusters, p.adjust.method = "fdr")
#	Pairwise comparisons using t tests with pooled SD 
#	  0       1       2       3       4      
#	1 9.5e-10 -       -       -       -      
#	2 2.3e-15 0.05445 -       -       -      
#	3 < 2e-16 < 2e-16 2.2e-14 -       -      
#	4 < 2e-16 0.00019 0.06276 1.1e-08 -      
#	5 1.0e-07 < 2e-16 < 2e-16 < 2e-16 < 2e-16

# P value adjustment method: fdr 

DefaultAssay(seu.sub.batch) <- "RNA"
for (i in unique((DEGs$cluster))){
  genes.oi <- DEGs[DEGs$cluster == i,]$gene
  seu.sub.batch <- AddModuleScore(seu.sub.batch, features = list(genes.oi), name = paste0("gEX_score_",i,"_"))
}

p1 <- FeaturePlot(seu.sub.batch, features = "Dose", pt.size = 2) + 
  scale_color_viridis_c(breaks = c(1, 8), labels = c(1, 8))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = NULL,  subtitle = "Dose")


cols.cluster <- c("#009E73", "#D55E00", "#56B4E9","lightgrey","#F0E442","#6A3D9A")
names(cols.cluster) <- c(0,1,2,3,4,5)

p2 <- DimPlot(seu.sub.batch, group.by = "seurat_clusters", cols = cols.cluster, label = T, label.size = 8, pt.size = 1)+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(subtitle = paste0("n = ", ncol(seu.sub.batch)), title = NULL)

p3 <- FeaturePlot(seu.sub.batch, features = "gEX_score_8_1", order = T, pt.size = 1)+
  scale_color_viridis_c(option = "B", breaks = c(min(seu.sub.batch$gEX_score_8_1), max(seu.sub.batch$gEX_score_8_1)), labels = c("min", "max"))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab("")+ylab("")+
  labs(title = "", subtitle = "Module 1")

p4 <- FeaturePlot(seu.sub.batch, features = "gEX_score_5_1", order = T, pt.size = 1)+
  scale_color_viridis_c(option = "B", breaks = c(min(seu.sub.batch$gEX_score_5_1), max(seu.sub.batch$gEX_score_5_1)), labels = c("min", "max"))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab("")+ylab("")+
  labs(title = "", subtitle = "Module 25")

p5 <- FeaturePlot(seu.sub.batch, features = "gEX_score_4_1", order = T, pt.size = 1)+
  scale_color_viridis_c(option = "B", breaks = c(min(seu.sub.batch$gEX_score_4_1), max(seu.sub.batch$gEX_score_4_1)), labels = c("min", "max"))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab("")+ylab("")+
  labs(title = "", subtitle = "Module 3")

p6 <- FeaturePlot(seu.sub.batch, features = "gEX_score_1_1", order = T, pt.size = 1)+
  scale_color_viridis_c(option = "B", breaks = c(min(seu.sub.batch$gEX_score_1_1), max(seu.sub.batch$gEX_score_1_1)), labels = c("min", "max"))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab("")+ylab("")+
  labs(title = "", subtitle = "Module 4")


p7 <- FeaturePlot(seu.sub.batch, features = "gEX_score_7_1", order = T, pt.size = 1)+
  scale_color_viridis_c(option = "B", breaks = c(min(seu.sub.batch$gEX_score_7_1), max(seu.sub.batch$gEX_score_7_1)), labels = c("min", "max"))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab("")+ylab("")+
  labs(title = "", subtitle = "Module 5")

pdf("figures/ED_Fig7j-Meis2_exp9-Dose-and-cluster.pdf", width = 6, height = 3)
marrangeGrob(list(p1,p2), ncol = 2, nrow = 1)
dev.off()


pdf("figures/ED_Fig7l-Meis2_exp9-modulescores.pdf", width = 9, height = 6)
marrangeGrob(list(p3,p4,p5,p6,p7), ncol = 3, nrow = 2)
dev.off()

p <- VlnPlot(seu.sub.batch, features = "Dose", cols = cols.cluster, pt.size = 0)+
  theme(legend.position = "none")+
  xlab("")+
  ylab("Meis2 dose")
ggsave(p, filename = "figures/ED_Fig7k_Meis2_exp9-vlnplot-dose.pdf", width = 4, height = 3)


## -----------------------------------------------Dissecting reprogramming heterogeneity (Myog) 
TFoi <- "Myog"
seu.sub <- Run_UMAP_on_subset(tfoi = TFoi,
                              data.seu = seu,
                              N.pcs = 30,
                              assay.oi = "corrected_integrated",
                              Res = 0.6,
                              fract_TF_in_clusters.cutoff = 0.6)

seu.sub$meta.new <- factor(seu.sub$meta.new, levels =  c("Controls", 3, 5, 4))
# check number of cells per cluster
table(seu.sub$meta.new)
# Controls        3        5        4 
#     1087      225       70      139 
# check number of TF cells per cluster
table(seu.sub$meta.new[seu.sub$TF == TFoi])
# Controls        3        5        4 
#       18      223       70      139 

p1 <- DimPlot(seu.sub, group.by = "batch", cols = cols.batch, shuffle = T)+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by rearranged batch"))
ggsave(p1, filename = paste0("figures/ED_Fig7a-",TFoi,"-batch.pdf"), width = 5, height = 4)

cols.tf <- c("gray85","gray70","#cc6699")
names(cols.tf) <- c("D0","D0_confluent",TFoi)
p2 <- DimPlot(seu.sub, group.by = "TF", cols = cols.tf)+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by TF/control"))
ggsave(p2, filename = paste0("figures/ED_Fig7b-",TFoi,"-TF.pdf"), width = 5.5, height = 4)

p3 <- FeaturePlot(seu.sub, features = "Dose", order = T) + 
  scale_color_viridis_c(breaks = c(0, 6), labels = c(0, 6))+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by Dose"))
ggsave(p3, filename = paste0("figures/ED_Fig7b-",TFoi,"-dose.pdf"), width = 4.5, height = 4)

cols.cluster <- c("lightgrey","#cab2d6","#cc99ff", "#984ea3")
names(cols.cluster) <- c("Controls", 3, 5, 4)
p4 <- DimPlot(seu.sub, group.by = "meta.new", cols = cols.cluster)+
  theme(axis.line = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = "")+
  labs(title = paste0(TFoi, " colored by clusters"))
ggsave(p4, filename = paste0("figures/ED_Fig7b-",TFoi,"-cluster.pdf"), width = 5, height = 4)

p5 <- VlnPlot(subset(seu.sub, TF == TFoi), group.by = "meta.new", features = "Dose", cols = cols.cluster, pt.size = 0)+
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        legend.position = "none",
        plot.background=element_blank(), 
        panel.background=element_blank())+
  xlab(label = "") + ylab(label = paste0(TFoi, " dose"))+
  labs(title = NULL)
ggsave(p5, filename = paste0("figures/ED_Fig7c-",TFoi,"-dose-vlnplot.pdf"), width = 3.5, height = 5)

pairwise.t.test(seu.sub$Dose[seu.sub$TF == TFoi], g = seu.sub$meta.new[seu.sub$TF == TFoi], p.adjust.method = "fdr")
# 	Pairwise comparisons using t tests with pooled SD 
# 
# data:  seu.sub$Dose[seu.sub$TF == TFoi] and seu.sub$meta.new[seu.sub$TF == TFoi] 
# 
#   Controls 3       5      
# 3 5.7e-05  -       -      
# 5 2.5e-08  0.00011 -      
# 4 4.2e-08  0.00019 0.34938
# 
# P value adjustment method: fdr 

