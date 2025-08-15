calculate_diffexp_discrete <- function(
  seu_diffexp, 
  col = "TF",
  model_cellcycle = FALSE,
  model_batch = FALSE
) {
  requireNamespace("edgeR")
  
  # subset the dataset
  counts_diffexp <- seu_diffexp@assays$RNA@counts
  metadata_diffexp <- seu_diffexp@meta.data %>% 
    dplyr::select(cell)
  metadata_diffexp[["x"]] <- factor(seu_diffexp@meta.data[[col]])
  
  # create base formula
  # https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_edgeRQLF.R
  dge <- edgeR::DGEList(counts_diffexp)
  dge <- edgeR::calcNormFactors(dge)
  metadata_diffexp$cdr <- scale(Matrix::colMeans(counts_diffexp > 0))
  formula <- "~ cdr + x"
  
  # create design
  design <- model.matrix(as.formula(formula), metadata_diffexp)
  colnames(design)[startsWith(colnames(design), "x")] <- "x"
  dge <- edgeR::estimateDisp(dge, design = design)
  
  # fit the model
  fit <- edgeR::glmQLFit(dge, design = design, coef = "x")
  
  qlf <- edgeR::glmQLFTest(fit, coef = "x")
  tt <- edgeR::topTags(qlf, n = Inf)
  tt
  
  # put into tidy format
  scores <- process_edger(tt)
  scores$coef <- qlf$coefficients[scores$gene, "x"]
  
  scores
}

process_edger <- function(tt) {
  tt$table %>% 
    rownames_to_column("gene") %>% 
    rename(
      lfc = logFC,
      pval = PValue,
      qval = FDR
    ) %>% 
    arrange(-lfc)
}