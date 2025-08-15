get_diffexp <- function(TFsoi, seu_all) {
  seu_diffexp <- seu_all[, seu_all@meta.data %>% filter(TF %in% c(TFsoi[[1]], TFsoi[[2]], "D0", paste0(TFsoi[[1]], "-", TFsoi[[2]]), paste0(TFsoi[[2]], "-", TFsoi[[1]]))) %>% pull(cell)]

  seu_diffexp$vector1 <- seu_all[["oeumi"]]@counts[TFsoi[[1]], ]
  seu_diffexp$vector2 <- seu_all[["oeumi"]]@counts[TFsoi[[2]], ]

  seu_diffexp$logvector1 <- log1p(seu_diffexp$vector1)
  seu_diffexp$logvector2 <- log1p(seu_diffexp$vector2)
  
  # ggplot(seu_diffexp@meta.data) + geom_point(aes(logvector1, logvector2, color = TF))
  
  ##
  
  # subset the dataset
  counts_diffexp <- seu_diffexp@assays$RNA@counts
  metadata_diffexp <- seu_diffexp@meta.data %>% 
    dplyr::select(c(logvector1, logvector2)) 
  
  # create base formula
  # dosage + s_score + dosage * s_score
  # then, compare to model without interaction dosage + s_score
  # https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_edgeRQLF.R
  dge <- edgeR::DGEList(counts = counts_diffexp)
  dge <- edgeR::calcNormFactors(dge) # calculate scaling factors to convert raw library size to effective library size
  metadata_diffexp$cdr <- scale(Matrix::colMeans(counts_diffexp > 0)) # center or scale the columns of matrix
  # cdr -> scaled gene expression mean
  formula <- "~cdr + logvector1 + logvector2 + logvector1*logvector2"
  
  # create design
  design <- model.matrix(as.formula(formula), metadata_diffexp)
  colnames(design)
  ## notice that the first level of TF group like Ctr is taking as base automatically.
  #colnames(design)[startsWith(colnames(design), "x")] <- "x"
  dge <- edgeR::estimateDisp(dge, design = design)# maximize negative binomial likelihood 
  
  
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
  
  # build the glm fit function input into mclapply
  fit_glm_coef <- function(Coef, Dge, Design){
    fit <- edgeR::glmQLFit(Dge, design = Design)
    # for a given coefficient, conduct genewise statistical test
    # glmQLFTest is similar to glmLRT except that it replaces likelihood ratio tests with empirical Bayes quasi-likelihood F-tests. 
    # use QL F-tests instead of the more usual likelihood ratio tests (LRT) as they give stricter error rate control by accounting for the uncertainty in dispersion estimation:
    # The p-values from glmQLFTest are always greater than or equal to those that would be obtained from glmLRT using the same negative binomial dispersions.
    # argument coef: integer or character index vector indicating which coefficients of the linear model are to be tested equal to zero.
    qlf <- edgeR::glmQLFTest(fit, coef = Coef)
    # use topTags to extract the most differentially expressed genes/tags 
    # ranked either by p-value or absolute log foldchange of coefficients 
    tt <- edgeR::topTags(qlf, n = Inf) 
    score <- as.data.frame(tt$table)
    score$gene <- rownames(score)
    # score <- cbind(score, qlf$coefficients[rownames(score), ])
    return(score)
  }
  
  score1 <- fit_glm_coef("logvector1", dge, design)
  score2 <- fit_glm_coef("logvector2", dge, design)
  score12 <- fit_glm_coef("logvector1:logvector2", dge, design)
  
  score1$significant <- score1$FDR < 0.05
  score2$significant <- score2$FDR < 0.05
  score12$significant <- score12$FDR < 0.05
  
  scores <- score1 %>% left_join(score2, by = "gene", suffix = c("_1", "_2")) %>% left_join(score12 %>% rename_at(vars(-one_of("gene")), ~paste0(., "_12")), by = "gene")
  # scores$symbol <- data.annot[scores$gene, "gene_short_name"]
  scores
}

get_diffexp_discrete <- function(TFsoi, seu_all) {
  seu_diffexp <- seu_all[, seu_all@meta.data %>% filter(TF %in% c(TFsoi[[1]], TFsoi[[2]], "D0", paste0(TFsoi[[1]], "-", TFsoi[[2]]), paste0(TFsoi[[2]], "-", TFsoi[[1]]))) %>% pull(cell)]

  seu_diffexp$vector1 <- seu_all[["oeumi"]]@counts[TFsoi[[1]], ]
  seu_diffexp$vector2 <- seu_all[["oeumi"]]@counts[TFsoi[[2]], ]

  seu_diffexp$logvector1 <- log1p(seu_diffexp$vector1)
  seu_diffexp$logvector2 <- log1p(seu_diffexp$vector2)
  
  seu_diffexp$high1 <- seu_diffexp$vector1 > 10
  seu_diffexp$high2 <- seu_diffexp$vector2 > 10
  
  ##
  
  # subset the dataset
  counts_diffexp <- seu_diffexp@assays$RNA@counts
  metadata_diffexp <- seu_diffexp@meta.data %>% 
    dplyr::select(c(high1, high2)) 
  
  # create base formula
  # dosage + s_score + dosage * s_score
  # then, compare to model without interaction dosage + s_score
  # https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_edgeRQLF.R
  dge <- edgeR::DGEList(counts = counts_diffexp)
  dge <- edgeR::calcNormFactors(dge) # calculate scaling factors to convert raw library size to effective library size
  metadata_diffexp$cdr <- scale(Matrix::colMeans(counts_diffexp > 0)) # center or scale the columns of matrix
  # cdr -> scaled gene expression mean
  formula <- "~cdr + high1 + high2 + high1*high2"
  
  # create design
  design <- model.matrix(as.formula(formula), metadata_diffexp)
  colnames(design)
  ## notice that the first level of TF group like Ctr is taking as base automatically.
  #colnames(design)[startsWith(colnames(design), "x")] <- "x"
  dge <- edgeR::estimateDisp(dge, design = design)# maximize negative binomial likelihood 
  
  
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
  
  # build the glm fit function input into mclapply
  fit_glm_coef <- function(Coef, Dge, Design){
    fit <- edgeR::glmQLFit(Dge, design = Design)
    # for a given coefficient, conduct genewise statistical test
    # glmQLFTest is similar to glmLRT except that it replaces likelihood ratio tests with empirical Bayes quasi-likelihood F-tests. 
    # use QL F-tests instead of the more usual likelihood ratio tests (LRT) as they give stricter error rate control by accounting for the uncertainty in dispersion estimation:
    # The p-values from glmQLFTest are always greater than or equal to those that would be obtained from glmLRT using the same negative binomial dispersions.
    # argument coef: integer or character index vector indicating which coefficients of the linear model are to be tested equal to zero.
    qlf <- edgeR::glmQLFTest(fit, coef = Coef)
    # use topTags to extract the most differentially expressed genes/tags 
    # ranked either by p-value or absolute log foldchange of coefficients 
    tt <- edgeR::topTags(qlf, n = Inf) 
    score <- as.data.frame(tt$table)
    score$gene <- rownames(score)
    # score <- cbind(score, qlf$coefficients[rownames(score), ])
    return(score)
  }
  
  score1 <- fit_glm_coef("high1TRUE", dge, design)
  score2 <- fit_glm_coef("high2TRUE", dge, design)
  score12 <- fit_glm_coef("high1TRUE:high2TRUE", dge, design)
  
  score1$significant <- score1$FDR < 0.05
  score2$significant <- score2$FDR < 0.05
  score12$significant <- score12$FDR < 0.05
  
  scores <- score1 %>% left_join(score2, by = "gene", suffix = c("_1", "_2")) %>% left_join(score12 %>% rename_at(vars(-one_of("gene")), ~paste0(., "_12")), by = "gene")
  # scores$symbol <- data.annot[scores$gene, "gene_short_name"]
  scores
}
