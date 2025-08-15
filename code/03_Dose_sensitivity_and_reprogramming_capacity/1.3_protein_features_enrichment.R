# title: "Features enriched for TF capacity"
# author: "Wangjie Liu"
# date: "2024/10/23"

setwd("./")

library(tidyverse)
library(ggplot2)
library(dplyr)


# function
convertToLower_except_firstLetter <- function(string) {
  first_letter <- substr(string, 1, 1)  # Extract the first letter
  rest_of_string <- substr(string, 2, nchar(string))  # Extract the rest of the string
  converted_string <- paste0(toupper(first_letter), tolower(rest_of_string))  # Convert the rest of the string to lowercase and concatenate with the uppercase first letter
  return(converted_string)
}



convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  return(genesV2)
}


##------------------------------------- enrichment for protein feature (Figure 4)
# load and preprocess TF categories
TF.all <- read.csv("data/TF_categories_on_potency_capacity_dosealigned.csv", header = T)
rownames(TF.all) <- TF.all$X
TF.all$X <- NULL
colnames(TF.all)

TF.all$TF[TF.all$TF == "T"] <- "Tbxt" # update the TF name
rownames(TF.all)[rownames(TF.all) == "T"] <- "Tbxt"

TFs.potent <- TF.all[TF.all$category %in% c("high-capacity & high-sensitivity", "high-capacity & low-sensitivity"),]$TF
TFs.rest <- TF.all$TF[!TF.all$TF %in% TFs.potent] 

# prepare mouse protein features
mouse_protein.final <- read.csv("data/mouse_protein_mapped_by_uniprotIDmapping.csv")
mouse_protein.final$X <- NULL
mouse_protein.final <- mouse_protein.final %>% filter(!is.na(gene_symbol))


# add TF category 
mouse_protein.final$TF <- mouse_protein.final$gene_symbol
mouse_protein.final <- left_join(mouse_protein.final, TF.all[, c("TF","category")], by="TF")
mouse_protein.final <- mouse_protein.final[!is.na(mouse_protein.final$category),] 

# map by uniprot ID mapping:
# 61 TFs are missing. Need to add manually from mouse_protein.csv
TF.miss <- TF.all$TF[! TF.all$TF %in% unique(mouse_protein.final$TF)]
TF.miss
TF.addit <- mouse_protein$protein_name[mouse_protein$protein_name %in% paste0(toupper(TF.miss),"_MOUSE")]
TF.addit
TF.addit <- c(TF.addit, "TBXT_MOUSE")

mouse_protein.addit <- mouse_protein[mouse_protein$protein_name %in% TF.addit,]
mouse_protein.addit$TF <- sub("_.*","",mouse_protein.addit$protein_name) %>% convertToLower_except_firstLetter()
mouse_protein.addit$category <- TF.all[mouse_protein.addit$TF, "category"] 

mouse_protein.final$gene_symbol <- NULL
colnames(mouse_protein.final)[97:ncol(mouse_protein.final)]
colnames(mouse_protein.addit)[97:ncol(mouse_protein.addit)]

mouse_protein.final <- rbind(mouse_protein.final, mouse_protein.addit)
dim(mouse_protein.final)

mouse_protein.final$Capacity <- "low"
mouse_protein.final$Capacity[mouse_protein.final$TF %in% TFs.potent] <- "high"


index<-  sapply(mouse_protein.final, function(x) !(is.numeric(x)))
index[which(index == TRUE)]

names(index) [index]


## remove the conditions that are not to be tested. 
condition.all <- names(mouse_protein.final)
condition.all <- condition.all[ !(condition.all %in%  names(index) [index]  )] 

# check if there is na in each condition
lapply(condition.all, function(x){
  value.na <- is.na(mouse_protein.final[,x]) %>% table()
  value.na <- value.na["TRUE"]
  if (isTRUE(value.na > 0)){
    print(x)
    # print(is.na(mouse_protein.final[,x]) %>% table())
  }
})


plist <- lapply(condition.all, function(x) {
  wilcox_test <- wilcox.test(mouse_protein.final[mouse_protein.final$Capacity == "high",x] ,   mouse_protein.final[ mouse_protein.final$Capacity == "low",x ] )
  print(formatC(wilcox_test$p.value, format = "e", digits = 3)  )
  p <- ggplot(mouse_protein.final, aes_string("Capacity", x, color = "Capacity")) + 
    geom_jitter(width = 0.2, size =0.3)
  # use geom_crossbar()
  return(p)
})
names(plist) <- condition.all


# get p value
list.test <- lapply(condition.all, function(x) {
  wilcox_test <- wilcox.test(mouse_protein.final[mouse_protein.final$Capacity == "high",x] ,   mouse_protein.final[ mouse_protein.final$Capacity == "low",x ] )
  pValue <- wilcox_test$p.value
  names(pValue) <- x
  return(pValue)
})
pValue_all <- unlist(list.test)
# FDR correction
p.adjust_all <- p.adjust(pValue_all, method = "fdr")


p.cutoff <- 5e-2

plist1<-  lapply(condition.all, function(x)  {
  if (p.adjust_all[x] < p.cutoff) {
    print(x)
    
    # print the details of  Wilcoxon Rank Sum test
    wilcox_test <- wilcox.test(mouse_protein.final[mouse_protein.final$Capacity == "high",x] ,   mouse_protein.final[ mouse_protein.final$Capacity == "low",x ] )
    print(wilcox_test)
    
    # prepare plot
    p <- ggplot(mouse_protein.final, aes_string("Capacity", x, color = "Capacity")) + 
      geom_jitter(width = 0.2, size =0.3)+
      cowplot::theme_cowplot()
    # use geom_crossbar()
    p <- p + 
      stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="crossbar", width=0.5, linewidth = 0.3) +
      ggtitle(paste("p.adj =", formatC(p.adjust_all[x], format = "e", digits = 3))) 
    return(p)
  }
} )

length(plist1)
length(plist1[!sapply(plist1, is.null)])
plist1 <-  plist1[!sapply(plist1, is.null)]
p <- plot_grid(plotlist = plist1)
p
ggsave("results/Mouse_Capable_TFs_high_VS_low_capacity_cutoff0.05_mapped_by_uniprot_IDmapping.pdf", 
       plot = p, 
       height = 14, width = 17, useDingbats =F)


##------------------------------------- enrichment for genetic constraints
# load and preprocess TF categories
TF_cate <- read.csv("data/TF_categories_on_potency_capacity_dosealigned.csv")
rownames(TF_cate) <- TF_cate$X
TF_cate$X <- NULL
TF_cate <- TF_cate[, c("TF","category")]
table(TF_cate$category)

# load genetic constraint data
constraints <- read.table("data/gnomad.v4.0.constraint_metrics.tsv", header = F, quote = "")
colnames(constraints) <- constraints[1,]
constraints <-  constraints[-1,]


# humanx <- unique(genesV2[, 2])
convert_result <- convertMouseGeneList(TF_cate$TF)
# ‘Atf7’, ‘E430018J23Rik’, ‘Pou5f1’, ‘Tgif2lx1’, ‘Zfp78’, ‘Zfp811’, ‘Zfp87’ 
# checked manually from NCBI gene:
# Atf7 -> ATF7
#  E430018J23Rik -> ZNF764 
# Pou5f1 -> POU5F1
# Tgif2lx1 -> TGIF2LX;TGIF2LY (choose one)
# Zfp78 -> ZNF337 (ensembl); ZNF471, ZNF470 (orthoDB) 
# Zfp811 -> ZNF44; ZNF442;ZNF563 (choose one)
# Zfp87 -> ZNF429 
convert_result[convert_result$MGI.symbol == "E430018J23Rik",] 


convert_result <- convert_result[-c(which(convert_result$MGI.symbol %in% c("Atf7",
                                                                           "E430018J23Rik",
                                                                           "Pou5f1",
                                                                           "Tgif2lx1",
                                                                           "Zfp78",
                                                                           "Zfp811",
                                                                           "Zfp87"))),]
convert_result <- rbind(convert_result, data.frame(MGI.symbol = c("Atf7",
                                                                  "E430018J23Rik",
                                                                  "Pou5f1",
                                                                  "Tgif2lx1",
                                                                  "Zfp78",
                                                                  "Zfp811",
                                                                  "Zfp87"),
                                                   HGNC.symbol =c("ATF7",
                                                                  "ZNF764",
                                                                  "POU5F1",
                                                                  "TGIF2LX",
                                                                  "ZNF471",
                                                                  "ZNF44",
                                                                  "ZNF429"))
)
# collect manually from NCBI/orthoDB 
rownames(convert_result) <- convert_result$MGI.symbol

TF_cate$TF_NewName <- convert_result[TF_cate$TF,"HGNC.symbol"]

# manually correct these based on NCBI gene:
TF_cate[TF_cate$TF == "Zfp93","TF_NewName"] <- "ZNF235"
TF_cate[TF_cate$TF == "Zfp729a","TF_NewName"] <- "ZNF91"
TF_cate[TF_cate$TF == "Tbx10","TF_NewName"] <- "TBX10"
TF_cate[TF_cate$TF == "Zfp747","TF_NewName"] <- "ZNF764"
TF_cate[TF_cate$TF == "Mycn","TF_NewName"] <- "MYCN"
TF_cate[TF_cate$TF == "Zfp707","TF_NewName"] <- "ZNF707"
TF_cate[TF_cate$TF == "Rhox12","TF_NewName"] <- "RHOXF1"
TF_cate[TF_cate$TF == "Pitx2","TF_NewName"] <- "PITX2"
TF_cate[TF_cate$TF == "Hmgn1","TF_NewName"] <- "HMGN1"
TF_cate[TF_cate$TF == "Spert","TF_NewName"] <- "CBY2"
TF_cate[TF_cate$TF == "A430033K04Rik","TF_NewName"] <- "ZNF487"
TF_cate[TF_cate$TF == "Hmgb1","TF_NewName"] <- "HMGB1"
TF_cate[TF_cate$TF == "Bmyc","TF_NewName"] <- "MYCL"

TF_cate <- TF_cate[!is.na(TF_cate$TF_NewName),]
rownames(TF_cate) <- TF_cate$TF_NewName
constraints.TFs <- subset(constraints, gene %in% c(TF_cate$TF_NewName))
table(constraints.TFs$gene) %>% length() 

# take only MANE selected transcripts
constraints.TFs <- constraints.TFs[constraints.TFs$mane_select == "true",]

# Data preparation 
constraints.TFs <- constraints.TFs[, c(
  "gene", "lof.pLI", "lof.z_score", "lof.oe_ci.upper",
  "mis.oe", "mis.z_score", "mis.oe_ci.upper"
)]
constraints.TFs <- constraints.TFs[!duplicated(constraints.TFs), ]
rownames(constraints.TFs) <- constraints.TFs$gene

# Add TF category and label
constraints.TFs$category <- TF_cate[constraints.TFs$gene, "category"]
constraints.TFs$tf       <- TF_cate[constraints.TFs$gene, "TF"]

# Identify numeric score columns
scores.oi <- setdiff(colnames(constraints.TFs), c("gene", "category", "tf"))

# Ensure numeric columns
for (i in scores.oi) {
  constraints.TFs[[i]] <- as.numeric(constraints.TFs[[i]])
}

# Plotting checks (by category) 
lapply(scores.oi, function(x) {
  ggplot(constraints.TFs, aes_string("category", x, color = "category")) +
    geom_jitter(width = 0.2, size = 0.3) +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1),
                 geom = "crossbar", width = 0.5, linewidth = 0.3) +
    scale_color_brewer(palette = "Dark2") +
    cowplot::theme_cowplot() +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle(x)
})

# Categorize TFs by capacity 
constraints.TFs$capacity <- "high"
constraints.TFs$capacity[constraints.TFs$category %in% "low-capacity"] <- "low"

#----------- LOEUF intolerance analysis 
constraints.TFs$LOF_intolerant <- ifelse(constraints.TFs$lof.oe_ci.upper < 0.6, "yes", "no")

# Contingency table
table_high <- table(constraints.TFs$LOF_intolerant[constraints.TFs$capacity == "high"])
table_low  <- table(constraints.TFs$LOF_intolerant[constraints.TFs$capacity == "low"])

data <- matrix(c(35, 36, 46, 94), nrow = 2,
               dimnames = list(c("intolerant", "tolerant"), c("high", "low")))
print(data)

# Fisher's Exact Test
fisher_res_LOEUF <- fisher.test(data)
print(fisher_res_LOEUF)

#----------- pLI intolerance analysis 
constraints.TFs$LOF_intolerant <- ifelse(constraints.TFs$lof.pLI > 0.9, "yes", "no")

# Contingency table
table_high <- table(constraints.TFs$LOF_intolerant[constraints.TFs$capacity == "high"])
table_low  <- table(constraints.TFs$LOF_intolerant[constraints.TFs$capacity == "low"])

data <- matrix(c(32, 39, 41, 99), nrow = 2,
               dimnames = list(c("intolerant", "tolerant"), c("high", "low")))
print(data)

# Fisher's Exact Test
fisher_res_pLI <- fisher.test(data)
print(fisher_res_pLI)

# Plotting by capacity --------------------------------------------------
lapply(scores.oi, function(x) {
  ggplot(constraints.TFs, aes_string("capacity", x, color = "capacity")) +
    geom_jitter(width = 0.2, size = 0.3) +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1),
                 geom = "crossbar", width = 0.5, linewidth = 0.3) +
    scale_color_brewer(palette = "Dark2") +
    mashaGgplot2Theme +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black")) +
    ggtitle(x)
})

# Specific plot: LOEUF vs capacity
p1 <- ggplot(constraints.TFs, aes(x = capacity, y = lof.oe_ci.upper, color = capacity)) +
  geom_jitter(width = 0.2, size = 1) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1),
               geom = "crossbar", width = 0.5, linewidth = 0.3) +
  mashaGgplot2Theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black")) +
  labs(subtitle = sprintf("Fisher's p-value = %.5f; OR = %.6f",
                          fisher_res_LOEUF$p.value,
                          fisher_res_LOEUF$estimate))
ggsave(p1, filename = "figure/lof.oe.ci.upper_capacity_high_VS_low.pdf",
       height = 5, width = 5)

# Specific plot: pLI vs capacity
p2 <- ggplot(constraints.TFs, aes(x = capacity, y = lof.pLI, color = capacity)) +
  geom_jitter(width = 0.2, size = 1) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1),
               geom = "crossbar", width = 0.5, linewidth = 0.3) +
  mashaGgplot2Theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black")) +
  labs(subtitle = sprintf("Fisher's p-value = %.5f; OR = %.6f",
                          fisher_res_pLI$p.value,
                          fisher_res_pLI$estimate))
ggsave(p2, filename = "figure/pLI_capacity_high_VS_low.pdf",
       height = 5, width = 5)