############################################################################
#                                                                          #
#         This script aims to :                                            #
#         1. DEG and functional genomics analyses for AT2                  #
#         2. DEG and functional genomics analyses for fibroblasts          #
#                                                                          #
############################################################################


# Libraries
library(Seurat)
library(enrichR)
library(dplyr)


# Source codes
source("/fs/ess/PCON0022/liyang/r_utilities/functions/combinatorial_tools.R")
source("/fs/ess/PCON0022/liyang/r_utilities/functions/transcriptome_tools.R")
source("/fs/ess/PCON0022/liyang/r_utilities/functions/visual_tools.R")


# Global parameters
R.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Rfiles/"
image.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Images/"
table.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Tables/"
n.AT2.DEGs <- 150
dir.create(paste0(table.dir, "AT2"))
AT2.dir <- paste0(table.dir, "AT2/")


############################################################################
#                                                                          #
#   Analyses of AT2 :                                                      #
#   1. Compare the top-150 DEGs between AT2 v.s. iAT2                      #
#   2. GO and KEGG enrichment analyses for AT2 and iAT2, respectively      #
#   3. Find overlapped GO terms and KEGG pathways and highlight the        #
#    overlapped ones                                                       #
#                                                                          #
############################################################################


# Local parameters
iAT2.DEGFile <- "COPD_healthy_DEGs_iAT2.csv"
AT2.DEGFile <- "COPD_healthy_DEGs.csv"


# Load the DEGs in iAT2
iAT2.DEGs <- read.csv(paste0(table.dir, iAT2.DEGFile))
dim(iAT2.DEGs)
colnames(iAT2.DEGs)
head(iAT2.DEGs[, 1 : 6])
nrow(iAT2.DEGs[iAT2.DEGs$cluster_name == "copd",])
nrow(iAT2.DEGs[iAT2.DEGs$cluster_name == "healthy",])
# Question : iAT2 have two cluster names, i.e., copd and healthy
# Shall we only use the copd DEGs?


# Sort the iAT2 DEGs and select the top-150
range(iAT2.DEGs$logfc)
iAT2.sortedDEGs <- iAT2.DEGs %>% dplyr::arrange(adj_pval, -logfc)
# iAT2.sortedDEGs <- iAT2.DEGs[iAT2.DEGs$adj_pval < 0.05,]
# iAT2.sortedDEGs <- iAT2.DEGs[iAT2.DEGs$pval < 0.05,]
# iAT2.sortedDEGs <- iAT2.sortedDEGs %>% dplyr::arrange(-logfc, adj_pval)
head(iAT2.sortedDEGs[, 1 : 6])
qs::qsave(iAT2.sortedDEGs, paste0(R.dir, "iAT2_DEGs_sorted.qsave"))
iAT2.top.COPD.DEGs <- iAT2.sortedDEGs[iAT2.sortedDEGs$cluster_name == "copd",] %>% 
  head(n = n.AT2.DEGs)
dim(iAT2.top.COPD.DEGs)
head(iAT2.top.COPD.DEGs[, 1 : 6])
qs::qsave(iAT2.top.COPD.DEGs, paste0(R.dir, "iAT2_DEGs_sorted_COPD.qsave"))


# Load the DEGs in AT2
AT2.DEGs <- read.csv(paste0(table.dir, AT2.DEGFile))
dim(AT2.DEGs)
colnames(AT2.DEGs)
colnames(AT2.DEGs)[1] <- "gene"
head(AT2.DEGs)
AT2.DEGs.sorted <- AT2.DEGs %>% dplyr::arrange(p_val_adj, -avg_log2FC) %>%
  head(n = n.AT2.DEGs)
# AT2.DEGs.sorted <- AT2.DEGs[AT2.DEGs$p_val_adj < 0.05,]
# AT2.DEGs.sorted <- AT2.DEGs[AT2.DEGs$p_val < 0.05,]
# AT2.DEGs.sorted <- AT2.DEGs.sorted %>% dplyr::arrange(-avg_log2FC, p_val_adj) %>% 
#   head(n = n.AT2.DEGs)
identical(AT2.DEGs, AT2.DEGs.sorted)
dim(AT2.DEGs.sorted)
head(AT2.DEGs.sorted)
range(AT2.DEGs.sorted$p_val_adj)
qs::qsave(AT2.DEGs.sorted, paste0(R.dir, "AT2_DEGs_sorted_COPD.qsave"))


# Find overlapped DEGs between AT2 and iAT2
AT2.DEGs.overlapped <- intersect(iAT2.top.COPD.DEGs$gene, 
                                 AT2.DEGs.sorted$gene)
length(AT2.DEGs.overlapped)
write.csv(iAT2.top.COPD.DEGs[which(iAT2.top.COPD.DEGs$gene %in% AT2.DEGs.overlapped),], 
          paste0(AT2.dir, "iAT2_DEGs_overlapped.csv"), 
          quote = F, row.names = F)
write.csv(AT2.DEGs.sorted[which(AT2.DEGs.sorted$gene %in% AT2.DEGs.overlapped),], 
          paste0(AT2.dir, "AT2_DEGs_overlapped.csv"), 
          quote = F, row.names = F)


# Since there are only sixteen overlapped DEGs, I will calculate the ranks of each 
# overlapped DEG among iAT2 and AT2 DEGs


# Get the vector of ranks of overlapped DEGs regarding iAT2 and AT2, respectively
# First, adjusted p-values

# iAT2
iAT2.sortedDEGs <- iAT2.DEGs %>% dplyr::arrange(adj_pval, -logfc)
iAT2.sorted.padj <- iAT2.sortedDEGs[iAT2.sortedDEGs$cluster_name == "copd", "gene"]
head(iAT2.sortedDEGs[iAT2.sortedDEGs$cluster_name == "copd", 1 : 6])
head(iAT2.sorted.padj)
range(iAT2.sortedDEGs$adj_pval)
range(iAT2.sortedDEGs$pval)

# AT2
AT2.sorted.padj <- AT2.DEGs %>% dplyr::arrange(p_val_adj, -avg_log2FC) %>% pull(gene)
head(AT2.DEGs %>% dplyr::arrange(p_val_adj, -avg_log2FC))
head(AT2.sorted.padj)
range(AT2.DEGs %>% dplyr::arrange(p_val_adj, -avg_log2FC) %>% pull(p_val_adj))
range(AT2.DEGs %>% dplyr::arrange(p_val_adj, -avg_log2FC) %>% pull(p_val))

# Calculate the overlap ratios
AT2.padj.overlaps <- get_overlap_ratio(x = iAT2.sorted.padj, y = AT2.sorted.padj)
# The maximum overlap ratio, 0.108695652173913, is achieved at 138.

get_lines(data = data.frame(seq_along(AT2.padj.overlaps), AT2.padj.overlaps), 
          path = paste0(image.dir, "AT2_padj_overlap_ratios.png"))

iAT2.overlapGenes <- iAT2.sorted.padj[1 : 138]
AT2.overlapGenes <- AT2.sorted.padj[1 : 138]
length(intersect(iAT2.overlapGenes, AT2.overlapGenes)) / 138


# Second, average log2 fold-change (log2-FC)

# iAT2
iAT2.sortedDEGs <- iAT2.DEGs %>% dplyr::arrange(-logfc, adj_pval)
iAT2.sorted.logFC <- iAT2.sortedDEGs[iAT2.sortedDEGs$cluster_name == "copd", "gene"]
head(iAT2.sortedDEGs[iAT2.sortedDEGs$cluster_name == "copd", 1 : 6])
length(iAT2.sorted.logFC)
head(iAT2.sorted.logFC)

# AT2
AT2.sorted.logFC <- AT2.DEGs %>% dplyr::arrange(-avg_log2FC, p_val_adj) %>% pull(gene)
head(AT2.DEGs %>% dplyr::arrange(-avg_log2FC, p_val_adj))
length(AT2.sorted.logFC)
head(AT2.sorted.logFC)

# Calculate the overlap ratios
AT2.logFC.overlaps <- get_overlap_ratio(x = iAT2.sorted.logFC, 
                                        y = AT2.sorted.logFC)
# The maximum overlap ratio, 0, is achieved at 1.


# Third, filter DEGs using adjusted p-values and then sort them via logFC

# iAT2
iAT2.filteredDEGs <- iAT2.DEGs[iAT2.DEGs$adj_pval < 0.05,]
iAT2.sortedDEGs <- iAT2.filteredDEGs %>% dplyr::arrange(-logfc, adj_pval)
iAT2.sorted.logFC <- iAT2.sortedDEGs[iAT2.sortedDEGs$cluster_name == "copd", "gene"]
head(iAT2.sortedDEGs[iAT2.sortedDEGs$cluster_name == "copd", 1 : 6])
length(iAT2.sorted.logFC)
head(iAT2.sorted.logFC)

# AT2
AT2.filteredDEGs <- AT2.DEGs[AT2.DEGs$p_val_adj < 0.05,]
AT2.sorted.logFC <- AT2.filteredDEGs %>% dplyr::arrange(-avg_log2FC, p_val_adj) %>% pull(gene)
head(AT2.filteredDEGs %>% dplyr::arrange(-avg_log2FC, p_val_adj))
length(AT2.sorted.logFC)
head(AT2.sorted.logFC)

# Calculate the overlap ratios
AT2.filtered.logFC.overlaps <- get_overlap_ratio(x = iAT2.sorted.logFC, 
                                        y = AT2.sorted.logFC)
# The maximum overlap ratio, 0.0552147239263804, is achieved at 163.

get_lines(data = data.frame(seq_along(AT2.filtered.logFC.overlaps), 
                            AT2.filtered.logFC.overlaps), 
          path = paste0(image.dir, "AT2_filter_padj_overlap_ratios.png"))


# Fourth, filter DEGs using p-values and then sort them via logFC

# iAT2
iAT2.filteredDEGs <- iAT2.DEGs[iAT2.DEGs$pval < 0.05,]
iAT2.sortedDEGs <- iAT2.filteredDEGs %>% dplyr::arrange(-logfc, adj_pval)
iAT2.sorted.logFC <- iAT2.sortedDEGs[iAT2.sortedDEGs$cluster_name == "copd", "gene"]
head(iAT2.sortedDEGs[iAT2.sortedDEGs$cluster_name == "copd", 1 : 6])
length(iAT2.sorted.logFC)
head(iAT2.sorted.logFC)

# AT2
AT2.filteredDEGs <- AT2.DEGs[AT2.DEGs$p_val < 0.05,]
AT2.sorted.logFC <- AT2.filteredDEGs %>% dplyr::arrange(-avg_log2FC, p_val_adj) %>% pull(gene)
head(AT2.filteredDEGs %>% dplyr::arrange(-avg_log2FC, p_val_adj))
length(AT2.sorted.logFC)
head(AT2.sorted.logFC)

# Calculate the overlap ratios
AT2.pval.logFC.overlaps <- get_overlap_ratio(x = iAT2.sorted.logFC, 
                                                 y = AT2.sorted.logFC)
# The maximum overlap ratio, 0.0416666666666667, is achieved at 168.
# Hence, we use the first plan, i.e., sort DEGs in increasing order of adjusted p-values

get_lines(data = data.frame(seq_along(AT2.pval.logFC.overlaps), 
                            AT2.pval.logFC.overlaps), 
          path = paste0(image.dir, "AT2_filter_pval_overlap_ratios.png"))


# Functional genomics analyses
iAT2.top.COPD.DEGs <- qs::qread(paste0(R.dir, "iAT2_DEGs_sorted_COPD.qsave"))
colnames(iAT2.top.COPD.DEGs)
head(iAT2.top.COPD.DEGs)
dim(iAT2.top.COPD.DEGs)
range(iAT2.top.COPD.DEGs$adj_pval) # 8.310000e-88 2.361047e-03
range(iAT2.top.COPD.DEGs$pval) # 3.28000e-91 2.51317e-04

AT2.top.COPD.DEGs <- qs::qread(paste0(R.dir, "AT2_DEGs_sorted_COPD.qsave"))
colnames(AT2.top.COPD.DEGs)
dim(AT2.top.COPD.DEGs)
range(AT2.top.COPD.DEGs$p_val_adj) # 2.828897e-284  5.335740e-11
range(AT2.top.COPD.DEGs$p_val) # 1.414449e-287  2.667870e-14

length(intersect(iAT2.top.COPD.DEGs$gene, AT2.top.COPD.DEGs$gene)) / 150

# Enrichment of iAT2
iAT2.enriched <- run_GO_and_KEGG(genes.ll = iAT2.top.COPD.DEGs$gene) # 250, 145, 1251, 147
AT2.enriched <- run_GO_and_KEGG(genes.ll = AT2.top.COPD.DEGs$gene) # 190, 139, 1025, 169
qs::qsave(iAT2.enriched, paste0(R.dir, "iAT2_enriched.qsave"))
qs::qsave(AT2.enriched, paste0(R.dir, "AT2_enriched.qsave"))

lapply(iAT2.enriched, function(x) {
  range(x$Adjusted.P.value)
})
# $GO_Molecular_Function_2018
# [1] 1.831069e-14 9.800072e-01
# 
# $GO_Cellular_Component_2018
# [1] 9.298868e-18 9.033578e-01
# 
# $GO_Biological_Process_2018
# [1] 7.176100e-21 9.721389e-01
# 
# $KEGG_2019_Human
# [1] 3.024941e-17 8.273958e-01

lapply(iAT2.enriched, function(x) {
  range(x$P.value)
})
# $GO_Molecular_Function_2018
# [1] 7.324276e-17 9.800072e-01
# 
# $GO_Cellular_Component_2018
# [1] 7.692639e-20 9.033578e-01
# 
# $GO_Biological_Process_2018
# [1] 5.736291e-24 9.721389e-01
# 
# $KEGG_2019_Human
# [1] 2.057783e-19 8.273958e-01

lapply(AT2.enriched, function(x) {
  range(x$Adjusted.P.value)
})
# $GO_Molecular_Function_2018
# [1] 6.856352e-23 9.602962e-01
# 
# $GO_Cellular_Component_2018
# [1] 1.346003e-69 9.790588e-01
# 
# $GO_Biological_Process_2018
# [1] 7.287636e-78 9.919876e-01
# 
# $KEGG_2019_Human
# [1] 3.949115e-70 8.273958e-01

lapply(AT2.enriched, function(x) {
  range(x$P.value)
})
# $GO_Molecular_Function_2018
# [1] 3.608606e-25 9.602962e-01
# 
# $GO_Cellular_Component_2018
# [1] 9.683474e-72 9.790588e-01
# 
# $GO_Biological_Process_2018
# [1] 7.109888e-81 9.919876e-01
# 
# $KEGG_2019_Human
# [1] 2.336754e-72 8.273958e-01


# GO molecular functions using all terms (107)
dim(iAT2.enriched$GO_Molecular_Function_2018)
dim(AT2.enriched$GO_Molecular_Function_2018)
AT2.overlapped.MF <- intersect(unique(iAT2.enriched$GO_Molecular_Function_2018$Term), 
                               unique(AT2.enriched$GO_Molecular_Function_2018$Term))
length(AT2.overlapped.MF) # 107


# GO molecular functions using significant terms (padj) : 4 (use this)
iAT2.sig.MF <- iAT2.enriched$GO_Molecular_Function_2018
iAT2.sig.MF <- iAT2.sig.MF[iAT2.sig.MF$Adjusted.P.value < 0.05,]
dim(iAT2.sig.MF)

AT2.sig.MF <- AT2.enriched$GO_Molecular_Function_2018
AT2.sig.MF <- AT2.sig.MF[AT2.sig.MF$Adjusted.P.value < 0.05,]
dim(AT2.sig.MF)

AT2.padj.overlapMF <- intersect(unique(iAT2.sig.MF$Term), unique(AT2.sig.MF$Term))
length(AT2.padj.overlapMF) # 4
length(AT2.padj.overlapMF) / 
  nrow(iAT2.sig.MF) # 0.2105263
length(AT2.padj.overlapMF) / 
  nrow(AT2.sig.MF) # 0.5714286

iAT2.padj.MF <- cbind(Overlap = iAT2.sig.MF$Term %in% AT2.padj.overlapMF, 
                      iAT2.sig.MF)
head(iAT2.padj.MF[, 1 : 5])
write.csv(iAT2.padj.MF, paste0(AT2.dir, "iAT2_padj_MF.csv"), 
          quote = F, row.names = F)

AT2.padj.MF <- cbind(Overlap = AT2.sig.MF$Term %in% AT2.padj.overlapMF, 
                      AT2.sig.MF)
head(AT2.padj.MF[, 1 : 5])
write.csv(AT2.padj.MF, paste0(AT2.dir, "AT2_padj_MF.csv"), 
          quote = F, row.names = F)


# GO molecular functions using significant terms (pval) : 6
iAT2.sig.MF <- iAT2.enriched$GO_Molecular_Function_2018
iAT2.sig.MF <- iAT2.sig.MF[iAT2.sig.MF$P.value < 0.05,]
dim(iAT2.sig.MF)

AT2.sig.MF <- AT2.enriched$GO_Molecular_Function_2018
AT2.sig.MF <- AT2.sig.MF[AT2.sig.MF$P.value < 0.05,]
dim(AT2.sig.MF)

AT2.overlapped.MF <- intersect(unique(iAT2.sig.MF$Term), unique(AT2.sig.MF$Term))
length(AT2.overlapped.MF) # 6
length(AT2.overlapped.MF) / nrow(iAT2.sig.MF) # 0.1132075
length(AT2.overlapped.MF) / nrow(AT2.sig.MF) # 0.1428571

iAT2.pval.MF <- cbind(Overlap = iAT2.sig.MF$Term %in% AT2.overlapped.MF, 
                      iAT2.sig.MF)
head(iAT2.pval.MF[, 1 : 5])
write.csv(iAT2.pval.MF, paste0(AT2.dir, "iAT2_pval_MF.csv"), 
          quote = F, row.names = F)


# GO cellular components using all terms (100)
dim(iAT2.enriched$GO_Cellular_Component_2018)
dim(AT2.enriched$GO_Cellular_Component_2018)
AT2.overlapped.CC <- intersect(unique(iAT2.enriched$GO_Cellular_Component_2018$Term), 
                               unique(AT2.enriched$GO_Cellular_Component_2018$Term))
length(AT2.overlapped.CC) # 100


# GO cellular component using significant terms (padj) : 4 (use this)
iAT2.sig.CC <- iAT2.enriched$GO_Cellular_Component_2018
iAT2.sig.CC <- iAT2.sig.CC[iAT2.sig.CC$Adjusted.P.value < 0.05,]
dim(iAT2.sig.CC)

AT2.sig.CC <- AT2.enriched$GO_Cellular_Component_2018
AT2.sig.CC <- AT2.sig.CC[AT2.sig.CC$Adjusted.P.value < 0.05,]
dim(AT2.sig.CC)

AT2.padj.overlapCC <- intersect(unique(iAT2.sig.CC$Term), unique(AT2.sig.CC$Term))
length(AT2.padj.overlapCC) # 11
length(AT2.padj.overlapCC) / 
  nrow(iAT2.sig.CC) # 0.2820513
length(AT2.padj.overlapCC) / 
  nrow(AT2.sig.CC) # 0.4782609

iAT2.padj.CC <- cbind(Overlap = iAT2.sig.CC$Term %in% AT2.padj.overlapCC, 
                      iAT2.sig.CC)
head(iAT2.padj.CC[, 1 : 5])
write.csv(iAT2.padj.CC, paste0(AT2.dir, "iAT2_padj_CC.csv"), 
          quote = F, row.names = F)

AT2.padj.CC <- cbind(Overlap = AT2.sig.CC$Term %in% AT2.padj.overlapCC, 
                     AT2.sig.CC)
head(AT2.padj.CC[, 1 : 5])
write.csv(AT2.padj.CC, paste0(AT2.dir, "AT2_padj_CC.csv"), 
          quote = F, row.names = F)


# GO cellular component using significant terms (pval) : 55
iAT2.sig.CC <- iAT2.enriched$GO_Cellular_Component_2018
iAT2.sig.CC <- iAT2.sig.CC[iAT2.sig.CC$P.value < 0.05,]
dim(iAT2.sig.CC) # 55

AT2.sig.CC <- AT2.enriched$GO_Cellular_Component_2018
AT2.sig.CC <- AT2.sig.CC[AT2.sig.CC$P.value < 0.05,]
dim(AT2.sig.CC) # 30

AT2.overlapped.CC <- intersect(unique(iAT2.sig.CC$Term), unique(AT2.sig.CC$Term))
length(AT2.overlapped.CC) # 17
length(AT2.overlapped.CC) / nrow(iAT2.sig.CC) # 0.3090909
length(AT2.overlapped.CC) / nrow(AT2.sig.CC) # 0.5666667

iAT2.pval.CC <- cbind(Overlap = iAT2.sig.CC$Term %in% AT2.overlapped.CC, 
                      iAT2.sig.CC)
head(iAT2.pval.CC[, 1 : 5])
write.csv(iAT2.pval.CC, paste0(AT2.dir, "iAT2_pval_CC.csv"), 
          quote = F, row.names = F)
AT2.pval.CC <- cbind(Overlap = AT2.sig.CC$Term %in% AT2.overlapped.CC, 
                      AT2.sig.CC)
head(AT2.pval.CC[, 1 : 5])
write.csv(AT2.pval.CC, paste0(AT2.dir, "AT2_pval_CC.csv"), 
          quote = F, row.names = F)


# GO Biological Process using significant terms (padj)
iAT2.sig.BP <- iAT2.enriched$GO_Biological_Process_2018
iAT2.sig.BP <- iAT2.sig.BP[iAT2.sig.BP$Adjusted.P.value < 0.05,]
dim(iAT2.sig.BP) # 84

AT2.sig.BP <- AT2.enriched$GO_Biological_Process_2018
AT2.sig.BP <- AT2.sig.BP[AT2.sig.BP$Adjusted.P.value < 0.05,]
dim(AT2.sig.BP) # 165

AT2.padj.overlapBP <- intersect(unique(iAT2.sig.BP$Term), unique(AT2.sig.BP$Term))
length(AT2.padj.overlapBP) # 28
length(AT2.padj.overlapBP) / 
  nrow(iAT2.sig.BP) # 0.3333333
length(AT2.padj.overlapBP) / 
  nrow(AT2.sig.BP) # 0.169697

iAT2.padj.BP <- cbind(Overlap = iAT2.sig.BP$Term %in% AT2.padj.overlapBP, 
                      iAT2.sig.BP)
head(iAT2.padj.BP[, 1 : 5])
write.csv(iAT2.padj.BP, paste0(AT2.dir, "iAT2_padj_BP.csv"), 
          quote = F, row.names = F)

AT2.padj.BP <- cbind(Overlap = AT2.sig.BP$Term %in% AT2.padj.overlapBP, 
                     AT2.sig.BP)
head(AT2.padj.BP[, 1 : 5])
write.csv(AT2.padj.BP, paste0(AT2.dir, "AT2_padj_BP.csv"), 
          quote = F, row.names = F)


# GO Biological Process using significant terms (pval)
iAT2.sig.BP <- iAT2.enriched$GO_Biological_Process_2018
iAT2.sig.BP <- iAT2.sig.BP[iAT2.sig.BP$P.value < 0.05,]
dim(iAT2.sig.BP) # 267

AT2.sig.BP <- AT2.enriched$GO_Biological_Process_2018
AT2.sig.BP <- AT2.sig.BP[AT2.sig.BP$P.value < 0.05,]
dim(AT2.sig.BP) # 305

AT2.pval.overlapBP <- intersect(unique(iAT2.sig.BP$Term), unique(AT2.sig.BP$Term))
length(AT2.pval.overlapBP) # 64
length(AT2.pval.overlapBP) / 
  nrow(iAT2.sig.BP) # 0.2397004
length(AT2.pval.overlapBP) / 
  nrow(AT2.sig.BP) # 0.2098361

iAT2.pval.BP <- cbind(Overlap = iAT2.sig.BP$Term %in% AT2.pval.overlapBP, 
                      iAT2.sig.BP)
head(iAT2.pval.BP[, 1 : 5])
write.csv(iAT2.pval.BP, paste0(AT2.dir, "iAT2_pval_BP.csv"), 
          quote = F, row.names = F)

AT2.pval.BP <- cbind(Overlap = AT2.sig.BP$Term %in% AT2.pval.overlapBP, 
                     AT2.sig.BP)
head(AT2.pval.BP[, 1 : 5])
write.csv(AT2.pval.BP, paste0(AT2.dir, "AT2_pval_BP.csv"), 
          quote = F, row.names = F)


# KEGG using significant pathways (padj)
iAT2.sig.KEGG <- iAT2.enriched$KEGG_2019_Human
iAT2.sig.KEGG <- iAT2.sig.KEGG[iAT2.sig.KEGG$Adjusted.P.value < 0.05,]
dim(iAT2.sig.KEGG) # 15

AT2.sig.KEGG <- AT2.enriched$KEGG_2019_Human
AT2.sig.KEGG <- AT2.sig.KEGG[AT2.sig.KEGG$Adjusted.P.value < 0.05,]
dim(AT2.sig.KEGG) # 21

AT2.padj.overlapKEGG <- intersect(unique(iAT2.sig.KEGG$Term), unique(AT2.sig.KEGG$Term))
length(AT2.padj.overlapKEGG) # 2
length(AT2.padj.overlapKEGG) / 
  nrow(iAT2.sig.KEGG) # 0.1333333
length(AT2.padj.overlapKEGG) / 
  nrow(AT2.sig.KEGG) # 0.0952381

iAT2.padj.KEGG <- cbind(Overlap = iAT2.sig.KEGG$Term %in% AT2.padj.overlapKEGG, 
                      iAT2.sig.KEGG)
head(iAT2.padj.KEGG[, 1 : 5])
write.csv(iAT2.padj.KEGG, paste0(AT2.dir, "iAT2_padj_KEGG.csv"), 
          quote = F, row.names = F)

AT2.padj.KEGG <- cbind(Overlap = AT2.sig.KEGG$Term %in% AT2.padj.overlapKEGG, 
                     AT2.sig.KEGG)
head(AT2.padj.KEGG[, 1 : 5])
write.csv(AT2.padj.KEGG, paste0(AT2.dir, "AT2_padj_KEGG.csv"), 
          quote = F, row.names = F)


# KEGG using significant pathways (pval)
iAT2.sig.KEGG <- iAT2.enriched$KEGG_2019_Human
iAT2.sig.KEGG <- iAT2.sig.KEGG[iAT2.sig.KEGG$P.value < 0.05,]
dim(iAT2.sig.KEGG) # 35

AT2.sig.KEGG <- AT2.enriched$KEGG_2019_Human
AT2.sig.KEGG <- AT2.sig.KEGG[AT2.sig.KEGG$P.value < 0.05,]
dim(AT2.sig.KEGG) # 48

AT2.pval.overlapKEGG <- intersect(unique(iAT2.sig.KEGG$Term), unique(AT2.sig.KEGG$Term))
length(AT2.pval.overlapKEGG) # 8
length(AT2.pval.overlapKEGG) / 
  nrow(iAT2.sig.KEGG) # 0.2285714
length(AT2.pval.overlapKEGG) / 
  nrow(AT2.sig.KEGG) # 0.1666667

iAT2.pval.KEGG <- cbind(Overlap = iAT2.sig.KEGG$Term %in% AT2.pval.overlapKEGG, 
                        iAT2.sig.KEGG)
head(iAT2.pval.KEGG[, 1 : 5])
write.csv(iAT2.pval.KEGG, paste0(AT2.dir, "iAT2_pval_KEGG.csv"), 
          quote = F, row.names = F)

AT2.pval.KEGG <- cbind(Overlap = AT2.sig.KEGG$Term %in% AT2.pval.overlapKEGG, 
                       AT2.sig.KEGG)
head(AT2.pval.KEGG[, 1 : 5])
write.csv(AT2.pval.KEGG, paste0(AT2.dir, "AT2_pval_KEGG.csv"), 
          quote = F, row.names = F)
