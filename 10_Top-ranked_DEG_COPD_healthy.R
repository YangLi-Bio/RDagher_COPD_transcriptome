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


# Global parameters
R.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Rfiles/"
image.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Images/"
table.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Tables/"
n.AT2.DEGs <- 150


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
          paste0(table.dir, "iAT2_DEGs_overlapped.csv"), 
          quote = F, row.names = F)
write.csv(AT2.DEGs.sorted[which(AT2.DEGs.sorted$gene %in% AT2.DEGs.overlapped),], 
          paste0(table.dir, "AT2_DEGs_overlapped.csv"), 
          quote = F, row.names = F)


# Since there are only sixteen overlapped DEGs, I will calculate the ranks of each 
# overlapped DEG among iAT2 and AT2 DEGs


# Get the vector of ranks of overlapped DEGs regarding iAT2 and AT2, respectively

