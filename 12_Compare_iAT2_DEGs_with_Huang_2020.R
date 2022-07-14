###################################################################
#                                                                 #
# The following tasks will be conducted:                          #
# 1. Extract DEGs in mock (Figs. 3E and 3F) from reference and    #
#    compare them with our iAT2 DEGs (control)                    #
# 2. Extract DEGs in 1 and 4 dpi (Figs. 3E amd 3F) from           #
#    reference and compare them with our iAT2 DEGs (1, and 3 dpi) #
#                                                                 #
###################################################################


# Libraries
library(Seurat)
library(dplyr)
library(readxl)
library(ggvenn)
library(fgsea)
library(qusage)
library(org.Hs.eg.db)
library(msigdbr)


# Global variables
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
R.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Rfiles/"
image.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Images/"
table.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Tables/"
iAT2.file <- "Supplemental_7._Figure_5E_DEG_results.xlsx"
pval.cutoff <- 0.05
logFC.cutoff <- 0.25
FC.cutoff <- 0.25


# Source codes
source(paste0(tool.dir, "transcriptome_tools.R"))
source(paste0(tool.dir, "visual_tools.R"))
source(paste0(tool.dir, "combinatorial_tools.R"))


# Identify DEGs between 1dpi v.s. Mock in iAT2 of our paper
DEGs.iAT2.1dpi <- read_excel(paste0(table.dir, 
                                    "iAT2_DEGs/Supplemental_7._Figure_5E_DEG_results.xlsx"), 
                            sheet = 2)
dim(DEGs.iAT2.1dpi)
head(DEGs.iAT2.1dpi)
length(unique(DEGs.iAT2.1dpi$hgnc_symbols))
range(DEGs.iAT2.1dpi$PValue)
range(DEGs.iAT2.1dpi$FC)
DEGs.iAT2.1dpi <- DEGs.iAT2.1dpi[DEGs.iAT2.1dpi$FC > logFC.cutoff & 
                                   DEGs.iAT2.1dpi$PValue < pval.cutoff & 
                                   !is.na(DEGs.iAT2.1dpi$hgnc_symbols),]
dim(DEGs.iAT2.1dpi)


# Identify DEGs between 1dpi v.s. Mock in AT2 of Kotton paper
DEGs.AT2 <- read_excel(paste0(table.dir, "[AT2_DEGs]_Jessie_Huang_et_al_2020.xlsx"), sheet = 1)
dim(DEGs.AT2)
range(DEGs.AT2$LogFC.dpi1_vs_ctr)
range(DEGs.AT2$FDR.dpi1_vs_ctr)
DEGs.AT2.1dpi <- DEGs.AT2[DEGs.AT2$LogFC.dpi1_vs_ctr > 
                            FC.cutoff & DEGs.AT2$FDR.dpi1_vs_ctr < pval.cutoff,]
dim(DEGs.AT2.1dpi)


# Calculate overlaps of DEGs between iAT2 1dpi v.s. AT2 1dpi
symbols.iAT2.1dpi <- unique(DEGs.iAT2.1dpi$hgnc_symbols)
length(symbols.iAT2.1dpi)
symbols.AT2.1dpi <- unique(DEGs.AT2.1dpi$Genes.SYMBOL)
length(symbols.AT2.1dpi)
symbols.overlap1 <- intersect(symbols.iAT2.1dpi, symbols.AT2.1dpi)
length(symbols.overlap1) # 2.9% (1/35)
# Only one intersected gene is discovered!


# Identify DEGs between 3dpi v.s. Mock in iAT2 of our paper
DEGs.iAT2.3dpi <- read_excel(paste0(table.dir, 
                                    "iAT2_DEGs/Supplemental_7._Figure_5E_DEG_results.xlsx"), 
                             sheet = 4)
dim(DEGs.iAT2.3dpi)
head(DEGs.iAT2.3dpi)
length(unique(DEGs.iAT2.3dpi$hgnc_symbols))
range(DEGs.iAT2.3dpi$PValue)
range(DEGs.iAT2.3dpi$FC)
DEGs.iAT2.3dpi <- DEGs.iAT2.3dpi[DEGs.iAT2.3dpi$FC > logFC.cutoff & 
                                   DEGs.iAT2.3dpi$PValue < pval.cutoff & 
                                   !is.na(DEGs.iAT2.3dpi$hgnc_symbols),]
dim(DEGs.iAT2.3dpi)


# Identify DEGs between 4dpi v.s. Mock in AT2 of Kotton paper
# This DEG list was calculated using iDEP based on the expression matrix in GSE153277
DEGs.AT2.4dpi <- read.csv(paste0(table.dir, "DEGs_AT2_4dpi_vs_1dpi.csv"))
dim(DEGs.AT2.4dpi)
head(DEGs.AT2.4dpi)
DEGs.AT2.4dpi <- DEGs.AT2.4dpi[DEGs.AT2.4dpi$logFC > FC.cutoff & DEGs.AT2.4dpi$adj.P.Val < pval.cutoff,]
dim(DEGs.AT2.4dpi)


# Calculate overlaps of DEGs between iAT2 3dpi v.s. AT2 4dpi
symbols.iAT2.3dpi <- unique(DEGs.iAT2.3dpi$hgnc_symbols)
length(symbols.iAT2.3dpi)
symbols.AT2.4dpi <- unique(DEGs.AT2.4dpi$Symbol)
length(symbols.AT2.4dpi)
symbols.overlap3 <- intersect(symbols.iAT2.3dpi, symbols.AT2.4dpi)
length(symbols.overlap3) # 56.4% (75)
length(symbols.iAT2.3dpi) # 133
length(symbols.AT2.4dpi) # 4697


# Calculate the intersection between combination of DEGs in 1dpi and 3/4 dpi
symbols.iAT2 <- union(symbols.iAT2.1dpi, symbols.iAT2.3dpi)
symbols.AT2 <- union(symbols.AT2.1dpi, symbols.AT2.4dpi)
symbols.overlap <- intersect(symbols.iAT2, symbols.AT2)
length(symbols.overlap)
length(symbols.iAT2)
length(symbols.AT2)


# # Calculate the intersection between combination of DEGs in 1dpi and 3 dpi and infected v.s. ctrl
# DEGs.AT2.infected <- read_excel(paste0(table.dir, "[AT2_DEGs]_Jessie_Huang_et_al_2020.xlsx"), sheet = 1)
# dim(DEGs.AT2.infected)
# range(DEGs.AT2.infected$LogFC.allinfected_vs_ctr)
# range(DEGs.AT2.infected$FDR.allinfected_vs_ctr)
# DEGs.AT2.infected <- DEGs.AT2.infected[DEGs.AT2.infected$LogFC.dpi1_vs_ctr > 
#                             FC.cutoff & DEGs.AT2.infected$FDR.dpi1_vs_ctr < pval.cutoff,]
# dim(DEGs.AT2.infected)


# # Calculate overlaps of DEGs between iAT2 1, 3dpi v.s. AT2 infected
# symbols.AT2.infected <- unique(DEGs.AT2.infected$Genes.SYMBOL)
# length(symbols.AT2.infected)
# symbols.infected.overlap <- intersect(symbols.AT2.infected, symbols.iAT2)
# length(symbols.infected.overlap) # 57
# length(symbols.iAT2) # 160
# length(symbols.AT2.infected) # 2393


# Get Venn diagrams between iAT2 DEGs (1 and 3 dpi) and AT2 DEGs (1 and 4 dpi)
p.Venn.DEGs <- ggvenn(list(iAT2 = symbols.iAT2, AT2 = symbols.AT2))
save_image(p = p.Venn.DEGs, path = paste0(image.dir, 
                                          "Venn_diagram_DEGs_iAT2[1,3dpi]_AT2[1,4dpi].tiff"))

# Save the overlapped DEGs in AT2 1dpi
flag <- rep(F, nrow(DEGs.AT2.1dpi))
flag[which(DEGs.AT2.1dpi$Genes.SYMBOL %in% symbols.overlap)] <- T
write.csv(cbind(Overlap = flag, DEGs.AT2.1dpi), 
          quote = F, paste0(table.dir, "DEGs_AT2_1dpi.csv"))


# Save the overlapped DEGs in AT2 4dpi
flag <- rep(F, nrow(DEGs.AT2.4dpi))
flag[which(DEGs.AT2.4dpi$Genes.SYMBOL %in% symbols.overlap)] <- T
write.csv(cbind(Overlap = flag, DEGs.AT2.4dpi), 
          quote = F, paste0(table.dir, "DEGs_AT2_4dpi.csv"))


# Save the overlapped DEGs in iAT2 1dpi
flag <- rep(F, nrow(DEGs.iAT2.1dpi))
flag[which(DEGs.iAT2.1dpi$hgnc_symbols %in% symbols.overlap)] <- T
write.csv(cbind(Overlap = flag, DEGs.iAT2.1dpi), 
          quote = F, paste0(table.dir, "DEGs_iAT2_1dpi.csv"))


# Save the overlapped DEGs in iAT2 3dpi
flag <- rep(F, nrow(DEGs.iAT2.3dpi))
flag[which(DEGs.iAT2.3dpi$hgnc_symbols %in% symbols.overlap)] <- T
write.csv(cbind(Overlap = flag, DEGs.iAT2.3dpi), 
          quote = F, paste0(table.dir, "DEGs_iAT2_3dpi.csv"))


# Load the GSEA result of iAT2
gsea.iAT2 <- read_excel(paste0(table.dir, "iAT2_DEGs/Supplemental_8._Figure_5F_GSEA_results.xlsx"),
                            sheet = 2)
dim(gsea.iAT2)
head(gsea.iAT2)
unique(gsea.iAT2$Database)
unique(gsea.iAT2$Cluster)
gsea.iAT2 <- gsea.iAT2[gsea.iAT2$Cluster != "2 dpi vs mock",]
gsea.iAT2 <- gsea.iAT2[order(gsea.iAT2$pvalue),]
dim(gsea.iAT2)


# GSEA analysis for DEG lists of AT2 1dpi
sorted.AT2.1dpi <- DEGs.AT2.1dpi[order(DEGs.AT2.1dpi$LogFC.dpi1_vs_ctr, decreasing = T), 
                                        "LogFC.dpi1_vs_ctr"]
class(sorted.AT2.1dpi)
dim(sorted.AT2.1dpi)
sorted.AT2.1dpi <- cbind(symbols = DEGs.AT2.1dpi[order(DEGs.AT2.1dpi$LogFC.dpi1_vs_ctr, decreasing = T), 
                                                 "Genes.SYMBOL"], 
                         sorted.AT2.1dpi)
sorted.AT2.1dpi <- cbind(Entrez = symbol_to_ENTREZ(gene.list = sorted.AT2.1dpi$Genes.SYMBOL, 
                                                   org.db = org.Hs.eg.db), 
                         sorted.AT2.1dpi)
dim(sorted.AT2.1dpi)
tail(sorted.AT2.1dpi)
sorted.AT2.1dpi <- sorted.AT2.1dpi[!is.na(sorted.AT2.1dpi$Entrez),]
dim(sorted.AT2.1dpi)
sorted.AT2.1dpi
Entrez.AT2.1dpi <- sorted.AT2.1dpi$LogFC.dpi1_vs_ctr
names(Entrez.AT2.1dpi) <- sorted.AT2.1dpi$Entrez
head(Entrez.AT2.1dpi)
all_gene_sets <- msigdbr(species = "Homo sapiens")
dim(all_gene_sets)
dim(all_gene_sets)
unique(all_gene_sets$gs_cat)
unique(all_gene_sets$gs_subcat)
msigdbr_list <- split(x = all_gene_sets$entrez_gene, f = all_gene_sets$gs_name)
gsea.AT2.1dpi <- fgsea(pathways = msigdbr_list, stats = Entrez.AT2.1dpi)
class(gsea.AT2.1dpi)
dim(gsea.AT2.1dpi)
gsea.AT2.1dpi <- gsea.AT2.1dpi[!is.na(gsea.AT2.1dpi$padj),]
dim(gsea.AT2.1dpi)
range(gsea.AT2.1dpi$padj)
gsea.AT2.1dpi <- gsea.AT2.1dpi[gsea.AT2.1dpi$padj < pval.cutoff,]
dim(gsea.AT2.1dpi)


# GSEA analysis for DEG lists of AT2 4dpi
sorted.AT2.4dpi <- DEGs.AT2.4dpi[order(DEGs.AT2.4dpi$logFC, decreasing = T), 
                                 "logFC"]
class(sorted.AT2.4dpi)
length(sorted.AT2.4dpi)
sorted.AT2.4dpi <- data.frame(symbols = DEGs.AT2.4dpi[order(DEGs.AT2.4dpi$logFC, decreasing = T), 
                                                 "Symbol"], 
                         logFC = sorted.AT2.4dpi)
sorted.AT2.4dpi <- sorted.AT2.4dpi[sorted.AT2.4dpi$symbols != "",]
Entrez.AT2.4dpi <- symbol_to_ENTREZ(gene.list = sorted.AT2.4dpi$symbols, 
                          org.db = org.Hs.eg.db)
sorted.AT2.4dpi <- sorted.AT2.4dpi[sorted.AT2.4dpi$symbols %in% names(Entrez.AT2.4dpi),]
dim(sorted.AT2.4dpi)
sorted.AT2.4dpi <- cbind(Entrez = symbol_to_ENTREZ(gene.list = sorted.AT2.4dpi$symbols, 
                                org.db = org.Hs.eg.db), 
      sorted.AT2.4dpi)
tail(sorted.AT2.4dpi)
sorted.AT2.4dpi <- sorted.AT2.4dpi[!is.na(sorted.AT2.4dpi$Entrez),]
dim(sorted.AT2.4dpi)
sorted.AT2.4dpi
Entrez.AT2.4dpi <- sorted.AT2.4dpi$logFC
names(Entrez.AT2.4dpi) <- sorted.AT2.4dpi$Entrez
head(Entrez.AT2.4dpi)
gsea.AT2.4dpi <- fgsea(pathways = msigdbr_list, stats = Entrez.AT2.4dpi)
class(gsea.AT2.4dpi)
dim(gsea.AT2.4dpi)
gsea.AT2.4dpi <- gsea.AT2.4dpi[!is.na(gsea.AT2.4dpi$padj),]
dim(gsea.AT2.4dpi)
range(gsea.AT2.4dpi$padj)
gsea.AT2.4dpi <- gsea.AT2.4dpi[gsea.AT2.4dpi$padj < pval.cutoff,]
dim(gsea.AT2.4dpi)


# Merge the 1dpi and 4dpi pathways
gsea.AT2.merged <- rbind(gsea.AT2.1dpi, gsea.AT2.4dpi)
dim(gsea.AT2.merged)
colnames(gsea.AT2.merged)
gsea.AT2.merged <- gsea.AT2.merged[order(gsea.AT2.merged$pval),]
head(gsea.AT2.merged$padj)
pathways.AT2 <- unique(gsea.AT2.merged$pathway)
length(pathways.AT2) # 1837
pathways.iAT2 <- unique(gsea.iAT2$ID)
length(pathways.iAT2) # 450
pathways.AT2.overlapped <- intersect(pathways.AT2, pathways.iAT2)
length(pathways.AT2.overlapped) # 115
pathway.overlap.ratios <- get_overlap_ratio(x = pathways.iAT2, pathways.AT2)
# The maximum overlap ratio, 0.416666666666667, is achieved at 12.
