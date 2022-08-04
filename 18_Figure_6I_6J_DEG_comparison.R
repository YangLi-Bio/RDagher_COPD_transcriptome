####################################################
#                                                  #
# Comparison analyses between co-cultured iAT2 and #
# Fibroblast with public data                      #
#                                                  #
####################################################


# 1.	Please compare Top 200 DEGs between iAT2-healthy fibroblasts (mock) with public data on healthy lung
# 2.	Please compare Top 200 DEGs between iAT2-COPD fibroblasts (mock) with public data on healthy lung
# 3.	Please compare Top 200 DEGs between mock v.s 1dpi of iAT2-healthy fibroblasts with public data on infected healthy lung 


# Comparison	Condition 1	Condition 2
# COPD fibro 1dpi vs COPD Fibro mock	Co-culture of iAT2 and fibroblasts from healthy lung, one day post infection	Co-culture of iAT2 and fibroblasts from healthy lung, without infection
# Healthy Fibro 1dpi vs Healthy Fibro mock	Co-culture of iAT2 and fibroblasts from COPD lung, one day post infection	Co-culture of iAT2 and fibroblasts from COPD lung, without infection
# COPD Fibro 1dpi vs Healthy Fibro 1dpi	Co-culture of iAT2 and fibroblasts from COPD lung, one day post infection	Co-culture of iAT2 and fibroblasts from healthy lung, one day post infection


# Libraries
library(Seurat)
library(dplyr)
library(parallel)
library(pbmcapply)
library(future)
library(readxl)


# Global variables
R.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Rfiles/"
table.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Tables/Figures_6I-6J_co-culture_0802200022/"
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
coculture.file <- "Supplemental_11_Figure_6J_DEG_results.xlsx"
padj.cutoff <- 0.05
FC.cutoff <- 0.25


# Source codes
source(paste0(tool.dir, "transcriptome_tools.R"))


# Load Excel files
table1 <- read_excel(paste0(table.dir, coculture.file), 
                     sheet = 2)
table2 <- read_excel(paste0(table.dir, coculture.file), 
                     sheet = 3)
table3 <- read_excel(paste0(table.dir, coculture.file), 
                     sheet = 4)
dim(table1) # 19148    
head(table1)
range(table1$PValue)
range(table1$FC)
dim(table1[table1$PValue < padj.cutoff,]) # 661
dim(table1[table1$FC > FC.cutoff,]) # 17787
dim(table1[!is.na(table1$hgnc_symbols),]) # 16450
# dim(table1[table1$PValue < padj.cutoff & table1$FC > FC.cutoff & 
#              !is.na(table1$hgnc_symbols),]) # 49
dim(table1[table1$PValue < padj.cutoff & 
             !is.na(table1$hgnc_symbols),]) # 593


dim(table2) # 19148
head(table2)
range(table2$PValue)
range(table2$FC)
dim(table2[table2$PValue < padj.cutoff,]) # 246
dim(table2[table2$FC > FC.cutoff,]) # 18197    
dim(table2[!is.na(table2$hgnc_symbols),]) # 16450
# dim(table2[table2$PValue < padj.cutoff & table2$FC > FC.cutoff & 
#              !is.na(table2$hgnc_symbols),]) # 48
dim(table2[table2$PValue < padj.cutoff & 
             !is.na(table2$hgnc_symbols),]) # 189


dim(table3) # 19148
head(table3)
range(table3$PValue)
range(table3$FC)
dim(table3[table3$PValue < padj.cutoff,]) # 68
dim(table3[table3$FC > FC.cutoff,]) # 18957        
dim(table3[!is.na(table3$hgnc_symbols),]) # 16450
# dim(table3[table3$PValue < padj.cutoff & table3$FC > FC.cutoff & 
#              !is.na(table3$hgnc_symbols),]) # 26
dim(table3[table3$PValue < padj.cutoff & 
             !is.na(table3$hgnc_symbols),]) # 49


# Load the public data
degs.public <- read.csv("/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/covid-19_vs_healthy/DEGs_Fibroblasts_COVID-19_vs_health.csv", row.names = 1)
dim(degs.public) # 1517
head(degs.public)
range(degs.public$avg_log2FC)
range(degs.public$p_val_adj)
dim(degs.public[degs.public$p_val_adj < padj.cutoff,]) # 1326
dim(degs.public[degs.public$avg_log2FC > FC.cutoff,]) # 853
dim(degs.public[degs.public$p_val_adj < padj.cutoff & 
                  degs.public$avg_log2FC > FC.cutoff,]) # 853
dim(degs.public[degs.public$p_val_adj < padj.cutoff & 
                  degs.public$avg_log2FC > 1,]) # 87
dim(degs.public[degs.public$p_val_adj < padj.cutoff & 
                  degs.public$avg_log2FC < -FC.cutoff,]) # 473


dim(degs.public[degs.public$p_val_adj < padj.cutoff & 
                  degs.public$avg_log2FC > 0,]) # 853
dim(degs.public[degs.public$p_val_adj < padj.cutoff & 
                  degs.public$avg_log2FC < 0,]) # 473


# Get common genes
# common.genes <- intersect(unique(unlist(table2[!is.na(table2$hgnc_symbols), "hgnc_symbols"])), 
#                           unique(rownames(degs.public)))
# length(common.genes) # 1428


# Classify DEGs from public 
# healthy.degs.public <- degs.public[degs.public$p_val_adj < padj.cutoff & 
#                                      degs.public$avg_log2FC < -FC.cutoff & 
#                                      rownames(degs.public) %in% common.genes,]
# copd.degs.public <- degs.public[degs.public$p_val_adj < padj.cutoff & 
#                                   degs.public$avg_log2FC > FC.cutoff & 
#                                   rownames(degs.public) %in% common.genes,]
# dim(healthy.degs.public) # 417
# dim(copd.degs.public) # 838


# Classify DEGs for our data
deg1 <- table1[table1$PValue < padj.cutoff & 
                 !is.na(table1$hgnc_symbols),]
deg1 <- deg1[order(deg1$PValue),]
deg1 <- deg1[!duplicated(deg1$hgnc_symbols),]
head(deg1)
dim(deg1) # 593
range(deg1$FC)


deg2 <- table2[table2$PValue < padj.cutoff  & 
                 !is.na(table2$hgnc_symbols),]
deg2 <- deg2[order(deg2$PValue),]
head(deg2)
dim(deg2) # 189


deg3 <- table3[table3$PValue < padj.cutoff & 
                 !is.na(table3$hgnc_symbols),]
deg3 <- deg3[order(deg3$PValue),]
head(deg3)
dim(deg3) # 49


# Parameter tuning for DEG analyses between COVID-19 v.s. ctrl
obj <- qs::qread("/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Rfiles/covid_and_health.qsave")
dim(obj)


# DEG analysis for Fibroblasts
Idents(obj) <- obj$cell.type
levels(Idents(obj))
dim(obj) # 34546
colnames(table1)
gene.list <- list(unique(unlist(table1[!is.na(table1$hgnc_symbols), "hgnc_symbols"])), 
                  unique(unlist(table2[!is.na(table2$hgnc_symbols), "hgnc_symbols"])), 
                  unique(unlist(table3[!is.na(table3$hgnc_symbols), "hgnc_symbols"])))
sapply(gene.list, length)
DefaultAssay(obj)
gene.union <- intersect(Reduce("union", gene.list), rownames(obj))
length(gene.union) # 15375
# plan("multicore", workers = 2)
# options(future.globals.maxSize = 3 * 1000 * 1024^2)
degs.public <- FindMarkers(obj, ident.1 = "COVID-19", group.by = "groups", 
                         subset.ident = "Fibroblasts", 
                         features = gene.union, logfc.threshold = 0.0, min.pct = 0.0, 
                         test.use = "bimod")
dim(degs.public) # 15375
head(degs.public)
qs::qsave(degs.public, paste0(R.dir, "DEGs_COVID_ctr_bimod.qsave"))
stop("Here!\n")


# Check overlaps for DEGs1
# dim(degs.public[degs.public$p_val_adj < padj.cutoff,]) # 9945
# length(intersect(rownames(degs.public), deg1$hgnc_symbols)) # 557
# length(intersect(rownames(degs.public), deg1$hgnc_symbols)) / 
#   dim(deg1) # 93.9%
# length(intersect(rownames(degs.public[degs.public$p_val_adj < padj.cutoff & 
#                                         degs.public$avg_log2FC > 0,]), deg1$hgnc_symbols)) # 246
# length(intersect(rownames(degs.public[degs.public$p_val_adj < padj.cutoff & 
#                                         degs.public$avg_log2FC > 0,]), deg1$hgnc_symbols)) / 
#   dim(deg1) # 41.5%
# common1 <- intersect(unique(unlist(table1[!is.na(table1$hgnc_symbols), "hgnc_symbols"])), 
#                      rownames(degs.public))
# length(common1)
# length(intersect(rownames(degs.public[degs.public$p_val_adj < padj.cutoff & 
#                                         degs.public$avg_log2FC > 0 & 
#                                         rownames(degs.public) %in% common1,]), 
#                  deg1[deg1$hgnc_symbols %in% common1,]$hgnc_symbols)) / 
#   dim(deg1[deg1$hgnc_symbols %in% common1,]) # 44.2%
# overlap1 <- intersect(rownames(degs.public[degs.public$p_val_adj < padj.cutoff & 
#                                              degs.public$avg_log2FC > 0 & 
#                                              rownames(degs.public) %in% common1,]), 
#                       deg1[deg1$hgnc_symbols %in% common1,]$hgnc_symbols)
# length(overlap1) # 246


degs.public <- degs.public[order(degs.public$p_val),]
common1 <- intersect(unique(unlist(table1[!is.na(table1$hgnc_symbols), "hgnc_symbols"])),
                     rownames(degs.public))
deg1.common <- deg1[deg1$hgnc_symbols %in% common1,]
dim(deg1.common) # 557
deg.copd.common <- degs.public[degs.public$p_val_adj < padj.cutoff & 
                                 degs.public$avg_log2FC > 0 & 
                                 rownames(degs.public) %in% common1,]
dim(deg.copd.common) # 6233
overlap1 <- intersect(rownames(deg.copd.common),
                      deg1.common$hgnc_symbols)
length(overlap1) # 246
length(overlap1) / nrow(deg1.common) # 44.2%
length(rownames(deg.copd.common)) - length(overlap1) # 5987
length(deg1.common$hgnc_symbols) - length(overlap1) # 311
write.csv(cbind(Overlap = deg1.common$hgnc_symbols %in% overlap1, 
                      deg1.common), 
          quote = F, paste0(table.dir, "Our_DEGs_comp_1.csv"))
write.csv(cbind(Overlap = rownames(deg.copd.common) %in% overlap1, 
                deg.copd.common), 
          quote = F, paste0(table.dir, "Public_DEGs_comp_1.csv"))


top.pub.deg1 <- deg.copd.common[1:200,]
top.our.deg1 <- deg1.common[1:200,]
dim(top.pub.deg1)
dim(top.our.deg1)
top.overlap1 <- intersect(rownames(top.pub.deg1),
                          top.our.deg1$hgnc_symbols)
length(top.overlap1) # 21
length(top.overlap1) / nrow(top.our.deg1) # 10.5%
length(rownames(top.pub.deg1)) - length(top.overlap1) # 179
length(top.our.deg1$hgnc_symbols) - length(top.overlap1) # 179
write.csv(cbind(Overlap = top.our.deg1$hgnc_symbols %in% top.overlap1, 
                top.our.deg1), 
          quote = F, paste0(table.dir, "Top200_Our_DEGs_comp_1.csv"))
write.csv(cbind(Overlap = rownames(top.pub.deg1) %in% top.overlap1, 
                top.pub.deg1), 
          quote = F, paste0(table.dir, "Top200_Public_DEGs_comp_1.csv"))



# GO enrichment analyses for our data
our.enriched1 <- run_GO_and_KEGG(genes.ll = deg1.common$hgnc_symbols, org = "Human")
our.mf1 <- our.enriched1$GO_Molecular_Function_2018
our.cc1 <- our.enriched1$GO_Cellular_Component_2018
our.bp1 <- our.enriched1$GO_Biological_Process_2018
our.kegg1 <- our.enriched1$KEGG_2019_Human


our.mf1 <- our.mf1[our.mf1$Adjusted.P.value < 0.05,]
our.cc1 <- our.cc1[our.cc1$Adjusted.P.value < 0.05,]
our.bp1 <- our.bp1[our.bp1$Adjusted.P.value < 0.05,]
our.kegg1 <- our.kegg1[our.kegg1$Adjusted.P.value < 0.05,]




# Enrichment for public data
pub.copd.enriched <- run_GO_and_KEGG(genes.ll = rownames(deg.copd.common), org = "Human")
pub.copd.mf1 <- pub.copd.enriched$GO_Molecular_Function_2018
pub.copd.cc1 <- pub.copd.enriched$GO_Cellular_Component_2018
pub.copd.bp1 <- pub.copd.enriched$GO_Biological_Process_2018
pub.copd.kegg1 <- pub.copd.enriched$KEGG_2019_Human


pub.copd.mf1 <- pub.copd.mf1[pub.copd.mf1$Adjusted.P.value < 0.05,]
pub.copd.cc1 <- pub.copd.cc1[pub.copd.cc1$Adjusted.P.value < 0.05,]
pub.copd.bp1 <- pub.copd.bp1[pub.copd.bp1$Adjusted.P.value < 0.05,]
pub.copd.kegg1 <- pub.copd.kegg1[pub.copd.kegg1$Adjusted.P.value < 0.05,]



# MF
overlap.mf <- intersect(our.mf1$Term, pub.copd.mf1$Term)
length(overlap.mf) / nrow(our.mf1) # 7.4%
length(overlap.mf) # 2
nrow(our.mf1) - length(overlap.mf) # 25
nrow(pub.copd.mf1) - length(overlap.mf) # 126
write.csv(cbind(Overlap = our.mf1$Term %in% overlap.mf, 
                our.mf1), quote = F, 
          paste0(table.dir, "Our_Molecular_Function_comp_1.csv"))
write.csv(cbind(Overlap = pub.copd.mf1$Term %in% overlap.mf, 
                pub.copd.mf1), quote = F, 
          paste0(table.dir, "Public_Molecular_Function_comp_1.csv"))




# CC
overlap.cc <- intersect(our.cc1$Term, pub.copd.cc1$Term)
length(overlap.cc) / nrow(our.cc1) # 71.4%
length(overlap.cc) # 5
nrow(our.cc1) - length(overlap.cc) # 2
nrow(pub.copd.cc1) - length(overlap.cc) # 154
write.csv(cbind(Overlap = our.cc1$Term %in% overlap.cc, 
                our.cc1), quote = F, 
          paste0(table.dir, "Our_Cellular_Component_comp_1.csv"))
write.csv(cbind(Overlap = pub.copd.cc1$Term %in% overlap.cc, 
                pub.copd.cc1), quote = F, 
          paste0(table.dir, "Public_Cellular_Component_comp_1.csv"))


# BP
overlap.bp <- intersect(our.bp1$Term, pub.copd.bp1$Term)
length(overlap.bp) / nrow(our.bp1) # 16.5%
length(overlap.bp) # 40
nrow(our.bp1) - length(overlap.bp) # 202
nrow(pub.copd.bp1) - length(overlap.bp) # 632
write.csv(cbind(Overlap = our.bp1$Term %in% overlap.bp, 
                our.bp1), quote = F, 
          paste0(table.dir, "Our_Biological_Process_comp_1.csv"))
write.csv(cbind(Overlap = pub.copd.bp1$Term %in% overlap.bp, 
                pub.copd.bp1), quote = F, 
          paste0(table.dir, "Public_Biological_Process_comp_1.csv"))


# kegg
overlap.kegg <- intersect(our.kegg1$Term, pub.copd.kegg1$Term)
length(overlap.kegg) / nrow(our.kegg1) # 33.3%
length(overlap.kegg) # 8
nrow(our.kegg1) - length(overlap.kegg) # 16
nrow(pub.copd.kegg1) - length(overlap.kegg) # 91
write.csv(cbind(Overlap = our.kegg1$Term %in% overlap.kegg, 
                our.kegg1), quote = F, 
          paste0(table.dir, "Our_KEGG_comp_1.csv"))
write.csv(cbind(Overlap = pub.copd.kegg1$Term %in% overlap.kegg, 
                pub.copd.kegg1), quote = F, 
          paste0(table.dir, "Public_KEGG_comp_1.csv"))


# Top-200


# Our Top-200
top.our.enriched1 <- run_GO_and_KEGG(genes.ll = top.our.deg1$hgnc_symbols, org = "Human")
top.our.mf1 <- top.our.enriched1$GO_Molecular_Function_2018
top.our.cc1 <- top.our.enriched1$GO_Cellular_Component_2018
top.our.bp1 <- top.our.enriched1$GO_Biological_Process_2018
top.our.kegg1 <- top.our.enriched1$KEGG_2019_Human


top.our.mf1 <- top.our.mf1[top.our.mf1$Adjusted.P.value < 0.05,]
top.our.cc1 <- top.our.cc1[top.our.cc1$Adjusted.P.value < 0.05,]
top.our.bp1 <- top.our.bp1[top.our.bp1$Adjusted.P.value < 0.05,]
top.our.kegg1 <- top.our.kegg1[top.our.kegg1$Adjusted.P.value < 0.05,]


# Public Top-200
top.pub.copd.enriched <- run_GO_and_KEGG(genes.ll = rownames(top.pub.deg1), org = "Human")
top.pub.copd.mf1 <- top.pub.copd.enriched$GO_Molecular_Function_2018
top.pub.copd.cc1 <- top.pub.copd.enriched$GO_Cellular_Component_2018
top.pub.copd.bp1 <- top.pub.copd.enriched$GO_Biological_Process_2018
top.pub.copd.kegg1 <- top.pub.copd.enriched$KEGG_2019_Human


top.pub.copd.mf1 <- top.pub.copd.mf1[top.pub.copd.mf1$Adjusted.P.value < 0.05,]
top.pub.copd.cc1 <- top.pub.copd.cc1[top.pub.copd.cc1$Adjusted.P.value < 0.05,]
top.pub.copd.bp1 <- top.pub.copd.bp1[top.pub.copd.bp1$Adjusted.P.value < 0.05,]
top.pub.copd.kegg1 <- top.pub.copd.kegg1[top.pub.copd.kegg1$Adjusted.P.value < 0.05,]



# MF top-200
top.overlap.mf <- intersect(top.our.mf1$Term, top.pub.copd.mf1$Term)
length(top.overlap.mf) / nrow(top.our.mf1) # 26.1%
length(top.overlap.mf) # 6
nrow(top.our.mf1) - length(top.overlap.mf) # 17
nrow(top.pub.copd.mf1) - length(top.overlap.mf) # 31
write.csv(cbind(top.overlap = top.our.mf1$Term %in% top.overlap.mf, 
                top.our.mf1), quote = F, 
          paste0(table.dir, "Top200_Our_Molecular_Function_comp_1.csv"))
write.csv(cbind(top.overlap = top.pub.copd.mf1$Term %in% top.overlap.mf, 
                top.pub.copd.mf1), quote = F, 
          paste0(table.dir, "Top200_Public_Molecular_Function_comp_1.csv"))



# cc top-200
top.overlap.cc <- intersect(top.our.cc1$Term, top.pub.copd.cc1$Term)
length(top.overlap.cc) / nrow(top.our.cc1) # 42.9%
length(top.overlap.cc) # 3
nrow(top.our.cc1) - length(top.overlap.cc) # 4
nrow(top.pub.copd.cc1) - length(top.overlap.cc) # 37
write.csv(cbind(top.overlap = top.our.cc1$Term %in% top.overlap.cc, 
                top.our.cc1), quote = F, 
          paste0(table.dir, "Top200_Our_Cellular_Component_comp_1.csv"))
write.csv(cbind(top.overlap = top.pub.copd.cc1$Term %in% top.overlap.cc, 
                top.pub.copd.cc1), quote = F, 
          paste0(table.dir, "Top200_Public_Cellular_Component_comp_1.csv"))




# bp top-200
top.overlap.bp <- intersect(top.our.bp1$Term, top.pub.copd.bp1$Term)
length(top.overlap.bp) / nrow(top.our.bp1) # 28.2%
length(top.overlap.bp) # 20
nrow(top.our.bp1) - length(top.overlap.bp) # 51
nrow(top.pub.copd.bp1) - length(top.overlap.bp) # 257
write.csv(cbind(top.overlap = top.our.bp1$Term %in% top.overlap.bp, 
                top.our.bp1), quote = F, 
          paste0(table.dir, "Top200_Our_Biological_Process_comp_1.csv"))
write.csv(cbind(top.overlap = top.pub.copd.bp1$Term %in% top.overlap.bp, 
                top.pub.copd.bp1), quote = F, 
          paste0(table.dir, "Top200_Public_Biological_Process_comp_1.csv"))





# kegg top-200
top.overlap.kegg <- intersect(top.our.kegg1$Term, top.pub.copd.kegg1$Term)
length(top.overlap.kegg) / nrow(top.our.kegg1) # 66.7%
length(top.overlap.kegg) # 8
nrow(top.our.kegg1) - length(top.overlap.kegg) # 4
nrow(top.pub.copd.kegg1) - length(top.overlap.kegg) # 27
write.csv(cbind(top.overlap = top.our.kegg1$Term %in% top.overlap.kegg, 
                top.our.kegg1), quote = F, 
          paste0(table.dir, "Top200_Our_KEGG_comp_1.csv"))
write.csv(cbind(top.overlap = top.pub.copd.kegg1$Term %in% top.overlap.kegg, 
                top.pub.copd.kegg1), quote = F, 
          paste0(table.dir, "Top200_Public_KEGG_comp_1.csv"))
