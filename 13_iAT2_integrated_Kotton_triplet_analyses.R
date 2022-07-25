###########################################################################
#                                                                         #
# This script includes the following analyses:                            #
# 1. Compare the top-200 DEGs of iAT2 cluster and iAT2 proliferative      #
#    cluster, respectively, with DEGs between healthy v.s.                #
#    infected from public data                                            #
#                                                                         #
###########################################################################


# Libraries
library(Seurat)
library(future)
library(dplyr)


# Global parameters
R.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Rfiles/"
table.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Tables/"
copd.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/2_copd_Kaminski_Sci_Adv_2020/10x/"
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"


# Source scripts
source(paste0(tool.dir, "combinatorial_tools.R"))
source(paste0(tool.dir, "transcriptome_tools.R"))


# Local parameters
iAT2.DEG.file <- "COPD_healthy_DEGs_iAT2.csv"
copd.meta.file <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/2_copd_Kaminski_Sci_Adv_2020/GSE136831_AllCells.Samples.CellType.MetadataTable.txt"
AT2.dir <- "AT2_DEG_overlap/"
padj.cutoff <- 0.05
logFC.cutoff <- 0
top.cutoff <- 200


# Identify the top-200 DEGs between healthy v.s. COPD in iAT2 & iAT2 proliferative, i.e., 
# clusters 1 & 2
DEGs.iAT2.healthy <- read.csv(paste0(table.dir, iAT2.DEG.file))
dim(DEGs.iAT2.healthy)
head(DEGs.iAT2.healthy[, 1:5])
range(DEGs.iAT2.healthy$adj_pval)
unique(DEGs.iAT2.healthy$cluster_name)
range(DEGs.iAT2.healthy$logfc)
DEGs.iAT2.healthy <- DEGs.iAT2.healthy[DEGs.iAT2.healthy$cluster_name == "healthy",]
dim(DEGs.iAT2.healthy)
range(DEGs.iAT2.healthy$rank_in_cluster)
head(DEGs.iAT2.healthy$rank_in_cluster)
range(DEGs.iAT2.healthy$adj_pval)
range(DEGs.iAT2.healthy$logfc)
DEGs.iAT2.healthy.top200 <- head(DEGs.iAT2.healthy, n = top.cutoff)
dim(DEGs.iAT2.healthy.top200)
range(DEGs.iAT2.healthy.top200$adj_pval)
range(DEGs.iAT2.healthy.top200$logfc)
write.csv(DEGs.iAT2.healthy.top200, paste0(table.dir, 
                                           "DEGs_healthy_iAT2_clusters_1_2_top200.csv"), 
          row.names = F, quote = F)


# Discover the the top-200 DEGs between healthy v.s. COPD in iAT2 from public data
copd.data <- Read10X(data.dir = copd.dir, gene.column = 2)
obj.copd <- CreateSeuratObject(counts = copd.data, project = "copd", 
                           min.cells = 0, min.features = 0)
ncol(obj.copd)
colnames(obj.copd)
obj.copd@meta.data
colnames(obj.copd@meta.data)
unique(obj.copd$cell.type)
copd.meta <- read.csv(copd.meta.file, header = T, sep = "\t")
dim(copd.meta)
head(copd.meta)
unique(copd.meta$Disease_Identity)
Idents(obj.copd) <- obj.copd$cell.type

# Add disease status
copd.disease <- copd.meta[, c(1, 7)]
copd.disease
copd.add.meta <- copd.disease[, 2]
names(copd.add.meta) <- copd.disease[, 1]
copd.add.meta <- copd.add.meta[colnames(obj.copd)]
length(copd.add.meta)
obj.copd <- AddMetaData(obj.copd, metadata = copd.add.meta, col.name = "Disease")
unique(obj.copd$Disease)


# Add cell cluster information
copd.ct <- copd.meta[, c(1, 5)]
copd.ct
copd.add.ct <- copd.ct[, 2]
names(copd.add.ct) <- copd.ct[, 1]
copd.add.ct <- copd.add.ct[colnames(obj.copd)]
length(copd.add.ct)
obj.copd <- AddMetaData(obj.copd, metadata = copd.add.ct, col.name = "cell.type")
levels(obj.copd$cell.type)
qs::qsave(obj.copd, paste0(R.dir, "Obj_COPD_healthy.qsave"))


# Identify DEGs in ATII (AT2) in control (healthy)
Idents(obj.copd) <- obj.copd$Disease
unique(Idents(obj.copd))
obj.ctrl <- subset(obj.copd, idents = "Control")
ncol(obj.ctrl)
colnames(obj.ctrl)
qs::qsave(obj.ctrl, paste0(R.dir, "Obj_AT2_ctrl.qsave"))

obj.ctrl.copd <- subset(obj.copd, idents = c("Control", "COPD"))
ncol(obj.ctrl.copd)
colnames(obj.ctrl.copd)
qs::qsave(obj.ctrl.copd, paste0(R.dir, "Obj_AT2_ctrl_COPD.qsave"))


# Identify DEGs in ATII (AT2) of healthy (control)
Idents(obj.ctrl) <- obj.ctrl$cell.type
unique(obj.ctrl$Disease)
unique(Idents(obj.ctrl))
plan("multicore", workers = 2)
options(future.globals.maxSize = 3 * 1000 * 1024^2)
DEGs.AT2.ctrl <- FindMarkers(obj.ctrl, ident.1 = "ATII")
qs::qsave(DEGs.AT2.ctrl, paste0(R.dir, "DEGs_AT2_ctrl.qsave"))


# # A small test for DEG discovery between ATII v.s. other cell types
# obj.sampled <- subset(obj.ctrl, cells = colnames(obj.ctrl)[sample(1:ncol(obj.ctrl), 1000)])
# Idents(obj.sampled) <- obj.sampled$cell.type
# DEGs.AT2.sampled <- FindMarkers(obj.sampled, ident.1 = "ATII")
# dim(DEGs.AT2.sampled)
# qs::qsave(DEGs.AT2.sampled, paste0(R.dir, "DEGs_AT2_ctrl_sampled.qsave"))
# DEGs.AT2.sampled <- DEGs.AT2.sampled[DEGs.AT2.sampled$p_val_adj < padj.cutoff & 
#                                        DEGs.AT2.sampled$avg_log2FC > logFC.cutoff,]
# dim(DEGs.AT2.sampled)


# Find DEGs between control v.s. COPD in ATII (AT2)
# we use this plan
Idents(obj.ctrl.copd) <- obj.ctrl.copd$cell.type
DEGs.AT2.ctrl.copd <- FindMarkers(obj.ctrl.copd, ident.1 = "Control", group.by = "Disease", 
                                  subset.ident = "ATII")
dim(DEGs.AT2.ctrl.copd)
head(DEGs.AT2.ctrl.copd)
colnames(DEGs.AT2.ctrl.copd)
range(DEGs.AT2.ctrl.copd$avg_log2FC)
range(DEGs.AT2.ctrl.copd$p_val_adj)
qs::qsave(DEGs.AT2.ctrl.copd, paste0(R.dir, "DEGs_ctrl_COPD_in_AT2.qsave"))

DEGs.AT2.ctrl.copd.top200 <- DEGs.AT2.ctrl.copd[DEGs.AT2.ctrl.copd$p_val_adj < padj.cutoff & 
                                           DEGs.AT2.ctrl.copd$avg_log2FC > logFC.cutoff,]
dim(DEGs.AT2.ctrl.copd.top200)
DEGs.AT2.ctrl.copd.top200 <- DEGs.AT2.ctrl.copd.top200[order(DEGs.AT2.ctrl.copd.top200$p_val),] %>% 
  head(n = top.cutoff)
dim(DEGs.AT2.ctrl.copd.top200)
DEGs.AT2.ctrl.copd.top200
write.csv(DEGs.AT2.ctrl.copd.top200, paste0(table.dir, "DEGs_ctrl_COPD_in_AT2_top200.csv"))


# Filter out the genes not included in iAT2 DEGs
DEGs.AT2.ctrl.copd.filtered <- DEGs.AT2.ctrl.copd[intersect(rownames(DEGs.AT2.ctrl.copd), 
                                                            DEGs.iAT2.healthy$gene),]
dim(DEGs.AT2.ctrl.copd.filtered)
DEGs.AT2.ctrl.copd.filtered <- DEGs.AT2.ctrl.copd.filtered[DEGs.AT2.ctrl.copd.filtered$p_val_adj < padj.cutoff & 
                                                             DEGs.AT2.ctrl.copd.filtered$avg_log2FC > logFC.cutoff,]
dim(DEGs.AT2.ctrl.copd.filtered)
# DEGs.AT2.ctrl.copd.filtered <- DEGs.AT2.ctrl.copd.filtered[rownames(DEGs.AT2.ctrl.copd.filtered) %in% 
#                                                              DEGs.iAT2.healthy$gene,]
# dim(DEGs.AT2.ctrl.copd.filtered)

DEGs.iAT2.healthy.filtered <- DEGs.iAT2.healthy[DEGs.iAT2.healthy$gene %in% rownames(DEGs.AT2.ctrl.copd),]
DEGs.AT2.ctrl.copd.filtered <- DEGs.AT2.ctrl.copd.filtered[order(DEGs.AT2.ctrl.copd.filtered$p_val),]
symbols.AT2.ctrl.copd.filtered <- rownames(DEGs.AT2.ctrl.copd.filtered) %>% head(n = top.cutoff)


# # Get top-200 DEGs in iAT2 ranked by logFC
# DEGs.iAT2.healthy
# range(DEGs.iAT2.healthy$adj_pval)
# range(DEGs.iAT2.healthy$logfc)
# DEGs.iAT2.healthy.logFC <- DEGs.iAT2.healthy[DEGs.iAT2.healthy$adj_pval < padj.cutoff,]
# DEGs.iAT2.healthy.logFC <- DEGs.iAT2.healthy.logFC[order(DEGs.iAT2.healthy.logFC$logfc,
#                                                          decreasing = T),]
# DEGs.iAT2.healthy.logFC <- head(DEGs.iAT2.healthy.logFC, n = top.cutoff)
# symbols.iAT2.healthy.logFC <- DEGs.iAT2.healthy.logFC$gene
# length(symbols.iAT2.healthy.logFC)
# intersect(symbols.iAT2.healthy.logFC, symbols.AT2.ctrl)
# # Only 29 overlapped genes


# Compare top-200 DEGs in iAT2 & iAT2 proliferative (Clusters 1 & 2) and AT2 in ctrl
symbols.iAT2.healthy <- DEGs.iAT2.healthy.top200$gene
symbols.AT2.ctrl <- rownames(DEGs.AT2.ctrl.copd.top200)
overlaps.healthy.iAT2.AT2 <- intersect(symbols.iAT2.healthy, symbols.AT2.ctrl)
qs::qsave(overlaps.healthy.iAT2.AT2, paste0(R.dir, "Overlaped_DEGs_AT2_top200.qsave"))

length(overlaps.healthy.iAT2.AT2)
length(overlaps.healthy.iAT2.AT2) / length(symbols.AT2.ctrl)

intersect(symbols.AT2.ctrl.copd.filtered, DEGs.iAT2.healthy$gene)


# Conclusion: if we selected the top-200 DEGs in iAT2 (Clusters 1 & 2) and AT2 (ctrl v.s. COPD)
# we only have 37 (18.5%) overlapped DEGs


# Generate plot of overlap ratios alongside the DEG ranks
DEGs.iAT2.healthy.filtered <- DEGs.iAT2.healthy[DEGs.iAT2.healthy$adj_pval < padj.cutoff & 
                                                  DEGs.iAT2.healthy$logfc > logFC.cutoff,]
DEGs.iAT2.healthy.filtered <- DEGs.iAT2.healthy.filtered[DEGs.iAT2.healthy.filtered$gene %in% 
                                                           rownames(DEGs.AT2.ctrl.copd),]
dim(DEGs.iAT2.healthy.filtered)
range(DEGs.iAT2.healthy.filtered$adj_pval)
range(DEGs.iAT2.healthy.filtered$logfc)
DEGs.iAT2.healthy.filtered <- DEGs.iAT2.healthy.filtered[order(DEGs.iAT2.healthy.filtered$pval),]
DEGs.iAT2.healthy.top <- head(DEGs.iAT2.healthy.filtered, n = top.cutoff)[, -1]
dim(DEGs.iAT2.healthy.top)
colnames(DEGs.iAT2.healthy.top)

DEGs.AT2.ctrl.copd <- qs::qread(paste0(R.dir, "DEGs_ctrl_COPD_in_AT2.qsave"))
DEGs.AT2.ctrl.copd <- DEGs.AT2.ctrl.copd[rownames(DEGs.AT2.ctrl.copd) %in% 
                                           DEGs.iAT2.healthy$gene,]
dim(DEGs.AT2.ctrl.copd)
DEGs.AT2.ctrl.copd.filtered <- DEGs.AT2.ctrl.copd[DEGs.AT2.ctrl.copd$p_val_adj < padj.cutoff &
                                                    DEGs.AT2.ctrl.copd$avg_log2FC > logFC.cutoff,]
dim(DEGs.AT2.ctrl.copd.filtered)
range(DEGs.AT2.ctrl.copd.filtered$p_val_adj)
range(DEGs.AT2.ctrl.copd.filtered$avg_log2FC)
DEGs.AT2.ctrl.copd.filtered <- DEGs.AT2.ctrl.copd.filtered[order(DEGs.AT2.ctrl.copd.filtered$p_val),]
DEGs.AT2.ctrl.copd.top <- head(DEGs.AT2.ctrl.copd.filtered, n = top.cutoff)[, -1]
dim(DEGs.AT2.ctrl.copd.top)
colnames(DEGs.AT2.ctrl.copd.top)
rownames(DEGs.AT2.ctrl.copd.top)

overlap.ratio <- get_overlap_ratio(x = rownames(DEGs.AT2.ctrl.copd.filtered), 
                                   DEGs.iAT2.healthy.filtered$gene)
plot(overlap.ratio, xlab = "Top-ranked DEGs", ylab = "Overlap ratio")

top.overlap <- intersect(rownames(DEGs.AT2.ctrl.copd.top), 
                                       DEGs.iAT2.healthy.top$gene)
length(top.overlap)
nrow(DEGs.AT2.ctrl.copd.top)
length(DEGs.iAT2.healthy.top$gene)
length(top.overlap) / 200


# Inspect overlap ratios
overlap.ratio[1000]
overlap.ratio[1500]


# Enrichment analyses (n= ALL)
enriched.iAT2 <- run_GO_and_KEGG(genes.ll = DEGs.iAT2.healthy.filtered$gene, org = "human")
enriched.AT2 <- run_GO_and_KEGG(genes.ll = rownames(DEGs.AT2.ctrl.copd.filtered), org = "human")
qs::qsave(enriched.iAT2, paste0(R.dir, "iAT2_enriched_all.qsave"))
qs::qsave(enriched.AT2, paste0(R.dir, "AT2_enriched_all.qsave"))


# Enrichment analyses (n = 200)
enriched.iAT2.top <- run_GO_and_KEGG(genes.ll = DEGs.iAT2.healthy.top$gene, org = "human")
enriched.AT2.top <- run_GO_and_KEGG(genes.ll = rownames(DEGs.AT2.ctrl.copd.top), org = "human")
qs::qsave(enriched.iAT2.top, paste0(R.dir, "iAT2_enriched_top200.qsave"))
qs::qsave(enriched.AT2.top, paste0(R.dir, "AT2_enriched_top200.qsave"))


# Overlapped terms in biological process
iAT2.BP <- enriched.iAT2$GO_Biological_Process_2018
iAT2.BP <- iAT2.BP[iAT2.BP$Adjusted.P.value < padj.cutoff,]
dim(iAT2.BP)
AT2.BP <- enriched.AT2$GO_Biological_Process_2018
AT2.BP <- AT2.BP[AT2.BP$Adjusted.P.value < padj.cutoff,]
dim(AT2.BP)
BP.enriched <- intersect(iAT2.BP$Term, AT2.BP$Term)
length(BP.enriched)
length(BP.enriched) / nrow(iAT2.BP) # 0.5301837

iAT2.BP.top <- enriched.iAT2.top$GO_Biological_Process_2018
iAT2.BP.top <- iAT2.BP.top[iAT2.BP$Adjusted.P.value < padj.cutoff,]
dim(iAT2.BP.top)
AT2.BP.top <- enriched.AT2.top$GO_Biological_Process_2018
AT2.BP.top <- AT2.BP.top[AT2.BP.top$Adjusted.P.value < padj.cutoff,]
dim(AT2.BP.top)
BP.enriched.top <- intersect(iAT2.BP.top$Term, AT2.BP.top$Term)
length(BP.enriched.top)
length(BP.enriched.top) / nrow(iAT2.BP.top) # 0.04104478


iAT2.BP <- cbind(Overlap = iAT2.BP$Term %in% BP.enriched.top, 
                 iAT2.BP)
dim(iAT2.BP)
colnames(iAT2.BP)
head(iAT2.BP$Overlap)
write.csv(iAT2.BP, paste0(table.dir, AT2.dir, "Biological_process_iAT2_clusters_1_2_iPSC.csv"), 
          quote = F)

iAT2.BP.top <- cbind(Overlap = iAT2.BP.top$Term %in% BP.enriched, 
                 iAT2.BP.top)
dim(iAT2.BP.top)
colnames(iAT2.BP.top)
dim(iAT2.BP.top)
head(iAT2.BP.top$Overlap)
write.csv(iAT2.BP.top, paste0(table.dir, AT2.dir, "Top200_Biological_process_iAT2_clusters_1_2_iPSC.csv"), 
          quote = F)

AT2.BP <- cbind(Overlap = AT2.BP.top$Term %in% BP.enriched.top, 
                AT2.BP)
dim(AT2.BP)
colnames(AT2.BP)
head(AT2.BP$Overlap)
write.csv(AT2.BP, paste0(table.dir, AT2.dir, "Biological_process_AT2_ctrl_COPD_public.csv"), 
          quote = F)

AT2.BP.top <- cbind(Overlap = AT2.BP.top$Term %in% BP.enriched.top, 
                    AT2.BP.top)
dim(AT2.BP.top)
colnames(AT2.BP.top)
head(AT2.BP.top$Overlap)
write.csv(AT2.BP.top, paste0(table.dir, AT2.dir, "Top200_Biological_process_AT2_ctrl_COPD_public.csv"), 
          quote = F)


# Overlapped terms in molecular function
iAT2.MF <- enriched.iAT2$GO_Molecular_Function_2018
iAT2.MF <- iAT2.MF[iAT2.MF$Adjusted.P.value < padj.cutoff,]
dim(iAT2.MF)
AT2.MF <- enriched.AT2$GO_Molecular_Function_2018
AT2.MF <- AT2.MF[AT2.MF$Adjusted.P.value < padj.cutoff,]
dim(AT2.MF)
MF.enriched <- intersect(iAT2.MF$Term, AT2.MF$Term)
length(MF.enriched)
length(MF.enriched) / nrow(iAT2.MF) # 0.4444444


iAT2.MF.top <- enriched.iAT2.top$GO_Molecular_Function_2018
iAT2.MF.top <- iAT2.MF.top[iAT2.MF.top$Adjusted.P.value < padj.cutoff,]
dim(iAT2.MF.top)
AT2.MF.top <- enriched.AT2.top$GO_Molecular_Function_2018
AT2.MF.top <- AT2.MF.top[AT2.MF.top$Adjusted.P.value < padj.cutoff,]
dim(AT2.MF.top)
MF.enriched.top <- intersect(iAT2.MF.top$Term, AT2.MF.top$Term)
length(MF.enriched.top)
length(MF.enriched.top) / nrow(iAT2.MF.top) # 0.4444444


iAT2.MF <- cbind(Overlap = iAT2.MF$Term %in% MF.enriched, 
                 iAT2.MF)
dim(iAT2.MF)
colnames(iAT2.MF)
head(iAT2.MF$Overlap)
write.csv(iAT2.MF, paste0(table.dir, AT2.dir, "Molecular_function_iAT2_clusters_1_2_iPSC.csv"), 
          quote = F)


iAT2.MF.top <- cbind(Overlap = iAT2.MF.top$Term %in% MF.enriched.top, 
                 iAT2.MF.top)
dim(iAT2.MF.top)
colnames(iAT2.MF.top)
head(iAT2.MF.top$Overlap)
write.csv(iAT2.MF.top, paste0(table.dir, AT2.dir, "Top200_Molecular_function_iAT2_clusters_1_2_iPSC.csv"), 
          quote = F)


AT2.MF <- cbind(Overlap = AT2.MF$Term %in% MF.enriched, 
                 AT2.MF)
dim(AT2.MF)
colnames(AT2.MF)
head(AT2.MF$Overlap)
write.csv(AT2.MF, paste0(table.dir, AT2.dir, "Molecular_function_AT2_ctrl_COPD_public.csv"), 
          quote = F)


AT2.MF.top <- cbind(Overlap = AT2.MF.top$Term %in% MF.enriched.top, 
                AT2.MF.top)
dim(AT2.MF.top)
colnames(AT2.MF.top)
head(AT2.MF.top$Overlap)
write.csv(AT2.MF.top, paste0(table.dir, AT2.dir, "Top200_Molecular_function_AT2_ctrl_COPD_public.csv"), 
          quote = F)


# Overlapped terms in cellular component
iAT2.CC <- enriched.iAT2$GO_Cellular_Component_2018
iAT2.CC <- iAT2.CC[iAT2.CC$Adjusted.P.value < padj.cutoff,]
dim(iAT2.CC)
AT2.CC <- enriched.AT2$GO_Cellular_Component_2018
AT2.CC <- AT2.CC[AT2.CC$Adjusted.P.value < padj.cutoff,]
dim(AT2.CC)
CC.enriched <- intersect(iAT2.CC$Term, AT2.CC$Term)
length(CC.enriched)
length(CC.enriched) / nrow(iAT2.CC) # 0.6893204


iAT2.CC.top <- enriched.iAT2.top$GO_Cellular_Component_2018
iAT2.CC.top <- iAT2.CC.top[iAT2.CC.top$Adjusted.P.value < padj.cutoff,]
dim(iAT2.CC.top)
AT2.CC.top <- enriched.AT2.top$GO_Cellular_Component_2018
AT2.CC.top <- AT2.CC.top[AT2.CC.top$Adjusted.P.value < padj.cutoff,]
dim(AT2.CC.top)
CC.enriched.top <- intersect(iAT2.CC.top$Term, AT2.CC.top$Term)
length(CC.enriched.top)
length(CC.enriched.top) / nrow(iAT2.CC.top) # 0.5789474


iAT2.CC <- cbind(Overlap = iAT2.CC$Term %in% CC.enriched, 
                 iAT2.CC)
dim(iAT2.CC)
colnames(iAT2.CC)
head(iAT2.CC$Overlap)
write.csv(iAT2.CC, paste0(table.dir, AT2.dir, "Cellular_component_iAT2_clusters_1_2_iPSC.csv"), 
          quote = F)


iAT2.CC.top <- cbind(Overlap = iAT2.CC.top$Term %in% CC.enriched.top, 
                iAT2.CC.top)
dim(iAT2.CC.top)
colnames(iAT2.CC.top)
head(iAT2.CC.top$Overlap)
write.csv(iAT2.CC.top, paste0(table.dir, AT2.dir, "Top200_Cellular_component_iAT2_clusters_1_2_iPSC.csv"), 
          quote = F)


AT2.CC <- cbind(Overlap = AT2.CC$Term %in% CC.enriched, 
                AT2.CC)
dim(AT2.CC)
colnames(AT2.CC)
head(AT2.CC$Overlap)
write.csv(AT2.CC, paste0(table.dir, AT2.dir, "Cellular_component_AT2_ctrl_COPD_public.csv"), 
          quote = F)


AT2.CC.top <- cbind(Overlap = AT2.CC.top$Term %in% CC.enriched.top, 
                AT2.CC.top)
dim(AT2.CC.top)
colnames(AT2.CC.top)
head(AT2.CC.top$Overlap)
write.csv(AT2.CC.top, paste0(table.dir, AT2.dir, "Top200_Cellular_component_AT2_ctrl_COPD_public.csv"), 
          quote = F)


# Overlapped terms in KEGG pathway
iAT2.KEGG <- enriched.iAT2$KEGG_2019_Human
iAT2.KEGG <- iAT2.KEGG[iAT2.KEGG$Adjusted.P.value < padj.cutoff,]
dim(iAT2.KEGG)
AT2.KEGG <- enriched.AT2$KEGG_2019_Human
AT2.KEGG <- AT2.KEGG[AT2.KEGG$Adjusted.P.value < padj.cutoff,]
dim(AT2.KEGG)
KEGG.enriched <- intersect(iAT2.KEGG$Term, AT2.KEGG$Term)
length(KEGG.enriched)
length(KEGG.enriched) / nrow(iAT2.KEGG) # 0.6
nrow(iAT2.KEGG) - length(KEGG.enriched)
nrow(AT2.KEGG) - length(KEGG.enriched)


iAT2.KEGG.top <- enriched.iAT2.top$KEGG_2019_Human
iAT2.KEGG.top <- iAT2.KEGG.top[iAT2.KEGG.top$Adjusted.P.value < padj.cutoff,]
dim(iAT2.KEGG.top)
AT2.KEGG.top <- enriched.AT2.top$KEGG_2019_Human
AT2.KEGG.top <- AT2.KEGG.top[AT2.KEGG.top$Adjusted.P.value < padj.cutoff,]
dim(AT2.KEGG.top)
KEGG.enriched.top <- intersect(iAT2.KEGG.top$Term, AT2.KEGG.top$Term)
length(KEGG.enriched.top)
length(KEGG.enriched.top) / nrow(iAT2.KEGG.top) # 0.5833333
nrow(iAT2.KEGG.top) - length(KEGG.enriched.top)
nrow(AT2.KEGG.top) - length(KEGG.enriched.top)


iAT2.KEGG <- cbind(Overlap = iAT2.KEGG$Term %in% KEGG.enriched, 
                 iAT2.KEGG)
dim(iAT2.KEGG)
colnames(iAT2.KEGG)
head(iAT2.KEGG$Overlap)
write.csv(iAT2.KEGG, paste0(table.dir, AT2.dir, "KEGG_iAT2_clusters_1_2_iPSC.csv"), 
          quote = F)


iAT2.KEGG.top <- cbind(Overlap = iAT2.KEGG.top$Term %in% KEGG.enriched.top, 
                   iAT2.KEGG.top)
dim(iAT2.KEGG.top)
colnames(iAT2.KEGG.top)
head(iAT2.KEGG.top$Overlap)
write.csv(iAT2.KEGG.top, paste0(table.dir, AT2.dir, "Top200_KEGG_iAT2_clusters_1_2_iPSC.csv"), 
          quote = F)


AT2.KEGG <- cbind(Overlap = AT2.KEGG$Term %in% KEGG.enriched, 
                AT2.KEGG)
dim(AT2.KEGG)
colnames(AT2.KEGG)
head(AT2.KEGG$Overlap)
write.csv(AT2.KEGG, paste0(table.dir, AT2.dir, "KEGG_AT2_ctrl_COPD_public.csv"), 
          quote = F)


AT2.KEGG.top <- cbind(Overlap = AT2.KEGG.top$Term %in% KEGG.enriched.top, 
                  AT2.KEGG.top)
dim(AT2.KEGG.top)
colnames(AT2.KEGG.top)
head(AT2.KEGG.top$Overlap)
write.csv(AT2.KEGG.top, paste0(table.dir, AT2.dir, "Top200_KEGG_AT2_ctrl_COPD_public.csv"), 
          quote = F)


# Save the DEGs of AT2
overlapped.DEGs.AT2 <- intersect(rownames(DEGs.AT2.ctrl.copd.filtered), 
                                 DEGs.iAT2.healthy.filtered$gene)
length(overlapped.DEGs.AT2)
DEGs.AT2.ctrl.copd.filtered <- cbind(Overlap = rownames(DEGs.AT2.ctrl.copd.filtered) %in% 
                                       overlapped.DEGs.AT2, 
                                     DEGs.AT2.ctrl.copd.filtered)
colnames(DEGs.AT2.ctrl.copd.filtered)
head(DEGs.AT2.ctrl.copd.filtered$Overlap)
dir.create(paste0(table.dir, AT2.dir))
write.csv(DEGs.AT2.ctrl.copd.filtered, paste0(table.dir, AT2.dir, "DEGs_AT2_ctrl_vs_COPD_public.csv"), 
          quote = F)
DEGs.iAT2.healthy.filtered <- cbind(Overlap = DEGs.iAT2.healthy.filtered$gene %in% overlapped.DEGs.AT2, 
                                    DEGs.iAT2.healthy.filtered)
colnames(DEGs.iAT2.healthy.filtered)
head(DEGs.iAT2.healthy.filtered$Overlap)
write.csv(DEGs.iAT2.healthy.filtered, paste0(table.dir, AT2.dir, "DEGs_iAT2_clusters_1_2_iPSC.csv"), 
          quote = F)


# Save the DEGs of AT2 (top 200)
overlapped.DEGs.AT2.top <- intersect(rownames(DEGs.AT2.ctrl.copd.top), 
                                 DEGs.iAT2.healthy.top$gene)
length(overlapped.DEGs.AT2.top)
DEGs.AT2.ctrl.copd.top <- cbind(Overlap = rownames(DEGs.AT2.ctrl.copd.top) %in% 
                                       overlapped.DEGs.AT2.top, 
                                     DEGs.AT2.ctrl.copd.top)
colnames(DEGs.AT2.ctrl.copd.top)
head(DEGs.AT2.ctrl.copd.top$Overlap)
dir.create(paste0(table.dir, AT2.dir))
write.csv(DEGs.AT2.ctrl.copd.top, paste0(table.dir, AT2.dir, "Top200_DEGs_AT2_ctrl_vs_COPD_public.csv"), 
          quote = F)


DEGs.iAT2.healthy.top <- cbind(Overlap = DEGs.iAT2.healthy.top$gene %in% overlapped.DEGs.AT2.top, 
                                    DEGs.iAT2.healthy.top)
colnames(DEGs.iAT2.healthy.top)
dim(DEGs.iAT2.healthy.top)
head(DEGs.iAT2.healthy.top$Overlap)
write.csv(DEGs.iAT2.healthy.top, paste0(table.dir, AT2.dir, "Top200_DEGs_iAT2_clusters_1_2_iPSC.csv"), 
          quote = F)
