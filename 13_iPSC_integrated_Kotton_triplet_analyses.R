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


# Local parameters
iAT2.DEG.file <- "COPD_healthy_DEGs_iAT2.csv"
copd.meta.file <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/2_copd_Kaminski_Sci_Adv_2020/GSE136831_AllCells.Samples.CellType.MetadataTable.txt"
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
DEGs.iAT2.healthy.top200 <- head(DEGs.iAT2.healthy, n = 200)
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


# A small test for DEG discovery between ATII v.s. other cell types
obj.sampled <- subset(obj.ctrl, cells = colnames(obj.ctrl)[sample(1:ncol(obj.ctrl), 1000)])
Idents(obj.sampled) <- obj.sampled$cell.type
DEGs.AT2.sampled <- FindMarkers(obj.sampled, ident.1 = "ATII")
dim(DEGs.AT2.sampled)
qs::qsave(DEGs.AT2.sampled, paste0(R.dir, "DEGs_AT2_ctrl_sampled.qsave"))
DEGs.AT2.sampled <- DEGs.AT2.sampled[DEGs.AT2.sampled$p_val_adj < padj.cutoff & 
                                       DEGs.AT2.sampled$avg_log2FC > logFC.cutoff,]
dim(DEGs.AT2.sampled)


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
DEGs.AT2.ctrl.copd.filtered <- DEGs.AT2.ctrl.copd.filtered[rownames(DEGs.AT2.ctrl.copd.filtered) %in% 
                                                             DEGs.iAT2.healthy$gene,]
dim(DEGs.AT2.ctrl.copd.filtered)

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
length(overlaps.healthy.iAT2.AT2)
intersect(symbols.AT2.ctrl.copd.filtered, symbols.iAT2.healthy)
# intersect(rownames(DEGs.AT2.sampled), symbols.iAT2.healthy)
