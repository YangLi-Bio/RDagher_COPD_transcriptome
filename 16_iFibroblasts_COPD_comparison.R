####################################################################
#                                                                  #
# Compare DEGs in iFibroblasts v.s. Fibroblasts in the public      #
#                                                                  #
####################################################################


# Libraries
library(Seurat)
library(dplyr)
library(future)
library(pbapply)


# Global parameters
R.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Rfiles/"
table.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Tables/"
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
padj.cutoff <- 0.05
logFC.cutoff <- 0


# Local parameters
iFibr.dir <- "iFibroblasts/Mahboobe_Box/"
obj.file <- "Obj_AT2_ctrl_COPD.qsave"
Fibr.dir <- "Fibr_COPD_cluster_specific/"


# Source scripts
source(paste0(tool.dir, "combinatorial_tools.R"))
source(paste0(tool.dir, "transcriptome_tools.R"))


DEGs.iFibr.list <- qs::qread(paste0(R.dir, "DEGs_iFibr_COPD_clusters_045.qsave"))
length(DEGs.iFibr.list)
colnames(DEGs.iFibr.list[[1]])
sapply(DEGs.iFibr.list, dim)


# Load data from the public
obj <- qs::qread(paste0(R.dir, obj.file))
dim(obj)
colnames(obj@meta.data)
unique(obj$Disease)
unique(obj$cell.type) # the 7-th cluster


# Identify DEGs between COPD v.s. ctrl in Fibroblasts
DEGs.Fibr.COPD <- qs::qread(paste0(R.dir, "DEGs_Fibr_COPD_public.qsave"))
dim(DEGs.Fibr.COPD)
colnames(DEGs.Fibr.COPD)


# Compare DEGs between each of clusters 0, 4, and 5 in iFibr v.s. Fibr
for (i in seq_along(DEGs.iFibr.list)) {
  DEGs.iFibr.COPD <- DEGs.iFibr.list[[i]]
  DEGs.commom.pri <- intersect(unique(DEGs.iFibr.COPD$gene), rownames(DEGs.Fibr.COPD))
  head(DEGs.iFibr.COPD)
  length(DEGs.commom.pri)
  length(DEGs.commom.pri) / length(unique(DEGs.iFibr.COPD$gene))
  
  
  DEGs.iFibr.COPD.filtered <- DEGs.iFibr.COPD[DEGs.iFibr.COPD$gene %in% DEGs.commom.pri,]
  dim(DEGs.iFibr.COPD.filtered)
  DEGs.Fibr.COPD.filtered <- DEGs.Fibr.COPD[rownames(DEGs.Fibr.COPD) %in% DEGs.commom.pri,]
  dim(DEGs.Fibr.COPD.filtered)
  
  
  DEGs.iFibr.COPD.filtered <- DEGs.iFibr.COPD.filtered[order(DEGs.iFibr.COPD.filtered$pval),]
  DEGs.Fibr.COPD.filtered <- DEGs.Fibr.COPD.filtered[order(DEGs.Fibr.COPD.filtered$p_val),]
  dim(DEGs.Fibr.COPD.filtered)
  dim(DEGs.iFibr.COPD.filtered)
  length(intersect(unique(DEGs.iFibr.COPD.filtered$gene), rownames(DEGs.Fibr.COPD.filtered)))
  
  
  DEGs.iFibr.COPD.sig <- DEGs.iFibr.COPD.filtered[DEGs.iFibr.COPD.filtered$adj_pval < padj.cutoff,]
  DEGs.iFibr.COPD.sig <- DEGs.iFibr.COPD.sig[!duplicated(DEGs.iFibr.COPD.sig$gene),]
  dim(DEGs.iFibr.COPD.sig)
  head(DEGs.iFibr.COPD.sig)
  
  
  DEGs.Fibr.COPD.sig <- DEGs.Fibr.COPD.filtered[DEGs.Fibr.COPD.filtered$p_val_adj < padj.cutoff,]
  dim(DEGs.Fibr.COPD.sig)
  
  
  DEGs.overlap <- intersect(DEGs.iFibr.COPD.sig$gene, 
                            rownames(DEGs.Fibr.COPD.sig)) # only seven intersected DEGs
  length(DEGs.overlap)
  nrow(DEGs.iFibr.COPD.sig)
  length(DEGs.overlap) / nrow(DEGs.iFibr.COPD.sig)
  overlap.ratio <- get_overlap_ratio(x = DEGs.iFibr.COPD.sig$gene,
                                     y = rownames(DEGs.Fibr.COPD.sig))
}





# Retain the common DEGs included in the DEG lists of iFibr and Fibr
DEGs.iFibr.COPD.filtered <- DEGs.iFibr.COPD[DEGs.iFibr.COPD$gene %in% DEGs.commom.pri,]
dim(DEGs.iFibr.COPD.filtered)
DEGs.Fibr.COPD.filtered <- DEGs.Fibr.COPD[rownames(DEGs.Fibr.COPD) %in% DEGs.commom.pri,]
dim(DEGs.Fibr.COPD.filtered)


# Rank the DEGs
DEGs.iFibr.COPD.filtered <- DEGs.iFibr.COPD.filtered[order(DEGs.iFibr.COPD.filtered$pval),]
DEGs.Fibr.COPD.filtered <- DEGs.Fibr.COPD.filtered[order(DEGs.Fibr.COPD.filtered$p_val),]
dim(DEGs.iFibr.COPD.filtered)


# Retain significant DEGs
# DEGs.iFibr.COPD.sig <- DEGs.iFibr.COPD.filtered[DEGs.iFibr.COPD.filtered$volcano_adjp < padj.cutoff & 
#                                                        DEGs.iFibr.COPD.filtered$volcano_logfc > logFC.cutoff,]
DEGs.iFibr.COPD.sig <- DEGs.iFibr.COPD.filtered[DEGs.iFibr.COPD.filtered$adj_pval < padj.cutoff,]
DEGs.iFibr.COPD.sig <- DEGs.iFibr.COPD.sig[!duplicated(DEGs.iFibr.COPD.sig$gene),]
dim(DEGs.iFibr.COPD.sig)
head(DEGs.iFibr.COPD.sig)
# DEGs.Fibr.COPD.sig <- DEGs.Fibr.COPD.filtered[DEGs.Fibr.COPD.filtered$p_val_adj < padj.cutoff & 
#                                                 DEGs.Fibr.COPD.filtered$avg_log2FC > logFC.cutoff,]
DEGs.Fibr.COPD.sig <- DEGs.Fibr.COPD.filtered[DEGs.Fibr.COPD.filtered$p_val_adj < padj.cutoff,]
dim(DEGs.Fibr.COPD.sig)
DEGs.overlap <- intersect(DEGs.iFibr.COPD.sig$gene, 
                                 rownames(DEGs.Fibr.COPD.sig)) # only seven intersected DEGs
overlap.ratio <- get_overlap_ratio(x = DEGs.iFibr.COPD.sig$gene, 
                                   y = rownames(DEGs.Fibr.COPD.sig))
plot(overlap.ratio, xlab = "Top-ranked DEGs", ylab = "Overlap ratio")
length(overlap.ratio)
length(DEGs.overlap)
length(DEGs.overlap) / nrow(DEGs.iFibr.COPD.sig)
nrow(DEGs.iFibr.COPD.sig) - length(DEGs.overlap)
nrow(DEGs.Fibr.COPD.sig) - length(DEGs.overlap)
# The maximum overlap ratio, 0.382352941176471, is achieved at 476.
overlap.ratio[476]
dir.create(paste0(table.dir, Fibr.dir))


overlap.ratio[200]# 0.36
length(DEGs.overlap) / nrow(DEGs.iFibr.COPD.sig)
nrow(DEGs.iFibr.COPD.sig)
nrow(DEGs.Fibr.COPD.sig)
write.csv(cbind(Overlap = DEGs.iFibr.COPD.sig$gene %in% DEGs.overlap, 
                DEGs.iFibr.COPD.sig), 
          paste0(table.dir, Fibr.dir, "DEGs_iFibr_COPD_clusters045.csv"), 
          quote = F)
write.csv(cbind(Overlap = rownames(DEGs.Fibr.COPD.sig) %in% DEGs.overlap, 
                DEGs.Fibr.COPD.sig), 
          paste0(table.dir, Fibr.dir, "DEGs_Fibr_COPD_public.csv"), 
          quote = F)


# Retain top-200 DEGs
DEGs.iFibr.COPD.top <- DEGs.iFibr.COPD.sig[1:200,]
DEGs.Fibr.COPD.top <- DEGs.Fibr.COPD.sig[1:200,]
overlap.top <- intersect(DEGs.iFibr.COPD.top$gene, rownames(DEGs.Fibr.COPD.top))
length(overlap.top)
length(overlap.top) / nrow(DEGs.iFibr.COPD.top)
nrow(DEGs.iFibr.COPD.top) - length(overlap.top)
nrow(DEGs.Fibr.COPD.top) - length(overlap.top)
write.csv(cbind(Overlap = DEGs.iFibr.COPD.top$gene %in% overlap.top, 
                DEGs.iFibr.COPD.top), 
          paste0(table.dir, Fibr.dir, "Top200_DEGs_iFibr_clusters_045.csv"), 
          quote = F)
write.csv(cbind(Overlap = rownames(DEGs.Fibr.COPD.top) %in% overlap.top, 
                DEGs.Fibr.COPD.top), 
          paste0(table.dir, Fibr.dir, "Top200_DEGs_Fibr_COPD_public.csv"), 
          quote = F)


# Enrichment analyses for all DEGs
iFibr.enriched <- run_GO_and_KEGG(genes.ll = DEGs.iFibr.COPD.top$gene, org = "human")
Fibr.enriched <- run_GO_and_KEGG(genes.ll = rownames(DEGs.Fibr.COPD.sig), org = "human")
qs::qsave(iFibr.enriched, paste0(R.dir, "Pathways_iFibr_COPD_clusters045.qsave"))
qs::qsave(Fibr.enriched, paste0(R.dir, "Pathways_Fibr_COPD_public.qsave"))


# Enrichment analyses for top-200 DEGs
iFibr.enriched.top <- run_GO_and_KEGG(genes.ll = DEGs.iFibr.COPD.top$gene, org = "human")
Fibr.enriched.top <- run_GO_and_KEGG(genes.ll = rownames(DEGs.Fibr.COPD.top), org = "human")
qs::qsave(iFibr.enriched, paste0(R.dir, "Pathways_iFibr_COPD_clusters045_top200.qsave"))
qs::qsave(Fibr.enriched, paste0(R.dir, "Pathways_Fibr_COPD_public_top200.qsave"))


# Overlaps of molecular function
iFibr.MF <- iFibr.enriched$GO_Molecular_Function_2018
iFibr.MF <- iFibr.MF[iFibr.MF$Adjusted.P.value < padj.cutoff,]
dim(iFibr.MF)
Fibr.MF <- Fibr.enriched$GO_Molecular_Function_2018
Fibr.MF <- Fibr.MF[Fibr.MF$Adjusted.P.value < padj.cutoff,]
dim(Fibr.MF)
overlap.MF <- intersect(iFibr.MF$Term, Fibr.MF$Term)
write.csv(cbind(Overlap = iFibr.MF$Term %in% overlap.MF, 
                iFibr.MF), 
          paste0(table.dir, Fibr.dir, "iFibr_molecular_function.csv"), 
          quote = F)
write.csv(cbind(Overlap = Fibr.MF$Term %in% overlap.MF, 
                Fibr.MF), 
          paste0(table.dir, Fibr.dir, "Fibr_molecular_function.csv"), 
          quote = F)
length(overlap.MF)
length(overlap.MF) / nrow(iFibr.MF)
nrow(iFibr.MF) - length(overlap.MF)
nrow(Fibr.MF) - length(overlap.MF)


# Overlaps of cellular component
iFibr.CC <- iFibr.enriched$GO_Cellular_Component_2018
iFibr.CC <- iFibr.CC[iFibr.CC$Adjusted.P.value < padj.cutoff,]
dim(iFibr.CC)
Fibr.CC <- Fibr.enriched$GO_Cellular_Component_2018
Fibr.CC <- Fibr.CC[Fibr.CC$Adjusted.P.value < padj.cutoff,]
dim(Fibr.CC)
overlap.CC <- intersect(iFibr.CC$Term, Fibr.CC$Term)
write.csv(cbind(Overlap = iFibr.CC$Term %in% overlap.CC, 
                iFibr.CC), 
          paste0(table.dir, Fibr.dir, "iFibr_cellular_component.csv"), 
          quote = F)
write.csv(cbind(Overlap = Fibr.CC$Term %in% overlap.CC, 
                Fibr.CC), 
          paste0(table.dir, Fibr.dir, "Fibr_cellular_component.csv"), 
          quote = F)
length(overlap.CC)
length(overlap.CC) / nrow(iFibr.CC)
nrow(iFibr.CC) - length(overlap.CC)
nrow(Fibr.CC) - length(overlap.CC)


# Overlaps of biological process
iFibr.BP <- iFibr.enriched$GO_Biological_Process_2018
iFibr.BP <- iFibr.BP[iFibr.BP$Adjusted.P.value < padj.cutoff,]
dim(iFibr.BP)
Fibr.BP <- Fibr.enriched$GO_Biological_Process_2018
Fibr.BP <- Fibr.BP[Fibr.BP$Adjusted.P.value < padj.cutoff,]
dim(Fibr.BP)
overlap.BP <- intersect(iFibr.BP$Term, Fibr.BP$Term)
write.csv(cbind(Overlap = iFibr.BP$Term %in% overlap.BP, 
                iFibr.BP), 
          paste0(table.dir, Fibr.dir, "iFibr_biological_process.csv"), 
          quote = F)
write.csv(cbind(Overlap = Fibr.BP$Term %in% overlap.BP, 
                Fibr.BP), 
          paste0(table.dir, Fibr.dir, "Fibr_biological_process.csv"), 
          quote = F)
length(overlap.BP)
length(overlap.BP) / nrow(iFibr.BP)
nrow(iFibr.BP) - length(overlap.BP)
nrow(Fibr.BP) - length(overlap.BP)


# Overlaps of KEGG pathways
iFibr.KEGG <- iFibr.enriched$KEGG_2019_Human
iFibr.KEGG <- iFibr.KEGG[iFibr.KEGG$Adjusted.P.value < padj.cutoff,]
dim(iFibr.KEGG)
Fibr.KEGG <- Fibr.enriched$KEGG_2019_Human
Fibr.KEGG <- Fibr.KEGG[Fibr.KEGG$Adjusted.P.value < padj.cutoff,]
dim(Fibr.KEGG)
overlap.KEGG <- intersect(iFibr.KEGG$Term, Fibr.KEGG$Term)
write.csv(cbind(Overlap = iFibr.KEGG$Term %in% overlap.KEGG, 
                iFibr.KEGG), 
          paste0(table.dir, Fibr.dir, "iFibr_KEGG.csv"), 
          quote = F)
write.csv(cbind(Overlap = Fibr.KEGG$Term %in% overlap.KEGG, 
                Fibr.KEGG), 
          paste0(table.dir, Fibr.dir, "Fibr_KEGG.csv"), 
          quote = F)
length(overlap.KEGG)
length(overlap.KEGG) / nrow(iFibr.KEGG)
nrow(iFibr.KEGG) - length(overlap.KEGG)
nrow(Fibr.KEGG) - length(overlap.KEGG)









# Top200 DEGs


# Overlaps of molecular function (top-200)
iFibr.MF.top <- iFibr.enriched.top$GO_Molecular_Function_2018
iFibr.MF.top <- iFibr.MF.top[iFibr.MF.top$Adjusted.P.value < padj.cutoff,]
dim(iFibr.MF.top)
Fibr.MF.top <- Fibr.enriched.top$GO_Molecular_Function_2018
Fibr.MF.top <- Fibr.MF.top[Fibr.MF.top$Adjusted.P.value < padj.cutoff,]
dim(Fibr.MF.top)
overlap.MF.top <- intersect(iFibr.MF.top$Term, Fibr.MF.top$Term)
length(overlap.MF.top)
length(overlap.MF.top) / nrow(iFibr.MF.top)
nrow(iFibr.MF.top) - length(overlap.MF.top)
nrow(Fibr.MF.top) - length(overlap.MF.top)
write.csv(cbind(Overlap = iFibr.MF.top$Term %in% overlap.MF.top, 
                iFibr.MF.top), 
          paste0(table.dir, Fibr.dir, "iFibr_molecular_function_top200.csv"), 
          quote = F)
write.csv(cbind(Overlap = Fibr.MF.top$Term %in% overlap.MF.top, 
                Fibr.MF.top), 
          paste0(table.dir, Fibr.dir, "Fibr_molecular_function_top200.csv"), 
          quote = F)


# Overlaps of cellular component (top-200)
iFibr.CC.top <- iFibr.enriched.top$GO_Cellular_Component_2018
iFibr.CC.top <- iFibr.CC.top[iFibr.CC.top$Adjusted.P.value < padj.cutoff,]
dim(iFibr.CC.top)
Fibr.CC.top <- Fibr.enriched.top$GO_Cellular_Component_2018
Fibr.CC.top <- Fibr.CC.top[Fibr.CC.top$Adjusted.P.value < padj.cutoff,]
dim(Fibr.CC.top)
overlap.CC.top <- intersect(iFibr.CC.top$Term, Fibr.CC.top$Term)
length(overlap.CC.top)
length(overlap.CC.top) / nrow(iFibr.CC.top)
nrow(iFibr.CC.top) - length(overlap.CC.top)
nrow(Fibr.CC.top) - length(overlap.CC.top)

write.csv(cbind(Overlap = iFibr.CC.top$Term %in% overlap.CC.top, 
                iFibr.CC.top), 
          paste0(table.dir, Fibr.dir, "iFibr_cellular_component_top200.csv"), 
          quote = F)
write.csv(cbind(Overlap = Fibr.CC.top$Term %in% overlap.CC.top, 
                Fibr.CC.top), 
          paste0(table.dir, Fibr.dir, "Fibr_cellular_component_top200.csv"), 
          quote = F)


# Overlaps of biological process
iFibr.BP.top <- iFibr.enriched.top$GO_Biological_Process_2018
iFibr.BP.top <- iFibr.BP.top[iFibr.BP.top$Adjusted.P.value < padj.cutoff,]
dim(iFibr.BP.top)
Fibr.BP.top <- Fibr.enriched.top$GO_Biological_Process_2018
Fibr.BP.top <- Fibr.BP.top[Fibr.BP.top$Adjusted.P.value < padj.cutoff,]
dim(Fibr.BP.top)
overlap.BP.top <- intersect(iFibr.BP.top$Term, Fibr.BP.top$Term)
write.csv(cbind(Overlap = iFibr.BP.top$Term %in% overlap.BP.top, 
                iFibr.BP.top), 
          paste0(table.dir, Fibr.dir, "iFibr_biological_process_top200.csv"), 
          quote = F)
write.csv(cbind(Overlap = Fibr.BP.top$Term %in% overlap.BP.top, 
                Fibr.BP.top), 
          paste0(table.dir, Fibr.dir, "Fibr_biological_process_top200.csv"), 
          quote = F)
length(overlap.BP.top)
length(overlap.BP.top) / nrow(iFibr.BP.top)
nrow(iFibr.BP.top) - length(overlap.BP.top)
nrow(Fibr.BP.top) - length(overlap.BP.top)



# Overlaps of KEGG pathways (top-200)
iFibr.KEGG.top <- iFibr.enriched.top$KEGG_2019_Human
iFibr.KEGG.top <- iFibr.KEGG.top[iFibr.KEGG.top$Adjusted.P.value < padj.cutoff,]
dim(iFibr.KEGG.top)
Fibr.KEGG.top <- Fibr.enriched.top$KEGG_2019_Human
Fibr.KEGG.top <- Fibr.KEGG.top[Fibr.KEGG.top$Adjusted.P.value < padj.cutoff,]
dim(Fibr.KEGG.top)
overlap.KEGG.top <- intersect(iFibr.KEGG.top$Term, Fibr.KEGG.top$Term)
write.csv(cbind(Overlap = iFibr.KEGG.top$Term %in% overlap.KEGG.top, 
                iFibr.KEGG.top), 
          paste0(table.dir, Fibr.dir, "iFibr_KEGG_top200.csv"), 
          quote = F)
write.csv(cbind(Overlap = Fibr.KEGG.top$Term %in% overlap.KEGG.top, 
                Fibr.KEGG.top), 
          paste0(table.dir, Fibr.dir, "Fibr_KEGG_top200.csv"), 
          quote = F)
length(overlap.KEGG.top)
length(overlap.KEGG.top) / nrow(iFibr.KEGG.top)
nrow(iFibr.KEGG.top) - length(overlap.KEGG.top)
nrow(Fibr.KEGG.top) - length(overlap.KEGG.top)
