#######################################################################
#                                                                     #
#         Analyses on scRNA-Seq data of COPD and COVID-19             #
#                                                                     #
#######################################################################


# Libraries
library(Seurat)
library(dplyr)
library(qs)
library(enrichR)
library(ggplot2)
library(RColorBrewer) # Provides color schemes for maps
library(Polychrome) # Qualitative Palettes with Many Colors
library(stringr)
library(tidyverse)
library(cowplot)
library(clusterProfiler)


# Parameters
setwd("/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19")
copd.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/2_copd_Kaminski_Sci_Adv_2020/10x"
copd.meta.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/2_copd_Kaminski_Sci_Adv_2020/COPD.MetadataTable.txt"
dbs <-
  c(
    "GO_Molecular_Function_2018",
    "GO_Cellular_Component_2018",
    "GO_Biological_Process_2018",
    "KEGG_2019_Human"
  )


# Functions
Plot.cluster2D <- function(object = combined, reduction.method = "umap", 
                           customized = T, pt_size = 1, txt = "Cell type", ...) {
  
  my.plot.all.source <- cbind.data.frame(Embeddings(object,reduction = reduction.method),
                                         Cell_type = Idents(object))
  
  tmp.celltype <- levels(unique(my.plot.all.source$Cell_type))
  p.cluster <- ggplot(my.plot.all.source,
                      aes(x = my.plot.all.source[,1], y = my.plot.all.source[,2])) + 
    xlab(colnames(my.plot.all.source)[1]) + ylab(colnames(my.plot.all.source)[2])
  p.cluster <- p.cluster + geom_point(stroke = pt_size,size=pt_size, aes(col = my.plot.all.source[, "Cell_type"])) 
  p.cluster <- p.cluster + guides(colour = guide_legend(override.aes = list(size = 5)))
  
  if (length(tmp.celltype) >= 5){
    p.cluster <- p.cluster + scale_colour_manual(name = paste(txt, ":(Cells)", sep = ""), 
                                                 values = as.character(palette36.colors(36)[-2][1:length(tmp.celltype)]),
                                                 breaks = tmp.celltype,
                                                 labels = paste0(tmp.celltype, 
                                                                 ":(",as.character(summary(my.plot.all.source$Cell_type)),")"))
  } else if (length(tmp.celltype) < 5) {
    p.cluster <- p.cluster + scale_colour_manual(name = paste(txt,":(Cells)",sep = ""), 
                                                 values = brewer.pal(4,"Spectral")[c(2,1,3,4)],
                                                 breaks = tmp.celltype,
                                                 labels = paste0(tmp.celltype, 
                                                                 ":(",as.character(summary(my.plot.all.source$Cell_type)),")"))
  } else {
    p.cluster <- p.cluster + scale_colour_manual(name = paste(txt,":(Cells)",sep = ""), 
                                                 values = brewer.pal(5,"Spectral")[c(1,5)],
                                                 breaks = tmp.celltype,
                                                 labels = paste0(tmp.celltype, 
                                                                 ":(",as.character(summary(my.plot.all.source$Cell_type)), ")"))
    
  }
  
  # + labs(col="cell type")           
  p.cluster <- p.cluster + theme_classic() 
  p.cluster <- p.cluster + coord_fixed(ratio = 1)
  p.cluster
}


# Quality control
copd.data <- Read10X(data.dir = copd.dir, gene.column = 2) # load data
copd.meta <- read.table(copd.meta.dir)
copd.data <- copd.data[, copd.meta[, 1]]
copd <- CreateSeuratObject(counts = copd.data, project = "copd", 
                            min.cells = 0, min.features = 0)
copd
copd[["percent.mt"]] <- PercentageFeatureSet(copd, pattern = "^MT-")
VlnPlot(copd, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(copd, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(copd, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


# Normalizing the data
copd <- NormalizeData(copd)


# Add meta data
copd.meta.df <- copd.meta[, c(1, 4, 5, 6)]
cell.types <- as.factor(copd.meta.df[, 3]) # extract the level 2 cell types
names(cell.types) <- copd.meta.df[, 1]
copd <- AddMetaData(object = copd, metadata = cell.types, col.name = "cell.type")
qs::qsave(copd, "copd.qsave")
covid <- qs::qread("covid.qsave")


# Unify cell types
levels(covid$cell.type)
levels(copd$cell.type)
levels(copd$cell.type) <- c(levels(copd$cell.type)[1], "AT1", "AT2", "B cells", levels(copd$cell.type)[5:13], 
                            "Fibroblasts", levels(copd$cell.type)[15:25], "NK cells", 
                            levels(copd$cell.type)[27:38])
intersect(levels(covid$cell.type), levels(copd$cell.type))
qs::qsave(copd, "copd_new_cell_types.qsave")


# Integrate two datasets
obj.ll <- c(covid, copd) # list of Seurat objects
obj.ll <- lapply(X = obj.ll, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = obj.ll)
obj.anchors <- FindIntegrationAnchors(object.list = obj.ll, anchor.features = features)

# this command creates an 'integrated' data assay
integrated <- IntegrateData(anchorset = obj.anchors)
qs::qsave(integrated, "integrated.qsave")
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
qs::qsave(integrated, "integrated.qsave")


Idents(integrated) <- integrated$cell.type
# p1 <- Plot.cluster2D(integrated, pt_size = 0.4,
#                      txt = "Cell types")
p1 <- DimPlot(integrated, reduction = "umap", group.by = "cell.type")
png(
  paste("combined_cell_types.png", sep = ""),
  width = 3000,
  height = 2000,
  res = 300
)
p1
dev.off()


# Idents(integrated) <- integrated$orig.ident
# covid.cells <- as.factor(integrated$orig.ident[integrated$orig.ident == "covid-19"])
# copd.cell.ids <- setdiff(names(integrated$orig.ident), names(covid.cells))
# copd.cells <- as.factor(rep("copd", length(copd.cell.ids)))
# names(copd.cells) <- copd.cell.ids
# disease.ids <- c(covid.cells, copd.cells)
# integrated <- AddMetaData(object = integrated, metadata = disease.ids, col.name = "diseases")
# qs::qsave(integrated, "integrated.qsave")

Idents(integrated) <- integrated$diseases
# p2 <- Plot.cluster2D(integrated, pt_size = 0.4,
#                      txt = "Diseases")
p2 <- DimPlot(integrated, reduction = "umap", group.by = "diseases")
png(
  paste("combined_diseases.png", sep = ""),
  width = 3000,
  height = 2000,
  res = 300
)
p2
dev.off()


# DEG analyses on integrated data
Idents(integrated) <- integrated$cell.type
DefaultAssay(integrated) <- "RNA"
integrated.at2.markers <- FindMarkers(integrated, ident.1 = "AT2") # find DEGs in AT2
integrated.fibro.markers <- FindMarkers(integrated, ident.1 = "Fibroblasts") # find DEGs in Fibroblasts
qs::qsave(integrated.at2.markers, "integrated_AT2_degs.qsave")
qs::qsave(integrated.fibro.markers, "integrated_Fibroblasts_degs.qsave")


# Pathway enrichment analysis
this_cts_genes <- rownames(integrated.at2.markers)
this_enriched <- enrichr(this_cts_genes, dbs)

write.csv(
  this_enriched$GO_Molecular_Function_2018[, c(-5, -6)],
  paste0("integrated_AT2_GO_MF.csv")
)
write.csv(
  this_enriched$GO_Cellular_Component_2018[, c(-5, -6)],
  paste0("integrated_AT2_GO_CC.csv")
)
write.csv(
  this_enriched$GO_Biological_Process_2018[, c(-5, -6)],
  paste0("integrated_AT2_GO_BP.csv")
)
write.csv(
  this_enriched$KEGG_2019_Human[, c(-5, -6)],
  paste0("integrated_AT2_KEGG.csv")
)


# Pathway enrichment analysis for Fibroblasts
this_cts_genes <- rownames(integrated.fibro.markers)
this_enriched <- enrichr(this_cts_genes, dbs)

write.csv(
  this_enriched$GO_Molecular_Function_2018[, c(-5, -6)],
  paste0("integrated_Fibroblasts_GO_MF.csv")
)
write.csv(
  this_enriched$GO_Cellular_Component_2018[, c(-5, -6)],
  paste0("integrated_Fibroblasts_GO_CC.csv")
)
write.csv(
  this_enriched$GO_Biological_Process_2018[, c(-5, -6)],
  paste0("integrated_Fibroblasts_GO_BP.csv")
)
write.csv(
  this_enriched$KEGG_2019_Human[, c(-5, -6)],
  paste0("integrated_Fibroblasts_KEGG.csv")
)
