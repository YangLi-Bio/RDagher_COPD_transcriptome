#######################################################################
#                                                                     #
#           Analyses on scRNA-Seq data of COVID-19                    #
#                                                                     #
#######################################################################


# Libraries
library(Seurat)
library(dplyr)
library(qs)
library(enrichR)


# Parameters
setwd("/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19")
data.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/1_covid-19_Izar_Nature_2021/10x_format"
cluster.info <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/1_covid-19_Izar_Nature_2021/lung_metaData.txt"
embedding.file <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/1_covid-19_Izar_Nature_2021/lung_clusterfile.txt"
dbs <-
  c(
    "GO_Molecular_Function_2018",
    "GO_Cellular_Component_2018",
    "GO_Biological_Process_2018",
    "KEGG_2019_Human"
  )


# Quality control
covid.data <- Read10X(data.dir = data.dir, gene.column = 1) # load data
cell.types.ll <- readLines(cluster.info) %>% strsplit(., split = "\t") %>% 
                  lapply(., "[", c(1, 13, 14, 15, 16)) # three levels of cell types
cell.types.ll <- cell.types.ll[-c(1, 2, 79639)]
covid <- CreateSeuratObject(counts = covid.data[, unlist(lapply(cell.types.ll, "[[", (1)))], 
                            project = "covid-19", 
                            min.cells = 0, min.features = 0)
covid
covid[["percent.mt"]] <- PercentageFeatureSet(covid, pattern = "^MT-")
VlnPlot(covid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(covid, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(covid, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


# Normalizing the data
# covid <- NormalizeData(covid, normalization.method = "LogNormalize", scale.factor = 10000)
covid <- NormalizeData(covid)


# Add meta data
cell.types <- as.factor(unlist(lapply(cell.types.ll, "[[", (3)))) # extract the level 2 cell types
names(cell.types) <- unlist(lapply(cell.types.ll, "[[", (1)))
# covid1 <- covid # back up
covid <- AddMetaData(object = covid, metadata = cell.types, col.name = "cell.type")
covid.embedding <- read.table(embedding.file)[-c(1, 2), ]
covid.embedding <- covid.embedding[covid.embedding$V1 %in% colnames(covid),]
numap <- as.matrix(covid.embedding[, 2:3])
rownames(numap) <- covid.embedding$V1
colnames(numap) <- paste0("umap_", 1:2)
# numap[, 1] <- as.numeric(numap[, 1])
# numap[, 2] <- as.numeric(numap[, 2])
covid@reductions[["umap"]] <- CreateDimReducObject(embeddings = numap, key = "umap_", assay = DefaultAssay(covid))
qs::qsave(covid, "covid.qsave")


# DEG analyses
Idents(covid) <- covid$cell.type
DefaultAssay(covid) <- "RNA"
covid.at2.markers <- FindMarkers(covid, ident.1 = "AT2") # find DEGs in AT2
covid.fibro.markers <- FindMarkers(covid, ident.1 = "Fibroblasts") # find DEGs in Fibroblasts
qs::qsave(covid.at2.markers, "covid_AT2_degs.qsave")
qs::qsave(covid.fibro.markers, "covid_Fibroblasts_degs.qsave")


# Pathway enrichment analysis for AT2
this_cts_genes <- rownames(covid.at2.markers)
this_enriched <- enrichr(this_cts_genes, dbs)

write.csv(
  this_enriched$GO_Molecular_Function_2018[, c(-5, -6)],
  paste0("covid_AT2_GO_MF.csv")
)
write.csv(
  this_enriched$GO_Cellular_Component_2018[, c(-5, -6)],
  paste0("covid_AT2_GO_CC.csv")
)
write.csv(
  this_enriched$GO_Biological_Process_2018[, c(-5, -6)],
  paste0("covid_AT2_GO_BP.csv")
)
write.csv(
  this_enriched$KEGG_2019_Human[, c(-5, -6)],
  paste0("covid_AT2_KEGG.csv")
)


# Pathway enrichment analysis for Fibroblasts
this_cts_genes <- rownames(covid.fibro.markers)
this_enriched <- enrichr(this_cts_genes, dbs)

write.csv(
  this_enriched$GO_Molecular_Function_2018[, c(-5, -6)],
  paste0("covid_Fibroblasts_GO_MF.csv")
)
write.csv(
  this_enriched$GO_Cellular_Component_2018[, c(-5, -6)],
  paste0("covid_Fibroblasts_GO_CC.csv")
)
write.csv(
  this_enriched$GO_Biological_Process_2018[, c(-5, -6)],
  paste0("covid_Fibroblasts_GO_BP.csv")
)
write.csv(
  this_enriched$KEGG_2019_Human[, c(-5, -6)],
  paste0("covid_Fibroblasts_KEGG.csv")
)
