################################################################
#                                                              #
#          DEG analysis and enrichment analysis in AT2         #
#                                                              #
################################################################


# Libraries
library(Seurat)
library(dplyr)
library(qs)
library(enrichR)
library(pbmcapply)


# Parameters
setwd("/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19")
data.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/1_covid-19_Izar_Nature_2021/10x_format"
cluster.info <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/1_covid-19_Izar_Nature_2021/lung_metaData_covid_control.txt"
embedding.file <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/1_covid-19_Izar_Nature_2021/lung_clusterfile.txt"
save.dir <- "covid-19_vs_healthy/"
dbs <-
  c(
    "GO_Molecular_Function_2018",
    "GO_Cellular_Component_2018",
    "GO_Biological_Process_2018",
    "KEGG_2019_Human"
  )


# Quality control
covid.data <- Read10X(data.dir = data.dir, gene.column = 1) # load data
dim(covid.data)
cell.types.ll <- readLines(cluster.info) %>% strsplit(., split = "\t") %>% 
  lapply(., "[", c(1, 12, 13, 14, 15, 16)) # three levels of cell types
length(cell.types.ll)
cell.types.ll <- cell.types.ll[-c(1, 2, 116316)]
length(cell.types.ll)
head(cell.types.ll)
tail(cell.types.ll)
covid <- CreateSeuratObject(counts = covid.data[, unlist(lapply(cell.types.ll, "[[", (1))), 
                                                drop = F], 
                            project = "covid-19", 
                            min.cells = 0, min.features = 0)
dim(covid)
covid[["percent.mt"]] <- PercentageFeatureSet(covid, pattern = "^MT-")
VlnPlot(covid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(covid, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(covid, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


# Normalizing the data
covid <- NormalizeData(covid)
dim(covid)


# Add cell types as meta data
cell.types <- as.factor(unlist(lapply(cell.types.ll, "[[", (4)))) # extract the level 2 cell types
names(cell.types) <- unlist(lapply(cell.types.ll, "[[", (1)))
cell.types
covid <- AddMetaData(object = covid, metadata = cell.types, col.name = "cell.type")


# Add groups as meta data
groups <- as.factor(unlist(lapply(cell.types.ll, "[[", (2)))) # extract the groups
names(groups) <- unlist(lapply(cell.types.ll, "[[", (1)))
groups
covid <- AddMetaData(object = covid, metadata = groups, col.name = "groups")


# Add UMAP embeddings
covid.embedding <- read.table(embedding.file)[-c(1, 2), ]
covid.embedding <- covid.embedding[covid.embedding$V1 %in% colnames(covid),]
numap <- as.matrix(covid.embedding[, 2:3])
rownames(numap) <- covid.embedding$V1
colnames(numap) <- paste0("umap_", 1:2)
covid@reductions[["umap"]] <- CreateDimReducObject(embeddings = numap, key = "umap_", assay = DefaultAssay(covid))
qs::qsave(covid, "covid_and_health.qsave")


# DEG analysis for AT2
Idents(covid) <- covid$cell.type
levels(Idents(covid))
dim(covid)
covid.at2 <- subset(covid, idents = "AT2")
dim(covid.at2)
DefaultAssay(covid.at2) <- "RNA"
Idents(covid.at2) <- covid.at2$groups
at2.deg <- FindMarkers(covid.at2, ident.1 = "COVID-19", ident.2 = "Control")
dim(at2.deg)
dir.create(save.dir)
write.csv(at2.deg, paste0(save.dir, "DEGs_AT2_COVID-19_vs_health.csv"))


# Pathway enrichment analyses for AT2
this_enriched <- enrichr(rownames(at2.deg), dbs)
write.csv(
  this_enriched$GO_Molecular_Function_2018[, c(-5, -6)],
  paste0(save.dir, "COVID-19_vs_health_AT2_GO_MF.csv")
)
write.csv(
  this_enriched$GO_Cellular_Component_2018[, c(-5, -6)],
  paste0(save.dir, "COVID-19_vs_health_AT2_GO_CC.csv")
)
write.csv(
  this_enriched$GO_Biological_Process_2018[, c(-5, -6)],
  paste0(save.dir, "COVID-19_vs_health_AT2_GO_BP.csv")
)
write.csv(
  this_enriched$KEGG_2019_Human[, c(-5, -6)],
  paste0(save.dir, "COVID-19_vs_health_AT2_KEGG.csv")
)


# DEG analysis for Fibroblasts
Idents(covid) <- covid$cell.type
levels(Idents(covid))
dim(covid)
covid.fibro <- subset(covid, idents = "Fibroblasts")
dim(covid.fibro)
DefaultAssay(covid.fibro) <- "RNA"
Idents(covid.fibro) <- covid.fibro$groups
fibro.deg <- FindMarkers(covid.fibro, ident.1 = "COVID-19", ident.2 = "Control")
dim(fibro.deg)
write.csv(fibro.deg, paste0(save.dir, "DEGs_Fibroblasts_COVID-19_vs_health.csv"))


# Pathway enrichment analyses for Fibroblasts
this_enriched <- enrichr(rownames(fibro.deg), dbs)
write.csv(
  this_enriched$GO_Molecular_Function_2018[, c(-5, -6)],
  paste0(save.dir, "COVID-19_vs_health_Fibroblasts_GO_MF.csv")
)
write.csv(
  this_enriched$GO_Cellular_Component_2018[, c(-5, -6)],
  paste0(save.dir, "COVID-19_vs_health_Fibroblasts_GO_CC.csv")
)
write.csv(
  this_enriched$GO_Biological_Process_2018[, c(-5, -6)],
  paste0(save.dir, "COVID-19_vs_health_Fibroblasts_GO_BP.csv")
)
write.csv(
  this_enriched$KEGG_2019_Human[, c(-5, -6)],
  paste0(save.dir, "COVID-19_vs_health_Fibroblasts_KEGG.csv")
)
