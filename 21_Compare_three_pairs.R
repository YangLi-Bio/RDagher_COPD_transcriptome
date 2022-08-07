############################################################
#                                                          #
# 1. Infected v.s. uninfected in COPD                      #
# 2. Infected v.s. uninfected in healthy                   #
# 3. COPD v.s. healthy in infected                         #
#                                                          #
############################################################


# Analysis plan:
# Public data: a COVID-19 & healthy dataset, a COPD & healthy dataset
# 1. Extract the AT2 and Fibroblast cells for both COVID-19 and COPD, respectively
# 2. Integrate the COVID-19 and COPD infected cells, leading to the COVID-19 & COPD cells
# 3. We have COVID-19 & COPD, COVID-19, COPD, and both control cells (four groups)
# 4. That is, we need to integrate COVID-19 and COPD cells


# Libraries
library(Seurat)
library(dplyr)
library(parallel)
library(pbmcapply)
library(future)
library(readxl)


# Global variables
R.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Rfiles/"
table.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Tables/Three_comparison_v2/"
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
coculture.file <- "Supplemental_11_Figure_6J_DEG_results.xlsx"
padj.cutoff <- 0.05
FC.cutoff <- 0.25


# Load the COVID-19 data
obj.covid <- qs::qread(paste0(R.dir, "covid_and_health.qsave"))
dim(obj.covid)
unique(obj.covid$cell.type)
unique(obj.covid$groups)
Idents(obj.covid) <- obj.covid$cell.type
obj.covid <- subset(obj.covid, idents = c("AT2", "Fibroblasts"))
table(obj.covid$cell.type)
table(obj.covid$groups)
Idents(obj.covid) <- obj.covid$groups
covid.meta.list <- SplitObject(obj.covid, split.by = "ident")
covid.meta.list
covid.infected <- covid.meta.list[[2]]
covid.ctrl <- covid.meta.list[[1]]
dim(covid.infected)
dim(covid.ctrl)


# Load the COPD data
obj.copd <- qs::qread(paste0(R.dir, "Obj_COPD_healthy.qsave"))
dim(obj.copd)
head(colnames(obj.copd))
table(obj.copd$cell.type)
table(obj.copd$Disease)
Idents(obj.copd) <- obj.copd$cell.type
obj.copd <- subset(obj.copd, idents = c("ATII", "Fibroblast"))
dim(obj.copd)
table(obj.copd$cell.type)
table(obj.copd$Disease)
Idents(obj.copd) <- obj.copd$Disease
obj.copd.list <- SplitObject(obj.copd, split.by = "Disease")
obj.copd.list
copd.infected <- obj.copd.list[[3]]
copd.ctrl <- obj.copd.list[[1]]
dim(copd.infected)
dim(copd.ctrl)


# Integrate the COVID-19 and COPD infected cells to construct COPD followed by COVID-19
infected.list <- list(covid.infected = covid.infected, 
                      copd.infected = copd.infected)
sapply(infected.list, dim)
infected.list <- lapply(X = infected.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = infected.list)
length(features)
infected.list <- lapply(X = infected.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = infected.list, reduction = "rpca",
                                  dims = 1:50)
infected.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
dim(infected.integrated)
sapply(infected.list, dim)
table(infected.integrated$orig.ident)
table(infected.integrated$cell.type)
table(infected.integrated$Disease)
Disease <- gsub("covid-19", "covid.infected", 
                infected.integrated$orig.ident)
Disease <- gsub("copd", "copd.infected", Disease)
table(Disease)
infected.integrated <- AddMetaData(infected.integrated, metadata = Disease, col.name = "Disease")
table(infected.integrated$Disease)
qs::qsave(infected.integrated, paste0(R.dir, "Obj_double_infected_COVID-19_COPD.qsave"))


# Integrate the two control datasets
ctrl.list <- list(covid.ctrl = covid.ctrl, 
                      copd.ctrl = copd.ctrl)
sapply(ctrl.list, dim)
ctrl.list <- lapply(X = ctrl.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = ctrl.list)
length(features)
ctrl.list <- lapply(X = ctrl.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = ctrl.list, reduction = "rpca",
                                  dims = 1:50)
ctrl.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
dim(ctrl.integrated)
sapply(ctrl.list, dim)
table(ctrl.integrated$orig.ident)
table(ctrl.integrated$cell.type)
table(ctrl.integrated$Disease)
Disease <- gsub("covid-19", "covid.ctrl", 
                ctrl.integrated$orig.ident)
Disease <- gsub("copd", "copd.ctrl", Disease)
table(Disease)
ctrl.integrated <- AddMetaData(ctrl.integrated, metadata = Disease, col.name = "Disease")
table(ctrl.integrated$Disease)
qs::qsave(ctrl.integrated, paste0(R.dir, "Obj_double_ctrl_COVID-19_COPD.qsave"))


# Integrate the double infected, COVID-19, COPD, and double control
covid.infected <- RenameCells(covid.infected, add.cell.id = "covid")
copd.infected <- RenameCells(copd.infected, add.cell.id = "copd")
head(colnames(covid.infected))
head(colnames(copd.infected))
infected.integrated <- RenameCells(infected.integrated, add.cell.id = "Infected")
colnames(infected.integrated)
ctrl.integrated <- RenameCells(ctrl.integrated, add.cell.id = "Ctrl")
colnames(ctrl.integrated)
DefaultAssay(infected.integrated) <- "RNA"
DefaultAssay(ctrl.integrated) <- "RNA"
DefaultAssay(covid.infected) <- "RNA"
DefaultAssay(copd.infected) <- "RNA"
quad.list <- list(Infected = infected.integrated, 
                  covid = covid.infected, 
                  copd = copd.infected, 
                  Ctrl = ctrl.integrated)
sapply(quad.list, dim)
quad.list <- lapply(X = quad.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = quad.list)
length(features)
quad.list <- lapply(X = quad.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = quad.list, reduction = "rpca",
                                  dims = 1:50)
quad.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
dim(quad.integrated)
sapply(quad.list, ncol) %>% sum
table(quad.integrated$orig.ident)
table(quad.integrated$cell.type)
table(quad.integrated$Disease)
colnames(quad.integrated) %>% head
colnames(quad.integrated)[30000]
quad.integrated <- AddMetaData(quad.integrated, metadata = Disease, col.name = "Disease")
table(quad.integrated$Disease)
Disease <- strsplit(colnames(quad.integrated), split = "_") %>% sapply(., "[[", 1)
length(Disease) == ncol(quad.integrated)
names(Disease) <- colnames(quad.integrated)
head(Disease)
qs::qsave(quad.integrated, paste0(R.dir, "Obj_infected_covid_copd_ctrl.qsave"))
