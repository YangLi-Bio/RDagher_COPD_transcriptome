#######################################################################
#                                                                     #
#           Find DEGs beteween COPD and healthy in AT2                #
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

dbs <-
  c(
    "GO_Molecular_Function_2018",
    "GO_Cellular_Component_2018",
    "GO_Biological_Process_2018",
    "KEGG_2019_Human"
  )


# Load data
setwd("/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19")
integrated <- qs::qread("integrated.qsave")
Idents(integrated) <- integrated$cell.type
integrated.AT2 <- subset(integrated, idents = "AT2")
ncol(integrated.AT2)
covid <- qs::qread("covid.qsave")
Idents(covid) <- covid$cell.type
covid.AT2 <- subset(covid, idents = "AT2")


# Set meta data
integrated.AT2.id <- rep("copd", ncol(integrated.AT2))
names(integrated.AT2.id) <- colnames(integrated.AT2)
integrated.AT2 <- AddMetaData(integrated.AT2, metadata = integrated.AT2.id, col.name = "new.id")
covid.AT2.id <- rep("healthy", ncol(covid.AT2))
names(covid.AT2.id) <- colnames(covid.AT2)
covid.AT2 <- AddMetaData(covid.AT2, metadata = covid.AT2.id, col.name = "new.id")
covid.AT2 <- RenameCells(covid.AT2, add.cell.id = "C")
colnames(covid.AT2)
qs::qsave(covid.AT2, "covid_AT2_renamed.qsave")


# Integration
DefaultAssay(integrated.AT2) <- "RNA"
DefaultAssay(covid.AT2) <- "RNA"
AT2.ll <- c(covid.AT2, integrated.AT2)
AT2.ll <- lapply(X = AT2.ll, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = AT2.ll)
obj.anchors <- FindIntegrationAnchors(object.list = AT2.ll, anchor.features = features)
integrated.again.AT2 <- IntegrateData(anchorset = obj.anchors)


# DEG analysis
DefaultAssay(integrated.again.AT2) <- "RNA"
Idents(integrated.again.AT2) <- integrated.again.AT2$new.id
copd.healthy.DEGs <- FindMarkers(integrated.again.AT2, ident.1 = "copd", ident.2 = "healthy", verbose = FALSE)
head(copd.healthy.DEGs)
selected.DEGs <- copd.healthy.DEGs[copd.healthy.DEGs$p_val_adj < 0.01,]
write.csv(copd.healthy.DEGs, "COPD_healthy_DEGs.csv")
iAT2.DEGs <- read.csv("COPD_healthy_DEGs_iAT2.csv")
head(iAT2.DEGs)
overlapped.genes <- intersect(rownames(copd.healthy.DEGs), unique(iAT2.DEGs$gene))
length(overlapped.genes)
DEG.pval <- phyper(length(overlapped.genes) - 1, length(unique(iAT2.DEGs$gene)), 
       nrow(integrated[['RNA']]) - length(unique(iAT2.DEGs$gene)), 
       nrow(copd.healthy.DEGs), lower.tail = F)


# Pathway enrichment analyses
this_enriched <- enrichr(rownames(copd.healthy.DEGs), dbs)
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
AT2.kegg <- this_enriched$KEGG_2019_Human[, c(-5, -6)]
iAT2.pathways <-read.csv("COPD_healthy_pathways_iAT2.csv")
head(iAT2.pathways)
overlap.pathways <- intersect(unique(AT2.kegg$Term), unique(iAT2.pathways$Term))
length(overlap.pathways)

pathway.pval <- phyper(length(overlap.pathways) - 1, length(unique(iAT2.pathways$Term)), 
                       308 - length(unique(iAT2.pathways$Term)), 
                       length(unique(AT2.kegg$Term)), lower.tail = F)


#######################################################################
#                                                                     #
#        Find DEGs beteween COPD and healthy in fibroblasts           #
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
library(readxl) # read xlsx files

dbs <-
  c(
    "GO_Molecular_Function_2018",
    "GO_Cellular_Component_2018",
    "GO_Biological_Process_2018",
    "KEGG_2019_Human"
  )


# Load data
setwd("/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19")
integrated <- qs::qread("integrated.qsave")
Idents(integrated) <- integrated$cell.type
unique(Idents(integrated))
integrated.fibroblast <- subset(integrated, idents = "Fibroblasts")
ncol(integrated.fibroblast)
covid <- qs::qread("covid.qsave")
Idents(covid) <- covid$cell.type
unique(Idents(covid))
covid.fibroblast <- subset(covid, idents = "Fibroblasts")
ncol(covid.fibroblast)


# Set meta data
integrated.fibroblast.id <- rep("copd", ncol(integrated.fibroblast))
names(integrated.fibroblast.id) <- colnames(integrated.fibroblast)
integrated.fibroblast <- AddMetaData(integrated.fibroblast, 
                                     metadata = integrated.fibroblast.id, 
                                     col.name = "new.id")
covid.fibroblast.id <- rep("healthy", ncol(covid.fibroblast))
names(covid.fibroblast.id) <- colnames(covid.fibroblast)
covid.fibroblast <- AddMetaData(covid.fibroblast, metadata = covid.fibroblast.id, 
                                col.name = "new.id")
covid.fibroblast <- RenameCells(covid.fibroblast, add.cell.id = "C")
colnames(covid.fibroblast)
qs::qsave(covid.fibroblast, 
          "/fs/ess/scratch/PCON0022/liyang/Astrazeneca/copd_covid-19/results/covid_fibroblast_renamed.qsave")


# Integration
DefaultAssay(integrated.fibroblast) <- "RNA"
DefaultAssay(covid.fibroblast) <- "RNA"
fibroblast.ll <- c(covid.fibroblast, integrated.fibroblast)
fibroblast.ll <- lapply(X = fibroblast.ll, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = fibroblast.ll)
obj.anchors <- FindIntegrationAnchors(object.list = fibroblast.ll, anchor.features = features)
integrated.again.fibroblast <- IntegrateData(anchorset = obj.anchors)


# DEG analysis
DefaultAssay(integrated.again.fibroblast) <- "RNA"
Idents(integrated.again.fibroblast) <- integrated.again.fibroblast$new.id
unique(Idents(integrated.again.fibroblast))
copd.healthy.DEGs <- FindMarkers(integrated.again.fibroblast, ident.1 = "copd", 
                                 ident.2 = "healthy", verbose = FALSE)
head(copd.healthy.DEGs)
dim(copd.healthy.DEGs)
# selected.DEGs <- copd.healthy.DEGs[copd.healthy.DEGs$p_val_adj < 0.01,]
write.csv(copd.healthy.DEGs, 
          "/fs/ess/scratch/PCON0022/liyang/Astrazeneca/copd_covid-19/results/COPD_healthy_DEGs_fibroblast.csv")
# iFibroblast.DEGs <- read.csv("/fs/ess/scratch/PCON0022/liyang/Astrazeneca/copd_covid-19/results/COPD_healthy_DEGs_iFibroblast.csv")
# head(iFibroblast.DEGs)
# overlapped.genes <- intersect(rownames(copd.healthy.DEGs), unique(iFibroblast.DEGs$gene))
# length(overlapped.genes)
# DEG.pval <- phyper(length(overlapped.genes) - 1, length(unique(iFibroblast.DEGs$gene)), 
#                    nrow(integrated[['RNA']]) - length(unique(iFibroblast.DEGs$gene)), 
#                    nrow(copd.healthy.DEGs), lower.tail = F)


# Pathway enrichment analyses
this_enriched <- enrichr(rownames(copd.healthy.DEGs), dbs)
dim(this_enriched$GO_Molecular_Function_2018[, c(-5, -6)])
write.csv(
  this_enriched$GO_Molecular_Function_2018[, c(-5, -6)],
  paste0("/fs/ess/scratch/PCON0022/liyang/Astrazeneca/copd_covid-19/results/integrated_fibroblast_GO_MF.csv")
)
dim(this_enriched$GO_Cellular_Component_2018[, c(-5, -6)])
write.csv(
  this_enriched$GO_Cellular_Component_2018[, c(-5, -6)],
  paste0("/fs/ess/scratch/PCON0022/liyang/Astrazeneca/copd_covid-19/results/integrated_fibroblast_GO_CC.csv")
)
dim(this_enriched$GO_Biological_Process_2018[, c(-5, -6)])
write.csv(
  this_enriched$GO_Biological_Process_2018[, c(-5, -6)],
  paste0("/fs/ess/scratch/PCON0022/liyang/Astrazeneca/copd_covid-19/results/integrated_fibroblast_GO_BP.csv")
)
dim(this_enriched$KEGG_2019_Human[, c(-5, -6)])
write.csv(
  this_enriched$KEGG_2019_Human[, c(-5, -6)],
  paste0("/fs/ess/scratch/PCON0022/liyang/Astrazeneca/copd_covid-19/results/integrated_fibroblast_KEGG.csv")
)


# Compare enriched GO terms and KEGG pathways
# fibroblast.kegg <- this_enriched$KEGG_2019_Human[, c(-5, -6)]
# cluster0.ll <- c("cluster0_BP", "cluster0_CC", "cluster0_MF", "cluster0_KEGG")
# cluster3.ll <- c("cluster3_BP", "cluster3_CC", "cluster3_MF", "cluster3_KEGG")
# cluster4.ll <- c("cluster4_BP", "cluster4_CC", "cluster4_MF", "cluster4_KEGG")


# GO BP: 30519
# GO CC: 4473
# GO MF: 12413
# Human KEGG: 308


# Cluster 0 BP
id <- 0


# GO BP: 30519
iFibroblast.pathways <-read_excel("/fs/ess/scratch/PCON0022/liyang/Astrazeneca/copd_covid-19/COPD_healthy_enrichment_iFibroblast.xlsx", sheet = id + 1)
dim(iFibroblast.pathways)
head(iFibroblast.pathways[, 1:5])
overlap.pathways <- intersect(unique(this_enriched$GO_Biological_Process_2018[, c(-5, -6)]$Term), 
                              unique(iFibroblast.pathways$Term))
length(overlap.pathways)
pathway.pval <- phyper(length(overlap.pathways) - 1, length(unique(iFibroblast.pathways$Term)), 
                       30519 - length(unique(iFibroblast.pathways$Term)), 
                       length(unique(this_enriched$GO_Biological_Process_2018[, c(-5, -6)]$Term)), 
                       lower.tail = F)
pathway.pval # 1.686129e-153


# GO CC: 4473
iFibroblast.pathways <-read_excel("/fs/ess/scratch/PCON0022/liyang/Astrazeneca/copd_covid-19/COPD_healthy_enrichment_iFibroblast.xlsx", sheet = id + 2)
dim(iFibroblast.pathways)
head(iFibroblast.pathways[, 1:5])
overlap.pathways <- intersect(unique(this_enriched$GO_Cellular_Component_2018[, c(-5, -6)]$Term), 
                              unique(iFibroblast.pathways$Term))
length(overlap.pathways)
pathway.pval <- phyper(length(overlap.pathways) - 1, length(unique(iFibroblast.pathways$Term)), 
                       4473 - length(unique(iFibroblast.pathways$Term)), 
                       length(unique(this_enriched$GO_Cellular_Component_2018[, c(-5, -6)]$Term)), 
                       lower.tail = F)
pathway.pval # 1.615417e-23


# GO MF: 12413
iFibroblast.pathways <-read_excel("/fs/ess/scratch/PCON0022/liyang/Astrazeneca/copd_covid-19/COPD_healthy_enrichment_iFibroblast.xlsx", sheet = id + 3)
dim(iFibroblast.pathways)
head(iFibroblast.pathways[, 1:5])
overlap.pathways <- intersect(unique(this_enriched$GO_Molecular_Function_2018[, c(-5, -6)]$Term), 
                              unique(iFibroblast.pathways$Term))
length(overlap.pathways)
pathway.pval <- phyper(length(overlap.pathways) - 1, length(unique(iFibroblast.pathways$Term)), 
                       12413 - length(unique(iFibroblast.pathways$Term)), 
                       length(unique(this_enriched$GO_Molecular_Function_2018[, c(-5, -6)]$Term)), 
                       lower.tail = F)
pathway.pval # 1.082409e-45


# Human KEGG: 308
iFibroblast.pathways <-read_excel("/fs/ess/scratch/PCON0022/liyang/Astrazeneca/copd_covid-19/COPD_healthy_enrichment_iFibroblast.xlsx", sheet = id + 4)
dim(iFibroblast.pathways)
head(iFibroblast.pathways[, 1:5])
overlap.pathways <- intersect(unique(this_enriched$KEGG_2019_Human[, c(-5, -6)]$Term), 
                              unique(iFibroblast.pathways$Term))
length(overlap.pathways)
pathway.pval <- phyper(length(overlap.pathways) - 1, length(unique(iFibroblast.pathways$Term)), 
                       308 - length(unique(iFibroblast.pathways$Term)), 
                       length(unique(this_enriched$KEGG_2019_Human[, c(-5, -6)]$Term)), 
                       lower.tail = F)
pathway.pval # 0.5889138


# Cluster 3 BP
id <- 4


# GO BP: 30519
iFibroblast.pathways <-read_excel("/fs/ess/scratch/PCON0022/liyang/Astrazeneca/copd_covid-19/COPD_healthy_enrichment_iFibroblast.xlsx", sheet = id + 1)
dim(iFibroblast.pathways)
head(iFibroblast.pathways[, 1:5])
overlap.pathways <- intersect(unique(this_enriched$GO_Biological_Process_2018[, c(-5, -6)]$Term), 
                              unique(iFibroblast.pathways$Term))
length(overlap.pathways)
pathway.pval <- phyper(length(overlap.pathways) - 1, length(unique(iFibroblast.pathways$Term)), 
                       30519 - length(unique(iFibroblast.pathways$Term)), 
                       length(unique(this_enriched$GO_Biological_Process_2018[, c(-5, -6)]$Term)), 
                       lower.tail = F)
pathway.pval # 7.485337e-161


# GO CC: 4473
iFibroblast.pathways <-read_excel("/fs/ess/scratch/PCON0022/liyang/Astrazeneca/copd_covid-19/COPD_healthy_enrichment_iFibroblast.xlsx", sheet = id + 2)
dim(iFibroblast.pathways)
head(iFibroblast.pathways[, 1:5])
overlap.pathways <- intersect(unique(this_enriched$GO_Cellular_Component_2018[, c(-5, -6)]$Term), 
                              unique(iFibroblast.pathways$Term))
length(overlap.pathways)
pathway.pval <- phyper(length(overlap.pathways) - 1, length(unique(iFibroblast.pathways$Term)), 
                       4473 - length(unique(iFibroblast.pathways$Term)), 
                       length(unique(this_enriched$GO_Cellular_Component_2018[, c(-5, -6)]$Term)), 
                       lower.tail = F)
pathway.pval # 9.608475e-27


# GO MF: 12413
iFibroblast.pathways <-read_excel("/fs/ess/scratch/PCON0022/liyang/Astrazeneca/copd_covid-19/COPD_healthy_enrichment_iFibroblast.xlsx", sheet = id + 3)
dim(iFibroblast.pathways)
head(iFibroblast.pathways[, 1:5])
overlap.pathways <- intersect(unique(this_enriched$GO_Molecular_Function_2018[, c(-5, -6)]$Term), 
                              unique(iFibroblast.pathways$Term))
length(overlap.pathways)
pathway.pval <- phyper(length(overlap.pathways) - 1, length(unique(iFibroblast.pathways$Term)), 
                       12413 - length(unique(iFibroblast.pathways$Term)), 
                       length(unique(this_enriched$GO_Molecular_Function_2018[, c(-5, -6)]$Term)), 
                       lower.tail = F)
pathway.pval # 2.862066e-44


# Human KEGG: 308
iFibroblast.pathways <-read_excel("/fs/ess/scratch/PCON0022/liyang/Astrazeneca/copd_covid-19/COPD_healthy_enrichment_iFibroblast.xlsx", sheet = id + 4)
dim(iFibroblast.pathways)
head(iFibroblast.pathways[, 1:5])
overlap.pathways <- intersect(unique(this_enriched$KEGG_2019_Human[, c(-5, -6)]$Term), 
                              unique(iFibroblast.pathways$Term))
length(overlap.pathways)
pathway.pval <- phyper(length(overlap.pathways) - 1, length(unique(iFibroblast.pathways$Term)), 
                       308 - length(unique(iFibroblast.pathways$Term)), 
                       length(unique(this_enriched$KEGG_2019_Human[, c(-5, -6)]$Term)), 
                       lower.tail = F)
pathway.pval # 0.5152996


# Cluster 4 BP
id <- 8


# GO BP: 30519
iFibroblast.pathways <-read_excel("/fs/ess/scratch/PCON0022/liyang/Astrazeneca/copd_covid-19/COPD_healthy_enrichment_iFibroblast.xlsx", sheet = id + 1)
dim(iFibroblast.pathways)
head(iFibroblast.pathways[, 1:5])
overlap.pathways <- intersect(unique(this_enriched$GO_Biological_Process_2018[, c(-5, -6)]$Term), 
                              unique(iFibroblast.pathways$Term))
length(overlap.pathways)
pathway.pval <- phyper(length(overlap.pathways) - 1, length(unique(iFibroblast.pathways$Term)), 
                       30519 - length(unique(iFibroblast.pathways$Term)), 
                       length(unique(this_enriched$GO_Biological_Process_2018[, c(-5, -6)]$Term)), 
                       lower.tail = F)
pathway.pval # 1.885395e-167


# GO CC: 4473
iFibroblast.pathways <-read_excel("/fs/ess/scratch/PCON0022/liyang/Astrazeneca/copd_covid-19/COPD_healthy_enrichment_iFibroblast.xlsx", sheet = id + 2)
dim(iFibroblast.pathways)
head(iFibroblast.pathways[, 1:5])
overlap.pathways <- intersect(unique(this_enriched$GO_Cellular_Component_2018[, c(-5, -6)]$Term), 
                              unique(iFibroblast.pathways$Term))
length(overlap.pathways)
pathway.pval <- phyper(length(overlap.pathways) - 1, length(unique(iFibroblast.pathways$Term)), 
                       4473 - length(unique(iFibroblast.pathways$Term)), 
                       length(unique(this_enriched$GO_Cellular_Component_2018[, c(-5, -6)]$Term)), 
                       lower.tail = F)
pathway.pval # 8.920108e-27


# GO MF: 12413
iFibroblast.pathways <-read_excel("/fs/ess/scratch/PCON0022/liyang/Astrazeneca/copd_covid-19/COPD_healthy_enrichment_iFibroblast.xlsx", sheet = id + 3)
dim(iFibroblast.pathways)
head(iFibroblast.pathways[, 1:5])
overlap.pathways <- intersect(unique(this_enriched$GO_Molecular_Function_2018[, c(-5, -6)]$Term), 
                              unique(iFibroblast.pathways$Term))
length(overlap.pathways)
pathway.pval <- phyper(length(overlap.pathways) - 1, length(unique(iFibroblast.pathways$Term)), 
                       12413 - length(unique(iFibroblast.pathways$Term)), 
                       length(unique(this_enriched$GO_Molecular_Function_2018[, c(-5, -6)]$Term)), 
                       lower.tail = F)
pathway.pval # 3.547779e-45


# Human KEGG: 308
iFibroblast.pathways <-read_excel("/fs/ess/scratch/PCON0022/liyang/Astrazeneca/copd_covid-19/COPD_healthy_enrichment_iFibroblast.xlsx", sheet = id + 4)
dim(iFibroblast.pathways)
head(iFibroblast.pathways[, 1:5])
overlap.pathways <- intersect(unique(this_enriched$KEGG_2019_Human[, c(-5, -6)]$Term), 
                              unique(iFibroblast.pathways$Term))
length(overlap.pathways)
pathway.pval <- phyper(length(overlap.pathways) - 1, length(unique(iFibroblast.pathways$Term)), 
                       308 - length(unique(iFibroblast.pathways$Term)), 
                       length(unique(this_enriched$KEGG_2019_Human[, c(-5, -6)]$Term)), 
                       lower.tail = F)
pathway.pval # 0.2289579
