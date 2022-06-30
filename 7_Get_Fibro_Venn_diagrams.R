#######################################################################
#                                                                     #
#          Generate 1 + 5 Venn diagrams for Fibroblast (Fig. S2J)     #
#                                                                     #
#######################################################################


# Parameters
work.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/"
scratch.dir <- "/fs/ess/scratch/PCON0022/liyang/Astrazeneca/copd_covid-19/results/"
out.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Rfiles/"
image.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Images/Fibroblast/"
table.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Tables/"
setwd(work.dir)
padj.cutoff <- 0.05
avgFC.cutoff <- 0.25
ngenes <- 59113


# Load the iFibroblast GO BP terms and select significant ones
source("/fs/ess/PCON0022/liyang/r_utilities/functions/visual_tools.R")
library(readxl) # read xlsx files
id <- 1
iFibroblast.BP1 <-read_excel(paste0(table.dir, "iFibroblast_terms_pathways.xlsx"), 
                           sheet = id + 1)
dim(iFibroblast.BP1)
head(iFibroblast.BP1)
range(iFibroblast.BP1$`Adjusted P-value`)
iFibroblast.sigBP1 <- iFibroblast.BP1[iFibroblast.BP1$`Adjusted P-value` < padj.cutoff,]
dim(iFibroblast.sigBP1)
range(iFibroblast.sigBP1$`Adjusted P-value`)


# Fibroblast GO BP terms
Fibroblast.BP <- read.csv("integrated_Fibroblasts_GO_BP.csv")
dim(Fibroblast.BP)
head(Fibroblast.BP[, 1:5])
range(Fibroblast.BP$Adjusted.P.value)
Fibroblast.sigBP <- Fibroblast.BP[Fibroblast.BP$Adjusted.P.value < padj.cutoff,]
dim(Fibroblast.sigBP)
range(Fibroblast.sigBP$Adjusted.P.value)


# Calculate hypergeometric p-value for Fibroblast BP Cluster 1
iFibroblast.sigBP1 <- iFibroblast.BP1
Fibroblast.sigBP <- Fibroblast.BP
Fibroblast.BP1Overlap <- length(intersect(unique(Fibroblast.sigBP$Term), unique(iFibroblast.sigBP1$Term)))
Fibroblast.BP1Pval <- phyper(Fibroblast.BP1Overlap - 1, length(unique(iFibroblast.sigBP1$Term)), 
                      ngenes - length(unique(Fibroblast.sigBP$Term)), 
                      length(unique(Fibroblast.sigBP$Term)), lower.tail = F)
Fibroblast.BP1Pval


# Generate the Venn diagram for GO BP terms between Fibroblast v.s. iFibroblast Cluster 1
get_Venn(set.list = list(unique(iFibroblast.sigBP1$Term), unique(Fibroblast.sigBP$Term)), 
         category.names = c("iPSC BP (C0)", "Ref BP"), 
         filename = paste0(image.dir, "Venn_GO_BP_C0.png"))


# Load the iFibroblast GO CC terms and select significant ones
library(readxl) # read xlsx files
id <- 2
iFibroblast.CC1 <-read_excel(paste0(table.dir, "iFibroblast_terms_pathways.xlsx"), 
                      sheet = id + 1)
dim(iFibroblast.CC1)
head(iFibroblast.CC1)
range(iFibroblast.CC1$`Adjusted P-value`)
iFibroblast.sigCC1 <- iFibroblast.CC1[iFibroblast.CC1$`Adjusted P-value` < padj.cutoff,]
dim(iFibroblast.sigCC1)
range(iFibroblast.sigCC1$`Adjusted P-value`)


# Fibroblast GO CC terms
Fibroblast.CC <- read.csv("integrated_Fibroblasts_GO_CC.csv")
dim(Fibroblast.CC)
head(Fibroblast.CC[, 1:5])
range(Fibroblast.CC$Adjusted.P.value)
Fibroblast.sigCC <- Fibroblast.CC[Fibroblast.CC$Adjusted.P.value < padj.cutoff,]
dim(Fibroblast.sigCC)
range(Fibroblast.sigCC$Adjusted.P.value)


# Calculate hypergeometric p-value for Fibroblast CC Cluster 1
iFibroblast.sigCC1 <- iFibroblast.CC1
Fibroblast.sigCC <- Fibroblast.CC
Fibroblast.CC1Overlap <- length(intersect(unique(Fibroblast.sigCC$Term), unique(iFibroblast.sigCC1$Term)))
Fibroblast.CC1Pval <- phyper(Fibroblast.CC1Overlap - 1, length(unique(iFibroblast.sigCC1$Term)), 
                      ngenes - length(unique(Fibroblast.sigCC$Term)), 
                      length(unique(Fibroblast.sigCC$Term)), lower.tail = F)
Fibroblast.CC1Pval


# Generate the Venn diagram for GO CC terms between Fibroblast v.s. iFibroblast Cluster 1
get_Venn(set.list = list(unique(iFibroblast.sigCC1$Term), unique(Fibroblast.sigCC$Term)), 
         category.names = c("iPSC CC (C0)", "Ref CC"), 
         filename = paste0(image.dir, "Venn_GO_CC_C0.png"))


# Load the iFibroblast GO MF terms and select significant ones
library(readxl) # read xlsx files
id <- 3
iFibroblast.MF1 <-read_excel(paste0(table.dir, "iFibroblast_terms_pathways.xlsx"), 
                      sheet = id + 1)
dim(iFibroblast.MF1)
head(iFibroblast.MF1)
range(iFibroblast.MF1$`Adjusted P-value`)
iFibroblast.sigMF1 <- iFibroblast.MF1[iFibroblast.MF1$`Adjusted P-value` < padj.cutoff,]
dim(iFibroblast.sigMF1)
range(iFibroblast.sigMF1$`Adjusted P-value`)


# Fibroblast GO MF terms
Fibroblast.MF <- read.csv("integrated_Fibroblasts_GO_MF.csv")
dim(Fibroblast.MF)
head(Fibroblast.MF[, 1:5])
range(Fibroblast.MF$Adjusted.P.value)
Fibroblast.sigMF <- Fibroblast.MF[Fibroblast.MF$Adjusted.P.value < padj.cutoff,]
dim(Fibroblast.sigMF)
range(Fibroblast.sigMF$Adjusted.P.value)


# Calculate hypergeometric p-value for Fibroblast MF Cluster 0
iFibroblast.sigMF1 <- iFibroblast.MF1
Fibroblast.sigMF <- Fibroblast.MF
Fibroblast.MF1Overlap <- length(intersect(unique(Fibroblast.sigMF$Term), unique(iFibroblast.sigMF1$Term)))
Fibroblast.MF1Pval <- phyper(Fibroblast.MF1Overlap - 1, length(unique(iFibroblast.sigMF1$Term)), 
                      ngenes - length(unique(Fibroblast.sigMF$Term)), 
                      length(unique(Fibroblast.sigMF$Term)), lower.tail = F)
Fibroblast.MF1Pval


# Generate the Venn diagram for GO MF terms between Fibroblast v.s. iFibroblast Cluster 1
get_Venn(set.list = list(unique(iFibroblast.sigMF1$Term), unique(Fibroblast.sigMF$Term)), 
         category.names = c("iPSC MF (C0)", "Ref MF"), 
         filename = paste0(image.dir, "Venn_GO_MF_C0.png"))


# Load the iFibroblast GO KEGG terms and select significant ones
library(readxl) # read xlsx files
id <- 4
iFibroblast.KEGG1 <-read_excel(paste0(table.dir, "iFibroblast_terms_pathways.xlsx"), 
                      sheet = id + 1)
dim(iFibroblast.KEGG1)
head(iFibroblast.KEGG1)
range(iFibroblast.KEGG1$`Adjusted P-value`)
iFibroblast.sigKEGG1 <- iFibroblast.KEGG1[iFibroblast.KEGG1$`Adjusted P-value` < padj.cutoff,]
dim(iFibroblast.sigKEGG1)
range(iFibroblast.sigKEGG1$`Adjusted P-value`)


# Fibroblast GO KEGG terms
Fibroblast.KEGG <- read.csv("integrated_Fibroblasts_KEGG.csv")
dim(Fibroblast.KEGG)
head(Fibroblast.KEGG[, 1:5])
range(Fibroblast.KEGG$Adjusted.P.value)
Fibroblast.sigKEGG <- Fibroblast.KEGG[Fibroblast.KEGG$Adjusted.P.value < padj.cutoff,]
dim(Fibroblast.sigKEGG)
range(Fibroblast.sigKEGG$Adjusted.P.value)


# Calculate hypergeometric p-value for Fibroblast KEGG Cluster 1
iFibroblast.sigKEGG1 <- iFibroblast.KEGG1
Fibroblast.sigKEGG <- Fibroblast.KEGG
Fibroblast.KEGG1Overlap <- length(intersect(unique(Fibroblast.sigKEGG$Term), unique(iFibroblast.sigKEGG1$Term)))
Fibroblast.KEGG1Pval <- phyper(Fibroblast.KEGG1Overlap - 1, length(unique(iFibroblast.sigKEGG1$Term)), 
                      ngenes - length(unique(Fibroblast.sigKEGG$Term)), 
                      length(unique(Fibroblast.sigKEGG$Term)), lower.tail = F)
Fibroblast.KEGG1Pval


# Generate the Venn diagram for GO KEGG terms between Fibroblast v.s. iFibroblast Cluster 1
get_Venn(set.list = list(unique(iFibroblast.sigKEGG1$Term), unique(Fibroblast.sigKEGG$Term)), 
         category.names = c("iPSC KEGG (C0)", "Ref KEGG"), 
         filename = paste0(image.dir, "Venn_KEGG_C0.png"))


# Load the iFibroblast GO BP terms and select significant ones (Cluster 3)
library(readxl) # read xlsx files
id <- 5
iFibroblast.BP1 <-read_excel(paste0(table.dir, "iFibroblast_terms_pathways.xlsx"), 
                      sheet = id + 1)
dim(iFibroblast.BP1)
head(iFibroblast.BP1)
range(iFibroblast.BP1$`Adjusted P-value`)
iFibroblast.sigBP1 <- iFibroblast.BP1[iFibroblast.BP1$`Adjusted P-value` < padj.cutoff,]
dim(iFibroblast.sigBP1)
range(iFibroblast.sigBP1$`Adjusted P-value`)


# Fibroblast GO BP terms
Fibroblast.BP <- read.csv("integrated_Fibroblasts_GO_BP.csv")
dim(Fibroblast.BP)
head(Fibroblast.BP[, 1:5])
range(Fibroblast.BP$Adjusted.P.value)
Fibroblast.sigBP <- Fibroblast.BP[Fibroblast.BP$Adjusted.P.value < padj.cutoff,]
dim(Fibroblast.sigBP)
range(Fibroblast.sigBP$Adjusted.P.value)


# Calculate hypergeometric p-value for Fibroblast BP Cluster 1
iFibroblast.sigBP1 <- iFibroblast.BP1
Fibroblast.sigBP <- Fibroblast.BP
Fibroblast.BP1Overlap <- length(intersect(unique(Fibroblast.sigBP$Term), unique(iFibroblast.sigBP1$Term)))
Fibroblast.BP1Pval <- phyper(Fibroblast.BP1Overlap - 1, length(unique(iFibroblast.sigBP1$Term)), 
                      ngenes - length(unique(Fibroblast.sigBP$Term)), 
                      length(unique(Fibroblast.sigBP$Term)), lower.tail = F)
Fibroblast.BP1Pval


# Generate the Venn diagram for GO BP terms between Fibroblast v.s. iFibroblast Cluster 1
get_Venn(set.list = list(unique(iFibroblast.sigBP1$Term), unique(Fibroblast.sigBP$Term)), 
         category.names = c("iPSC BP (C3)", "Ref BP"), 
         filename = paste0(image.dir, "Venn_GO_BP_C3.png"))


# Load the iFibroblast GO CC terms and select significant ones (Cluster 2)
library(readxl) # read xlsx files
id <- 6
iFibroblast.CC1 <-read_excel(paste0(table.dir, "iFibroblast_terms_pathways.xlsx"), 
                      sheet = id + 1)
dim(iFibroblast.CC1)
head(iFibroblast.CC1)
range(iFibroblast.CC1$`Adjusted P-value`)
iFibroblast.sigCC1 <- iFibroblast.CC1[iFibroblast.CC1$`Adjusted P-value` < padj.cutoff,]
dim(iFibroblast.sigCC1)
range(iFibroblast.sigCC1$`Adjusted P-value`)


# Fibroblast GO CC terms
Fibroblast.CC <- read.csv("integrated_Fibroblasts_GO_CC.csv")
dim(Fibroblast.CC)
head(Fibroblast.CC[, 1:5])
range(Fibroblast.CC$Adjusted.P.value)
Fibroblast.sigCC <- Fibroblast.CC[Fibroblast.CC$Adjusted.P.value < padj.cutoff,]
dim(Fibroblast.sigCC)
range(Fibroblast.sigCC$Adjusted.P.value)


# Calculate hypergeometric p-value for Fibroblast CC Cluster 1
iFibroblast.sigCC1 <- iFibroblast.CC1
Fibroblast.sigCC <- Fibroblast.CC
Fibroblast.CC1Overlap <- length(intersect(unique(Fibroblast.sigCC$Term), unique(iFibroblast.sigCC1$Term)))
Fibroblast.CC1Pval <- phyper(Fibroblast.CC1Overlap - 1, length(unique(iFibroblast.sigCC1$Term)), 
                      ngenes - length(unique(Fibroblast.sigCC$Term)), 
                      length(unique(Fibroblast.sigCC$Term)), lower.tail = F)
Fibroblast.CC1Pval


# Generate the Venn diagram for GO CC terms between Fibroblast v.s. iFibroblast Cluster 1
get_Venn(set.list = list(unique(iFibroblast.sigCC1$Term), unique(Fibroblast.sigCC$Term)), 
         category.names = c("iPSC CC (C3)", "Ref CC"), 
         filename = paste0(image.dir, "Venn_GO_CC_C3.png"))


# Load the iFibroblast GO MF terms and select significant ones (Cluster 2)
library(readxl) # read xlsx files
id <- 7
iFibroblast.MF1 <-read_excel(paste0(table.dir, "iFibroblast_terms_pathways.xlsx"), 
                      sheet = id + 1)
dim(iFibroblast.MF1)
head(iFibroblast.MF1)
range(iFibroblast.MF1$`Adjusted P-value`)
iFibroblast.sigMF1 <- iFibroblast.MF1[iFibroblast.MF1$`Adjusted P-value` < padj.cutoff,]
dim(iFibroblast.sigMF1)
range(iFibroblast.sigMF1$`Adjusted P-value`)


# Fibroblast GO MF terms
Fibroblast.MF <- read.csv("integrated_Fibroblasts_GO_MF.csv")
dim(Fibroblast.MF)
head(Fibroblast.MF[, 1:5])
range(Fibroblast.MF$Adjusted.P.value)
Fibroblast.sigMF <- Fibroblast.MF[Fibroblast.MF$Adjusted.P.value < padj.cutoff,]
dim(Fibroblast.sigMF)
range(Fibroblast.sigMF$Adjusted.P.value)


# Calculate hypergeometric p-value for Fibroblast MF Cluster 1
iFibroblast.sigMF1 <- iFibroblast.MF1
Fibroblast.sigMF <- Fibroblast.MF
Fibroblast.MF1Overlap <- length(intersect(unique(Fibroblast.sigMF$Term), unique(iFibroblast.sigMF1$Term)))
Fibroblast.MF1Pval <- phyper(Fibroblast.MF1Overlap - 1, length(unique(iFibroblast.sigMF1$Term)), 
                      ngenes - length(unique(Fibroblast.sigMF$Term)), 
                      length(unique(Fibroblast.sigMF$Term)), lower.tail = F)
Fibroblast.MF1Pval


# Generate the Venn diagram for GO MF terms between Fibroblast v.s. iFibroblast Cluster 1
get_Venn(set.list = list(unique(iFibroblast.sigMF1$Term), unique(Fibroblast.sigMF$Term)), 
         category.names = c("iPSC MF (C3)", "Ref MF"), 
         filename = paste0(image.dir, "Venn_GO_MF_C3.png"))


# Load the iFibroblast KEGG terms and select significant ones (Cluster 2)
library(readxl) # read xlsx files
id <- 8
iFibroblast.MF1 <-read_excel(paste0(table.dir, "iFibroblast_terms_pathways.xlsx"), 
                      sheet = id + 1)
dim(iFibroblast.MF1)
head(iFibroblast.MF1)
range(iFibroblast.MF1$`Adjusted P-value`)
iFibroblast.sigMF1 <- iFibroblast.MF1[iFibroblast.MF1$`Adjusted P-value` < padj.cutoff,]
dim(iFibroblast.sigMF1)
range(iFibroblast.sigMF1$`Adjusted P-value`)


# Fibroblast KEGG terms
Fibroblast.MF <- read.csv("integrated_Fibroblasts_KEGG.csv")
dim(Fibroblast.MF)
head(Fibroblast.MF[, 1:5])
range(Fibroblast.MF$Adjusted.P.value)
Fibroblast.sigMF <- Fibroblast.MF[Fibroblast.MF$Adjusted.P.value < padj.cutoff,]
dim(Fibroblast.sigMF)
range(Fibroblast.sigMF$Adjusted.P.value)


# Calculate hypergeometric p-value for Fibroblast MF Cluster 1
iFibroblast.sigMF1 <- iFibroblast.MF1
Fibroblast.sigMF <- Fibroblast.MF
Fibroblast.MF1Overlap <- length(intersect(unique(Fibroblast.sigMF$Term), unique(iFibroblast.sigMF1$Term)))
Fibroblast.MF1Pval <- phyper(Fibroblast.MF1Overlap - 1, length(unique(iFibroblast.sigMF1$Term)), 
                      ngenes - length(unique(Fibroblast.sigMF$Term)), 
                      length(unique(Fibroblast.sigMF$Term)), lower.tail = F)
Fibroblast.MF1Pval


# Generate the Venn diagram for KEGG terms between Fibroblast v.s. iFibroblast Cluster 1
get_Venn(set.list = list(unique(iFibroblast.sigMF1$Term), unique(Fibroblast.sigMF$Term)), 
         category.names = c("iPSC KEGG (C3)", "Ref KEGG"), 
         filename = paste0(image.dir, "Venn_KEGG_C3.png"))


# Load the iFibroblast GO_BP terms and select significant ones (Cluster 2)
library(readxl) # read xlsx files
id <- 9
iFibroblast.MF1 <-read_excel(paste0(table.dir, "iFibroblast_terms_pathways.xlsx"), 
                             sheet = id + 1)
dim(iFibroblast.MF1)
head(iFibroblast.MF1)
range(iFibroblast.MF1$`Adjusted P-value`)
iFibroblast.sigMF1 <- iFibroblast.MF1[iFibroblast.MF1$`Adjusted P-value` < padj.cutoff,]
dim(iFibroblast.sigMF1)
range(iFibroblast.sigMF1$`Adjusted P-value`)


# Fibroblast GO_BP terms
Fibroblast.MF <- read.csv("integrated_Fibroblasts_GO_BP.csv")
dim(Fibroblast.MF)
head(Fibroblast.MF[, 1:5])
range(Fibroblast.MF$Adjusted.P.value)
Fibroblast.sigMF <- Fibroblast.MF[Fibroblast.MF$Adjusted.P.value < padj.cutoff,]
dim(Fibroblast.sigMF)
range(Fibroblast.sigMF$Adjusted.P.value)


# Calculate hypergeometric p-value for Fibroblast MF Cluster 1
iFibroblast.sigMF1 <- iFibroblast.MF1
Fibroblast.sigMF <- Fibroblast.MF
Fibroblast.MF1Overlap <- length(intersect(unique(Fibroblast.sigMF$Term), unique(iFibroblast.sigMF1$Term)))
Fibroblast.MF1Pval <- phyper(Fibroblast.MF1Overlap - 1, length(unique(iFibroblast.sigMF1$Term)), 
                             ngenes - length(unique(Fibroblast.sigMF$Term)), 
                             length(unique(Fibroblast.sigMF$Term)), lower.tail = F)
Fibroblast.MF1Pval


# Generate the Venn diagram for GO_BP terms between Fibroblast v.s. iFibroblast Cluster 1
get_Venn(set.list = list(unique(iFibroblast.sigMF1$Term), unique(Fibroblast.sigMF$Term)), 
         category.names = c("iPSC BP (C4)", "Ref BP"), 
         filename = paste0(image.dir, "Venn_GO_BP_C4.png"))


# Load the iFibroblast GO_CC terms and select significant ones (Cluster 2)
library(readxl) # read xlsx files
id <- 10
iFibroblast.MF1 <-read_excel(paste0(table.dir, "iFibroblast_terms_pathways.xlsx"), 
                             sheet = id + 1)
dim(iFibroblast.MF1)
head(iFibroblast.MF1)
range(iFibroblast.MF1$`Adjusted P-value`)
iFibroblast.sigMF1 <- iFibroblast.MF1[iFibroblast.MF1$`Adjusted P-value` < padj.cutoff,]
dim(iFibroblast.sigMF1)
range(iFibroblast.sigMF1$`Adjusted P-value`)


# Fibroblast GO_CC terms
Fibroblast.MF <- read.csv("integrated_Fibroblasts_GO_CC.csv")
dim(Fibroblast.MF)
head(Fibroblast.MF[, 1:5])
range(Fibroblast.MF$Adjusted.P.value)
Fibroblast.sigMF <- Fibroblast.MF[Fibroblast.MF$Adjusted.P.value < padj.cutoff,]
dim(Fibroblast.sigMF)
range(Fibroblast.sigMF$Adjusted.P.value)


# Calculate hypergeometric p-value for Fibroblast MF Cluster 1
iFibroblast.sigMF1 <- iFibroblast.MF1
Fibroblast.sigMF <- Fibroblast.MF
Fibroblast.MF1Overlap <- length(intersect(unique(Fibroblast.sigMF$Term), unique(iFibroblast.sigMF1$Term)))
Fibroblast.MF1Pval <- phyper(Fibroblast.MF1Overlap - 1, length(unique(iFibroblast.sigMF1$Term)), 
                             ngenes - length(unique(Fibroblast.sigMF$Term)), 
                             length(unique(Fibroblast.sigMF$Term)), lower.tail = F)
Fibroblast.MF1Pval


# Generate the Venn diagram for GO_CC terms between Fibroblast v.s. iFibroblast Cluster 1
get_Venn(set.list = list(unique(iFibroblast.sigMF1$Term), unique(Fibroblast.sigMF$Term)), 
         category.names = c("iPSC CC (C4)", "Ref CC"), 
         filename = paste0(image.dir, "Venn_GO_CC_C4.png"))


# Load the iFibroblast GO_MF terms and select significant ones (Cluster 2)
library(readxl) # read xlsx files
id <- 11
iFibroblast.MF1 <-read_excel(paste0(table.dir, "iFibroblast_terms_pathways.xlsx"), 
                             sheet = id + 1)
dim(iFibroblast.MF1)
head(iFibroblast.MF1)
range(iFibroblast.MF1$`Adjusted P-value`)
iFibroblast.sigMF1 <- iFibroblast.MF1[iFibroblast.MF1$`Adjusted P-value` < padj.cutoff,]
dim(iFibroblast.sigMF1)
range(iFibroblast.sigMF1$`Adjusted P-value`)


# Fibroblast GO_MF terms
Fibroblast.MF <- read.csv("integrated_Fibroblasts_GO_MF.csv")
dim(Fibroblast.MF)
head(Fibroblast.MF[, 1:5])
range(Fibroblast.MF$Adjusted.P.value)
Fibroblast.sigMF <- Fibroblast.MF[Fibroblast.MF$Adjusted.P.value < padj.cutoff,]
dim(Fibroblast.sigMF)
range(Fibroblast.sigMF$Adjusted.P.value)


# Calculate hypergeometric p-value for Fibroblast MF Cluster 1
iFibroblast.sigMF1 <- iFibroblast.MF1
Fibroblast.sigMF <- Fibroblast.MF
Fibroblast.MF1Overlap <- length(intersect(unique(Fibroblast.sigMF$Term), unique(iFibroblast.sigMF1$Term)))
Fibroblast.MF1Pval <- phyper(Fibroblast.MF1Overlap - 1, length(unique(iFibroblast.sigMF1$Term)), 
                             ngenes - length(unique(Fibroblast.sigMF$Term)), 
                             length(unique(Fibroblast.sigMF$Term)), lower.tail = F)
Fibroblast.MF1Pval


# Generate the Venn diagram for GO_MF terms between Fibroblast v.s. iFibroblast Cluster 1
get_Venn(set.list = list(unique(iFibroblast.sigMF1$Term), unique(Fibroblast.sigMF$Term)), 
         category.names = c("iPSC MF (C4)", "Ref MF"), 
         filename = paste0(image.dir, "Venn_GO_MF_C4.png"))


# Load the iFibroblast KEGG terms and select significant ones (Cluster 2)
library(readxl) # read xlsx files
id <- 12
iFibroblast.MF1 <-read_excel(paste0(table.dir, "iFibroblast_terms_pathways.xlsx"), 
                             sheet = id + 1)
dim(iFibroblast.MF1)
head(iFibroblast.MF1)
range(iFibroblast.MF1$`Adjusted P-value`)
iFibroblast.sigMF1 <- iFibroblast.MF1[iFibroblast.MF1$`Adjusted P-value` < padj.cutoff,]
dim(iFibroblast.sigMF1)
range(iFibroblast.sigMF1$`Adjusted P-value`)


# Fibroblast KEGG terms
Fibroblast.MF <- read.csv("integrated_Fibroblasts_KEGG.csv")
dim(Fibroblast.MF)
head(Fibroblast.MF[, 1:5])
range(Fibroblast.MF$Adjusted.P.value)
Fibroblast.sigMF <- Fibroblast.MF[Fibroblast.MF$Adjusted.P.value < padj.cutoff,]
dim(Fibroblast.sigMF)
range(Fibroblast.sigMF$Adjusted.P.value)


# Calculate hypergeometric p-value for Fibroblast MF Cluster 1
iFibroblast.sigMF1 <- iFibroblast.MF1
Fibroblast.sigMF <- Fibroblast.MF
Fibroblast.MF1Overlap <- length(intersect(unique(Fibroblast.sigMF$Term), unique(iFibroblast.sigMF1$Term)))
Fibroblast.MF1Pval <- phyper(Fibroblast.MF1Overlap - 1, length(unique(iFibroblast.sigMF1$Term)), 
                             ngenes - length(unique(Fibroblast.sigMF$Term)), 
                             length(unique(Fibroblast.sigMF$Term)), lower.tail = F)
Fibroblast.MF1Pval


# Generate the Venn diagram for KEGG terms between Fibroblast v.s. iFibroblast Cluster 1
get_Venn(set.list = list(unique(iFibroblast.sigMF1$Term), unique(Fibroblast.sigMF$Term)), 
         category.names = c("iPSC KEGG (C4)", "Ref KEGG"), 
         filename = paste0(image.dir, "Venn_KEGG_C4.png"))
