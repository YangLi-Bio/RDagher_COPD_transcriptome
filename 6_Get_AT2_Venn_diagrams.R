#######################################################################
#                                                                     #
#          Generate 1 + 5 Venn diagrams for AT2 (Fig. S2J)            #
#                                                                     #
#######################################################################


# Parameters
source("/fs/ess/PCON0022/liyang/r_utilities/functions/visual_tools.R")
work.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/"
scratch.dir <- "/fs/ess/scratch/PCON0022/liyang/Astrazeneca/copd_covid-19/results/"
out.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Rfiles/"
image.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Images/AT2/"
table.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Tables/"
setwd(work.dir)
padj.cutoff <- 0.05
avgFC.cutoff <- 0.25
ngenes <- 59113


# Load iAT2 DEGs and select the significant ones
iAT2.DEGs <- read.csv("COPD_healthy_DEGs_iAT2.csv", header = T)
dim(iAT2.DEGs)
head(iAT2.DEGs)
range(iAT2.DEGs$adj_pval)
range(abs(iAT2.DEGs$logfc))
iAT2.sigDEGs <- iAT2.DEGs[iAT2.DEGs$adj_pval < padj.cutoff & 
                            abs(iAT2.DEGs$logfc) > avgFC.cutoff,]
dim(iAT2.sigDEGs)
range(iAT2.sigDEGs$adj_pval)
range(abs(iAT2.sigDEGs$logfc))
qs::qsave(iAT2.sigDEGs, paste0(out.dir, "DEGs_iAT2_significant.qsave"))


# Load the AT2 DEGs and select significant ones
AT2.DEGs <- read.csv("COPD_healthy_DEGs.csv", header = T)
dim(AT2.DEGs)
head(AT2.DEGs)
colnames(AT2.DEGs)[1] <- "gene"
head(AT2.DEGs)
AT2.sigDEGs <- AT2.DEGs[AT2.DEGs$p_val_adj < padj.cutoff & 
                          abs(AT2.DEGs$avg_log2FC) > avgFC.cutoff,]
dim(AT2.sigDEGs)
range(AT2.sigDEGs$p_val_adj)
range(abs(AT2.sigDEGs$avg_log2FC))
qs::qsave(AT2.sigDEGs, paste0(out.dir, "DEGs_AT2_significant.qsave"))


# Calculate hypergeometric p-value
AT2.sigDEGs <- AT2.DEGs
iAT2.sigDEGs <- iAT2.DEGs
AT2.overlap <- length(intersect(AT2.sigDEGs$gene, iAT2.sigDEGs$gene))
AT2.DEGPval <- phyper(AT2.overlap - 1, length(unique(iAT2.sigDEGs$gene)), 
                   ngenes - length(unique(iAT2.sigDEGs$gene)), 
                   length(unique(AT2.sigDEGs$gene)), lower.tail = F)
AT2.DEGPval


# Generate the Venn diagram for DEGs between AT2 v.s. iAT2
get_Venn(set.list = list(unique(iAT2.sigDEGs$gene), unique(AT2.sigDEGs$gene)), 
         category.names = c("iAT2", "AT2"), 
         filename = paste0(image.dir, "Venn_DEGs_AT2.png"))


# Load the iAT2 GO BP terms and select significant ones
library(readxl) # read xlsx files
id <- 1
iAT2.BP1 <-read_excel(paste0(table.dir, "iAT2_terms_pathways.xlsx"), 
                           sheet = id + 1)
dim(iAT2.BP1)
head(iAT2.BP1)
range(iAT2.BP1$`Adjusted P-value`)
iAT2.sigBP1 <- iAT2.BP1[iAT2.BP1$`Adjusted P-value` < padj.cutoff,]
dim(iAT2.sigBP1)
range(iAT2.sigBP1$`Adjusted P-value`)


# AT2 GO BP terms
AT2.BP <- read.csv("integrated_AT2_GO_BP.csv")
dim(AT2.BP)
head(AT2.BP[, 1:5])
range(AT2.BP$Adjusted.P.value)
AT2.sigBP <- AT2.BP[AT2.BP$Adjusted.P.value < padj.cutoff,]
dim(AT2.sigBP)
range(AT2.sigBP$Adjusted.P.value)


# Calculate hypergeometric p-value for AT2 BP Cluster 1
iAT2.sigBP1 <- iAT2.BP1
AT2.sigBP <- AT2.BP
AT2.BP1Overlap <- length(intersect(unique(AT2.sigBP$Term), unique(iAT2.sigBP1$Term)))
AT2.BP1Pval <- phyper(AT2.BP1Overlap - 1, length(unique(iAT2.sigBP1$Term)), 
                      ngenes - length(unique(AT2.sigBP$Term)), 
                      length(unique(AT2.sigBP$Term)), lower.tail = F)
AT2.BP1Pval


# Generate the Venn diagram for GO BP terms between AT2 v.s. iAT2 Cluster 1
get_Venn(set.list = list(unique(iAT2.sigBP1$Term), unique(AT2.sigBP$Term)), 
         category.names = c("iPSC BP (C1)", "Ref BP"), 
         filename = paste0(image.dir, "Venn_GO_BP_C1.png"))


# Load the iAT2 GO CC terms and select significant ones
library(readxl) # read xlsx files
id <- 2
iAT2.CC1 <-read_excel(paste0(table.dir, "iAT2_terms_pathways.xlsx"), 
                      sheet = id + 1)
dim(iAT2.CC1)
head(iAT2.CC1)
range(iAT2.CC1$`Adjusted P-value`)
iAT2.sigCC1 <- iAT2.CC1[iAT2.CC1$`Adjusted P-value` < padj.cutoff,]
dim(iAT2.sigCC1)
range(iAT2.sigCC1$`Adjusted P-value`)


# AT2 GO CC terms
AT2.CC <- read.csv("integrated_AT2_GO_CC.csv")
dim(AT2.CC)
head(AT2.CC[, 1:5])
range(AT2.CC$Adjusted.P.value)
AT2.sigCC <- AT2.CC[AT2.CC$Adjusted.P.value < padj.cutoff,]
dim(AT2.sigCC)
range(AT2.sigCC$Adjusted.P.value)


# Calculate hypergeometric p-value for AT2 CC Cluster 1
iAT2.sigCC1 <- iAT2.CC1
AT2.sigCC <- AT2.CC
AT2.CC1Overlap <- length(intersect(unique(AT2.sigCC$Term), unique(iAT2.sigCC1$Term)))
AT2.CC1Pval <- phyper(AT2.CC1Overlap - 1, length(unique(iAT2.sigCC1$Term)), 
                      ngenes - length(unique(AT2.sigCC$Term)), 
                      length(unique(AT2.sigCC$Term)), lower.tail = F)
AT2.CC1Pval


# Generate the Venn diagram for GO CC terms between AT2 v.s. iAT2 Cluster 1
get_Venn(set.list = list(unique(iAT2.sigCC1$Term), unique(AT2.sigCC$Term)), 
         category.names = c("iPSC CC (C1)", "Ref CC"), 
         filename = paste0(image.dir, "Venn_GO_CC_C1.png"))


# Load the iAT2 GO MF terms and select significant ones
library(readxl) # read xlsx files
id <- 3
iAT2.MF1 <-read_excel(paste0(table.dir, "iAT2_terms_pathways.xlsx"), 
                      sheet = id + 1)
dim(iAT2.MF1)
head(iAT2.MF1)
range(iAT2.MF1$`Adjusted P-value`)
iAT2.sigMF1 <- iAT2.MF1[iAT2.MF1$`Adjusted P-value` < padj.cutoff,]
dim(iAT2.sigMF1)
range(iAT2.sigMF1$`Adjusted P-value`)


# AT2 GO MF terms
AT2.MF <- read.csv("integrated_AT2_GO_MF.csv")
dim(AT2.MF)
head(AT2.MF[, 1:5])
range(AT2.MF$Adjusted.P.value)
AT2.sigMF <- AT2.MF[AT2.MF$Adjusted.P.value < padj.cutoff,]
dim(AT2.sigMF)
range(AT2.sigMF$Adjusted.P.value)


# Calculate hypergeometric p-value for AT2 MF Cluster 1
iAT2.sigMF1 <- iAT2.MF1
AT2.sigMF <- AT2.MF
AT2.MF1Overlap <- length(intersect(unique(AT2.sigMF$Term), unique(iAT2.sigMF1$Term)))
AT2.MF1Pval <- phyper(AT2.MF1Overlap - 1, length(unique(iAT2.sigMF1$Term)), 
                      ngenes - length(unique(AT2.sigMF$Term)), 
                      length(unique(AT2.sigMF$Term)), lower.tail = F)
AT2.MF1Pval


# Generate the Venn diagram for GO MF terms between AT2 v.s. iAT2 Cluster 1
get_Venn(set.list = list(unique(iAT2.sigMF1$Term), unique(AT2.sigMF$Term)), 
         category.names = c("iPSC MF (C1)", "Ref MF"), 
         filename = paste0(image.dir, "Venn_GO_MF_C1.png"))


# Load the iAT2 GO KEGG terms and select significant ones
library(readxl) # read xlsx files
id <- 4
iAT2.KEGG1 <-read_excel(paste0(table.dir, "iAT2_terms_pathways.xlsx"), 
                      sheet = id + 1)
dim(iAT2.KEGG1)
head(iAT2.KEGG1)
range(iAT2.KEGG1$`Adjusted P-value`)
iAT2.sigKEGG1 <- iAT2.KEGG1[iAT2.KEGG1$`Adjusted P-value` < padj.cutoff,]
dim(iAT2.sigKEGG1)
range(iAT2.sigKEGG1$`Adjusted P-value`)


# AT2 GO KEGG terms
AT2.KEGG <- read.csv("integrated_AT2_KEGG.csv")
dim(AT2.KEGG)
head(AT2.KEGG[, 1:5])
range(AT2.KEGG$Adjusted.P.value)
AT2.sigKEGG <- AT2.KEGG[AT2.KEGG$Adjusted.P.value < padj.cutoff,]
dim(AT2.sigKEGG)
range(AT2.sigKEGG$Adjusted.P.value)


# Calculate hypergeometric p-value for AT2 KEGG Cluster 1
iAT2.sigKEGG1 <- iAT2.KEGG1
AT2.sigKEGG <- AT2.KEGG
AT2.KEGG1Overlap <- length(intersect(unique(AT2.sigKEGG$Term), unique(iAT2.sigKEGG1$Term)))
AT2.KEGG1Pval <- phyper(AT2.KEGG1Overlap - 1, length(unique(iAT2.sigKEGG1$Term)), 
                      ngenes - length(unique(AT2.sigKEGG$Term)), 
                      length(unique(AT2.sigKEGG$Term)), lower.tail = F)
AT2.KEGG1Pval


# Generate the Venn diagram for GO KEGG terms between AT2 v.s. iAT2 Cluster 1
get_Venn(set.list = list(unique(iAT2.sigKEGG1$Term), unique(AT2.sigKEGG$Term)), 
         category.names = c("iPSC KEGG (C1)", "Ref KEGG"), 
         filename = paste0(image.dir, "Venn_KEGG_C1.png"))


# Load the iAT2 GO BP terms and select significant ones (Cluster 2)
library(readxl) # read xlsx files
id <- 5
iAT2.BP1 <-read_excel(paste0(table.dir, "iAT2_terms_pathways.xlsx"), 
                      sheet = id + 1)
dim(iAT2.BP1)
head(iAT2.BP1)
range(iAT2.BP1$`Adjusted P-value`)
iAT2.sigBP1 <- iAT2.BP1[iAT2.BP1$`Adjusted P-value` < padj.cutoff,]
dim(iAT2.sigBP1)
range(iAT2.sigBP1$`Adjusted P-value`)


# AT2 GO BP terms
AT2.BP <- read.csv("integrated_AT2_GO_BP.csv")
dim(AT2.BP)
head(AT2.BP[, 1:5])
range(AT2.BP$Adjusted.P.value)
AT2.sigBP <- AT2.BP[AT2.BP$Adjusted.P.value < padj.cutoff,]
dim(AT2.sigBP)
range(AT2.sigBP$Adjusted.P.value)


# Calculate hypergeometric p-value for AT2 BP Cluster 1
iAT2.sigBP1 <- iAT2.BP1
AT2.sigBP <- AT2.BP
AT2.BP1Overlap <- length(intersect(unique(AT2.sigBP$Term), unique(iAT2.sigBP1$Term)))
AT2.BP1Pval <- phyper(AT2.BP1Overlap - 1, length(unique(iAT2.sigBP1$Term)), 
                      ngenes - length(unique(AT2.sigBP$Term)), 
                      length(unique(AT2.sigBP$Term)), lower.tail = F)
AT2.BP1Pval


# Generate the Venn diagram for GO BP terms between AT2 v.s. iAT2 Cluster 1
get_Venn(set.list = list(unique(iAT2.sigBP1$Term), unique(AT2.sigBP$Term)), 
         category.names = c("iPSC BP (C2)", "Ref BP"), 
         filename = paste0(image.dir, "Venn_GO_BP_C2.png"))


# Load the iAT2 GO CC terms and select significant ones (Cluster 2)
library(readxl) # read xlsx files
id <- 6
iAT2.CC1 <-read_excel(paste0(table.dir, "iAT2_terms_pathways.xlsx"), 
                      sheet = id + 1)
dim(iAT2.CC1)
head(iAT2.CC1)
range(iAT2.CC1$`Adjusted P-value`)
iAT2.sigCC1 <- iAT2.CC1[iAT2.CC1$`Adjusted P-value` < padj.cutoff,]
dim(iAT2.sigCC1)
range(iAT2.sigCC1$`Adjusted P-value`)


# AT2 GO CC terms
AT2.CC <- read.csv("integrated_AT2_GO_CC.csv")
dim(AT2.CC)
head(AT2.CC[, 1:5])
range(AT2.CC$Adjusted.P.value)
AT2.sigCC <- AT2.CC[AT2.CC$Adjusted.P.value < padj.cutoff,]
dim(AT2.sigCC)
range(AT2.sigCC$Adjusted.P.value)


# Calculate hypergeometric p-value for AT2 CC Cluster 1
iAT2.sigCC1 <- iAT2.CC1
AT2.sigCC <- AT2.CC
AT2.CC1Overlap <- length(intersect(unique(AT2.sigCC$Term), unique(iAT2.sigCC1$Term)))
AT2.CC1Pval <- phyper(AT2.CC1Overlap - 1, length(unique(iAT2.sigCC1$Term)), 
                      ngenes - length(unique(AT2.sigCC$Term)), 
                      length(unique(AT2.sigCC$Term)), lower.tail = F)
AT2.CC1Pval


# Generate the Venn diagram for GO CC terms between AT2 v.s. iAT2 Cluster 1
get_Venn(set.list = list(unique(iAT2.sigCC1$Term), unique(AT2.sigCC$Term)), 
         category.names = c("iPSC CC (C2)", "Ref CC"), 
         filename = paste0(image.dir, "Venn_GO_CC_C2.png"))


# Load the iAT2 GO MF terms and select significant ones (Cluster 2)
library(readxl) # read xlsx files
id <- 7
iAT2.MF1 <-read_excel(paste0(table.dir, "iAT2_terms_pathways.xlsx"), 
                      sheet = id + 1)
dim(iAT2.MF1)
head(iAT2.MF1)
range(iAT2.MF1$`Adjusted P-value`)
iAT2.sigMF1 <- iAT2.MF1[iAT2.MF1$`Adjusted P-value` < padj.cutoff,]
dim(iAT2.sigMF1)
range(iAT2.sigMF1$`Adjusted P-value`)


# AT2 GO MF terms
AT2.MF <- read.csv("integrated_AT2_GO_MF.csv")
dim(AT2.MF)
head(AT2.MF[, 1:5])
range(AT2.MF$Adjusted.P.value)
AT2.sigMF <- AT2.MF[AT2.MF$Adjusted.P.value < padj.cutoff,]
dim(AT2.sigMF)
range(AT2.sigMF$Adjusted.P.value)


# Calculate hypergeometric p-value for AT2 MF Cluster 1
iAT2.sigMF1 <- iAT2.MF1
AT2.sigMF <- AT2.MF
AT2.MF1Overlap <- length(intersect(unique(AT2.sigMF$Term), unique(iAT2.sigMF1$Term)))
AT2.MF1Pval <- phyper(AT2.MF1Overlap - 1, length(unique(iAT2.sigMF1$Term)), 
                      ngenes - length(unique(AT2.sigMF$Term)), 
                      length(unique(AT2.sigMF$Term)), lower.tail = F)
AT2.MF1Pval


# Generate the Venn diagram for GO MF terms between AT2 v.s. iAT2 Cluster 1
get_Venn(set.list = list(unique(iAT2.sigMF1$Term), unique(AT2.sigMF$Term)), 
         category.names = c("iPSC MF (C2)", "Ref MF"), 
         filename = paste0(image.dir, "Venn_GO_MF_C2.png"))


# Load the iAT2 KEGG terms and select significant ones (Cluster 2)
library(readxl) # read xlsx files
id <- 8
iAT2.MF1 <-read_excel(paste0(table.dir, "iAT2_terms_pathways.xlsx"), 
                      sheet = id + 1)
dim(iAT2.MF1)
head(iAT2.MF1)
range(iAT2.MF1$`Adjusted P-value`)
iAT2.sigMF1 <- iAT2.MF1[iAT2.MF1$`Adjusted P-value` < padj.cutoff,]
dim(iAT2.sigMF1)
range(iAT2.sigMF1$`Adjusted P-value`)


# AT2 KEGG terms
AT2.MF <- read.csv("integrated_AT2_KEGG.csv")
dim(AT2.MF)
head(AT2.MF[, 1:5])
range(AT2.MF$Adjusted.P.value)
AT2.sigMF <- AT2.MF[AT2.MF$Adjusted.P.value < padj.cutoff,]
dim(AT2.sigMF)
range(AT2.sigMF$Adjusted.P.value)


# Calculate hypergeometric p-value for AT2 MF Cluster 1
iAT2.sigMF1 <- iAT2.MF1
AT2.sigMF <- AT2.MF
AT2.MF1Overlap <- length(intersect(unique(AT2.sigMF$Term), unique(iAT2.sigMF1$Term)))
AT2.MF1Pval <- phyper(AT2.MF1Overlap - 1, length(unique(iAT2.sigMF1$Term)), 
                      ngenes - length(unique(AT2.sigMF$Term)), 
                      length(unique(AT2.sigMF$Term)), lower.tail = F)
AT2.MF1Pval


# Generate the Venn diagram for KEGG terms between AT2 v.s. iAT2 Cluster 1
get_Venn(set.list = list(unique(iAT2.sigMF1$Term), unique(AT2.sigMF$Term)), 
         category.names = c("iPSC KEGG (C2)", "Ref KEGG"), 
         filename = paste0(image.dir, "Venn_KEGG_C2.png"))
