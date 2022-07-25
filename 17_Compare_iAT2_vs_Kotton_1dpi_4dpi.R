###################################################################
#                                                                 #
# Compare top-ranked DEGs between 1dpi v.s. mock (our data)       #
# and 4dpi v.s. mock (Kotton paper)                               #
#                                                                 #
###################################################################


# Check line 233



# Libraries
library(Seurat)
library(dplyr)
library(readxl)
library(ggvenn)
library(fgsea)
library(qusage)
library(org.Hs.eg.db)
library(msigdbr)


# Global variables
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
R.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Rfiles/"
image.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Images/"
table.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Tables/"
iAT2.file <- "Supplemental_7._Figure_5E_DEG_results.xlsx"
pval.cutoff <- 0.05
logFC.cutoff <- 0.25


# Source codes
source(paste0(tool.dir, "transcriptome_tools.R"))
source(paste0(tool.dir, "visual_tools.R"))
source(paste0(tool.dir, "combinatorial_tools.R"))


# Identify DEGs between 1dpi v.s. Mock in iAT2 of our paper
DEGs.iAT2.1dpi <- read_excel(paste0(table.dir, 
                                    "iAT2_DEGs/Supplemental_7._Figure_5E_DEG_results.xlsx"), 
                             sheet = 2)
dim(DEGs.iAT2.1dpi)
head(DEGs.iAT2.1dpi)
length(unique(DEGs.iAT2.1dpi$hgnc_symbols))
range(DEGs.iAT2.1dpi$PValue)
range(DEGs.iAT2.1dpi$FC)
dim(DEGs.iAT2.1dpi[DEGs.iAT2.1dpi$FC > logFC.cutoff,]) # 18802    
dim(DEGs.iAT2.1dpi[DEGs.iAT2.1dpi$PValue < pval.cutoff,]) # 111  
dim(DEGs.iAT2.1dpi[!is.na(DEGs.iAT2.1dpi$hgnc_symbols),]) # 16450    
DEGs.iAT2.1dpi <- DEGs.iAT2.1dpi[DEGs.iAT2.1dpi$FC > logFC.cutoff & 
                                   DEGs.iAT2.1dpi$PValue < pval.cutoff & 
                                   !is.na(DEGs.iAT2.1dpi$hgnc_symbols),]
dim(DEGs.iAT2.1dpi) # Only 35 DEGs
head(DEGs.iAT2.1dpi)
DEGs.iAT2.1dpi <- DEGs.iAT2.1dpi[order(DEGs.iAT2.1dpi$PValue),] # only 35 DEGs


# Identify DEGs between 1dpi v.s. Mock in AT2 of Kotton paper
DEGs.AT2 <- read_excel(paste0(table.dir, "[AT2_DEGs]_Jessie_Huang_et_al_2020.xlsx"), 
                       sheet = 1)
dim(DEGs.AT2) # 17612    
range(DEGs.AT2$LogFC.dpi1_vs_ctr)
range(DEGs.AT2$FDR.dpi1_vs_ctr)
dim(DEGs.AT2[DEGs.AT2$LogFC.dpi1_vs_ctr > 
               logFC.cutoff,]) # 3965   
dim(DEGs.AT2[DEGs.AT2$FDR.dpi1_vs_ctr < pval.cutoff,]) # 4519   
DEGs.AT2.1dpi <- DEGs.AT2[DEGs.AT2$LogFC.dpi1_vs_ctr > 
                            logFC.cutoff & DEGs.AT2$FDR.dpi1_vs_ctr < pval.cutoff,]
dim(DEGs.AT2.1dpi) # 2393   
colnames(DEGs.AT2.1dpi)
head(DEGs.AT2.1dpi)
DEGs.AT2.1dpi$Genes.SYMBOL
DEGs.AT2.1dpi <- DEGs.AT2.1dpi[order(DEGs.AT2.1dpi$FDR.dpi1_vs_ctr),]
DEGs.AT2.1dpi <- DEGs.AT2.1dpi[!duplicated(DEGs.AT2.1dpi$Genes.SYMBOL),]
AT2.1dpi.overlaps <- get_overlap_ratio(x = DEGs.iAT2.1dpi$hgnc_symbols, y = DEGs.AT2.1dpi$Genes.SYMBOL)
# No overlapped DEG exists!!!


# Try 4dpi v.s. 1dpi
dim(DEGs.AT2) # 17612    
range(DEGs.AT2$LogFC.dpi4_vs_dpi1)
range(DEGs.AT2$FDR.dpi4_vs_dpi1)
dim(DEGs.AT2[DEGs.AT2$LogFC.dpi4_vs_dpi1 > 
               logFC.cutoff,]) # 6021      
dim(DEGs.AT2[DEGs.AT2$FDR.dpi4_vs_dpi1 < pval.cutoff,]) # 10725       
DEGs.AT2.4dpi <- DEGs.AT2[DEGs.AT2$LogFC.dpi4_vs_dpi1 > 
                            logFC.cutoff & DEGs.AT2$FDR.dpi4_vs_dpi1 < pval.cutoff,]
dim(DEGs.AT2.4dpi) # 4849      
colnames(DEGs.AT2.4dpi)
head(DEGs.AT2.4dpi)
DEGs.AT2.4dpi$Genes.SYMBOL
length(DEGs.AT2.4dpi$Genes.SYMBOL)
length(unique(DEGs.AT2.4dpi$Genes.SYMBOL))
DEGs.AT2.4dpi <- DEGs.AT2.4dpi[order(DEGs.AT2.4dpi$FDR.dpi4_vs_dpi1),]
DEGs.AT2.4dpi <- DEGs.AT2.4dpi[!duplicated(DEGs.AT2.4dpi$Genes.SYMBOL),]
dim(DEGs.AT2.4dpi) # 4848
AT2.4dpi.overlaps <- get_overlap_ratio(x = DEGs.iAT2.1dpi$hgnc_symbols, 
                                       y = DEGs.AT2.4dpi$Genes.SYMBOL)


# Try all infection v.s. ctrl
dim(DEGs.AT2) # 17612    
range(DEGs.AT2$LogFC.allinfected_vs_ctr)
range(DEGs.AT2$FDR.allinfected_vs_ctr)
dim(DEGs.AT2[DEGs.AT2$LogFC.allinfected_vs_ctr > 
               logFC.cutoff,]) # 5292         
dim(DEGs.AT2[DEGs.AT2$FDR.allinfected_vs_ctr < pval.cutoff,]) # 10462           
DEGs.AT2.infected <- DEGs.AT2[DEGs.AT2$LogFC.allinfected_vs_ctr > 
                            logFC.cutoff & DEGs.AT2$FDR.allinfected_vs_ctr < pval.cutoff,]
dim(DEGs.AT2.infected) # 4288         
colnames(DEGs.AT2.infected)
head(DEGs.AT2.infected)
DEGs.AT2.infected$Genes.SYMBOL
length(DEGs.AT2.infected$Genes.SYMBOL)
length(unique(DEGs.AT2.infected$Genes.SYMBOL))
DEGs.AT2.infected <- DEGs.AT2.infected[order(DEGs.AT2.infected$FDR.allinfected_vs_ctr),]
DEGs.AT2.infected <- DEGs.AT2.infected[!duplicated(DEGs.AT2.infected$Genes.SYMBOL),]
dim(DEGs.AT2.infected) # 4288
AT2.infected.overlaps <- get_overlap_ratio(x = DEGs.iAT2.1dpi$hgnc_symbols, 
                                       y = DEGs.AT2.infected$Genes.SYMBOL)
# No overlapped DEG exists!
# The possible reason may be the small number of DEGs in our data!!!


# Increase the number of DEGs in iAT2 by loosing the p-value cutoff
DEGs.iAT2.1dpi <- read_excel(paste0(table.dir, 
                                    "iAT2_DEGs/Supplemental_7._Figure_5E_DEG_results.xlsx"), 
                             sheet = 2)
dim(DEGs.iAT2.1dpi) # 19148    
DEGs.iAT2.1dpi <- DEGs.iAT2.1dpi[!is.na(DEGs.iAT2.1dpi$hgnc_symbols),]
dim(DEGs.iAT2.1dpi) # 16450    
range(DEGs.iAT2.1dpi$PValue)
dim(DEGs.iAT2.1dpi[DEGs.iAT2.1dpi$PValue <= pval.cutoff,]) # 81
dim(DEGs.iAT2.1dpi[DEGs.iAT2.1dpi$PValue < pval.cutoff,]) # 81
DEGs.iAT2.1dpi <- DEGs.iAT2.1dpi[DEGs.iAT2.1dpi$PValue < pval.cutoff,]
dim(DEGs.iAT2.1dpi) # 81
length(intersect(DEGs.iAT2.1dpi$hgnc_symbols, DEGs.AT2.1dpi$Genes.SYMBOL)) # only 6!!!
length(intersect(DEGs.iAT2.1dpi$hgnc_symbols, DEGs.AT2.4dpi$Genes.SYMBOL)) # only 16!!!
length(intersect(DEGs.iAT2.1dpi$hgnc_symbols, DEGs.AT2.infected$Genes.SYMBOL)) # only 17!!!


# Check the distribution of negative logarithm of  p-values
DEGs.iAT2.1dpi <- read_excel(paste0(table.dir, 
                                    "iAT2_DEGs/Supplemental_7._Figure_5E_DEG_results.xlsx"), 
                             sheet = 2)
dim(DEGs.iAT2.1dpi) # 19148    
boxplot(-log(DEGs.iAT2.1dpi$PValue))


DEGs.AT2 <- read_excel(paste0(table.dir, "[AT2_DEGs]_Jessie_Huang_et_al_2020.xlsx"), 
                       sheet = 1)
dim(DEGs.AT2) # 17612    
head(DEGs.AT2$Genes.SYMBOL)
overlap.1dpi <- intersect(unique(DEGs.iAT2.1dpi$hgnc_symbols), unique(DEGs.AT2$Genes.SYMBOL))
length(overlap.1dpi) # 13705
overlapped.iAT2.1dpi <- DEGs.iAT2.1dpi[DEGs.iAT2.1dpi$hgnc_symbols %in% overlap.1dpi,]
dim(overlapped.iAT2.1dpi) # 13710    
overlapped.AT2.1dpi <- DEGs.AT2[DEGs.AT2$Genes.SYMBOL %in% overlap.1dpi,]
dim(overlapped.AT2.1dpi) # 13713    


pval.AT2.1dpi <- overlapped.AT2.1dpi[overlapped.AT2.1dpi$FDR.dpi1_vs_ctr < pval.cutoff,]


pval.iAT2.1dpi <- overlapped.iAT2.1dpi[overlapped.iAT2.1dpi$PValue < pval.cutoff,]
dim(pval.AT2.1dpi) # 3863
dim(pval.iAT2.1dpi) # 40
length(intersect(unique(pval.AT2.1dpi$Genes.SYMBOL), 
                 unique(pval.iAT2.1dpi$hgnc_symbols))) # 7


length(intersect(unique(pval.AT2.1dpi$Genes.SYMBOL), 
                 unique(DEGs.iAT2.2dpi[DEGs.iAT2.2dpi$PValue < pval.cutoff,]$hgnc_symbols))) # 28
dim(DEGs.iAT2.3dpi[DEGs.iAT2.3dpi$PValue < pval.cutoff,]) # 151
length(unique(DEGs.iAT2.3dpi[DEGs.iAT2.3dpi$PValue < pval.cutoff,]$hgnc_symbols)) # 105


length(intersect(unique(pval.AT2.1dpi$Genes.SYMBOL), 
                 unique(DEGs.iAT2.3dpi[DEGs.iAT2.3dpi$PValue < pval.cutoff,]$hgnc_symbols))) # 58
dim(DEGs.iAT2.3dpi[DEGs.iAT2.3dpi$PValue < pval.cutoff,]) # 217
length(unique(DEGs.iAT2.3dpi[DEGs.iAT2.3dpi$PValue < pval.cutoff,]$hgnc_symbols)) # 170; 34.1%


sig.AT2.1dpi <- overlapped.AT2.1dpi[overlapped.AT2.1dpi$LogFC.dpi1_vs_ctr > logFC.cutoff & 
                                      overlapped.AT2.1dpi$FDR.dpi1_vs_ctr < pval.cutoff, ]
dim(sig.AT2.1dpi) # 1837   
sorted.AT2.1dpi <- sig.AT2.1dpi[order(sig.AT2.1dpi$FDR.dpi1_vs_ctr),]
sorted.AT2.1dpi <- sorted.AT2.1dpi[!duplicated(sorted.AT2.1dpi$Genes.SYMBOL),]
dim(sorted.AT2.1dpi) # 1837   


length(intersect(sorted.AT2.1dpi$Genes.SYMBOL, overlapped.iAT2.1dpi$hgnc_symbols)) # 1837
high.iAT2.1dpi <- overlapped.iAT2.1dpi[overlapped.iAT2.1dpi$FC > logFC.cutoff,]
dim(high.iAT2.1dpi) # 13564    
sig.iAT2.1dpi <- high.iAT2.1dpi[high.iAT2.1dpi$PValue < pval.cutoff,]
dim(sig.iAT2.1dpi) # only 18!!!


sig.iAT2.1dpi <- sig.iAT2.1dpi[order(sig.iAT2.1dpi$PValue),]
sig.iAT2.1dpi <- sig.iAT2.1dpi[!duplicated(sig.iAT2.1dpi$hgnc_symbols),]
dim(sig.iAT2.1dpi)
get_overlap_ratio(x = sorted.AT2.1dpi$Genes.SYMBOL, y = sig.iAT2.1dpi$hgnc_symbols)
# No overlapped DEG!!!
get_overlap_ratio(x = sorted.AT2.1dpi$Genes.SYMBOL, 
                  y = overlapped.iAT2.1dpi[overlapped.iAT2.1dpi$PValue < pval.cutoff,]$hgnc_symbols)
# No overlapped DEG!!!


# Use DEGs from 2dpi in our data
DEGs.iAT2.2dpi <- read_excel(paste0(table.dir, 
                                    "iAT2_DEGs/Supplemental_7._Figure_5E_DEG_results.xlsx"), 
                             sheet = 3)
dim(DEGs.iAT2.2dpi) # 19148        
boxplot(-log(DEGs.iAT2.2dpi$PValue))
dim(DEGs.iAT2.2dpi[DEGs.iAT2.2dpi$PValue < pval.cutoff,]) #151
dim(DEGs.iAT2.2dpi[DEGs.iAT2.2dpi$FC > logFC.cutoff,]) # 18659
length(intersect(DEGs.iAT2.2dpi[DEGs.iAT2.2dpi$PValue < pval.cutoff,]$hgnc_symbols, 
                 sorted.AT2.1dpi$Genes.SYMBOL)) # 28 / 151 = 18.5%


# Use DEGs from 3dpi in our data
DEGs.iAT2.3dpi <- read_excel(paste0(table.dir, 
                                    "iAT2_DEGs/Supplemental_7._Figure_5E_DEG_results.xlsx"), 
                             sheet = 4)
dim(DEGs.iAT2.3dpi) # 19148
boxplot(-log(DEGs.iAT2.3dpi$PValue))
dim(DEGs.iAT2.3dpi[DEGs.iAT2.3dpi$PValue < pval.cutoff,]) # 217
dim(DEGs.iAT2.3dpi[DEGs.iAT2.3dpi$FC > logFC.cutoff,]) # 18733
length(intersect(DEGs.iAT2.3dpi[DEGs.iAT2.3dpi$PValue < pval.cutoff,]$hgnc_symbols, 
                 sorted.AT2.1dpi$Genes.SYMBOL)) # 58 / 151 = 26.7%



# Comparison of DEGs between 3dpi (our data) v.s. 4dpi (Kotton paper)
DEGs.iAT2.3dpi <- read_excel(paste0(table.dir, 
                                    "iAT2_DEGs/Supplemental_7._Figure_5E_DEG_results.xlsx"), 
                             sheet = 4)
dim(DEGs.iAT2.3dpi) # 19148    
dim(DEGs.AT2) # 17612   
omitted.AT2.4dpi <- DEGs.AT2[!is.na(DEGs.AT2$Genes.SYMBOL),]
dim(omitted.AT2.4dpi) # 17612    
overlap.4dpi <- intersect(unique(omitted.AT2.4dpi$Genes.SYMBOL), 
                          unique(DEGs.iAT2.3dpi$hgnc_symbols))
length(overlap.4dpi) # 13705
overlap.AT2.4dpi <- omitted.AT2.4dpi[omitted.AT2.4dpi$Genes.SYMBOL %in% overlap.4dpi, ]
dim(overlap.AT2.4dpi) # 13713    
overlap.iAT2.3dpi <- DEGs.iAT2.3dpi[DEGs.iAT2.3dpi$hgnc_symbols %in% overlap.4dpi,]
dim(overlap.iAT2.3dpi) # 13710    


sig.AT2.4dpi <- overlap.AT2.4dpi[overlap.AT2.4dpi$FDR.dpi1_vs_ctr < pval.cutoff,]
dim(sig.AT2.4dpi) # 3863   
sig.iAT2.3dpi <- overlap.iAT2.3dpi[overlap.iAT2.3dpi$PValue < pval.cutoff,]
dim(sig.iAT2.3dpi) # 115
length(unique(sig.iAT2.3dpi$hgnc_symbols)) # 115


uniq.AT2.4dpi <- sig.AT2.4dpi[!duplicated(sig.AT2.4dpi$Genes.SYMBOL),]
dim(uniq.AT2.4dpi) # 3863   
uniq.iAT2.3dpi <- sig.iAT2.3dpi[!duplicated(sig.iAT2.3dpi$hgnc_symbols),]
dim(uniq.iAT2.3dpi) # 115
get_overlap_ratio(x = uniq.AT2.4dpi$Genes.SYMBOL, 
                  y = uniq.iAT2.3dpi$hgnc_symbols)
# The maximum overlap ratio, 0.25, is achieved at 4.

length(intersect(x = uniq.AT2.4dpi$Genes.SYMBOL, 
                 y = uniq.iAT2.3dpi$hgnc_symbols)) / 
  length(uniq.iAT2.3dpi$hgnc_symbols) # 50.4%
