###################################################################
#                                                                 #
# Compare top-ranked DEGs between 1dpi v.s. mock (our data)       #
# and 4dpi v.s. mock (Kotton paper)                               #
#                                                                 #
###################################################################


# Libraries
library(Seurat)
library(dplyr)
library(readxl)
library(ggvenn)
library(fgsea)
library(qusage)
library(org.Hs.eg.db)
library(msigdbr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)


# Global variables
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
R.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Rfiles/"
image.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Images/"
table.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Tables/"
iAT2.file <- "Supplemental_7._Figure_5E_DEG_results.xlsx"
orig.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Orig_data/"
kotton.file <- "GSE153277_cpm.csv"
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
dim(DEGs.iAT2.2dpi[DEGs.iAT2.2dpi$PValue < pval.cutoff & 
                     DEGs.iAT2.2dpi$FC > logFC.cutoff,]) # 80
length(intersect(unique(DEGs.iAT2.2dpi[DEGs.iAT2.2dpi$PValue < pval.cutoff & 
                                  DEGs.iAT2.2dpi$FC > logFC.cutoff,]$hgnc_symbols), 
                 unique(omitted.AT2.4dpi[omitted.AT2.4dpi$LogFC.dpi4_vs_dpi1 > logFC.cutoff & 
                                           omitted.AT2.4dpi$FDR.dpi4_vs_dpi1 < pval.cutoff,]$Genes.SYMBOL))) / 
  length(unique(DEGs.iAT2.2dpi[DEGs.iAT2.2dpi$PValue < pval.cutoff & 
                          DEGs.iAT2.2dpi$FC > logFC.cutoff,]$hgnc_symbols)) # 47.5% ()


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
DEGs.AT2 <- read_excel(paste0(table.dir, "[AT2_DEGs]_Jessie_Huang_et_al_2020.xlsx"), 
                       sheet = 1)
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
length(intersect(x = uniq.AT2.4dpi$Genes.SYMBOL, 
                 y = uniq.iAT2.3dpi$hgnc_symbols)) # 58


uniq.AT2.4dpi <- uniq.AT2.4dpi[order(uniq.AT2.4dpi$FDR.dpi4_vs_dpi1),]
dpi4_1.ratio <- get_overlap_ratio(x = uniq.iAT2.3dpi$hgnc_symbols, y = uniq.AT2.4dpi$Genes.SYMBOL)
# The maximum overlap ratio, 0.0571428571428571, is achieved at 105.
plot(dpi4_1.ratio, xlab = "Top-ranked DEGs", ylab = "Ratio of overlapped DEGs")


# The analyses above focus on the comparison between 
# iAT2 DEGs (3dpi v.s. mock) and AT2 DEGs (4dpi v.s. 1dpi)
# We might try comparing iAT2 DEGs with AT2 DEGs from 4dpi v.s. mock
# identified by Yang


kotton.df <- read.csv(paste0(orig.dir, kotton.file))
dim(kotton.df)
head(kotton.df)
rownames(kotton.df)
colnames(kotton.df)[1] <- "Symbol"
head(kotton.df)
colnames(kotton.df)


# Build the metadata
meta.df <- t(data.frame(Study_design = colnames(kotton.df)[-1], 
                      Treatment = c(rep("ctr", 3), 
                                    rep("dpi1", 3), 
                                    rep("dpi4", 3))))
dim(meta.df)
head(meta.df)
write.csv(meta.df, paste0(orig.dir, "metadata.csv"), 
          row.names = T, quote = F)
all.genes <- readLines(paste0(orig.dir, "All_gene_lists_GMT.txt")) %>% 
  strsplit(., split = "\t") %>% `[[` (1)
all.genes <- all.genes[-c(1, 2)]
head(all.genes)
tail(all.genes)
length(which(grepl("^ENSG", all.genes)))
length(all.genes)
ensembl.ids <- grep("^ENSG", all.genes, value = T)




# Build correspondence between Ensembl IDs with symbols
annotation.info <- AnnotationDbi:::select(EnsDb.Hsapiens.v86, keys = ensembl.ids, keytype = "GENEID",
                                          columns = c("GENEID", "GENENAME", 
                                                      "SYMBOL", "GENEBIOTYPE", "ENTREZID",
                                                      "SEQNAME", "SEQSTRAND", "PROTEINID", 
                                                      "UNIPROTID",  "UNIPROTDB"))
head(annotation.info)
colnames(annotation.info)
ensembl.symbol <- annotation.info$SYMBOL
names(ensembl.symbol) <- annotation.info$GENEID
head(ensembl.symbol)
setdiff(ensembl.ids, names(ensembl.symbol))
intersect(ensembl.ids, names(ensembl.symbol))
setdiff(ensembl.ids, names(ensembl.symbol)) %>% length
my.DEGs.AT2.4dpi <- ensembl.symbol[ensembl.ids] %>% na.omit
length(my.DEGs.AT2.4dpi)
dim(uniq.iAT2.3dpi) # 115
length(intersect(my.DEGs.AT2.4dpi, unique(uniq.iAT2.3dpi$hgnc_symbols))) / 
  length(unique(uniq.iAT2.3dpi$hgnc_symbols)) # 47.0%; 53
length(intersect(my.DEGs.AT2.4dpi, unique(sig.iAT2.1dpi$hgnc_symbols))) / 
  length(unique(sig.iAT2.1dpi$hgnc_symbols)) # 5.6%; 1 (No!!!)


# Load all diff genes from AT2 between 4dpi v.s. mock
diff.AT2.4dpi.vs.mock <- read.csv(paste0(orig.dir, "Diff_expression_all_comparisons.csv"))
dim(diff.AT2.4dpi.vs.mock) # 28435    
head(diff.AT2.4dpi.vs.mock)
DEGs.AT2.diff.4dpi.vs.mock <- diff.AT2.4dpi.vs.mock[diff.AT2.4dpi.vs.mock$logFC > logFC.cutoff & 
                                                      diff.AT2.4dpi.vs.mock$adj.P.Val < pval.cutoff,]
DEGs.AT2.diff.4dpi.vs.mock <- DEGs.AT2.diff.4dpi.vs.mock[order(DEGs.AT2.diff.4dpi.vs.mock$adj.P.Val),]
dim(DEGs.AT2.diff.4dpi.vs.mock) # 5250   
qs::qsave(DEGs.AT2.diff.4dpi.vs.mock, paste0(R.dir, "Kotton_iAT2/DEGs_AT2_4dpi_vs_mock_Yang.qsave"))
dim(DEGs.AT2.diff.4dpi.vs.mock) # 5250


uniq.iAT2.3dpi <- uniq.iAT2.3dpi[order(uniq.iAT2.3dpi$PValue),]
mock.ratio <- get_overlap_ratio(x = uniq.iAT2.3dpi$hgnc_symbols, y = DEGs.AT2.diff.4dpi.vs.mock$Symbol)
# The maximum overlap ratio, 0.0571428571428571, is achieved at 105.
plot(mock.ratio, xlab = "Top-ranked DEGs", ylab = "Ratio of overlapped DEGs")


length(intersect(unique(DEGs.AT2.diff.4dpi.vs.mock$Symbol), unique(uniq.iAT2.3dpi$hgnc_symbols))) # 83
length(intersect(unique(DEGs.AT2.diff.4dpi.vs.mock$Symbol), unique(uniq.iAT2.3dpi$hgnc_symbols))) / 
  length(unique(uniq.iAT2.3dpi$hgnc_symbols)) # 72.2%
length(unique(uniq.iAT2.3dpi$hgnc_symbols)) - 
  length(intersect(unique(DEGs.AT2.diff.4dpi.vs.mock$Symbol), unique(uniq.iAT2.3dpi$hgnc_symbols))) # 32
length(unique(DEGs.AT2.diff.4dpi.vs.mock$Symbol)) - 
  length(intersect(unique(DEGs.AT2.diff.4dpi.vs.mock$Symbol), unique(uniq.iAT2.3dpi$hgnc_symbols))) # 4535


uniq.iAT2.1dpi <- sig.iAT2.1dpi[!duplicated(sig.iAT2.1dpi$hgnc_symbols),]
length(intersect(unique(DEGs.AT2.diff.4dpi.vs.mock$Symbol), unique(uniq.iAT2.1dpi$hgnc_symbols))) # 7
length(intersect(unique(DEGs.AT2.diff.4dpi.vs.mock$Symbol), unique(uniq.iAT2.1dpi$hgnc_symbols))) / 
  length(unique(uniq.iAT2.1dpi$hgnc_symbols)) # 38.9%


# Enrichment analyses for DEGs (iAT2 between 3dpi v.s. mock) and (AT2 between 4dpi v.s. mock)
overlap.3dpi.4dpi.mock <- intersect(unique(DEGs.AT2.diff.4dpi.vs.mock$Symbol), 
                                    unique(uniq.iAT2.3dpi$hgnc_symbols))
length(overlap.3dpi.4dpi.mock) # 83
write.csv(cbind(Overlap = uniq.iAT2.3dpi$hgnc_symbols %in% overlap.3dpi.4dpi.mock, 
                uniq.iAT2.3dpi), quote = F, row.names = F, 
          paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_mock/", 
                 "iAT2_DEGs_3dpi_vs_mock.csv"))
write.csv(cbind(Overlap = DEGs.AT2.diff.4dpi.vs.mock$Symbol %in% overlap.3dpi.4dpi.mock, 
                DEGs.AT2.diff.4dpi.vs.mock), quote = F, row.names = F, 
          paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_mock/", 
                 "AT2_DEGs_4dpi_vs_mock.csv"))



iAT2.enriched.3dpi.4dpi.mock <- run_GO_and_KEGG(genes.ll = unique(uniq.iAT2.3dpi$hgnc_symbols), 
                                               org = "human")
AT2.enriched.3dpi.4dpi.mock <- run_GO_and_KEGG(genes.ll = unique(DEGs.AT2.diff.4dpi.vs.mock$Symbol), 
                                                org = "human")


# MF
iAT2.3dpi.4dpi.MF <- iAT2.enriched.3dpi.4dpi.mock$GO_Molecular_Function_2018
AT2.3dpi.4dpi.MF <- AT2.enriched.3dpi.4dpi.mock$GO_Molecular_Function_2018
iAT2.3dpi.4dpi.MF <- iAT2.3dpi.4dpi.MF[iAT2.3dpi.4dpi.MF$Adjusted.P.value < pval.cutoff,]
dim(iAT2.3dpi.4dpi.MF) # 2
AT2.3dpi.4dpi.MF <- AT2.3dpi.4dpi.MF[AT2.3dpi.4dpi.MF$Adjusted.P.value < pval.cutoff,]
dim(AT2.3dpi.4dpi.MF) # 75
length(intersect(unique(iAT2.3dpi.4dpi.MF$Term), unique(AT2.3dpi.4dpi.MF$Term))) # 1
length(intersect(unique(iAT2.3dpi.4dpi.MF$Term), unique(AT2.3dpi.4dpi.MF$Term))) / 
  length(unique(iAT2.3dpi.4dpi.MF$Term)) # 50.0%
length(unique(iAT2.3dpi.4dpi.MF$Term)) - 
  length(intersect(unique(iAT2.3dpi.4dpi.MF$Term), unique(AT2.3dpi.4dpi.MF$Term))) # 1
length(unique(AT2.3dpi.4dpi.MF$Term)) - 
  length(intersect(unique(iAT2.3dpi.4dpi.MF$Term), unique(AT2.3dpi.4dpi.MF$Term))) # 86
overlap.MF <- intersect(unique(iAT2.3dpi.4dpi.MF$Term), unique(AT2.3dpi.4dpi.MF$Term))
overlap.MF
dir.create(paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_mock/"))
write.csv(cbind(Overlap = iAT2.3dpi.4dpi.MF$Term %in% overlap.MF, iAT2.3dpi.4dpi.MF), 
          quote = F, row.names = F, 
          paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_mock/", 
                 "iAT2_molecular_function.csv"))
write.csv(cbind(Overlap = AT2.3dpi.4dpi.MF$Term %in% overlap.MF, AT2.3dpi.4dpi.MF), 
          quote = F, row.names = F, 
          paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_mock/", 
                 "AT2_molecular_function.csv"))


# CC
iAT2.3dpi.4dpi.CC <- iAT2.enriched.3dpi.4dpi.mock$GO_Cellular_Component_2018
AT2.3dpi.4dpi.CC <- AT2.enriched.3dpi.4dpi.mock$GO_Cellular_Component_2018
iAT2.3dpi.4dpi.CC <- iAT2.3dpi.4dpi.CC[iAT2.3dpi.4dpi.CC$Adjusted.P.value < pval.cutoff,]
dim(iAT2.3dpi.4dpi.CC) # 0
AT2.3dpi.4dpi.CC <- AT2.3dpi.4dpi.CC[AT2.3dpi.4dpi.CC$Adjusted.P.value < pval.cutoff,]
dim(AT2.3dpi.4dpi.CC) # 7
length(intersect(unique(iAT2.3dpi.4dpi.CC$Term), unique(AT2.3dpi.4dpi.CC$Term))) # 0
length(intersect(unique(iAT2.3dpi.4dpi.CC$Term), unique(AT2.3dpi.4dpi.CC$Term))) / 
  length(unique(iAT2.3dpi.4dpi.CC$Term)) # NA
length(unique(iAT2.3dpi.4dpi.CC$Term)) - 
  length(intersect(unique(iAT2.3dpi.4dpi.CC$Term), unique(AT2.3dpi.4dpi.CC$Term))) # 0
length(unique(AT2.3dpi.4dpi.CC$Term)) - 
  length(intersect(unique(iAT2.3dpi.4dpi.CC$Term), unique(AT2.3dpi.4dpi.CC$Term))) # 7
overlap.CC <- intersect(unique(iAT2.3dpi.4dpi.CC$Term), unique(AT2.3dpi.4dpi.CC$Term))
overlap.CC
write.csv(cbind(Overlap = iAT2.3dpi.4dpi.CC$Term %in% overlap.CC, iAT2.3dpi.4dpi.CC), 
          quote = F, row.names = F, 
          paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_mock/", 
                 "iAT2_cellular_component.csv"))
write.csv(cbind(Overlap = AT2.3dpi.4dpi.CC$Term %in% overlap.CC, AT2.3dpi.4dpi.CC), 
          quote = F, row.names = F, 
          paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_mock/", 
                 "AT2_cellular_component.csv"))


# BP
iAT2.3dpi.4dpi.BP <- iAT2.enriched.3dpi.4dpi.mock$GO_Biological_Process_2018
AT2.3dpi.4dpi.BP <- AT2.enriched.3dpi.4dpi.mock$GO_Biological_Process_2018
iAT2.3dpi.4dpi.BP <- iAT2.3dpi.4dpi.BP[iAT2.3dpi.4dpi.BP$Adjusted.P.value < pval.cutoff,]
dim(iAT2.3dpi.4dpi.BP) # 31
AT2.3dpi.4dpi.BP <- AT2.3dpi.4dpi.BP[AT2.3dpi.4dpi.BP$Adjusted.P.value < pval.cutoff,]
dim(AT2.3dpi.4dpi.BP) # 183
length(intersect(unique(iAT2.3dpi.4dpi.BP$Term), unique(AT2.3dpi.4dpi.BP$Term))) # 5
length(intersect(unique(iAT2.3dpi.4dpi.BP$Term), unique(AT2.3dpi.4dpi.BP$Term))) / 
  length(unique(iAT2.3dpi.4dpi.BP$Term)) # 0.1612903
length(unique(iAT2.3dpi.4dpi.BP$Term)) - 
  length(intersect(unique(iAT2.3dpi.4dpi.BP$Term), unique(AT2.3dpi.4dpi.BP$Term))) # 26
length(unique(AT2.3dpi.4dpi.BP$Term)) - 
  length(intersect(unique(iAT2.3dpi.4dpi.BP$Term), unique(AT2.3dpi.4dpi.BP$Term))) # 178
overlap.BP <- intersect(unique(iAT2.3dpi.4dpi.BP$Term), unique(AT2.3dpi.4dpi.BP$Term))
overlap.BP
write.csv(cbind(Overlap = iAT2.3dpi.4dpi.BP$Term %in% overlap.BP, iAT2.3dpi.4dpi.BP), 
          quote = F, row.names = F, 
          paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_mock/", 
                 "iAT2_biological_process.csv"))
write.csv(cbind(Overlap = AT2.3dpi.4dpi.BP$Term %in% overlap.BP, AT2.3dpi.4dpi.BP), 
          quote = F, row.names = F, 
          paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_mock/", 
                 "AT2_biological_process.csv"))



# KEGG
iAT2.3dpi.4dpi.KEGG <- iAT2.enriched.3dpi.4dpi.mock$KEGG_2019_Human
AT2.3dpi.4dpi.KEGG <- AT2.enriched.3dpi.4dpi.mock$KEGG_2019_Human
iAT2.3dpi.4dpi.KEGG <- iAT2.3dpi.4dpi.KEGG[iAT2.3dpi.4dpi.KEGG$Adjusted.P.value < pval.cutoff,]
dim(iAT2.3dpi.4dpi.KEGG) # 4
AT2.3dpi.4dpi.KEGG <- AT2.3dpi.4dpi.KEGG[AT2.3dpi.4dpi.KEGG$Adjusted.P.value < pval.cutoff,]
dim(AT2.3dpi.4dpi.KEGG) # 25
length(intersect(unique(iAT2.3dpi.4dpi.KEGG$Term), unique(AT2.3dpi.4dpi.KEGG$Term))) # 0
length(intersect(unique(iAT2.3dpi.4dpi.KEGG$Term), unique(AT2.3dpi.4dpi.KEGG$Term))) / 
  length(unique(iAT2.3dpi.4dpi.KEGG$Term)) # 0
length(unique(iAT2.3dpi.4dpi.KEGG$Term)) - 
  length(intersect(unique(iAT2.3dpi.4dpi.KEGG$Term), unique(AT2.3dpi.4dpi.KEGG$Term))) # 4
length(unique(AT2.3dpi.4dpi.KEGG$Term)) - 
  length(intersect(unique(iAT2.3dpi.4dpi.KEGG$Term), unique(AT2.3dpi.4dpi.KEGG$Term))) # 25
overlap.KEGG <- intersect(unique(iAT2.3dpi.4dpi.KEGG$Term), unique(AT2.3dpi.4dpi.KEGG$Term))
overlap.KEGG
write.csv(cbind(Overlap = iAT2.3dpi.4dpi.KEGG$Term %in% overlap.KEGG, iAT2.3dpi.4dpi.KEGG), 
          quote = F, row.names = F, 
          paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_mock/", 
                 "iAT2_KEGG.csv"))
write.csv(cbind(Overlap = AT2.3dpi.4dpi.KEGG$Term %in% overlap.KEGG, AT2.3dpi.4dpi.KEGG), 
          quote = F, row.names = F, 
          paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_mock/", 
                 "AT2_KEGG.csv"))








# iAT2 DEGs between 3dpi v.s. mock
# AT2 DEGs between 4dpi v.s. 1dpi
length(intersect(x = uniq.AT2.4dpi$Genes.SYMBOL, 
                 y = uniq.iAT2.3dpi$hgnc_symbols)) / 
  length(uniq.iAT2.3dpi$hgnc_symbols) # 50.4%
length(intersect(x = uniq.AT2.4dpi$Genes.SYMBOL, 
                 y = uniq.iAT2.3dpi$hgnc_symbols)) # 58
length(uniq.iAT2.3dpi$hgnc_symbols) - 
  length(intersect(x = uniq.AT2.4dpi$Genes.SYMBOL, 
                   y = uniq.iAT2.3dpi$hgnc_symbols)) # 57
length(uniq.AT2.4dpi$Genes.SYMBOL) - 
  length(intersect(x = uniq.AT2.4dpi$Genes.SYMBOL, 
                   y = uniq.iAT2.3dpi$hgnc_symbols)) # 3805
dir.create(paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_1dpi/"))
overlap.3dpi.4dpi <- intersect(x = uniq.AT2.4dpi$Genes.SYMBOL, 
                               y = uniq.iAT2.3dpi$hgnc_symbols)
length(overlap.3dpi.4dpi)
write.csv(cbind(Overlap = uniq.iAT2.3dpi$hgnc_symbols %in% overlap.3dpi.4dpi, 
                uniq.iAT2.3dpi), quote = F, row.names = F, 
          paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_1dpi/", 
                 "iAT2_DEGs_3dpi_vs_mock.csv"))
write.csv(cbind(Overlap = uniq.AT2.4dpi$Genes.SYMBOL %in% overlap.3dpi.4dpi, 
                uniq.AT2.4dpi), quote = F, row.names = F, 
          paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_1dpi/", 
                 "AT2_DEGs_4dpi_vs_1dpi.csv"))


iAT2.enriched.3dpi.4dpi <- run_GO_and_KEGG(genes.ll = unique(uniq.iAT2.3dpi$hgnc_symbols), 
                                                org = "human")
AT2.enriched.3dpi.4dpi <- run_GO_and_KEGG(genes.ll = unique(uniq.AT2.4dpi$Genes.SYMBOL ), 
                                               org = "human")


# MF
iAT2.3dpi.4dpi.MF <- iAT2.enriched.3dpi.4dpi$GO_Molecular_Function_2018
AT2.3dpi.4dpi.MF <- AT2.enriched.3dpi.4dpi$GO_Molecular_Function_2018
iAT2.3dpi.4dpi.MF <- iAT2.3dpi.4dpi.MF[iAT2.3dpi.4dpi.MF$Adjusted.P.value < pval.cutoff,]
dim(iAT2.3dpi.4dpi.MF) # 2
AT2.3dpi.4dpi.MF <- AT2.3dpi.4dpi.MF[AT2.3dpi.4dpi.MF$Adjusted.P.value < pval.cutoff,]
dim(AT2.3dpi.4dpi.MF) # 87
length(intersect(unique(iAT2.3dpi.4dpi.MF$Term), unique(AT2.3dpi.4dpi.MF$Term))) # 1
length(intersect(unique(iAT2.3dpi.4dpi.MF$Term), unique(AT2.3dpi.4dpi.MF$Term))) / 
  length(unique(iAT2.3dpi.4dpi.MF$Term)) # 50.0%
length(unique(iAT2.3dpi.4dpi.MF$Term)) - 
  length(intersect(unique(iAT2.3dpi.4dpi.MF$Term), unique(AT2.3dpi.4dpi.MF$Term))) # 1
length(unique(AT2.3dpi.4dpi.MF$Term)) - 
  length(intersect(unique(iAT2.3dpi.4dpi.MF$Term), unique(AT2.3dpi.4dpi.MF$Term))) # 86
overlap.MF <- intersect(unique(iAT2.3dpi.4dpi.MF$Term), unique(AT2.3dpi.4dpi.MF$Term))
overlap.MF
dir.create(paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_1dpi/"))
write.csv(cbind(Overlap = iAT2.3dpi.4dpi.MF$Term %in% overlap.MF, iAT2.3dpi.4dpi.MF), 
          quote = F, row.names = F, 
          paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_1dpi/", 
                 "iAT2_molecular_function.csv"))
write.csv(cbind(Overlap = AT2.3dpi.4dpi.MF$Term %in% overlap.MF, AT2.3dpi.4dpi.MF), 
          quote = F, row.names = F, 
          paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_1dpi/", 
                 "AT2_molecular_function.csv"))


# CC
iAT2.3dpi.4dpi.CC <- iAT2.enriched.3dpi.4dpi.mock$GO_Cellular_Component_2018
AT2.3dpi.4dpi.CC <- AT2.enriched.3dpi.4dpi.mock$GO_Cellular_Component_2018
iAT2.3dpi.4dpi.CC <- iAT2.3dpi.4dpi.CC[iAT2.3dpi.4dpi.CC$Adjusted.P.value < pval.cutoff,]
dim(iAT2.3dpi.4dpi.CC)
AT2.3dpi.4dpi.CC <- AT2.3dpi.4dpi.CC[AT2.3dpi.4dpi.CC$Adjusted.P.value < pval.cutoff,]
dim(AT2.3dpi.4dpi.CC)
length(intersect(unique(iAT2.3dpi.4dpi.CC$Term), unique(AT2.3dpi.4dpi.CC$Term))) # 0
length(intersect(unique(iAT2.3dpi.4dpi.CC$Term), unique(AT2.3dpi.4dpi.CC$Term))) / 
  length(unique(iAT2.3dpi.4dpi.CC$Term)) # NA
length(unique(iAT2.3dpi.4dpi.CC$Term)) - 
  length(intersect(unique(iAT2.3dpi.4dpi.CC$Term), unique(AT2.3dpi.4dpi.CC$Term))) # 0
length(unique(AT2.3dpi.4dpi.CC$Term)) - 
  length(intersect(unique(iAT2.3dpi.4dpi.CC$Term), unique(AT2.3dpi.4dpi.CC$Term))) # 54
overlap.CC <- intersect(unique(iAT2.3dpi.4dpi.CC$Term), unique(AT2.3dpi.4dpi.CC$Term))
overlap.CC
write.csv(cbind(Overlap = iAT2.3dpi.4dpi.CC$Term %in% overlap.CC, iAT2.3dpi.4dpi.CC), 
          quote = F, row.names = F, 
          paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_1dpi/", 
                 "iAT2_cellular_component.csv"))
write.csv(cbind(Overlap = AT2.3dpi.4dpi.CC$Term %in% overlap.CC, AT2.3dpi.4dpi.CC), 
          quote = F, row.names = F, 
          paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_1dpi/", 
                 "AT2_cellular_component.csv"))


# BP
iAT2.3dpi.4dpi.BP <- iAT2.enriched.3dpi.4dpi.mock$GO_Biological_Process_2018
AT2.3dpi.4dpi.BP <- AT2.enriched.3dpi.4dpi.mock$GO_Biological_Process_2018
iAT2.3dpi.4dpi.BP <- iAT2.3dpi.4dpi.BP[iAT2.3dpi.4dpi.BP$Adjusted.P.value < pval.cutoff,]
dim(iAT2.3dpi.4dpi.BP)
AT2.3dpi.4dpi.BP <- AT2.3dpi.4dpi.BP[AT2.3dpi.4dpi.BP$Adjusted.P.value < pval.cutoff,]
dim(AT2.3dpi.4dpi.BP)
length(intersect(unique(iAT2.3dpi.4dpi.BP$Term), unique(AT2.3dpi.4dpi.BP$Term))) # 1
length(intersect(unique(iAT2.3dpi.4dpi.BP$Term), unique(AT2.3dpi.4dpi.BP$Term))) / 
  length(unique(iAT2.3dpi.4dpi.BP$Term)) # 0.03225806
length(unique(iAT2.3dpi.4dpi.BP$Term)) - 
  length(intersect(unique(iAT2.3dpi.4dpi.BP$Term), unique(AT2.3dpi.4dpi.BP$Term))) # 30
length(unique(AT2.3dpi.4dpi.BP$Term)) - 
  length(intersect(unique(iAT2.3dpi.4dpi.BP$Term), unique(AT2.3dpi.4dpi.BP$Term))) # 228
overlap.BP <- intersect(unique(iAT2.3dpi.4dpi.BP$Term), unique(AT2.3dpi.4dpi.BP$Term))
overlap.BP
write.csv(cbind(Overlap = iAT2.3dpi.4dpi.BP$Term %in% overlap.BP, iAT2.3dpi.4dpi.BP), 
          quote = F, row.names = F, 
          paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_1dpi/", 
                 "iAT2_biological_process.csv"))
write.csv(cbind(Overlap = AT2.3dpi.4dpi.BP$Term %in% overlap.BP, AT2.3dpi.4dpi.BP), 
          quote = F, row.names = F, 
          paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_1dpi/", 
                 "AT2_biological_process.csv"))



# KEGG
iAT2.3dpi.4dpi.KEGG <- iAT2.enriched.3dpi.4dpi.mock$KEGG_2019_Human
AT2.3dpi.4dpi.KEGG <- AT2.enriched.3dpi.4dpi.mock$KEGG_2019_Human
iAT2.3dpi.4dpi.KEGG <- iAT2.3dpi.4dpi.KEGG[iAT2.3dpi.4dpi.KEGG$Adjusted.P.value < pval.cutoff,]
dim(iAT2.3dpi.4dpi.KEGG) # 4
AT2.3dpi.4dpi.KEGG <- AT2.3dpi.4dpi.KEGG[AT2.3dpi.4dpi.KEGG$Adjusted.P.value < pval.cutoff,]
dim(AT2.3dpi.4dpi.KEGG) # 42
length(intersect(unique(iAT2.3dpi.4dpi.KEGG$Term), unique(AT2.3dpi.4dpi.KEGG$Term))) # 2
length(intersect(unique(iAT2.3dpi.4dpi.KEGG$Term), unique(AT2.3dpi.4dpi.KEGG$Term))) / 
  length(unique(iAT2.3dpi.4dpi.KEGG$Term)) # 0.50
length(unique(iAT2.3dpi.4dpi.KEGG$Term)) - 
  length(intersect(unique(iAT2.3dpi.4dpi.KEGG$Term), unique(AT2.3dpi.4dpi.KEGG$Term))) # 2
length(unique(AT2.3dpi.4dpi.KEGG$Term)) - 
  length(intersect(unique(iAT2.3dpi.4dpi.KEGG$Term), unique(AT2.3dpi.4dpi.KEGG$Term))) # 40
overlap.KEGG <- intersect(unique(iAT2.3dpi.4dpi.KEGG$Term), unique(AT2.3dpi.4dpi.KEGG$Term))
overlap.KEGG
write.csv(cbind(Overlap = iAT2.3dpi.4dpi.KEGG$Term %in% overlap.KEGG, iAT2.3dpi.4dpi.KEGG), 
          quote = F, row.names = F, 
          paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_1dpi/", 
                 "iAT2_KEGG.csv"))
write.csv(cbind(Overlap = AT2.3dpi.4dpi.KEGG$Term %in% overlap.KEGG, AT2.3dpi.4dpi.KEGG), 
          quote = F, row.names = F, 
          paste0(table.dir, "iAT2_3dpi_vs_mock_and_AT2_4dpi_vs_1dpi/", 
                 "AT2_KEGG.csv"))