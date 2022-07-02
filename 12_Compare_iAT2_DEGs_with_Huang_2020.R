###################################################################
#                                                                 #
# The following tasks will be conducted:                          #
# 1. Extract DEGs in mock (Figs. 3E and 3F) from reference and    #
#    compare them with our iAT2 DEGs (control)                    #
# 2. Extract DEGs in 1 and 4 dpi (Figs. 3E amd 3F) from           #
#    reference and compare them with our iAT2 DEGs (1, 2, 3 dpi)  #
# Question: now that only two pathways are selected in the        #
# reference, how about ours?                                      #
#                                                                 #
###################################################################


# Libraries
library(Seurat)
library(dplyr)
library(readxl)


# Global variables
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
R.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Rfiles/"
image.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Images/"
table.dir <- "/fs/ess/PCON0022/liyang/astrazeneca/copd_covid-19/Tables/"
iAT2.file <- "Supplemental_7._Figure_5E_DEG_results.xlsx"


# Select the DEGs of control and 1, 2, 3, dpi in iAT2 of our data
