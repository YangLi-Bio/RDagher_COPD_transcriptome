common2 <- intersect(unique(unlist(table2[!is.na(table2$hgnc_symbols), "hgnc_symbols"])),
                     rownames(degs.public))
deg2.common <- deg2[deg2$hgnc_symbols %in% common2,]
dim(deg2.common) # 165
deg.copd.common <- degs.public[degs.public$p_val_adj < padj.cutoff & 
                                 degs.public$avg_log2FC > 0 & 
                                 rownames(degs.public) %in% common2,]
dim(deg.copd.common) # 6233
overlap2 <- intersect(rownames(deg.copd.common),
                      deg2.common$hgnc_symbols)
length(overlap2) # 43
length(overlap2) / nrow(deg2.common) # 26.1%
length(rownames(deg.copd.common)) - length(overlap2) # 6190
length(deg2.common$hgnc_symbols) - length(overlap2) # 122
write.csv(cbind(Overlap = deg2.common$hgnc_symbols %in% overlap2, 
                deg2.common), 
          quote = F, paste0(table.dir, "Our_DEGs_comp_2.csv"))
write.csv(cbind(Overlap = rownames(deg.copd.common) %in% overlap2, 
                deg.copd.common), 
          quote = F, paste0(table.dir, "Public_DEGs_comp_2.csv"))


top.pub.deg2 <- deg.copd.common[1:200,]
top.our.deg2 <- deg2.common[1:200,]
dim(top.pub.deg2)
dim(top.our.deg2)
top.overlap2 <- intersect(rownames(top.pub.deg2),
                          top.our.deg2$hgnc_symbols)
length(top.overlap2) # 5
length(top.overlap2) / nrow(top.our.deg2) # 2.5%
length(rownames(top.pub.deg2)) - length(top.overlap2) # 195
length(top.our.deg2$hgnc_symbols) - length(top.overlap2) # 195
write.csv(cbind(Overlap = top.our.deg2$hgnc_symbols %in% top.overlap2, 
                top.our.deg2), 
          quote = F, paste0(table.dir, "Top200_Our_DEGs_comp_2.csv"))
write.csv(cbind(Overlap = rownames(top.pub.deg2) %in% top.overlap2, 
                top.pub.deg2), 
          quote = F, paste0(table.dir, "Top200_Public_DEGs_comp_2.csv"))



# GO enrichment analyses for our data
our.enriched2 <- run_GO_and_KEGG(genes.ll = deg2.common$hgnc_symbols, org = "Human")
our.mf2 <- our.enriched2$GO_Molecular_Function_2018
our.cc2 <- our.enriched2$GO_Cellular_Component_2018
our.bp2 <- our.enriched2$GO_Biological_Process_2018
our.kegg2 <- our.enriched2$KEGG_2019_Human


our.mf2 <- our.mf2[our.mf2$Adjusted.P.value < 0.05,]
our.cc2 <- our.cc2[our.cc2$Adjusted.P.value < 0.05,]
our.bp2 <- our.bp2[our.bp2$Adjusted.P.value < 0.05,]
our.kegg2 <- our.kegg2[our.kegg2$Adjusted.P.value < 0.05,]




# Enrichment for public data
pub.copd.enriched <- run_GO_and_KEGG(genes.ll = rownames(deg.copd.common), org = "Human")
pub.copd.mf2 <- pub.copd.enriched$GO_Molecular_Function_2018
pub.copd.cc2 <- pub.copd.enriched$GO_Cellular_Component_2018
pub.copd.bp2 <- pub.copd.enriched$GO_Biological_Process_2018
pub.copd.kegg2 <- pub.copd.enriched$KEGG_2019_Human


pub.copd.mf2 <- pub.copd.mf2[pub.copd.mf2$Adjusted.P.value < 0.05,]
pub.copd.cc2 <- pub.copd.cc2[pub.copd.cc2$Adjusted.P.value < 0.05,]
pub.copd.bp2 <- pub.copd.bp2[pub.copd.bp2$Adjusted.P.value < 0.05,]
pub.copd.kegg2 <- pub.copd.kegg2[pub.copd.kegg2$Adjusted.P.value < 0.05,]



# MF
overlap.mf <- intersect(our.mf2$Term, pub.copd.mf2$Term)
length(overlap.mf) / nrow(our.mf2) # 7.1%
length(overlap.mf) # 1
nrow(our.mf2) - length(overlap.mf) # 13
nrow(pub.copd.mf2) - length(overlap.mf) # 127
write.csv(cbind(Overlap = our.mf2$Term %in% overlap.mf, 
                our.mf2), quote = F, 
          paste0(table.dir, "Our_Molecular_Function_comp_2.csv"))
write.csv(cbind(Overlap = pub.copd.mf2$Term %in% overlap.mf, 
                pub.copd.mf2), quote = F, 
          paste0(table.dir, "Public_Molecular_Function_comp_2.csv"))




# CC
overlap.cc <- intersect(our.cc2$Term, pub.copd.cc2$Term)
length(overlap.cc) / nrow(our.cc2) # 0%
length(overlap.cc) # 0
nrow(our.cc2) - length(overlap.cc) # 0
nrow(pub.copd.cc2) - length(overlap.cc) # 159
write.csv(cbind(Overlap = our.cc2$Term %in% overlap.cc, 
                our.cc2), quote = F, 
          paste0(table.dir, "Our_Cellular_Component_comp_2.csv"))
write.csv(cbind(Overlap = pub.copd.cc2$Term %in% overlap.cc, 
                pub.copd.cc2), quote = F, 
          paste0(table.dir, "Public_Cellular_Component_comp_2.csv"))


# BP
overlap.bp <- intersect(our.bp2$Term, pub.copd.bp2$Term)
length(overlap.bp) / nrow(our.bp2) # 11.1%
length(overlap.bp) # 1
nrow(our.bp2) - length(overlap.bp) # 8
nrow(pub.copd.bp2) - length(overlap.bp) # 671
write.csv(cbind(Overlap = our.bp2$Term %in% overlap.bp, 
                our.bp2), quote = F, 
          paste0(table.dir, "Our_Biological_Process_comp_2.csv"))
write.csv(cbind(Overlap = pub.copd.bp2$Term %in% overlap.bp, 
                pub.copd.bp2), quote = F, 
          paste0(table.dir, "Public_Biological_Process_comp_2.csv"))


# kegg
overlap.kegg <- intersect(our.kegg2$Term, pub.copd.kegg2$Term)
length(overlap.kegg) / nrow(our.kegg2) # 0%
length(overlap.kegg) # 0
nrow(our.kegg2) - length(overlap.kegg) # 0
nrow(pub.copd.kegg2) - length(overlap.kegg) # 99
write.csv(cbind(Overlap = our.kegg2$Term %in% overlap.kegg, 
                our.kegg2), quote = F, 
          paste0(table.dir, "Our_KEGG_comp_2.csv"))
write.csv(cbind(Overlap = pub.copd.kegg2$Term %in% overlap.kegg, 
                pub.copd.kegg2), quote = F, 
          paste0(table.dir, "Public_KEGG_comp_2.csv"))


# Top-200


# Our Top-200
top.our.enriched2 <- run_GO_and_KEGG(genes.ll = top.our.deg2$hgnc_symbols, org = "Human")
top.our.mf2 <- top.our.enriched2$GO_Molecular_Function_2018
top.our.cc2 <- top.our.enriched2$GO_Cellular_Component_2018
top.our.bp2 <- top.our.enriched2$GO_Biological_Process_2018
top.our.kegg2 <- top.our.enriched2$KEGG_2019_Human


top.our.mf2 <- top.our.mf2[top.our.mf2$Adjusted.P.value < 0.05,]
top.our.cc2 <- top.our.cc2[top.our.cc2$Adjusted.P.value < 0.05,]
top.our.bp2 <- top.our.bp2[top.our.bp2$Adjusted.P.value < 0.05,]
top.our.kegg2 <- top.our.kegg2[top.our.kegg2$Adjusted.P.value < 0.05,]


# Public Top-200
top.pub.copd.enriched <- run_GO_and_KEGG(genes.ll = rownames(top.pub.deg2), org = "Human")
top.pub.copd.mf2 <- top.pub.copd.enriched$GO_Molecular_Function_2018
top.pub.copd.cc2 <- top.pub.copd.enriched$GO_Cellular_Component_2018
top.pub.copd.bp2 <- top.pub.copd.enriched$GO_Biological_Process_2018
top.pub.copd.kegg2 <- top.pub.copd.enriched$KEGG_2019_Human


top.pub.copd.mf2 <- top.pub.copd.mf2[top.pub.copd.mf2$Adjusted.P.value < 0.05,]
top.pub.copd.cc2 <- top.pub.copd.cc2[top.pub.copd.cc2$Adjusted.P.value < 0.05,]
top.pub.copd.bp2 <- top.pub.copd.bp2[top.pub.copd.bp2$Adjusted.P.value < 0.05,]
top.pub.copd.kegg2 <- top.pub.copd.kegg2[top.pub.copd.kegg2$Adjusted.P.value < 0.05,]



# MF top-200
top.overlap.mf <- intersect(top.our.mf2$Term, top.pub.copd.mf2$Term)
length(top.overlap.mf) / nrow(top.our.mf2) # 18.8%
length(top.overlap.mf) # 3
nrow(top.our.mf2) - length(top.overlap.mf) # 14
nrow(top.pub.copd.mf2) - length(top.overlap.mf) # 37
write.csv(cbind(top.overlap = top.our.mf2$Term %in% top.overlap.mf, 
                top.our.mf2), quote = F, 
          paste0(table.dir, "Top200_Our_Molecular_Function_comp_2.csv"))
write.csv(cbind(top.overlap = top.pub.copd.mf2$Term %in% top.overlap.mf, 
                top.pub.copd.mf2), quote = F, 
          paste0(table.dir, "Top200_Public_Molecular_Function_comp_2.csv"))



# cc top-200
top.overlap.cc <- intersect(top.our.cc2$Term, top.pub.copd.cc2$Term)
length(top.overlap.cc) / nrow(top.our.cc2) # 0%
length(top.overlap.cc) # 0
nrow(top.our.cc2) - length(top.overlap.cc) # 0
nrow(top.pub.copd.cc2) - length(top.overlap.cc) # 40
write.csv(cbind(top.overlap = top.our.cc2$Term %in% top.overlap.cc, 
                top.our.cc2), quote = F, 
          paste0(table.dir, "Top200_Our_Cellular_Component_comp_2.csv"))
write.csv(cbind(top.overlap = top.pub.copd.cc2$Term %in% top.overlap.cc, 
                top.pub.copd.cc2), quote = F, 
          paste0(table.dir, "Top200_Public_Cellular_Component_comp_2.csv"))




# bp top-200
top.overlap.bp <- intersect(top.our.bp2$Term, top.pub.copd.bp2$Term)
length(top.overlap.bp) / nrow(top.our.bp2) # 33.3%
length(top.overlap.bp) # 3
nrow(top.our.bp2) - length(top.overlap.bp) # 6
nrow(top.pub.copd.bp2) - length(top.overlap.bp) # 274
write.csv(cbind(top.overlap = top.our.bp2$Term %in% top.overlap.bp, 
                top.our.bp2), quote = F, 
          paste0(table.dir, "Top200_Our_Biological_Process_comp_2.csv"))
write.csv(cbind(top.overlap = top.pub.copd.bp2$Term %in% top.overlap.bp, 
                top.pub.copd.bp2), quote = F, 
          paste0(table.dir, "Top200_Public_Biological_Process_comp_2.csv"))





# kegg top-200
top.overlap.kegg <- intersect(top.our.kegg2$Term, top.pub.copd.kegg2$Term)
length(top.overlap.kegg) / nrow(top.our.kegg2) # 0%
length(top.overlap.kegg) # 0
nrow(top.our.kegg2) - length(top.overlap.kegg) # 0
nrow(top.pub.copd.kegg2) - length(top.overlap.kegg) # 35
write.csv(cbind(top.overlap = top.our.kegg2$Term %in% top.overlap.kegg, 
                top.our.kegg2), quote = F, 
          paste0(table.dir, "Top200_Our_KEGG_comp_2.csv"))
write.csv(cbind(top.overlap = top.pub.copd.kegg2$Term %in% top.overlap.kegg, 
                top.pub.copd.kegg2), quote = F, 
          paste0(table.dir, "Top200_Public_KEGG_comp_2.csv"))
