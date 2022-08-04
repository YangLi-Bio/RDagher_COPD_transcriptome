degs.public <- qs::qread(paste0(R.dir, "DEGs_bimod_Fibr_COPD_public.qsave"))
range(degs.public$avg_log2FC)
range(degs.public$p_val_adj)
common3 <- intersect(unique(unlist(table3[!is.na(table3$hgnc_symbols), "hgnc_symbols"])),
                     rownames(degs.public))
length(common3) # 5618
deg3.common <- deg3[deg3$hgnc_symbols %in% common3,]
dim(deg3.common) # 5
deg.copd.common <- degs.public[degs.public$p_val_adj < padj.cutoff & 
                                 # degs.public$avg_log2FC > 0 & 
                                 rownames(degs.public) %in% common3,]
dim(deg.copd.common) # 2779
overlap3 <- intersect(rownames(deg.copd.common),
                      deg3.common$hgnc_symbols)
length(overlap3) # 1
length(overlap3) / nrow(deg3.common) # 20%
length(rownames(deg.copd.common)) - length(overlap3) # 2778
length(deg3.common$hgnc_symbols) - length(overlap3) # 4
write.csv(cbind(Overlap = deg3.common$hgnc_symbols %in% overlap3, 
                deg3.common), 
          quote = F, paste0(table.dir, "Our_DEGs_comp_3.csv"))
write.csv(cbind(Overlap = rownames(deg.copd.common) %in% overlap3, 
                deg.copd.common), 
          quote = F, paste0(table.dir, "Public_DEGs_comp_3.csv"))


top.pub.deg3 <- deg.copd.common[1:200,]
top.our.deg3 <- deg3.common[1:200,]
dim(top.pub.deg3)
dim(top.our.deg3)
top.overlap3 <- intersect(rownames(top.pub.deg3),
                          top.our.deg3$hgnc_symbols)
length(top.overlap3) # 0
length(top.overlap3) / nrow(top.our.deg3) # 0%
length(rownames(top.pub.deg3)) - length(top.overlap3) # 200
length(top.our.deg3$hgnc_symbols) - length(top.overlap3) # 200
write.csv(cbind(Overlap = top.our.deg3$hgnc_symbols %in% top.overlap3, 
                top.our.deg3), 
          quote = F, paste0(table.dir, "Top200_Our_DEGs_comp_3.csv"))
write.csv(cbind(Overlap = rownames(top.pub.deg3) %in% top.overlap3, 
                top.pub.deg3), 
          quote = F, paste0(table.dir, "Top200_Public_DEGs_comp_3.csv"))



# GO enrichment analyses for our data
our.enriched3 <- run_GO_and_KEGG(genes.ll = deg3.common$hgnc_symbols, org = "Human")
our.mf3 <- our.enriched3$GO_Molecular_Function_2018
our.cc3 <- our.enriched3$GO_Cellular_Component_2018
our.bp3 <- our.enriched3$GO_Biological_Process_2018
our.kegg3 <- our.enriched3$KEGG_2019_Human


our.mf3 <- our.mf3[our.mf3$Adjusted.P.value < 0.05,]
our.cc3 <- our.cc3[our.cc3$Adjusted.P.value < 0.05,]
our.bp3 <- our.bp3[our.bp3$Adjusted.P.value < 0.05,]
our.kegg3 <- our.kegg3[our.kegg3$Adjusted.P.value < 0.05,]




# Enrichment for public data
pub.copd.enriched <- run_GO_and_KEGG(genes.ll = rownames(deg.copd.common), org = "Human")
pub.copd.mf3 <- pub.copd.enriched$GO_Molecular_Function_2018
pub.copd.cc3 <- pub.copd.enriched$GO_Cellular_Component_2018
pub.copd.bp3 <- pub.copd.enriched$GO_Biological_Process_2018
pub.copd.kegg3 <- pub.copd.enriched$KEGG_2019_Human


pub.copd.mf3 <- pub.copd.mf3[pub.copd.mf3$Adjusted.P.value < 0.05,]
pub.copd.cc3 <- pub.copd.cc3[pub.copd.cc3$Adjusted.P.value < 0.05,]
pub.copd.bp3 <- pub.copd.bp3[pub.copd.bp3$Adjusted.P.value < 0.05,]
pub.copd.kegg3 <- pub.copd.kegg3[pub.copd.kegg3$Adjusted.P.value < 0.05,]



# MF
overlap.mf <- intersect(our.mf3$Term, pub.copd.mf3$Term)
length(overlap.mf) / nrow(our.mf3) # 0%
length(overlap.mf) # 0
nrow(our.mf3) - length(overlap.mf) # 1
nrow(pub.copd.mf3) - length(overlap.mf) # 53
write.csv(cbind(Overlap = our.mf3$Term %in% overlap.mf, 
                our.mf3), quote = F, 
          paste0(table.dir, "Our_Molecular_Function_comp_3.csv"))
write.csv(cbind(Overlap = pub.copd.mf3$Term %in% overlap.mf, 
                pub.copd.mf3), quote = F, 
          paste0(table.dir, "Public_Molecular_Function_comp_3.csv"))




# CC
overlap.cc <- intersect(our.cc3$Term, pub.copd.cc3$Term)
length(overlap.cc) / nrow(our.cc3) # 0%
length(overlap.cc) # 0
nrow(our.cc3) - length(overlap.cc) # 0
nrow(pub.copd.cc3) - length(overlap.cc) # 112
write.csv(cbind(Overlap = our.cc3$Term %in% overlap.cc, 
                our.cc3), quote = F, 
          paste0(table.dir, "Our_Cellular_Component_comp_3.csv"))
write.csv(cbind(Overlap = pub.copd.cc3$Term %in% overlap.cc, 
                pub.copd.cc3), quote = F, 
          paste0(table.dir, "Public_Cellular_Component_comp_3.csv"))


# BP
overlap.bp <- intersect(our.bp3$Term, pub.copd.bp3$Term)
length(overlap.bp) / nrow(our.bp3) # 0%
length(overlap.bp) # 0
nrow(our.bp3) - length(overlap.bp) # 35
nrow(pub.copd.bp3) - length(overlap.bp) # 481
write.csv(cbind(Overlap = our.bp3$Term %in% overlap.bp, 
                our.bp3), quote = F, 
          paste0(table.dir, "Our_Biological_Process_comp_3.csv"))
write.csv(cbind(Overlap = pub.copd.bp3$Term %in% overlap.bp, 
                pub.copd.bp3), quote = F, 
          paste0(table.dir, "Public_Biological_Process_comp_3.csv"))


# kegg
overlap.kegg <- intersect(our.kegg3$Term, pub.copd.kegg3$Term)
length(overlap.kegg) / nrow(our.kegg3) # 25%
length(overlap.kegg) # 1
nrow(our.kegg3) - length(overlap.kegg) # 3
nrow(pub.copd.kegg3) - length(overlap.kegg) # 92
write.csv(cbind(Overlap = our.kegg3$Term %in% overlap.kegg, 
                our.kegg3), quote = F, 
          paste0(table.dir, "Our_KEGG_comp_3.csv"))
write.csv(cbind(Overlap = pub.copd.kegg3$Term %in% overlap.kegg, 
                pub.copd.kegg3), quote = F, 
          paste0(table.dir, "Public_KEGG_comp_3.csv"))


# Top-200


# Our Top-200
top.our.enriched3 <- run_GO_and_KEGG(genes.ll = top.our.deg3$hgnc_symbols, org = "Human")
top.our.mf3 <- top.our.enriched3$GO_Molecular_Function_2018
top.our.cc3 <- top.our.enriched3$GO_Cellular_Component_2018
top.our.bp3 <- top.our.enriched3$GO_Biological_Process_2018
top.our.kegg3 <- top.our.enriched3$KEGG_2019_Human


top.our.mf3 <- top.our.mf3[top.our.mf3$Adjusted.P.value < 0.05,]
top.our.cc3 <- top.our.cc3[top.our.cc3$Adjusted.P.value < 0.05,]
top.our.bp3 <- top.our.bp3[top.our.bp3$Adjusted.P.value < 0.05,]
top.our.kegg3 <- top.our.kegg3[top.our.kegg3$Adjusted.P.value < 0.05,]


# Public Top-200
top.pub.copd.enriched <- run_GO_and_KEGG(genes.ll = rownames(top.pub.deg3), org = "Human")
top.pub.copd.mf3 <- top.pub.copd.enriched$GO_Molecular_Function_2018
top.pub.copd.cc3 <- top.pub.copd.enriched$GO_Cellular_Component_2018
top.pub.copd.bp3 <- top.pub.copd.enriched$GO_Biological_Process_2018
top.pub.copd.kegg3 <- top.pub.copd.enriched$KEGG_2019_Human


top.pub.copd.mf3 <- top.pub.copd.mf3[top.pub.copd.mf3$Adjusted.P.value < 0.05,]
top.pub.copd.cc3 <- top.pub.copd.cc3[top.pub.copd.cc3$Adjusted.P.value < 0.05,]
top.pub.copd.bp3 <- top.pub.copd.bp3[top.pub.copd.bp3$Adjusted.P.value < 0.05,]
top.pub.copd.kegg3 <- top.pub.copd.kegg3[top.pub.copd.kegg3$Adjusted.P.value < 0.05,]



# MF top-200
top.overlap.mf <- intersect(top.our.mf3$Term, top.pub.copd.mf3$Term)
length(top.overlap.mf) / nrow(top.our.mf3) # 00%
length(top.overlap.mf) # 00
nrow(top.our.mf3) - length(top.overlap.mf) # 1
nrow(top.pub.copd.mf3) - length(top.overlap.mf) # 14
write.csv(cbind(top.overlap = top.our.mf3$Term %in% top.overlap.mf, 
                top.our.mf3), quote = F, 
          paste0(table.dir, "Top200_Our_Molecular_Function_comp_3.csv"))
write.csv(cbind(top.overlap = top.pub.copd.mf3$Term %in% top.overlap.mf, 
                top.pub.copd.mf3), quote = F, 
          paste0(table.dir, "Top200_Public_Molecular_Function_comp_3.csv"))



# cc top-200
top.overlap.cc <- intersect(top.our.cc3$Term, top.pub.copd.cc3$Term)
length(top.overlap.cc) / nrow(top.our.cc3) # 0%
length(top.overlap.cc) # 0
nrow(top.our.cc3) - length(top.overlap.cc) # 0
nrow(top.pub.copd.cc3) - length(top.overlap.cc) # 13
write.csv(cbind(top.overlap = top.our.cc3$Term %in% top.overlap.cc, 
                top.our.cc3), quote = F, 
          paste0(table.dir, "Top200_Our_Cellular_Component_comp_3.csv"))
write.csv(cbind(top.overlap = top.pub.copd.cc3$Term %in% top.overlap.cc, 
                top.pub.copd.cc3), quote = F, 
          paste0(table.dir, "Top200_Public_Cellular_Component_comp_3.csv"))




# bp top-200
top.overlap.bp <- intersect(top.our.bp3$Term, top.pub.copd.bp3$Term)
length(top.overlap.bp) / nrow(top.our.bp3) # 0%
length(top.overlap.bp) # 0
nrow(top.our.bp3) - length(top.overlap.bp) # 34
nrow(top.pub.copd.bp3) - length(top.overlap.bp) # 41
write.csv(cbind(top.overlap = top.our.bp3$Term %in% top.overlap.bp, 
                top.our.bp3), quote = F, 
          paste0(table.dir, "Top200_Our_Biological_Process_comp_3.csv"))
write.csv(cbind(top.overlap = top.pub.copd.bp3$Term %in% top.overlap.bp, 
                top.pub.copd.bp3), quote = F, 
          paste0(table.dir, "Top200_Public_Biological_Process_comp_3.csv"))





# kegg top-200
top.overlap.kegg <- intersect(top.our.kegg3$Term, top.pub.copd.kegg3$Term)
length(top.overlap.kegg) / nrow(top.our.kegg3) # 0%
length(top.overlap.kegg) # 0
nrow(top.our.kegg3) - length(top.overlap.kegg) # 2
nrow(top.pub.copd.kegg3) - length(top.overlap.kegg) # 1
write.csv(cbind(top.overlap = top.our.kegg3$Term %in% top.overlap.kegg, 
                top.our.kegg3), quote = F, 
          paste0(table.dir, "Top200_Our_KEGG_comp_3.csv"))
write.csv(cbind(top.overlap = top.pub.copd.kegg3$Term %in% top.overlap.kegg, 
                top.pub.copd.kegg3), quote = F, 
          paste0(table.dir, "Top200_Public_KEGG_comp_3.csv"))
