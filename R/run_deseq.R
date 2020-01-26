
## DEseq for T-LGLL STAT5b mut vs wt

library(DESeq)


build_deseq_objects <- function(group1, group2, group1_name = "A", group2_name = "B"){
  
  # group1 = stat5mt_cd8
  # group2= lgl_stat5b_wt
  # group1_name = "STAT5bmt_cd8"
  # group2_name = "STAT5bwt"
  
  keep <- c(group1, group2)
  
  design <- as.factor(rep(c(group1_name, group2_name), c(length(group1), length(group2))))
  
  cds_temp <- newCountDataSet(countData = lgl_counts[ ,keep], conditions = design)
  
  dim(lgl_counts[ ,keep]) == length(design)
  
  # Remove lowly expressed genees
  cds_temp <- cds_temp[rowSums(1e+06 * counts(cds_temp)/expandAsMatrix(colSums(counts(cds_temp)), dim(cds_temp)) > 1) >= 3, ]
  
  ## Estimate size factors and dispersions
  cds_temp <- estimateSizeFactors(cds_temp)
  cds_temp <- estimateDispersions(cds_temp, fitType = "local")
  
  plotDispEsts(cds_temp)
  
  return(cds_temp)
  
}




estimate_de_deseq <- function(cds_temp, group1_name = "A", group2_name = "B"){
  
  # Estimate de-genes
  # group1_name = "STAT5bmt_all"
  # group2_name = "Healthy_all"
  # cds_temp = stat5mt_all_healthy_all_cds
  
  deseq_res <- nbinomTest(cds_temp, group2_name, group1_name)
  
  deseq_res$qval <- p.adjust(deseq_res$pval, method = "BH")
  deseq_res$sigf <- ifelse(deseq_res$qval < 0.05 & abs(deseq_res$log2FoldChange) > 1, "Significant", "Not_significant")
  deseq_res$direction <- ifelse(deseq_res$log2FoldChange > 0, "Up", "Down")
  
  ## Add biomart annotation
  de_bm <- get_bm(deseq_res$id)
  
  #  de_df_full <- merge(de_df, de_bm, by.x = "Row.names", by.y = "ensembl_gene_id")
  de_df_full <- merge(deseq_res, de_bm, by.x = "id", by.y = "ensembl_gene_id")
  
  return(de_df_full)
  
}




plot_volcano_plots_deseq <- function(de_df, title = ""){
  
  a <- ggplot(de_df, aes(x = log2FoldChange, y = -log10(pval), color = sigf, size = log10(baseMean))) + geom_point(alpha = 0.5) + labs(color = "qval < 0.05 and |logFC| > 1", subtitle = "DEseq") + scale_color_manual(values = c("grey", "dodgerblue")) + theme(legend.position = "top") + 
    geom_text(data = subset(de_df, sigf == "Significant"), aes(x = log2FoldChange, y = -log10(pval), label = external_gene_name), nudge_y = 0.5, fontface = "italic")
  
  return(grid.arrange(a, top = textGrob(title,gp=gpar(fontsize=20,font=3))))
  
}






# Building deseq-objects
stat5mt_all_healthy_all_deseq      <- build_deseq_objects(stat5mt_all, healthy_all,   "STAT5bmt_all", "Healthy_all")
stat5mt_all_lgl_wtl_deseq          <- build_deseq_objects(stat5mt_all, lgl_stat5b_wt, "STAT5bmt_all", "STAT5bwt")

stat5mt_cd8_cd48_healthy_all_deseq <- build_deseq_objects(c(stat5mt_cd8, stat5mt_cd4cd8), healthy_cd8, group1_name = "STAT5bmt_cd8_cd48", group2_name = "Healthy_cd8")
stat5mt_cd8_cd48_lgl_wt_deseq      <- build_deseq_objects(c(stat5mt_cd8, lgl_stat5b_wt), healthy_cd8, "STAT5bmt_cd8_cd48", "STAT5bwt")

stat5mt_cd8_healthy_cd8_deseq      <- build_deseq_objects(stat5mt_cd8, healthy_cd8, "STAT5bmt_cd8", "Healthy_cd8")
stat5mt_cd8_lgl_wt_deseq           <- build_deseq_objects(stat5mt_cd8, lgl_stat5b_wt, "STAT5bmt_cd8", "STAT5bwt")

healthy_cd8_healthy_cd4_deseq      <- build_deseq_objects(healthy_cd8, healthy_cd4, "Healthy_cd8", "Healthy_cd4")

?estimateDispersions




## Estimate the fit of dispersion

# Plot mean-var plots
pdf("stat5b_results/deseq_dispersion_plot_local.pdf", width = 12, height = 20)

par(mfrow=c(5,2))
plotDispEsts(stat5mt_all_healthy_all_deseq, main = a_title)
plotDispEsts(stat5mt_all_lgl_wtl_deseq, main = b_title)
plotDispEsts(stat5mt_cd8_cd48_healthy_all_deseq, main = c_title)
plotDispEsts(stat5mt_cd8_cd48_lgl_wt_deseq, main = d_title)
plotDispEsts(stat5mt_cd8_healthy_cd8_deseq, main = e_title)
plotDispEsts(stat5mt_cd8_lgl_wt_deseq, main = f_title)
plotDispEsts(healthy_cd8_healthy_cd4_deseq, main = g_title)

par(mfrow=c(1,1))
dev.off()

# Estimate de_deseq from de_deseqseq-objects
stat5mt_all_healthy_all_deseq_df      <- estimate_de_deseq(stat5mt_all_healthy_all_deseq, "STAT5bmt_all", "Healthy_all")
stat5mt_all_lgl_wtl_deseq_df          <- estimate_de_deseq(stat5mt_all_lgl_wtl_deseq, "STAT5bmt_all", "STAT5bwt")

stat5mt_cd8_cd48_healthy_all_deseq_df <- estimate_de_deseq(stat5mt_cd8_cd48_healthy_all_deseq, "STAT5bmt_cd8_cd48", "Healthy_cd8")
stat5mt_cd8_cd48_lgl_wt_deseq_df      <- estimate_de_deseq(stat5mt_cd8_cd48_lgl_wt_deseq, "STAT5bmt_cd8_cd48", "STAT5bwt")

stat5mt_cd8_healthy_cd8_deseq_df      <- estimate_de_deseq(stat5mt_cd8_healthy_cd8_deseq, "STAT5bmt_cd8", "Healthy_cd8")
stat5mt_cd8_lgl_wt_deseq_df           <- estimate_de_deseq(stat5mt_cd8_lgl_wt_deseq, "STAT5bmt_cd8", "STAT5bwt")

healthy_cd8_healthy_cd4_deseq_df      <- estimate_de_deseq(healthy_cd8_healthy_cd4_deseq, "Healthy_cd8", "Healthy_cd4")



write.table(stat5mt_all_healthy_all_deseq_df, paste0("stat5b_results/DE_tables_deseq/DE_local", a_name, ".txt"), sep = "\t")
write.table(stat5mt_all_lgl_wtl_deseq_df, paste0("stat5b_results/DE_tables_deseq/DE_local", b_name, ".txt"), sep = "\t")

write.table(stat5mt_cd8_cd48_healthy_all_deseq_df, paste0("stat5b_results/DE_tables_deseq/DE_local", c_name, ".txt"), sep = "\t")
write.table(stat5mt_cd8_cd48_lgl_wt_deseq_df, paste0("stat5b_results/DE_tables_deseq/DE_local", d_name, ".txt"), sep = "\t")

write.table(stat5mt_cd8_healthy_cd8_deseq_df, paste0("stat5b_results/DE_tables_deseq/DE_local", e_name, ".txt"), sep = "\t")
write.table(stat5mt_cd8_lgl_wt_deseq_df, paste0("stat5b_results/DE_tables_deseq/DE_local", f_name, ".txt"), sep = "\t")

write.table(healthy_cd8_healthy_cd4_deseq_df, paste0("stat5b_results/DE_tables_deseq/DE_local", g_name, ".txt"), sep = "\t")





## Read the files, too long to calculate again
stat5mt_all_healthy_all_deseq_df      <- read.delim(paste0("stat5b_results/DE_tables_deseq/DE_", a_name, ".txt"))
stat5mt_all_lgl_wtl_deseq_df          <- read.delim(paste0("stat5b_results/DE_tables_deseq/DE_", b_name, ".txt"))

stat5mt_cd8_cd48_healthy_all_deseq_df <- read.delim(paste0("stat5b_results/DE_tables_deseq/DE_", c_name, ".txt"))
stat5mt_cd8_cd48_lgl_wt_deseq_df      <- read.delim(paste0("stat5b_results/DE_tables_deseq/DE_", d_name, ".txt"))

stat5mt_cd8_healthy_cd8_deseq_df      <- read.delim(paste0("stat5b_results/DE_tables_deseq/DE_", e_name, ".txt"))
stat5mt_cd8_lgl_wt_deseq_df           <- read.delim(paste0("stat5b_results/DE_tables_deseq/DE_", f_name, ".txt"))

healthy_cd8_healthy_cd4_deseq_df      <- read.delim(paste0("stat5b_results/DE_tables_deseq/DE_", g_name, ".txt"))




# Plot volcanoplots_deseq of DE-genes
a <- plot_volcano_plots_deseq(stat5mt_all_healthy_all_deseq_df, title = a_title)
b <- plot_volcano_plots_deseq(stat5mt_all_lgl_wtl_deseq_df, title = b_title)

c <- plot_volcano_plots_deseq(stat5mt_cd8_cd48_healthy_all_deseq_df, title = c_title)
d <- plot_volcano_plots_deseq(stat5mt_cd8_cd48_lgl_wt_deseq_df,  title = d_title)

e <- plot_volcano_plots_deseq(stat5mt_cd8_healthy_cd8_deseq_df, title = e_title)
f <- plot_volcano_plots_deseq(stat5mt_cd8_lgl_wt_deseq_df, title = f_title)

g <- plot_volcano_plots_deseq(healthy_cd8_healthy_cd4_deseq_df, title = g_title)

pdf("stat5b_results/volcanoplots_deseq_local_meta.pdf", width = 18, height = 38)
grid.arrange(a,b,c,d,e,f,g,ncol=1)
dev.off()



healthy_cd8_healthy_cd4_deseq_df
ggplot(df_temp, aes(x = log2FoldChange, y = -log10(pval), color = sigf)) + geom_point()
ggplot(healthy_cd8_healthy_cd4_deseq_df, aes(x = log2FoldChange, y = -log10(pval), color = sigf)) + geom_point()




# Save individual plots
grob_list <- list(a,b,c,d,e,f,g)
for(i in 1:7){
  grob_to_pdf(grob_list[[i]], title = "Volcanoplot", name = name_list[[i]], folder = "Volcanoplots_deseq/", width = 18, height = 12)
}




## Look at the differences between edgeR and DEseq


## Load edgeR results
stat5mt_all_healthy_all_df <- read.delim(paste0("stat5b_results/DE_tables/DE_", a_name, ".txt"))
stat5mt_all_lgl_wtl_df    <- read.delim(paste0("stat5b_results/DE_tables/DE_", b_name, ".txt"))

stat5mt_cd8_cd48_healthy_all_df <-read.delim(paste0("stat5b_results/DE_tables/DE_", c_name, ".txt"))
stat5mt_cd8_cd48_lgl_wt_df <- read.delim(paste0("stat5b_results/DE_tables/DE_", d_name, ".txt"))

stat5mt_cd8_healthy_cd8_df <- read.delim(paste0("stat5b_results/DE_tables/DE_", e_name, ".txt"))
stat5mt_cd8_lgl_wt_df <- read.delim(paste0("stat5b_results/DE_tables/DE_", f_name, ".txt"))

healthy_cd8_healthy_cd4_df <- read.delim(paste0("stat5b_results/DE_tables/DE_", g_name, ".txt"))


healthy_cd8_healthy_cd4_df

draw_venn <- function(genelist_a, genelist_b){
  
  require(VennDiagram)
  
  genelist_a = subset(genelist_a, sigf_cmn == "Significant")
  genelist_b = subset(genelist_b, sigf     == "Significant")
  
  area1 = nrow(genelist_a)
  area2 = nrow(genelist_b)
  cross.area = length(intersect(genelist_a$external_gene_name, genelist_b$external_gene_name))
  
  venn.plot <- draw.pairwise.venn(area1, area2, cross.area, fill = c("salmon", "dodgerblue"),  col = c("salmon", "dodgerblue"), alpha = c(0.1, 0.1), cat.pos = c(9, 3), category = c("edgeR, cmn", "DESeq"), main = "Test")
  return(venn.plot)
}

a <- draw_venn(stat5mt_all_healthy_all_df, stat5mt_all_healthy_all_deseq_df)
b <- draw_venn(stat5mt_all_lgl_wtl_df, stat5mt_all_lgl_wtl_deseq_df)
c <- draw_venn(stat5mt_cd8_cd48_healthy_all_df, stat5mt_cd8_cd48_healthy_all_deseq_df)
d <- draw_venn(stat5mt_cd8_cd48_lgl_wt_df, stat5mt_cd8_cd48_lgl_wt_deseq_df)

e <- draw_venn(stat5mt_cd8_healthy_cd8_df, stat5mt_cd8_healthy_cd8_deseq_df)
f <- draw_venn(stat5mt_cd8_lgl_wt_df, stat5mt_cd8_lgl_wt_deseq_df)
g <- draw_venn(healthy_cd8_healthy_cd4_df, healthy_cd8_healthy_cd4_deseq_df)


grob_list <- list(a,b,c,d,e,f,g)
for(i in 1:7){
  grob_to_pdf(grid.draw(grob_list[[i]]), title = "Venn_diagram", name = name_list[[i]], folder = "Venn_diagrams_edgeR_deseq/", width = 9, height = 9)
}






                   
                   

# ====== Enrichment analysis

library(GSEABase)
library(grid)


enrich_gsea_deseq <- function(df_temp){
  
  
  # df_temp = stat5mt_cd8_lgl_wt_deseq_df
  # hallmark <- read.gmt("h.all.v6.2.symbols.gmt")
  
  percentage <- length(intersect(df_temp$external_gene_name, hallmark$gene)) / nrow(hallmark)
  print(paste(round(percentage, 3) * 100, "% of HALLMARK genes found in df..."))
  
  ## Remove genes from Hallmark that are not found in df
  hallmark_temp <- hallmark[hallmark$gene %in% intersect(hallmark$gene, df_temp$external_gene_name), ]
  # df_temp       <- df_temp[df_temp$external_gene_name %in% intersect(hallmark$gene, df_temp$external_gene_name), ]
  
  nrow(hallmark_temp)
  nrow(df_temp)
  
  # GSEA
  df                <- df_temp[order(df_temp$log2FoldChange, decreasing = T), ]
  deseq_gsea        <- df$log2FoldChange
  names(deseq_gsea) <- df$external_gene_name
  deseq_gsea        <- deseq_gsea[!duplicated(names(deseq_gsea))]
  
  gsea        <- GSEA(deseq_gsea[1:length(deseq_gsea)], nPerm = 100, TERM2GENE = hallmark, verbose = FALSE, pvalueCutoff = 0.05)

  
  gsea        <- GSEA(deseq_gsea[1:10000], nPerm = 10, TERM2GENE = hallmark, verbose = FALSE, pvalueCutoff = 0.05)
  
  results     <- do.call(rbind, gsea@result)
  results     <- as.data.frame(t(results))

  gseaplot(gsea, geneSetID = "HALLMARK_APOPTOSIS")

  
  return(results)
  
  
} 

enrich_hypergeometric <- function(df_temp, dir_temp){
  
  # df_temp = stat5mt_cd8_lgl_wt_df
  # dir_temp = "Down"
  
  # Enrichment as in hypergeometric test
  df_cmn      <- subset(df_temp, sigf_cmn == "Significant" & direction_cmn == dir_temp)
  sigf_cmn    <- df_cmn[order(df_cmn$qval_cmn, decreasing = T), ]$external_gene_name
  enrich_cmn  <- enricher(sigf_cmn, universe = universe_df$external_gene_name, TERM2GENE = hallmark)
  results_cmn <- do.call(rbind, enrich_cmn@result)
  results_cmn <- as.data.frame(t(results_cmn))
  
  df_tgw      <- subset(df_temp, sigf_tgw == "Significant" & direction_tgw == dir_temp)
  sigf_tgw    <- df_tgw[order(df_tgw$qval_tgw, decreasing = T), ]$external_gene_name
  enrich_tgw  <- enricher(sigf_tgw, universe = universe_df$external_gene_name, TERM2GENE = hallmark)
  results_tgw <- do.call(rbind, enrich_tgw@result)
  results_tgw <- as.data.frame(t(results_tgw))
  
  results_cmn$dispersion_type <- "Common"
  results_tgw$dispersion_type <- "Tagwise"
  results                     <- rbind(results_cmn, results_tgw)
  
  colnames(results)
  results[,c(5:7, 9)] <- sapply(results[,c(5:7, 9)], function(x) {as.numeric(as.character(x))})
  return(results)
  
}

plot_enrichment <- function(enrichment_df, title = ""){
  
  # enrichment_df = stat5mt_cd8_healthy_cd8_enrich 
  
  enrichment_df$ID <- substr(enrichment_df$ID, 10, nchar(as.character(enrichment_df$ID)))
  interferon_df    <- enrichment_df[grep("*INTERFERON*", as.character(enrichment_df$ID)), ]
  enrichment_df    <- subset(enrichment_df, qvalue < 0.05)
  
  a <- ggplot(enrichment_df, aes(x = ID, y = Count, fill = dispersion_type)) + geom_bar(stat = "identity", position = "dodge") + coord_flip() + labs(x = "", title = "Significantly (qvalue < 0.05) enriched HALLMARK genesets", fill = "Dispersion") + theme(legend.position = "none")
  b <- ggplot(interferon_df, aes(x = ID, y = qvalue, fill = dispersion_type)) + geom_bar(stat = "identity", position = "dodge") + coord_flip()  + labs(x = "", title = "Interferon HALLMARK genesets", fill = "Dispersion") + geom_hline(yintercept = 0.05, linetype = "dotted")
  
  return(grid.arrange(a,b,ncol=1, top = textGrob(title,gp=gpar(fontsize=20,font=3))))
  
}


hallmark <- read.gmt("h.all.v6.2.symbols.gmt")


## Calculate enrichment as df:s
# DOWN in stat5b
stat5mt_all_healthy_all_enrich_dn      <- enrich_hypergeometric(stat5mt_all_healthy_all_df, dir_temp = "Down")
stat5mt_all_lgl_wtl_enrich_dn          <- enrich_hypergeometric(stat5mt_all_lgl_wtl_df, dir_temp = "Down")

stat5mt_cd8_cd48_healthy_all_enrich_dn <- enrich_hypergeometric(stat5mt_cd8_cd48_healthy_all_df, dir_temp = "Down")
stat5mt_cd8_cd48_lgl_wt_enrich_dn      <- enrich_hypergeometric(stat5mt_cd8_cd48_lgl_wt_df, dir_temp = "Down")

stat5mt_cd8_healthy_cd8_enrich_dn      <- enrich_hypergeometric(stat5mt_cd8_healthy_cd8_df, dir_temp = "Down")
stat5mt_cd8_lgl_wt_enrich_dn           <- enrich_hypergeometric(stat5mt_cd8_lgl_wt_df, dir_temp = "Down")

healthy_cd8_healthy_cd4_enrich_dn      <- enrich_hypergeometric(healthy_cd8_healthy_cd4_df, dir_temp = "Down")



enrich_gsea(subset(stat5mt_all_healthy_all_df, direction_cmn == "Down"))
enrich_gsea(subset(stat5mt_all_lgl_wtl_df, direction_cmn == "Down"))

enrich_gsea(subset(stat5mt_cd8_cd48_healthy_all_df, direction_cmn == "Down"))
enrich_gsea(subset(stat5mt_cd8_cd48_lgl_wt_df, direction_cmn == "Down"))

enrich_gsea(subset(stat5mt_cd8_healthy_cd8_df, direction_cmn == "Down"))
enrich_gsea(subset(stat5mt_cd8_lgl_wt_df, direction_cmn == "Down"))

enrich_gsea(subset(healthy_cd8_healthy_cd4_df, direction_cmn == "Down"))





## Write down
write.table(stat5mt_all_healthy_all_enrich_dn, paste0("stat5b_results/Enrichment/Down_enrich", a_name, ".txt"), sep = "\t")
write.table(stat5mt_all_lgl_wtl_enrich_dn, paste0("stat5b_results/Enrichment/Down_enrich", b_name, ".txt"), sep = "\t")

write.table(stat5mt_cd8_cd48_healthy_all_enrich_dn, paste0("stat5b_results/Enrichment/Down_enrich", c_name, ".txt"), sep = "\t")
write.table(stat5mt_cd8_cd48_lgl_wt_enrich_dn, paste0("stat5b_results/Enrichment/Down_enrich", d_name, ".txt"), sep = "\t")

write.table(stat5mt_cd8_healthy_cd8_enrich_dn, paste0("stat5b_results/Enrichment/Down_enrich", e_name, ".txt"), sep = "\t")
write.table(stat5mt_cd8_lgl_wt_enrich_dn, paste0("stat5b_results/Enrichment/Down_enrich", f_name, ".txt"), sep = "\t")

write.table(healthy_cd8_healthy_cd4_enrich_dn, paste0("stat5b_results/Enrichment/Down_enrich", g_name, ".txt"), sep = "\t")


