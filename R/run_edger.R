
# Datamanipulation and plotting
library(gridExtra)
library(ggplot2)
library(reshape2)
library(grid)


# Analyses
library(edgeR)
library(biomaRt)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DESeq)


theme_set(theme_bw())





## ====== SETTING UP
me     <- system("whoami", intern = T)
folder <- paste0("/Users/", me, "Dropbox/lgl_stat5b/")
setwd(folder)

## Slight preprocssing of data
lgl_counts           <- read.delim("raw_counts_lgl.txt")
rownames(lgl_counts) <- lgl_counts$X
lgl_counts           <- lgl_counts[,-1]
all_genes_ens        <- rownames(counts)

universe_df <- get_bm(all_genes_ens)


## For patient C_HRUH1254, there's 4 copies of RNAseq transcriptomics.
# Decided to sum these columns up
lgl_counts$C_HRUH_1254 <- rowSums(lgl_counts[,c(15:18)])
lgl_counts <- lgl_counts[,-c(16:18)]

## Get the patient data
patients <- read.csv("LGL_rnaseq_overview.csv")
patients <- patients[match(colnames(lgl_counts), patients$ID), ]
head(patients)
patients$raw_id <- patients$ID

## Pseudonymise
patients$ID <- paste0(as.character(patients$Celltype), "_", as.character(patients$Status), "_", 1:nrow(patients))
colnames(lgl_counts) <- patients$ID

## Add STAT5b-status
patients$stat5b_status <- ifelse(patients$STAT5b == "WT", "WT", "MT")
table(patients$stat5b_status)




stat5mt_cd8_lgl_wt_cds
stat5mt_all_lgl_wtl_cds
stat5mt_cd8_lgl_wt_cds




## Add different study designs

# 1) STAT5bmt all vs healthy
# 2) STAT5bmt all vs STAT5bwt LGL

# 3) STAT5bmt CD8 and CD4/CD8 vs healthy CD8
# 4) STAT5bmt CD8 and CD4/CD8 vs STAT5bwt LGL

# 5) STAT5bmt CD8 vs healthy CD8
# 6) STAT5bmt CD8 vs STAT5bwt LGL

# and bonus;
# 7) Healthy CD8 vs healthy CD4

stat5mt_all    <- which(patients$STAT5b != "WT")
stat5mt_cd8    <- which(patients$STAT5b != "WT" & patients$Celltype == "CD8")
stat5mt_cd4    <- which(patients$STAT5b != "WT" & patients$Celltype == "CD4")
stat5mt_cd4cd8 <- which(patients$STAT5b != "WT" & patients$Celltype == "CD4CD8")

healthy_all    <- which(patients$Status == "Healthy")
healthy_cd8    <- which(patients$Status == "Healthy" & patients$Celltype == "CD8")
healthy_cd4    <- which(patients$Status == "Healthy" & patients$Celltype == "CD4")

lgl_stat5b_wt  <-  which(patients$STAT5b == "WT" & patients$Status == "T-LGL" & patients$Celltype == "CD8")


# Graphic titles
a_title = "STAT5mt all vs healthy all"
b_title = "STAT5mt all vs STAT5wt LGL"

c_title = "STAT5mt CD8 and CD4/CD8 vs healthy CD8"
d_title = "STAT5mt CD8 and CD4/CD8 vs STAT5wt LGL"

e_title = "STAT5mt CD8 vs healthy CD8"
f_title = "STAT5mt CD8 vs STAT5wt LGL"

g_title = "Healthy CD8 vs healthy CD4"


a_name <- paste(c(strsplit(a_title, " ")[[1]]), collapse = "_")
b_name <- paste(c(strsplit(b_title, " ")[[1]]), collapse = "_")
c_name <- paste(c(strsplit(c_title, " ")[[1]]), collapse = "_")
d_name <- paste(c(strsplit(d_title, " ")[[1]]), collapse = "_")
e_name <- paste(c(strsplit(e_title, " ")[[1]]), collapse = "_")
f_name <- paste(c(strsplit(f_title, " ")[[1]]), collapse = "_")
g_name <- paste(c(strsplit(g_title, " ")[[1]]), collapse = "_")

c_name <- paste0(substr(c_name, 1, 19), substr(c_name, 21, nchar(c_name)))
d_name <- paste0(substr(c_name, 1, 19), substr(c_name, 21, nchar(c_name)))

name_list <- list(a_name, b_name, c_name, d_name, e_name, f_name, g_name)

grob_to_pdf <- function(grob_temp, title, name = "", folder = "", width = 9, height = 9){

  ggsave(grob_temp, file = paste0("stat5b_results/", folder, title, "_", name, ".pdf"),  width = width, height = height)

}






## Build diagnostic PCA:s
calc_and_plot_pca <- function(group1, group2, title = ""){

  # group1 = stat5mt_all
  # group2 = healthy_all

  pca_temp     <- prcomp(t(lgl_counts[,c(group1, group2)]))
  patient_temp <- patients[c(group1, group2),]

  a <- ggplot(as.data.frame(pca_temp$x), aes(x = PC1, y = PC2, label = patient_temp$ID, color = patient_temp$STAT5b)) +
    geom_point(size = 3) + labs(color = "STAT5b", title = "STAT5b") + geom_text(nudge_y = 100, size = 3)

  b <- ggplot(as.data.frame(pca_temp$x), aes(x = PC1, y = PC2, label = patient_temp$ID, color = patient_temp$Celltype)) +
    geom_point(size = 3) + labs(color = "Celltype", title = "Celltype") + geom_text(nudge_y = 101, size = 3)

  return(grid.arrange(a,b,top = textGrob(title,gp=gpar(fontsize=20,font=3)), ncol = 2))

}


a <- calc_and_plot_pca(stat5mt_all, healthy_all, title = a_title)
b <- calc_and_plot_pca(stat5mt_all, lgl_stat5b_wt, title = b_title)

c <- calc_and_plot_pca(c(stat5mt_cd8, stat5mt_cd4cd8), healthy_cd8, title = c_title)
d <- calc_and_plot_pca(c(stat5mt_cd8, stat5mt_cd4cd8), lgl_stat5b_wt, title = d_title)

e <- calc_and_plot_pca(stat5mt_cd8, healthy_cd8, title = e_title)
f <- calc_and_plot_pca(stat5mt_cd8, lgl_stat5b_wt, title = f_title)

g <- calc_and_plot_pca(healthy_cd8, healthy_cd4, title = g_title)

pdf("stat5b_results/pca_meta.pdf", width = 16, height = 38)
grid.arrange(a,b,c,d,e,f,g,ncol = 1)
dev.off()


# Save individual plots
grob_list <- list(a,b,c,d,e,f,g)
for(i in 1:7){
  grob_to_pdf(grob_list[[i]], title = "PCA", name = name_list[[i]] ,folder = "PCA/", width = 16, height = 8)
}



# ========== edgeR ===========


build_edgeR_objects <- function(group1, group2, group1_name = "A", group2_name = "B"){

  keep <- c(group1, group2)

  design <- rep(c(group1_name, group2_name), c(length(group1), length(group2)))

  cds_temp <- DGEList(counts = lgl_counts[ ,keep], group = design)

  dim(lgl_counts[ ,keep]) == length(design)

  # Remove lowly expressed genees
  cds_temp <- cds_temp[rowSums(1e+06 * cds_temp$counts/expandAsMatrix(cds_temp$samples$lib.size, dim(cds_temp)) > 1) >= 3, ]

  # Prior: by recommendation or by heuristic
  # prior <- floor(50 / length(keep) - 2)
  prior = 10

  # Count normalisation factors and cmn and tgw dispersions
  cds_temp <- calcNormFactors(cds_temp)
  cds_temp <- estimateCommonDisp(cds_temp)
  cds_temp <- estimateTagwiseDisp(cds_temp, prior.n = prior)

  return(cds_temp)

}

plot_mean_var <- function(cds_temp, title = ""){

  plotMeanVar(cds_temp ,show.raw.vars=TRUE , show.tagwise.vars=TRUE , show.binned.common.disp.vars=FALSE , show.ave.raw.vars=FALSE , dispersion.method = "qcml" , NBline = TRUE , nbins = 100 , #these are arguments about what is plotted
              pch = 16 , xlab ="Mean Expression (Log10 Scale)" , ylab = "Variance (Log10 Scale)" , main = paste0("Mean-Variance Plot \n", title))


}

estimate_de <- function(cds_temp, group1_name = "A", group2_name = "B"){

  # Estimate de-genes
  # group1_name = "STAT5bmt_all"
  # group2_name = "Healthy_all"
  # cds_temp = stat5mt_all_healthy_all_cds

  de.cmn <- exactTest(cds_temp, dispersion = "common",  pair = c(group2_name, group1_name))$table
  de.tgw <- exactTest(cds_temp, dispersion = "tagwise", pair = c(group2_name ,group1_name))$table

  de.cmn$qval <- p.adjust(de.cmn$PValue, method = "BH")
  de.tgw$qval <- p.adjust(de.tgw$PValue, method = "BH")

  de.cmn$sigf <- ifelse(de.cmn$qval < 0.05 & abs(de.cmn$logFC) > 1, "Significant", "Not_significant")
  de.tgw$sigf <- ifelse(de.tgw$qval < 0.05 & abs(de.tgw$logFC) > 1, "Significant", "Not_significant")

  de.cmn$direction <- ifelse(de.cmn$logFC > 0, "Up", "Down")
  de.tgw$direction <- ifelse(de.tgw$logFC > 0, "Up", "Down")

  print("De.cmn")
  print(table(de.cmn$sigf))

  print("De.tgq")
  print(table(de.tgw$sigf))


  colnames(de.cmn) <- paste0(colnames(de.cmn), "_cmn")
  colnames(de.tgw) <- paste0(colnames(de.tgw), "_tgw")

  de_df <- cbind(de.cmn, de.tgw)

  # xprs  <- data.frame(Mean_exprs = rowSums(cds_temp$counts))
  # de_df <- merge(de_df, xprs, by = 0)

  ## Add biomart annotation
  de_bm <- get_bm(rownames(de_df))

#  de_df_full <- merge(de_df, de_bm, by.x = "Row.names", by.y = "ensembl_gene_id")
  de_df_full <- merge(de_df, de_bm, by.x = 0, by.y = "ensembl_gene_id")

  print(paste("Lost information of", nrow(de_df_full) - nrow(de_df), "genes"))

  return(de_df_full)

}


get_bm <- function(genes){

  ensembl    = useMart("ENSEMBL_MART_ENSEMBL")
  ensembl    = useDataset("hsapiens_gene_ensembl", mart = ensembl)
  attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype")

  gene.df <- getBM(attributes = attributes,
                   filters = "ensembl_gene_id",
                   values = genes,
                   mart = ensembl)

}



plot_volcano_plots <- function(de_df, title = ""){

  a <- ggplot(de_df, aes(x = logFC_cmn, y = -log10(PValue_cmn), color = sigf_cmn, size = logCPM_cmn)) + geom_point(alpha = 0.5) + labs(color = "qval < 0.05 and |logFC| > 1", subtitle = "Common") + scale_color_manual(values = c("grey", "dodgerblue")) + theme(legend.position = "top") +
    geom_text(data = subset(de_df, sigf_cmn == "Significant"), aes(x = logFC_cmn, y = -log10(PValue_cmn), label = external_gene_name), nudge_y = 0.5, fontface = "italic")

  b <- ggplot(de_df, aes(x = logFC_tgw, y = -log10(PValue_tgw), color = sigf_tgw, size = logCPM_tgw)) + geom_point(alpha = 0.5) + labs(color = "qval < 0.05 and |logFC| > 1", subtitle = "Tagwise") + scale_color_manual(values = c("grey", "salmon")) + theme(legend.position = "top") +
    geom_text(data = subset(de_df, sigf_tgw == "Significant"), aes(x = logFC_tgw, y = -log10(PValue_tgw), label = external_gene_name), nudge_y = 0.5, fontface = "italic")

  return(grid.arrange(a,b,ncol=2, top = textGrob(title,gp=gpar(fontsize=20,font=3))))

}



# Building edgeR-objects
stat5mt_all_healthy_all_cds      <- build_edgeR_objects(stat5mt_all, healthy_all,   "STAT5bmt_all", "Healthy_all")
stat5mt_all_lgl_wtl_cds          <- build_edgeR_objects(stat5mt_all, lgl_stat5b_wt, "STAT5bmt_all", "STAT5bwt")

stat5mt_cd8_cd48_healthy_all_cds <- build_edgeR_objects(c(stat5mt_cd8, stat5mt_cd4cd8), healthy_cd8, group1_name = "STAT5bmt_cd8_cd48", group2_name = "Healthy_cd8")
stat5mt_cd8_cd48_lgl_wt_cds      <- build_edgeR_objects(c(stat5mt_cd8, lgl_stat5b_wt), healthy_cd8, "STAT5bmt_cd8_cd48", "STAT5bwt")

stat5mt_cd8_healthy_cd8_cds      <- build_edgeR_objects(stat5mt_cd8, healthy_cd8, "STAT5bmt_cd8", "Healthy_cd8")
stat5mt_cd8_lgl_wt_cds           <- build_edgeR_objects(stat5mt_cd8, lgl_stat5b_wt, "STAT5bmt_cd8", "STAT5bwt")

healthy_cd8_healthy_cd4_cds      <- build_edgeR_objects(healthy_cd8, healthy_cd4, "Healthy_cd8", "Healthy_cd4")



# Plot mean-var plots
pdf("stat5b_results/mean_var_plots.pdf", width = 12, height = 20)
par(mfrow=c(5,2))
plot_mean_var(stat5mt_all_healthy_all_cds, a_title)
plot_mean_var(stat5mt_all_lgl_wtl_cds, b_title)

plot_mean_var(stat5mt_cd8_cd48_healthy_all_cds, c_title)
plot_mean_var(stat5mt_cd8_cd48_lgl_wt_cds, d_title)

plot_mean_var(stat5mt_cd8_healthy_cd8_cds, e_title)
plot_mean_var(stat5mt_cd8_lgl_wt_cds, f_title)

plot_mean_var(healthy_cd8_healthy_cd4_cds, g_title)

par(mfrow=c(1,1))
dev.off()



# Estimate DE
stat5mt_all_healthy_all_df      <- estimate_de(stat5mt_all_healthy_all_cds, "STAT5bmt_all", "Healthy_all")
stat5mt_all_lgl_wtl_df          <- estimate_de(stat5mt_all_lgl_wtl_cds, "STAT5bmt_all", "STAT5bwt")

stat5mt_cd8_cd48_healthy_all_df <- estimate_de(stat5mt_cd8_cd48_healthy_all_cds, "STAT5bmt_cd8_cd48", "Healthy_cd8")
stat5mt_cd8_cd48_lgl_wt_df      <- estimate_de(stat5mt_cd8_cd48_lgl_wt_cds, "STAT5bmt_cd8_cd48", "STAT5bwt")

stat5mt_cd8_healthy_cd8_df      <- estimate_de(stat5mt_cd8_healthy_cd8_cds, "STAT5bmt_cd8", "Healthy_cd8")
stat5mt_cd8_lgl_wt_df           <- estimate_de(stat5mt_cd8_lgl_wt_cds, "STAT5bmt_cd8", "STAT5bwt")

healthy_cd8_healthy_cd4_df      <- estimate_de(healthy_cd8_healthy_cd4_cds, "Healthy_cd8", "Healthy_cd4")


write.table(stat5mt_all_healthy_all_df, paste0("stat5b_results/DE_tables/DE_", a_name, ".txt"), sep = "\t")
write.table(stat5mt_all_lgl_wtl_df, paste0("stat5b_results/DE_tables/DE_", b_name, ".txt"), sep = "\t")

write.table(stat5mt_cd8_cd48_healthy_all_df, paste0("stat5b_results/DE_tables/DE_", c_name, ".txt"), sep = "\t")
write.table(stat5mt_cd8_cd48_lgl_wt_df, paste0("stat5b_results/DE_tables/DE_", d_name, ".txt"), sep = "\t")

write.table(stat5mt_cd8_healthy_cd8_df, paste0("stat5b_results/DE_tables/DE_", e_name, ".txt"), sep = "\t")
write.table(stat5mt_cd8_lgl_wt_df, paste0("stat5b_results/DE_tables/DE_", f_name, ".txt"), sep = "\t")

write.table(healthy_cd8_healthy_cd4_df, paste0("stat5b_results/DE_tables/DE_", g_name, ".txt"), sep = "\t")




# Plot volcanoplots of DE-genes
a <- plot_volcano_plots(stat5mt_all_healthy_all_df, title = a_title)
b <- plot_volcano_plots(stat5mt_all_lgl_wtl_df, title = b_title)

c <- plot_volcano_plots(stat5mt_cd8_cd48_healthy_all_df, title = c_title)
d <- plot_volcano_plots(stat5mt_cd8_cd48_lgl_wt_df,  title = d_title)

e <- plot_volcano_plots(stat5mt_cd8_healthy_cd8_df, title = e_title)
f <- plot_volcano_plots(stat5mt_cd8_lgl_wt_df, title = f_title)

g <- plot_volcano_plots(healthy_cd8_healthy_cd4_df, title = g_title)

pdf("stat5b_results/volcanoplots_meta.pdf", width = 18, height = 38)
grid.arrange(a,b,c,d,e,f,g,ncol=1)
dev.off()



# Save individual plots
grob_list <- list(a,b,c,d,e,f,g)
for(i in 1:7){
  grob_to_pdf(grob_list[[i]], title = "Volcanoplot", name = name_list[[i]], folder = "Volcanoplots/", width = 18, height = 12)
}




# ====== Enrichment analysis

library(GSEABase)
library(grid)


enrich_gsea <- function(df_temp){

  percentage <- length(intersect(df_temp$external_gene_name, hallmark$gene)) / nrow(hallmark)
  print(paste(round(percentage, 3), "% of HALLMARK genes found in df..."))

  # GSEA
  df_cmn          <- df_temp[order(df_temp$PValue_cmn, decreasing = T), ]
  cmn_gsea        <- df_cmn$PValue_cmn
  names(cmn_gsea) <- df_cmn$external_gene_name
  cmn_gsea        <- cmn_gsea[!duplicated(names(cmn_gsea))]

  gsea_cmn        <- GSEA(cmn_gsea, TERM2GENE = hallmark, verbose = FALSE, pvalueCutoff = 0.05)
  results_cmn     <- do.call(rbind, gsea_cmn@result)
  results_cmn     <- as.data.frame(t(results_cmn))

  df_tgw          <- df_temp[order(df_temp$PValue_tgw, decreasing = T), ]
  tgw_gsea        <- df_tgw$PValue_tgw
  names(tgw_gsea) <- df_tgw$external_gene_name
  gsea_tgw        <- GSEA(tgw_gsea, TERM2GENE = hallmark, verbose = FALSE)
  results_tgw     <- do.call(rbind, gsea_tgw@result)
  results_tgw     <- as.data.frame(t(results_tgw))

  results <- rbind(results_cmn, results_tgw)

  return(results)


}

enrich_hypergeometric <- function(df_temp, dir_temp){

  # df_temp = stat5mt_cd8_lgl_wt_df
  # dir_temp = "Down"

  # Enrichment as in hypergeometric test
  df_cmn      <- subset(df_temp, sigf_cmn == "Significant" & direction_cmn == dir_temp)
  sigf_cmn    <- df_cmn[order(df_cmn$qval_cmn, decreasing = T), ]$external_gene_name
  enrich_cmn  <- enricher(sigf_cmn, universe = un iverse_df$external_gene_name, TERM2GENE = hallmark)
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




# UP in stat5b
stat5mt_all_healthy_all_enrich_up      <- enrich_hypergeometric(stat5mt_all_healthy_all_df, dir_temp = "Up")
stat5mt_all_lgl_wtl_enrich_up          <- enrich_hypergeometric(stat5mt_all_lgl_wtl_df, dir_temp = "Up")

stat5mt_cd8_cd48_healthy_all_enrich_up <- enrich_hypergeometric(stat5mt_cd8_cd48_healthy_all_df, dir_temp = "Up")
stat5mt_cd8_cd48_lgl_wt_enrich_up      <- enrich_hypergeometric(stat5mt_cd8_cd48_lgl_wt_df, dir_temp = "Up")

stat5mt_cd8_healthy_cd8_enrich_up      <- enrich_hypergeometric(stat5mt_cd8_healthy_cd8_df, dir_temp = "Up")
stat5mt_cd8_lgl_wt_enrich_up           <- enrich_hypergeometric(healthy_cd8_healthy_cd4_df, dir_temp = "Up")

healthy_cd8_healthy_cd4_enrich_up      <- enrich_hypergeometric(healthy_cd8_healthy_cd4_df, dir_temp = "Up")



enrich_gsea(subset(stat5mt_all_healthy_all_df, direction_cmn == "Up"))
df <- enrich_gsea(subset(stat5mt_all_lgl_wtl_df, direction_cmn == "Up"))

enrich_gsea(subset(stat5mt_cd8_cd48_healthy_all_df, direction_cmn == "Up"))
enrich_gsea(subset(stat5mt_cd8_cd48_lgl_wt_df, direction_cmn == "Up"))

enrich_gsea(subset(stat5mt_cd8_healthy_cd8_df, direction_cmn == "Up"))
enrich_gsea(subset(stat5mt_cd8_lgl_wt_df, direction_cmn == "Up"))

enrich_gsea(subset(healthy_cd8_healthy_cd4_df, direction_cmn == "Up"))


## Write down
write.table(stat5mt_all_healthy_all_enrich_up, paste0("stat5b_results/Enrichment/Up_enrich", a_name, ".txt"), sep = "\t")
write.table(stat5mt_all_lgl_wtl_enrich_up, paste0("stat5b_results/Enrichment/Up_enrich", b_name, ".txt"), sep = "\t")

write.table(stat5mt_cd8_cd48_healthy_all_enrich_up, paste0("stat5b_results/Enrichment/Up_enrich", c_name, ".txt"), sep = "\t")
write.table(stat5mt_cd8_cd48_lgl_wt_enrich_up, paste0("stat5b_results/Enrichment/Up_enrich", d_name, ".txt"), sep = "\t")

write.table(stat5mt_cd8_healthy_cd8_enrich_up, paste0("stat5b_results/Enrichment/Up_enrich", e_name, ".txt"), sep = "\t")
write.table(stat5mt_cd8_lgl_wt_enrich_up, paste0("stat5b_results/Enrichment/Up_enrich", f_name, ".txt"), sep = "\t")

write.table(healthy_cd8_healthy_cd4_enrich_up, paste0("stat5b_results/Enrichment/Up_enrich", g_name, ".txt"), sep = "\t")






## Plot enrichments as barplots
a <- plot_enrichment(stat5mt_all_healthy_all_enrich_dn, title = paste("Down in",  a_title))
b <- plot_enrichment(stat5mt_all_lgl_wtl_enrich_dn, title = paste("Down in", b_title))

c <- plot_enrichment(stat5mt_cd8_cd48_healthy_all_enrich_dn, title = paste("Down in", c_title))
d <- plot_enrichment(stat5mt_cd8_cd48_lgl_wt_enrich_dn, title = paste("Down in", d_title))

e <- plot_enrichment(stat5mt_cd8_healthy_cd8_enrich_dn, title = paste("Down in", e_title))
f <- plot_enrichment(stat5mt_cd8_lgl_wt_enrich_dn, title = paste("Down in", f_title))

g <- plot_enrichment(healthy_cd8_healthy_cd4_enrich_dn, title = paste("Down in", g_title))


pdf("stat5b_results/geneset_enrichment_meta_down.pdf", height = 32, width = 8)
grid.arrange(a,b,c,d,e,f,g,ncol=1)
dev.off()


# Save individual plots
grob_list <- list(a,b,c,d,e,f,g)
for(i in 1:7){
  grob_to_pdf(grob_list[[i]], title = "Down_geneset_enrichment", name = name_list[[i]], folder = "Enrichment/", width = 12, height = 12)
}




a <- plot_enrichment(stat5mt_all_healthy_all_enrich_up, title = paste("Up in",  a_title))
b <- plot_enrichment(stat5mt_all_lgl_wtl_enrich_up, title = paste("Up in", b_title))

c <- plot_enrichment(stat5mt_cd8_cd48_healthy_all_enrich_up, title = paste("Up in", c_title))
d <- plot_enrichment(stat5mt_cd8_cd48_lgl_wt_enrich_up, title = paste("Up in", d_title))

e <- plot_enrichment(stat5mt_cd8_healthy_cd8_enrich_up, title = paste("Up in", e_title))
f <- plot_enrichment(stat5mt_cd8_lgl_wt_enrich_up, title = paste("Up in", f_title))

g <- plot_enrichment(healthy_cd8_healthy_cd4_enrich_up, title = paste("Up in", g_title))


pdf("stat5b_results/geneset_enrichment_meta_Up.pdf", height = 32, width = 8)
grid.arrange(a,b,c,d,e,f,ncol=1)
dev.off()


# Save individual plots
grob_list <- list(a,b,c,d,e,f,g)
for(i in 1:7){
  grob_to_pdf(grob_list[[i]], title = "Up_geneset_enrichment", name = name_list[[i]], folder = "Enrichment/", width = 12, height = 12)
}




## => So statistically enriched in the IFNa and/or IFNg are

# 2) STAT5bmt all vs STAT5bwt LGL
# 6) STAT5bmt CD8 vs STAT5bwt LGL

# and in the unfavorable direction;
# 4) STAT5bmt CD8 and CD4/CD8 vs STAT5bwt LGL, where STAT5b mt shows upregulation


## Look at the volcanoplots of just IFNa/IFNg -genes



plot_geneset_volcano_plots <- function(de_df, title = "", geneset){

  de_df <- de_df[de_df$external_gene_name %in% geneset, ]

  a <- ggplot(de_df, aes(x = logFC_cmn, y = -log10(PValue_cmn), color = sigf_cmn, size = logCPM_cmn)) + geom_point(alpha = 0.5) + labs(color = "qval < 0.05 and |logFC| > 1", subtitle = "Common") + scale_color_manual(values = c("grey", "dodgerblue")) + theme(legend.position = "top") +
    geom_text(data = subset(de_df, sigf_cmn == "Significant"), aes(x = logFC_cmn, y = -log10(PValue_cmn), size = 3, label = external_gene_name), nudge_y = 0.2, fontface = "italic")

  if(length(unique(de_df$sigf_tgw)) == 2){
    b <- ggplot(de_df, aes(x = logFC_tgw, y = -log10(PValue_tgw), color = sigf_tgw, size = logCPM_tgw)) + geom_point(alpha = 0.5) + labs(color = "qval < 0.05 and |logFC| > 1", subtitle = "Tagwise") + scale_color_manual(values = c("grey", "salmon")) + theme(legend.position = "top") +
      geom_text(data = subset(de_df, sigf_tgw == "Significant"), aes(x = logFC_tgw, y = -log10(PValue_tgw), size = 3, label = external_gene_name), nudge_y = 0.2, fontface = "italic")
  }

  if(length(unique(de_df$sigf_tgw)) == 1){
    b <- ggplot(de_df, aes(x = logFC_tgw, y = -log10(PValue_tgw), color = sigf_tgw, size = logCPM_tgw)) + geom_point(alpha = 0.5, color = "grey") + labs(color = "qval < 0.05 and |logFC| > 1", subtitle = "Tagwise")  + theme(legend.position = "top")
      # geom_text(data = subset(de_df, sigf_tgw == "Significant"), aes(x = logFC_tgw, y = -log10(PValue_tgw), size = 3, label = external_gene_name), nudge_y = 0.2, fontface = "italic")
  }

  return(grid.arrange(a,b,ncol=2, top = textGrob(title,gp=gpar(fontsize=20,font=3))))

}

## IFNg
a <- plot_geneset_volcano_plots(stat5mt_all_healthy_all_df, geneset = interferon_gamma, title = paste("Interferon gamma response genes:\n", a_title))
b <- plot_geneset_volcano_plots(stat5mt_all_lgl_wtl_df, geneset = interferon_gamma, title = paste("Interferon gamma response genes:\n", b_title))

c <- plot_geneset_volcano_plots(stat5mt_cd8_cd48_healthy_all_df, geneset = interferon_gamma, title = paste("Interferon gamma response genes:\n", c_title))
d <- plot_geneset_volcano_plots(stat5mt_cd8_cd48_lgl_wt_df, geneset = interferon_gamma, title = paste("Interferon gamma response genes:\n", d_title))

e <- plot_geneset_volcano_plots(stat5mt_cd8_healthy_cd8_df, geneset = interferon_gamma, title = paste("Interferon gamma response genes:\n", e_title))
f <- plot_geneset_volcano_plots(stat5mt_cd8_lgl_wt_df,  geneset = interferon_gamma, title = paste("Interferon gamma response genes:\n", f_title))

g <- plot_geneset_volcano_plots(healthy_cd8_healthy_cd4_df,  geneset = interferon_gamma, title = paste("Interferon gamma response genes:\n", g_title))



pdf("stat5b_results/ifng_hallmark_volcanoplot_meta.pdf", height = 38, width = 18)
grid.arrange(a,b,c,d,e,f,g,ncol=1)
dev.off()


# Save individual plots
grob_list <- list(a,b,c,d,e,f,g)
for(i in 1:7){
  grob_to_pdf(grob_list[[i]], title = "Volcano_ifng", name = name_list[[i]], folder = "Volcanoplots/", width = 18, height = 12)
}





## IFNa
a <- plot_geneset_volcano_plots(stat5mt_all_healthy_all_df, geneset = interferon_alpha, title = paste("Interferon alpha response genes:\n", a_title))
b <- plot_geneset_volcano_plots(stat5mt_all_lgl_wtl_df, geneset = interferon_alpha, title = paste("Interferon alpha response genes:\n", b_title))

c <- plot_geneset_volcano_plots(stat5mt_cd8_cd48_healthy_all_df, geneset = interferon_alpha, title = paste("Interferon alpha response genes:\n", c_title))
d <- plot_geneset_volcano_plots(stat5mt_cd8_cd48_lgl_wt_df, geneset = interferon_alpha, title = paste("Interferon alpha response genes:\n", d_title))

e <- plot_geneset_volcano_plots(stat5mt_cd8_healthy_cd8_df, geneset = interferon_alpha, title = paste("Interferon alpha response genes:\n", e_title))
f <- plot_geneset_volcano_plots(stat5mt_cd8_lgl_wt_df,  geneset = interferon_alpha, title = paste("Interferon alpha response genes:\n", f_title))

g <- plot_geneset_volcano_plots(healthy_cd8_healthy_cd4_df,  geneset = interferon_alpha, title = paste("Interferon alpha response genes:\n", g_title))



pdf("stat5b_results/ifna_hallmark_volcanoplot_meta.pdf", height = 38, width = 18)
grid.arrange(a,b,c,d,e,f,g,ncol=1)
dev.off()


# Save individual plots
grob_list <- list(a,b,c,d,e,f,g)
for(i in 1:7){
  grob_to_pdf(grob_list[[i]], title = "Volcano_ifna", name = name_list[[i]], folder = "Volcanoplots/", width = 18, height = 12)
}





## Building a heatmap of interferon gamma and alpha responses
library(ComplexHeatmap)
library(circlize)


interferon_gamma <- hallmark[grep("*INTERFERON_GAMMA*", as.character(hallmark$ont)), ]$gene
interferon_alpha <- hallmark[grep("*INTERFERON_ALPHA*", as.character(hallmark$ont)), ]$gene

temp <- hallmark[grep("*INTERFERON_GAMMA*", as.character(hallmark$ont)), ]

## Count sytotoxicity scoe (Dufva, Pölönen, unpublished)
  cytolytic_genes = c("GZMA", "GZMM", "GZMH", "PRF1", "GNLY")

  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }


  cyto_bm <- universe_df[universe_df$external_gene_name %in% cytolytic_genes, ]
  cyto_df <- lgl_counts[cyto_bm$ensembl_gene_id, ]

  cytoscore <- apply(cyto_df, 2, gm_mean)


  patients$cytoscore <- cytoscore





plot_hm <- function(genes, patient_index, cluster_columns = F, title = ""){

  keep = patient_index

  temp_bm           <- universe_df[universe_df$external_gene_name %in% genes, ]

  temp_df           <- lgl_counts[temp_bm$ensembl_gene_id, ]
  rownames(temp_df) <- temp_bm$external_gene_name


  # temp_df <- temp_df[,keep]
  temp_df <- log(temp_df[,keep] + 1)
  patient_temp <- patients[keep,]

  ## Build annotation
  ha <- HeatmapAnnotation(STAT5b        = patient_temp$stat5b_status,
                          Celltype      = patient_temp$Celltype,
                          # Cytoscore     = anno_barplot(patient_temp$cytoscore),
                          col = list(STAT5b   = c("MT" =  "salmon", "WT" = "darkblue"),
                                     Celltype = c("CD4" = "deepskyblue3", "CD8" = "burlywood1", "CD4CD8" = "darkolivegreen4")),
                          show_annotation_name = T)


  ## Plot heatmap
  Heatmap(scale(temp_df),
          top_annotation = ha,
          # col = colorRamp2(c(-4,0,4), colors = c("dodgerblue", "white", "salmon")),
          cluster_rows = T,
          cluster_columns = cluster_columns,
          name = "Z-score \nof log(counts)",
          column_title = title,
          column_title_gp = gpar(fontsize = 20, font = 3))

}


## === Plot heatmaps

# STAT5mt all vs STAT5wt LGL; all genes
pdf("stat5b_results/STAT5mt_all_vs_STAT5wt_LGL_hallmark_ifna_hm.pdf", height = 24, width = 10)
plot_hm(interferon_alpha, patient_index = c(stat5mt_all, lgl_stat5b_wt), title = paste("IFNa HALLMARK-genes\n", b_title))
dev.off()

pdf("stat5b_results/STAT5mt_all_vs_STAT5wt_LGL_hallmark_ifng_hm.pdf", height = 28, width = 10)
plot_hm(interferon_gamma, patient_index = c(stat5mt_all, lgl_stat5b_wt), title = paste("IFNg HALLMARK-genes\n", b_title))
dev.off()


# STAT5mt CD8 vs STAT5wt LGL; all genes
pdf("stat5b_results/STAT5mt_cd8_vs_STAT5wt_LGL_hallmark_ifna_hm.pdf", height = 24, width = 10)
plot_hm(interferon_alpha, patient_index = c(stat5mt_cd8, lgl_stat5b_wt), title = paste("IFNa HALLMARK-genes\n", f_title))
dev.off()

pdf("stat5b_results/STAT5mt_cd8_vs_STAT5wt_LGL_hallmark_ifng_hm.pdf", height = 28, width = 10)
plot_hm(interferon_gamma, patient_index = c(stat5mt_cd8, lgl_stat5b_wt), title = paste("IFNg HALLMARK-genes\n", f_title))
dev.off()



# STAT5mt CD8 vs healthy; all genes
pdf("stat5b_results/STAT5mt_cd8_vs_healthy_cd8_hallmark_ifna_hm.pdf", height = 24, width = 10)
plot_hm(interferon_alpha, patient_index = c(stat5mt_cd8, healthy_cd8), title = paste("IFNa HALLMARK-genes\n", e_title))
dev.off()

pdf("stat5b_results/STAT5mt_cd8_vs_healthy_cd8_hallmark_ifng_hm.pdf", height = 28, width = 10)
plot_hm(interferon_gamma, patient_index = c(stat5mt_cd8, healthy_cd8), title = paste("IFNg HALLMARK-genes\n", e_title))
dev.off()



# ====== DE


# STAT5mt all vs STAT5wt LGL; DE-genes
subset_df <- function(de_df, geneset, sigf){

  # de_df = stat5mt_all_lgl_wtl_df
  de_df <- de_df[de_df$external_gene_name %in% geneset, ]
  de_df <- subset(de_df, sigf_cmn == "Significant")

  return(de_df)
}

stat5m_all_lgl_infg <- subset_df(stat5mt_all_lgl_wtl_df, interferon_gamma)
stat5m_all_lgl_infa <- subset_df(stat5mt_all_lgl_wtl_df, interferon_alpha)

pdf("stat5b_results/STAT5mt_all_vs_STAT5wt_LGL_DE_hallmark_ifna_hm.pdf", height = 6, width = 10)
plot_hm(genes = stat5m_all_lgl_infa$external_gene_name, patient_index =c(stat5mt_all, lgl_stat5b_wt), cluster_columns = F, title = paste("DE IFNa HALLMARK-genes\n", b_title))
dev.off()

pdf("stat5b_results/STAT5mt_all_vs_STAT5wt_LGL_DE_hallmark_ifng_hm.pdf", height = 12, width = 10)
plot_hm(genes = stat5m_all_lgl_infg$external_gene_name, patient_index =c(stat5mt_all, lgl_stat5b_wt), cluster_columns = F, title = paste("DE IFNg HALLMARK-genes\n", b_title))
dev.off()



# STAT5mt CD8 vs STAT5wt LGL; DE-genes
stat5m_cd8_lgl_infg <- subset_df(stat5mt_cd8_lgl_wt_df, interferon_gamma)
stat5m_cd8_lgl_infa <- subset_df(stat5mt_cd8_lgl_wt_df, interferon_alpha)


pdf("stat5b_results/STAT5mt_cd8_vs_STAT5wt_LGL_DE_hallmark_ifna_hm.pdf", height = 6, width = 8)
plot_hm(genes = stat5m_cd8_lgl_infa$external_gene_name, patient_index =c(stat5mt_cd8, lgl_stat5b_wt), cluster_columns = F, title = paste("DE IFNa HALLMARK-genes\n", f_title))
dev.off()

pdf("stat5b_results/STAT5mt_cd8_vs_STAT5wt_LGL_DE_hallmark_ifng_hm.pdf", height = 12, width = 10)
plot_hm(genes = stat5m_cd8_lgl_infg$external_gene_name, patient_index =c(stat5mt_cd8, lgl_stat5b_wt), cluster_columns = F, title = paste("DE IFNg HALLMARK-genes\n", f_title))
dev.off()




# STAT5mt CD8 vs healthy CD8; DE-genes
stat5m_cd8_hcd8_infg <- subset_df(stat5mt_cd8_healthy_cd8_df, interferon_gamma)
stat5m_cd8_hcd8_infa <- subset_df(stat5mt_cd8_healthy_cd8_df, interferon_alpha)


pdf("stat5b_results/STAT5mt_cd8_vs_healthy_cd8_DE_hallmark_ifna_hm.pdf", height = 6, width = 8)
plot_hm(genes = stat5m_cd8_hcd8_infa$external_gene_name, patient_index =c(stat5mt_cd8, healthy_cd8), cluster_columns = F, title = paste("DE IFNa HALLMARK-genes\n", e_title))
dev.off()

pdf("stat5b_results/STAT5mt_cd8_vs_healthy_cd8_DE_hallmark_ifng_hm.pdf", height = 12, width = 10)
plot_hm(genes = stat5m_cd8_hcd8_infa$external_gene_name, patient_index =c(stat5mt_cd8, healthy_cd8), cluster_columns = F, title = paste("DE IFNa HALLMARK-genes\n", e_title))
dev.off()







# === Boxplots of DE-genes

plot_bp <- function(genes, group1, group2, group1_name, group2_name){

  # group1 = stat5mt_all
  # group2 = lgl_stat5b_wt
  # group1_name = "stat5mt_all"
  # group2_name = "lgl_stat5b_wt"

  # genes = stat5m_cd8_hcd8_infg$external_gene_name; group1 = stat5mt_cd8; group2 = healthy_cd8; group1_name = "STAT5mt_CD8"; group2_name = "Healthy_CD8"
  keep = c(group1, group2)

  # Get the external genenames to ensembl
  temp_bm           <- universe_df[universe_df$external_gene_name %in% genes, ]
  temp_df           <- lgl_counts[temp_bm$ensembl_gene_id, ]
  rownames(temp_df) <- temp_bm$external_gene_name

  # Make a df with counts and patients mathcing request
  temp_df <- log(temp_df[,keep] + 1)
  patient_temp <- patients[keep,]

  # Make a grouping vector
  patient_vector <-  data.frame(grouping = rep(c(group1_name, group2_name), c(length(group1), length(group2))))
  rownames(patient_vector) <- patient_temp$ID


  # Make a melted df for boxplotting, with grouping vector as well
  temp_df$gene <- rownames(temp_df)
  temp_df_melt <- melt(temp_df)
  temp_df_melt <- merge(temp_df_melt, patient_vector, by.x = "variable", by.y = 0)

  a <- ggplot(temp_df_melt, aes(x = gene, y = value, color = grouping, fill = grouping)) + geom_boxplot(alpha = 0.3, position = "dodge", outlier.colour = NA) + labs(y = "log(counts)", x = "", fill = "") + geom_jitter(position = position_dodge(width = 0.5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


  return(a)
}



## MT vs WT LGL; all
a <- plot_bp(genes = stat5m_all_lgl_infa$external_gene_name, group1 = stat5mt_all, group2 = lgl_stat5b_wt, group1_name = "STAT5mt_all", group2_name = "STAT5wt_lgl") + ggtitle(paste("DE INFa HALLMARK-genes\n", b_title))
b <- plot_bp(genes = stat5m_all_lgl_infg$external_gene_name, group1 = stat5mt_all, group2 = lgl_stat5b_wt, group1_name = "STAT5mt_all", group2_name = "STAT5wt_lgl") + ggtitle(paste("DE INFg HALLMARK-genes\n", b_title))


pdf("stat5b_results/STAT5mt_all_vs_stat5wt_lgl_DE_infa_box.pdf", height = 12, width = 10)
a
dev.off()


pdf("stat5b_results/STAT5mt_all_vs_stat5wt_lgl_DE_infg_box.pdf", height = 12, width = 18)
b
dev.off()




## MT vs WT LGL; CD()
a <- plot_bp(genes = stat5m_cd8_lgl_infa$external_gene_name, group1 = stat5mt_cd8, group2 = lgl_stat5b_wt, group1_name = "STAT5mt_cd8", group2_name = "STAT5wt_lgl") + ggtitle(paste("DE INFa Hcd8MARK-genes\n", b_title))
b <- plot_bp(genes = stat5m_cd8_lgl_infg$external_gene_name, group1 = stat5mt_cd8, group2 = lgl_stat5b_wt, group1_name = "STAT5mt_cd8", group2_name = "STAT5wt_lgl") + ggtitle(paste("DE INFg Hcd8MARK-genes\n", b_title))


pdf("stat5b_results/STAT5mt_cd8_vs_stat5wt_lgl_DE_infa_box.pdf", height = 12, width = 10)
a
dev.off()


pdf("stat5b_results/STAT5mt_cd8_vs_stat5wt_lgl_DE_infg_box.pdf", height = 12, width = 18)
b
dev.off()




## Healthy
a <- plot_bp(genes = stat5m_cd8_hcd8_infa$external_gene_name, group1 = stat5mt_cd8, group2 = healthy_cd8, group1_name = "STAT5mt_CD8", group2_name = "Healthy_CD8") + ggtitle(paste("DE INFa HALLMARK-genes\n", e_title))
b <- plot_bp(genes = stat5m_cd8_hcd8_infg$external_gene_name, group1 = stat5mt_cd8, group2 = healthy_cd8, group1_name = "STAT5mt_CD8", group2_name = "Healthy_CD8") + ggtitle(paste("DE INFg HALLMARK-genes\n", e_title))


pdf("stat5b_results/STAT5mt_cd8_vs_healthy_cd8_DE_infa_box.pdf", height = 12, width = 10)
a
dev.off()



pdf("stat5b_results/STAT5mt_cd8_vs_healthy_cd8_DE_infg_box.pdf", height = 12, width = 18)
b
dev.off()






plot_jitter <- function(genes, group1, group2, group1_name, group2_name){

  # group1 = stat5mt_all
  # group2 = lgl_stat5b_wt
  # group1_name = "stat5mt_all"
  # group2_name = "lgl_stat5b_wt"

  # genes = stat5m_cd8_hcd8_infg$external_gene_name; group1 = stat5mt_cd8; group2 = healthy_cd8; group1_name = "STAT5mt_CD8"; group2_name = "Healthy_CD8"
  keep = c(group1, group2)

  # Get the external genenames to ensembl
  temp_bm           <- universe_df[universe_df$external_gene_name %in% genes, ]
  temp_df           <- lgl_counts[temp_bm$ensembl_gene_id, ]
  rownames(temp_df) <- temp_bm$external_gene_name

  # Make a df with counts and patients mathcing request
  temp_df <- log(temp_df[,keep] + 1)
  patient_temp <- patients[keep,]

  # Make a grouping vector
  patient_vector <-  data.frame(grouping = rep(c(group1_name, group2_name), c(length(group1), length(group2))))
  rownames(patient_vector) <- patient_temp$ID


  # Make a melted df for boxplotting, with grouping vector as well
  temp_df$gene <- rownames(temp_df)
  temp_df_melt <- melt(temp_df)
  temp_df_melt <- merge(temp_df_melt, patient_vector, by.x = "variable", by.y = 0)

  a <- ggplot(temp_df_melt, aes(x = variable, y = value, label = gene, fill = grouping)) + labs(y = "log(counts)", x = "", fill = "") + geom_boxplot(alpha = 0.5) +
    # geom_jitter(position = position_dodge(width = 0.5)) +
    # geom_jitter(aes(color = gene)) +
    geom_text(aes(color = gene), size = 3, position = position_jitter()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  a

  return(a)
}






## MT vs WT LGL; all
a <- plot_jitter(genes = stat5m_all_lgl_infa$external_gene_name, group1 = stat5mt_all, group2 = lgl_stat5b_wt, group1_name = "STAT5mt_all", group2_name = "STAT5wt_lgl") + ggtitle(paste("DE INFa HALLMARK-genes\n", b_title))
b <- plot_jitter(genes = stat5m_all_lgl_infg$external_gene_name, group1 = stat5mt_all, group2 = lgl_stat5b_wt, group1_name = "STAT5mt_all", group2_name = "STAT5wt_lgl") + ggtitle(paste("DE INFg HALLMARK-genes\n", b_title))


pdf("stat5b_results/STAT5mt_all_vs_stat5wt_lgl_DE_infa_patient_wise.pdf", height = 12, width = 10)
a
dev.off()


pdf("stat5b_results/STAT5mt_all_vs_stat5wt_lgl_DE_infg_patient_wise.pdf", height = 12, width = 18)
b
dev.off()




## MT vs WT LGL; CD()
a <- plot_jitter(genes = stat5m_cd8_lgl_infa$external_gene_name, group1 = stat5mt_cd8, group2 = lgl_stat5b_wt, group1_name = "STAT5mt_cd8", group2_name = "STAT5wt_lgl") + ggtitle(paste("DE INFa Hcd8MARK-genes\n", b_title))
b <- plot_jitter(genes = stat5m_cd8_lgl_infg$external_gene_name, group1 = stat5mt_cd8, group2 = lgl_stat5b_wt, group1_name = "STAT5mt_cd8", group2_name = "STAT5wt_lgl") + ggtitle(paste("DE INFg Hcd8MARK-genes\n", b_title))


pdf("stat5b_results/STAT5mt_cd8_vs_stat5wt_lgl_DE_infa_patient_wise.pdf", height = 12, width = 10)
a
dev.off()


pdf("stat5b_results/STAT5mt_cd8_vs_stat5wt_lgl_DE_infg_patient_wise.pdf", height = 12, width = 18)
b
dev.off()




## Healthy
a <- plot_jitter(genes = stat5m_cd8_hcd8_infa$external_gene_name, group1 = stat5mt_cd8, group2 = healthy_cd8, group1_name = "STAT5mt_CD8", group2_name = "Healthy_CD8") + ggtitle(paste("DE INFa HALLMARK-genes\n", e_title))
b <- plot_jitter(genes = stat5m_cd8_hcd8_infg$external_gene_name, group1 = stat5mt_cd8, group2 = healthy_cd8, group1_name = "STAT5mt_CD8", group2_name = "Healthy_CD8") + ggtitle(paste("DE INFg HALLMARK-genes\n", e_title))


pdf("stat5b_results/STAT5mt_cd8_vs_healthy_cd8_DE_infa_patient_wise.pdf", height = 12, width = 10)
a
dev.off()



pdf("stat5b_results/STAT5mt_cd8_vs_healthy_cd8_DE_infg_patient_wise.pdf", height = 12, width = 18)
b
dev.off()
