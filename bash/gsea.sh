#!/bin/bash



# RNKFILE=$RNKPATH/Healthy_CD8_vs_healthy_CD4.rnk # rnkfile.rnk



RNKPATH=/Users/hru/Documents/Laaketieteen_tohtori/LGL_RNAseq/stat5b_results/rnk/
GMT=/Users/hru/Documents/Laaketieteen_tohtori/Applications/h.all.v6.2.symbols.gmt # path to gmt file
OUTDIR=/Users/hru/Documents/Laaketieteen_tohtori/LGL_rnaseq/stat5b_results/GSEA #outdirectory



FILES=$(ls $RNKPATH)
    


for RNKNAME in $FILES

    do

    LABEL=$RNKNAME # prefix-label for subfolder of $OUTDIR
    RNKFILE=$RNKPATH/$RNKNAME

    java -cp /Users/hru/Documents/Laaketieteen_tohtori/Applications/gsea-3.0.jar \
        -Xmx5g xtools.gsea.GseaPreranked \
        -rpt_label $LABEL \
        -rnk $RNKFILE \
        -gmx $GMT \
        -out $OUTDIR \
        -plot_top_x 250 \
        -collapse false \
        -mode Max_probe \
        -norm meandiv \
        -scoring_scheme weighted \
        -include_only_symbols true \
        -make_sets true \
        -rnd_seed 149 \
        -zip_report false \
        -gui false \
        -nperm 1000 \
        -set_min 5 \
        -set_max 500
    
    done


# Healthy_CD8_vs_healthy_CD4.rnk			STAT5mt_CD8_vs_STAT5wt_LGL.rnk			STAT5mt_all_vs_healthy_all.rnk
# STAT5mt_CD8_and_CD4CD8_vs_healthy_CD8.rnk	STAT5mt_CD8_vs_healthy_CD8.rnk
# STAT5mt_CD8_and_CD4D8_vs_healthy_CD8.rnk	STAT5mt_all_vs_STAT5wt_LGL.rnk