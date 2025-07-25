library(ComplexHeatmap)
library(stringr)
library(data.table)
library(reshape2)
library(circlize)

source("PATHS.R")
source("Scripts/Oncoprint_options_v2.R")

all_target_genes = T
sort.oncoprint.byTC = F

ai <- fread(paste0("Analyses/", run_name, "/CN_SNVs_calls.csv"))

gl <- fread("AnnotationFiles/GeneList_PCF-SELECT_Annotation_Sept2020.tsv")
target_genes <- gl[PCF_SELECT == "yes"]$hgnc_symbol

#### Genes to plot
if (all_target_genes) {
  genes <- target_genes
}else{
  genes <- c("TP53", "RB1", "MYCN", "AURKA", "ATM", "ATR","BRCA1", "BRCA2", "FOXA1", "MYC", "CHD1", "SPOP", "PTEN", "AR", "NKX3-1", "ERG_TMPRSS2")
}

dt <- ai[gene %in% genes]
dt[cn.call.corr == "wt"]$cn.call.corr <- ""

correct_order_samples <- unique(dt$sample)

double_genes <- c("ERCC1", "ERCC2", "MSH2", "MSH6", "RNF43", "RAD51C")
mat <- data.frame(matrix(nrow = length(genes), ncol = length(correct_order_samples)))
colnames(mat) <- correct_order_samples
rownames(mat) <- genes
for (s in correct_order_samples) {
  for (g in genes) {
    if (g %in% double_genes) {
      if (grep(g, double_genes) %% 2 == 1) {
        tmp <- dt[sample == s & gene == g][1]
        # tmp_l <- dt_lowTC[sample == s & gene == g][1]
      }else{
        tmp <- dt[sample == s & gene == g][2]
        # tmp_l <- dt_lowTC[sample == s & gene == g][2]
      }
    }else{
      tmp <- dt[sample == s & gene == g]
      # tmp_l <- dt_lowTC[sample == s & gene == g]
    }
    is.snv <- tmp$non_syn_snv != 0
    i <- which(rownames(mat) == g)
    j <- which(colnames(mat) == s)
    if(is.snv){
      if (!is.na(tmp$cn.call.corr)) {
        mat[i,j] <- paste(tmp$cn.call.corr, "SNV", sep = ";")
        # mat_lowTC[i,j] <- paste(tmp_l$cn.call.corr, "SNV", sep = ";")
      }else{
        mat[i,j] <- "SNV"
        # mat_lowTC[i,j] <- "SNV"
      }
    }else{
      mat[i,j] <- tmp$cn.call.corr
      # mat_lowTC[i,j] <- tmp_l$cn.call.corr
    }
  }
}

s <- correct_order_samples

tc.anno <- dt$tc
n <- dt$sample
tc.anno <- tc.anno[!duplicated(n)]
names(tc.anno) <- n[!duplicated(n)]
tc.anno <- tc.anno[which(names(tc.anno) %in% correct_order_samples)]
tc.anno <- tc.anno[match(correct_order_samples, names(tc.anno))]

split <- duplicated(sort(str_split(correct_order_samples, "\\.", simplify = T)[,1]))
split_samples <- c()
n <- 0
for (i in 1:length(split)) {
  if (split[i]) {
    n <- n
  }else{
    n <- n + 1
  }
  split_samples <- c(split_samples, n)
}

### Put ERG_TMPRSS2 in TS
gl[hgnc_symbol == "ERG_TMPRSS2"]$`Oncogene/TS` <- "TS"
row_order <- c(which(sort(rownames(mat)) %in% gl[`Oncogene/TS` == "oncogene"]$hgnc_symbol), #oncogene
               which(sort(rownames(mat)) %in% gl[`Oncogene/TS` == "TS"]$hgnc_symbol), #TS
               which(sort(rownames(mat)) %in% gl[!`Oncogene/TS` %in% c("oncogene", "TS")]$hgnc_symbol)
)

if (sort.oncoprint.byTC) {
  ord_byTC <- order(tc.anno, decreasing = T)
}else{
  ord_byTC <- 1:length(tc.anno)
}
tc.anno <- tc.anno[ord_byTC]
# sample_anno <- sample_anno[ord_byTC]
mat <- mat[,ord_byTC]
# mat_lowTC <- mat_lowTC[,ord_byTC]

select_cols <- which(tc.anno > thr.tc & !is.na(tc.anno))
tc.anno1 <- tc.anno[select_cols]
# sample_anno1 <- sample_anno[select_cols]
mat1 <- mat[,select_cols]

select_cols <- setdiff(1:ncol(mat), select_cols)
tc.anno2 <- tc.anno[select_cols]
# sample_anno2 <- sample_anno[select_cols]
mat2 <- mat[,select_cols]

### Masking of Gains among TS and Losses among oncogenes
# for (gene in rownames(mat2)) {
#   if (gene %in% gl[`Oncogene/TS` == "oncogene"]$hgnc_symbol) {
#     mat2[which(rownames(mat2) == gene),] <- gsub("Deletion", "", mat2[which(rownames(mat2) == gene),])
#   } else if (gene %in% gl[`Oncogene/TS` == "TS"]$hgnc_symbol) {
#     mat2[which(rownames(mat2) == gene),] <- ifelse(mat2[which(rownames(mat2) == gene),] != "DelGain",
#                                                    gsub("Gain", "", mat2[which(rownames(mat2) == gene),]),
#                                                    mat2[which(rownames(mat2) == gene),])
#   }
# }

column_title = "OncoPrint for samples with TC > 15%"
onco_byTC_1 <- oncoPrint(mat1,
                       alter_fun = alter_fun, col = col, show_column_names = T, 
                       column_title = column_title, heatmap_legend_param = heatmap_legend_param, 
                       column_order = 1:ncol(mat1), #column_split = split_samples,
                       row_order = row_order,
                       bottom_annotation = HeatmapAnnotation("Tumor Content" = tc.anno1, 
                                                             # "Samples from" = sample_anno1, 
                                                             col = list("Tumor Content"= colorRamp2(c(0, 1), c("white", "firebrick4"))
                                                                        # , "Samples from"=c("Cornell" = pal_jco(alpha = 0.8)(3)[2], "Vancouver" = pal_jco(alpha = 0.8)(3)[1], "UCL" = pal_jco(alpha = 0.8)(3)[3])
                                                                        ))
)

column_title = "OncoPrint for samples with TC <= 15%"
onco_byTC_2 <- oncoPrint(mat2,
                         alter_fun = alter_fun_lowTC, col = col_lowTC, show_column_names = T, 
                         column_title = column_title, heatmap_legend_param = heatmap_legend_param_lowTC, 
                         column_order = 1:ncol(mat2), #column_split = split_samples,
                         row_order = row_order,
                         bottom_annotation = HeatmapAnnotation("Tumor Content" = tc.anno2, 
                                                               # "Samples from" = sample_anno2, 
                                                               col = list("Tumor Content"= colorRamp2(c(0, 1), c("white", "firebrick4"))
                                                                          # , "Samples from"=c("Cornell" = pal_jco(alpha = 0.8)(3)[2], "Vancouver" = pal_jco(alpha = 0.8)(3)[1], "UCL" = pal_jco(alpha = 0.8)(3)[3])
                                                                          ))
)


pdf(paste0("Analyses/", run_name, "/Plots/Oncoprint_highTC.pdf"), width = 13, height = 14)
onco_byTC_1
dev.off()

pdf(paste0("Analyses/", run_name, "/Plots/Oncoprint_lowTC.pdf"), width = 13, height = 14)
onco_byTC_2
dev.off()
