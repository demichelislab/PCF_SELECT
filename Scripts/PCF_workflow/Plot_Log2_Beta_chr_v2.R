library(data.table)
#library(plotrix)
library(ggrepel)
library(dplyr)

source("PATHS.R")

load(paste0("Analyses/", run_name, "/betaTable.RData"))
beta.table <- ConvertBetaTable(beta.table)
beta.table <- beta.table[,c("sample", "gene", "beta", "evidence", "evidence.n")]
seg <- fread(paste0("Analyses/", run_name, "/Segmentation_focal.seg"))

data <- merge(beta.table, seg, by = c("sample", "gene"), all = T)

data$ID <- data$sample
data <- data[order(data$ID)]

tc <- fread(paste0("Analyses/", run_name, "/tc_estimations_CLONETv2.tsv"))

asP <- compute_asP(fread(paste0("Analyses/", run_name, "/ai_log2_table.csv")))

gl <- fread(path_annotations)
control.genes <- gl[PCF_SELECT != "yes"]$hgnc_symbol
target.genes <- gl[PCF_SELECT == "yes"]$hgnc_symbol

dir.create(paste0("Analyses/", run_name, "/Plots"))
pdf(paste0("Analyses/", run_name, "/Plots/Log2_Beta_labels_chr_v2.pdf"), width = 12, height = 12)

for (id in unique(data$ID)) {
  x <- data[ID == id,]
  x$col <- chr.cols[match(x$chr, names(chr.cols))]
  x$imb <- ifelse(x$evidence < thr.imbT | is.na(x$evidence), "No", "Yes")
  x$gene.class <- ifelse(x$gene %in% control.genes, "control", x$chr)
  
  for.labels <- x[imb == "Yes" & gene.class != "control"] %>% dplyr::arrange(log2ratio)
  left.genes <- for.labels[1:round(nrow(for.labels)/2),]$gene
  right.genes <- setdiff(for.labels$gene, left.genes)
  
  log2_shift <- -log2(tc[sample == id]$ploidy/2)
  
  print((ggplot(x, aes(x = log2ratio, y = beta, color = chr, fill = gene.class, shape = imb, label = gene))
    + geom_point(size = 2)
    + geom_vline(xintercept = 0+log2_shift, lty = "dashed", col = "red")
    + scale_shape_manual(name = "Allelic Imbalance", values = c("Yes" = 24, "No" = 21))
    + scale_color_manual(name = "Chromosome", values = chr.cols)
    + scale_fill_manual(values = c(chr.cols, "control" = "white"))
    + xlim(c(-2,2))
    + ylim(c(0,1))
    + coord_fixed(ratio = 4)
    + geom_text_repel(data = x[gene %in% left.genes],
                      # colour = "black",
                      xlim = c(-2, -1.5),
                      ylim = c(1, 0),
                      direction = "y",
                      )
    + geom_text_repel(data = x[gene %in% right.genes],
                      # colour = "black",
                      xlim = c(1.5, 2),
                      ylim = c(1, 0),
                      direction = "y",
                      )
    + ylab("Beta")
    + xlab("Log2(R)")
    + labs(title = id,
	   caption = paste0("TC = ", tc[sample == id]$tc, "\t\tasP = ", asP[sample == id]$asP))
    + annotate("text",
               x = -2,
               y = 0.1,
               adj = 0,
               label = paste0("AR Log2R:\n",
                              "Uncorrected = ", round(data[gene == "AR" & sample == id]$log2ratio, digits = 3), "\n",
                              "Corrected = ", round(correct_log2(data[gene == "AR" & sample == id]$log2ratio, ploidy = tc[sample == id]$ploidy, tc = tc[sample == id]$tc_CLONETv2), digits = 3)),
               size = 15/.pt
               )
    + theme_bw(base_size = 15)
    + theme(plot.caption = element_text(size = 15, hjust = 0.5),
            axis.text = element_text(color = "black"))
    + guides(fill = F)
    ))
}
dev.off()

pdf(paste0("Analyses/", run_name, "/Plots/Log2_Beta_chr_v2.pdf"), width = 12, height = 12)

for (id in unique(data$ID)) {
  x <- data[ID == id,]
  x$col <- chr.cols[match(x$chr, names(chr.cols))]
  x$imb <- ifelse(x$evidence < thr.imbT | is.na(x$evidence), "No", "Yes")
  x$gene.class <- ifelse(x$gene %in% control.genes, "control", x$chr)
  
  # for.labels <- x[imb == "Yes" & gene.class != "control"] %>% dplyr::arrange(log2ratio)
  # left.genes <- for.labels[1:round(nrow(for.labels)/2),]$gene
  # right.genes <- setdiff(for.labels$gene, left.genes)
  
  log2_shift <- -log2(tc[sample == id]$ploidy/2)
  
  print((ggplot(x, aes(x = log2ratio, y = beta, color = chr, fill = gene.class, shape = imb, label = gene))
         + geom_point(size = 2)
         + geom_vline(xintercept = 0+log2_shift, lty = "dashed", col = "red")
         + scale_shape_manual(name = "Allelic Imbalance", values = c("Yes" = 24, "No" = 21))
         + scale_color_manual(name = "Chromosome", values = chr.cols)
         + scale_fill_manual(values = c(chr.cols, "control" = "white"))
         + xlim(c(-2,2))
         + ylim(c(0,1))
         + coord_fixed(ratio = 4)
         + ylab("Beta")
         + xlab("Log2(R)")
         + labs(title = id,
                caption = paste0("TC = ", tc[sample == id]$tc, "\t\tasP = ", asP[sample == id]$asP))
         + annotate("text",
                    x = -2,
                    y = 0.1,
                    adj = 0,
                    label = paste0("AR Log2R:\n",
                                   "Uncorrected = ", round(data[gene == "AR" & sample == id]$log2ratio, digits = 3), "\n",
                                   "Corrected = ", round(correct_log2(data[gene == "AR" & sample == id]$log2ratio, ploidy = tc[sample == id]$ploidy, tc = tc[sample == id]$tc_CLONETv2), digits = 3)),
                    size = 15/.pt
         )
         + theme_bw(base_size = 15)
         + theme(plot.caption = element_text(size = 15, hjust = 0.5),
                 axis.text = element_text(color = "black"))
         + guides(fill = F)
  ))
}
dev.off()
