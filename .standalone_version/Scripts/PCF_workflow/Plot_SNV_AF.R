library(data.table)
#library(plotrix)
library(patchwork)
library(ggrepel)
library(dplyr)

source("PATHS.R")

# load(paste0("Analyses/", run_name, "/betaTable.RData"))
# beta.table <- ConvertBetaTable(beta.table)
# beta.table <- beta.table[,c("sample", "gene", "beta", "evidence", "evidence.n")]
# seg <- fread(paste0("Analyses/", run_name, "/Segmentation_focal.seg"))
# 
# data <- merge(beta.table, seg, by = c("sample", "gene"), all = T)
# 
# data$ID <- data$sample
# data <- data[order(data$ID)]

snv <- fread(paste0("Analyses/", run_name, "/SNVs_calls_corrected.csv"))
snv$ID <- snv$sample
## Put to 1 all the SNVs with AF.corr > 1
snv$t_af_corr <- ifelse(snv$t_af_corr > 1, 1, snv$t_af_corr)

tc <- fread(paste0("Analyses/", run_name, "/tc_estimations_CLONETv2.tsv"))

asP <- compute_asP(fread(paste0("Analyses/", run_name, "/ai_log2_table.csv")))

gl <- fread(path_annotations)
control.genes <- gl[PCF_SELECT != "yes"]$hgnc_symbol
target.genes <- gl[PCF_SELECT == "yes"]$hgnc_symbol

dir.create(paste0("Analyses/", run_name, "/Plots"))
pdf(paste0("Analyses/", run_name, "/Plots/SNV_AF.pdf"))

for (id in unique(snv$ID)) {
  x <- snv[ID == id & gene %in% gl$hgnc_symbol]
  x$col <- chr.cols[match(x$Chromosome, names(chr.cols))]
  # x$imb <- ifelse(x$evidence < thr.imbT | is.na(x$evidence), "No", "Yes")
  # x$gene.class <- ifelse(x$gene %in% control.genes, "control", x$chr)
  
  # for.labels <- x[imb == "Yes" & gene.class != "control"] %>% dplyr::arrange(log2ratio)
  # left.genes <- for.labels[1:round(nrow(for.labels)/2),]$gene
  # right.genes <- setdiff(for.labels$gene, left.genes)
  
  # log2_shift <- -log2(tc[sample == id]$ploidy/2)
  
  p.unc <- (ggplot(x, aes(x = af_case, y = id, color = Chromosome, label = gene
                       # fill = gene.class, shape = imb, label = gene
                       ))
         + geom_point(size = 2)
         + scale_color_manual(name = "Chromosome", values = chr.cols)
         # + xlim(c(0,1))
         + scale_x_continuous(name = "Uncorrected AF", breaks = seq(0,1,0.1), limits = c(0,1))
         + coord_fixed(ratio = 0.2)
         + geom_text_repel(data = x,
                           size = 2,
                           # colour = "black",
                           # xlim = c(-2, -1.5),
                           ylim = c(1.1, 2),
                           # direction = "y"
         )
         + ylab("")
         # + xlab("Uncorrected AF")
         + labs(title = id,
                subtitle = paste0("TC = ", tc[sample == id]$tc, "\t\tasP = ", asP[sample == id]$asP)
                )
         + theme_bw(base_size = 15)
         + theme(plot.title = element_text(hjust = 0.5),
                 plot.subtitle = element_text(size = 15, hjust = 0.5),
                 legend.position="bottom",
                 axis.text = element_text(color = "black"),
                 axis.text.y = element_blank(),
                 axis.ticks.y = element_blank())
         + guides(color = F)
  )
  p.corr <- (ggplot(x, aes(x = t_af_corr, y = id, color = Chromosome, label = gene
                           # fill = gene.class, shape = imb, label = gene
                           ))
             + geom_point(size = 2)
             + scale_color_manual(name = "Chromosome", values = chr.cols)
             # + xlim(c(0,1))
             + scale_x_continuous(name = "Corrected AF", breaks = seq(0,1,0.1), limits = c(0,1))
             + coord_fixed(ratio = 0.2)
             + geom_text_repel(data = x,
                               size = 2,
                               # colour = "black",
                               # xlim = c(-2, -1.5),
                               ylim = c(1.1, 2),
                               # direction = "y"
             )
             + ylab("")
             # + xlab("Uncorrected AF")
             # + labs(title = id
             #        # caption = paste0("TC = ", tc[sample == id]$tc, "\t\tasP = ", asP[sample == id]$asP)
             # )
             + theme_bw(base_size = 15)
             + theme(legend.position="bottom",
                     axis.text = element_text(color = "black"),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank())
             + guides(color = F)
             )
  print(p.unc/p.corr)
}
dev.off()