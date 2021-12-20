library(data.table)
library(plotrix)

source("PATHS.R")

FOLDER_PCF_ANALYSES <- "/shares/CIBIO-Storage/CO/SPICE/downloads/PCF_SELECT/analyses/2020.02/"
FOLDER_PCF_DATA <- "/shares/CIBIO-Storage/CO/SPICE/downloads/PCF_SELECT/"
setwd(FOLDER_PCF_ANALYSES)

load(paste0("Analyses/", run_name, "/betaTable.RData"))
beta.table <- ConvertBetaTable(beta.table)
beta.table <- beta.table[,c("sample", "gene", "beta", "evidence", "evidence.n")]
seg <- fread(paste0("Analyses/", run_name, "/Segmentation_focal.seg"))

data <- merge(beta.table, seg, by = c("sample", "gene"), all = T)

data$ID <- data$sample
data <- data[order(data$ID)]

tc <- fread(paste0("Analyses/", run_name, "/tc_estimations_CLONETv2.tsv"))

dir.create(paste0("Analyses/", run_name, "/Plots"))
pdf(paste0("Analyses/", run_name, "/Plots/Log2_Beta_labels_chr.pdf"))

for (id in unique(data$ID)) {
  x <- data[ID == id,]
  x$col <- chr.cols[match(x$chr, names(chr.cols))]
  plot(x = x$log2ratio,
       y = x$beta, 
       xlim = c(-2,2), 
       ylim = c(0,1), 
       main = id,
       pch = ifelse(x$evidence >= thr.imbT, 17, 19),
       ylab = "Beta", 
       xlab = "Log2(R)", 
       # col = ifelse(x$evidence >= thr.imbT, "red", "black"),
       col = x$col,
       sub = paste0("TC = ", tc[sample == id]$tc, "\t\tPloidy = ", tc[sample == id]$ploidy)
       )
  log2_shift <- -log2(tc[sample == id]$ploidy/2)
  abline(v = 0+log2_shift, lty = "dashed", col = "red")
  legend("topright", pch = 19, legend = names(chr.cols), col = chr.cols, title = "Chromosome",
         cex=0.7, xjust=1, ncol=2, x.intersp=0.5, y.intersp=0.65)
  legend("bottomright", pch = c(17, 19),
         legend = c(as.expression(bquote(E(AI[T])>=.(thr.imbT))), bquote(E(AI[T])<.(thr.imbT))),
         title = "Allelic Imbalance",
         cex=0.7, xjust=1, ncol=1, x.intersp=0.5, y.intersp=0.65)
  i <- ifelse(x$evidence >= thr.imbT, T, F)
  if (TRUE %in% unique(i)) {
    spread.labels(na.omit(x$log2ratio[i]), na.omit(x$beta[i]), labels=na.omit(x$gene[i]), ony = T, cex = 0.55)
    #text(x$log2ratio[i], x$beta[i], labels=x$gene[i], cex= 0.55, font = 2)
  }
}
dev.off()

pdf(paste0("Analyses/", run_name, "/Plots/Log2_Beta_chr.pdf"))

for (id in unique(data$ID)) {
  x <- data[ID == id,]
  x$col <- chr.cols[match(x$chr, names(chr.cols))]
  plot(x = x$log2ratio,
       y = x$beta, 
       xlim = c(-2,2), 
       ylim = c(0,1), 
       main = id,
       pch = ifelse(x$evidence >= thr.imbT, 17, 19),
       ylab = "Beta", 
       xlab = "Log2(R)", 
       # col = ifelse(x$evidence >= thr.imbT, "red", "black"),
       col = x$col,
       sub = paste0("TC = ", tc[sample == id]$tc, "\t\tPloidy = ", tc[sample == id]$ploidy)
  )
  log2_shift <- -log2(tc[sample == id]$ploidy/2)
  abline(v = 0+log2_shift, lty = "dashed", col = "red")
  legend("topright", pch = 19, legend = names(chr.cols), col = chr.cols, title = "Chromosome",
         cex=0.7, xjust=1, ncol=2, x.intersp=0.5, y.intersp=0.65)
  legend("bottomright", pch = c(17, 19),
         legend = c(as.expression(bquote(E(AI[T])>=.(thr.imbT))), bquote(E(AI[T])<.(thr.imbT))),
         title = "Allelic Imbalance",
         cex=0.7, xjust=1, ncol=1, x.intersp=0.5, y.intersp=0.65)
  i <- ifelse(x$evidence >= thr.imbT, T, F)
}
dev.off()
