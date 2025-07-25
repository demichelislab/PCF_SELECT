library(data.table)
#library(basicPlotteR)

FOLDER_PCF_ANALYSES <- "/shares/CIBIO-Storage/CO/SPICE/downloads/PCF_SELECT/analyses/2020.02/"
FOLDER_PCF_DATA <- "/shares/CIBIO-Storage/CO/SPICE/downloads/PCF_SELECT/"
setwd(FOLDER_PCF_ANALYSES)
source("Scripts/utility_functions.R")
source("Scripts/parameters.R")

run_name = "run001"

plot_cnA_cnB <- function(AI, sample_name, yMax = NA, xMax = NA, add.labels = F){
  thisAI <- AI
  thisAI$isUnbalancedGain <- ifelse(thisAI$cnA.int > thisAI$cnB.int + 1 & thisAI$cnB.int >= 1, 1, 0)
  thisAI$isCNNL <- ifelse(thisAI$cnA.int == 2 & thisAI$cnB.int == 0, 1, 0)
  thisAI$isLOHg <- ifelse(thisAI$cnA.int > 2 & thisAI$cnB.int == 0, 1, 0)
  ## predefinite colors
  copyNeutralLOH.col<- rgb(101,86,67,maxColorValue = 255, alpha = 255)
  copyAberrantLOH.col <- rgb(128,188,163,maxColorValue = 255, alpha = 255)
  AIgain.col <- rgb(230,172,39,maxColorValue = 255, alpha = 255)
  
  ## set colors
  thisAI$color <- "gray60"
  thisAI$color[which(thisAI$isCNNL == 1)] <- copyNeutralLOH.col
  thisAI$color[which(thisAI$isLOHg == 1)] <- copyAberrantLOH.col
  thisAI$color[which(thisAI$isUnbalancedGain == 1)] <- AIgain.col
  
  xMin <- -0.1
  yMin <- -0.1
  
  if (is.na(yMax) & is.na(xMax)) {
    yMax <- max(ceiling(thisAI$cnB), na.rm = T)
    xMax <- max(ceiling(thisAI$cnA), na.rm = T) 
  }
  
  layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
  
  par(mar=c(5.1,4.1,2,1))
  plot(thisAI$cnA,thisAI$cnB, pch=20, xlim=c(xMin,xMax), ylim=c(yMin,yMax), col=thisAI$color,asp=1,axes=F,xlab="",ylab="")
  axis(1,at=seq(0,xMax,1))  
  axis(2,at=seq(0,yMax,1))  
  abline(v = seq(0,xMax,1), lty=3, col="gray60")
  abline(h = seq(0,yMax,1), lty=3, col="gray60")
  abline(coef = c(0,1), col="gray40")
  
  if (add.labels) {
    text_anno <- thisAI[(cnA.int == 0 & cnB.int == 0) | (cnA.int == 1 & cnB.int == 0) | (cnA.int == 2 & cnB.int == 0) | isUnbalancedGain == 1]
    text_anno <- thisAI
    
    text(x = text_anno$cnA, y = text_anno$cnB, labels = text_anno$gene, cex = .8)
    #basicPlotteR::addTextLabels(text_anno$cnA, text_anno$cnB, labels = text_anno$gene_symbol, col.label = "black")
  }
  title(xlab = "cnA", ylab = "cnB", main = sample_name)
  
  par(mar=c(0,0,2,0))
  plot(NA,xlim=c(xMin,xMax), ylim=c(yMin,yMax),axes=F,xlab="",ylab="")
  legend(0,yMax, legend = c("NA","Copy-neutral LOH","Copy-aberrant LOH", "Allele specific gain"),bty = "n",
         fill = c("gray60",copyNeutralLOH.col,copyAberrantLOH.col,AIgain.col))
}

gl <- fread("AnnotationFiles/GeneList_PCF-SELECT_Annotation_Sept2020.tsv")
control.genes = gl[PCF_SELECT == "control_gene"]$hgnc_symbol
other.genes = gl[PCF_SELECT == "other"]$hgnc_symbol
ai <- fread(paste0("Analyses/", run_name,"/ai_log2_table.csv"))
ai <- ai[!is.na(cnA)]

pdf(paste0("Analyses/", run_name, "/Plots/cnAcnB.pdf"))
for(s_name in sort(unique(ai$sample))){
  thisAI <- ai[sample == s_name]
  thisAI <- thisAI[!duplicated(gene)]
  thisAI <- thisAI[!gene %in% c(control.genes, other.genes)]
  plot_cnA_cnB(thisAI, s_name, add.labels = F)
}
dev.off()

pdf(paste0("Analyses/", run_name, "/Plots/cnAcnB_samescale.pdf"))
for(s_name in sort(unique(ai$sample))){
  thisAI <- ai[sample == s_name]
  thisAI <- thisAI[!duplicated(gene)]
  thisAI <- thisAI[!gene %in% c(control.genes, other.genes)]
  plot_cnA_cnB(thisAI, s_name, xMax = 4, yMax = 3, add.labels = F)
}
dev.off()

pdf(paste0("Analyses/", run_name, "/Plots/cnAcnB_labels.pdf"))
for(s_name in sort(unique(ai$sample))){
  thisAI <- ai[sample == s_name]
  thisAI <- thisAI[!duplicated(gene)]
  thisAI <- thisAI[!gene %in% c(control.genes, other.genes)]
  plot_cnA_cnB(thisAI, s_name, add.labels = T)
}
dev.off() 

pdf(paste0("Analyses/", run_name, "/Plots/cnAcnB_samescale_labels.pdf"))
for(s_name in sort(unique(ai$sample))){
  thisAI <- ai[sample == s_name]
  thisAI <- thisAI[!duplicated(gene)]
  thisAI <- thisAI[!gene %in% c(control.genes, other.genes)]
  plot_cnA_cnB(thisAI, s_name, xMax = 4, yMax = 3, add.labels = T)
}
dev.off()