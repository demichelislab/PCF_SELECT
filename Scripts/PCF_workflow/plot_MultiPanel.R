library(data.table)
library(parallel)
library(GenomicRanges)
library(ggsci)
#library(annotatr)

source("PATHS.R")

load(path_germ_distr)
tmp = germ.distr[[1]]
tmp = tmp[which(tmp$af.mean>0.35&tmp$af.mean<0.65),]
tmp = tmp[which((tmp$af.cv<=0.1&tmp$af.alt>1)),]
tmp = tmp[which(tmp$af.freq<0.8),]
germ.distr[[1]] = tmp
load(path_beta_table)
beta.table <- ConvertBetaTable(beta.table)

#### Correct cn.call.corr names if TC lower than TC_thr ####
TC_thr <- thr.tc

Calls <- fread(paste0("Analyses/", run_name, "/CN_SNVs_calls.csv"))
Calls <- correctDoubleGenes(Calls, unique(Calls$gene), unique(Calls$sample))
####

load(paste0("Analyses/", run_name, "/.focal_tables/bed_with_rc.RData"))
bed.rc$V4 <- gsub("\\|", ";", bed.rc$V4)
bed.rc$V4 <- paste0(bed.rc$V4, ";")
#annotations_data <- suppressMessages(annotatr::build_annotations(genome = "hg19", annotations = paste("hg19", "genes", c("exons", "introns", "intergenic"), sep = "_")))
load("AnnotationFiles/hg19_annotation.RData")
annotations_data <- hg19_annotations
annotate_regions_PCF <- function(d1, d2 = NULL, genome = "hg19", annotation_name = "type", simplify_annotations = F, NA_as_intergenic = T, remove_chr_label = T){
  # use default if annotatr::build_annotations() is used to generate annotations
  if (is.null(d2)) {
    d2 <- suppressMessages(annotatr::build_annotations(genome = genome, annotations = paste(genome, "genes", c("exons", "introns", "intergenic"), sep = "_")))
    simplify_annotations <- T
  }
  d1 <- as.data.table(d1)
  d2 <- as.data.table(d2)
  colnames(d1)[1:3] <- c("chr", "from", "to")
  colnames(d2)[1:3] <- c("chr", "from", "to")
  if (!unique(grepl("chr", d1$chr))) {
    d1$chr <- paste0("chr", d1$chr)
  }
  query <- GRanges(seqnames=d1$chr, IRanges(start=(d1$from),end=(d1$to)))
  bed.ranges <- GRanges(seqnames=d2$chr, IRanges(start=d2$from,end=d2$to))
  hits <- findOverlaps(query, bed.ranges)
  dt <- cbind(d1[queryHits(hits),], d2[subjectHits(hits), get(annotation_name)])
  dt <- dt[!duplicated.data.frame(dt),]
  colnames(dt)[ncol(dt)] <- annotation_name
  dt <- plyr::ddply(dt, .variables = setdiff(colnames(dt), annotation_name), 
              function(x) paste(sort(x[,which(colnames(x) == annotation_name)]), collapse = "-"))#, .progress = "text")
  colnames(dt)[ncol(dt)] <- annotation_name
  out <- merge(d1, dt, all.x = T, by = setdiff(colnames(d1), annotation_name))
  if (simplify_annotations) {
    out$type <- gsub("hg19_genes_", "", out$type)
    if (NA_as_intergenic) {
      out$type[is.na(out$type)] <- "intergenic" 
    }
  }
  if (remove_chr_label) {
    out$chr <- gsub("chr", "", out$chr) 
  }
  return(out)
}

sif = fread(path_sif ,data.table=F)
#sif$Normal <- gsub("/shares/CIBIO-Storage/CO", "", sif$Normal)
#sif$Plasma <- gsub("/shares/CIBIO-Storage/CO", "", sif$Plasma)
sif$Normal_pile <- gsub(".snps", ".pileup", sif$Normal)
sif$Plasma_pile <- gsub(".snps", ".pileup", sif$Plasma)
pileup.dir = path_dirRC
sif$Normal_rc <- paste0(pileup.dir, basename(gsub("snps", "rc", sif$Normal)))
sif$Plasma_rc <- paste0(pileup.dir, basename(gsub("snps", "rc", sif$Plasma)))
coords <- fread(path_panelgenescoords)
#### TC annotation with manual or CLONETv2 estimation ####
tc <- fread(paste0("Analyses/", run_name, "/tc_estimations_CLONETv2.tsv"), header = T)
tc$tc <- tc$tc*100
tc$is.manual <- ifelse(is.na(tc$tc_CLONETv2 & !is.na(tc$tc_manual)), "y", "n")

extract_snps <- function(gene, sample, sif, coords, germ.distr, sample_type){
  ### extract het snps limited to gene
  j <- which(coords$gene_symbol == gene)
  p = fread(sif$Normal[grep(paste0("/", sample), sif$Plasma)])
  p = p[which(p$chr==coords$chr[j]&p$pos>=coords$start[j]&p$pos<=coords$end[j])]
  p = p[which(p$af<0.8 & p$af>0.2),]
  het.snps = p$rsid
  
  if (sample_type == "Normal") {
    ### extract AF of hetsnps in germ sample
    p = p[which(p$rsid%in%germ.distr[[1]][,1]),] 
    return(p)
  }else{
    ### extract AF of hetsnps in tumor sample
    p = fread(sif$Plasma[grep(paste0("/", sample), sif$Plasma)])
    p = p[which(p$rsid%in%het.snps),]
    p = p[which(p$rsid%in%germ.distr[[1]][,1]),]
    return(p)
  }
}

extract_cov <- function(gene, sample, sif, coords, sample_type){
  j <- which(coords$gene_symbol == gene)
  if(sample_type == "Normal"){
    p = fread(sif$Normal_pile[grep(paste0("/", sample), sif$Plasma_pile)])
    p = p[which(p$chr==coords$chr[j]&p$pos>=coords$start[j]&p$pos<=coords$end[j])]
    return(p)
  }else{
    p = fread(sif$Plasma_pile[grep(paste0("/", sample), sif$Plasma_pile)])
    p = p[which(p$chr==coords$chr[j]&p$pos>=coords$start[j]&p$pos<=coords$end[j])]
    return(p)
  }
}

extract_cov2 <- function(gene, sample, sif, bed.rc){
  s <- sample
  tmp <- bed.rc[sample == s & grepl(paste0(gene,";"), V4)]
  zzz <- annotate_regions_PCF(tmp, annotations_data, simplify_annotations = T)
  zzz <- zzz[, c("chr", "from", "to", "log2", "type")]
  zzz$from <- zzz$from + 1
  p = fread(sif$Normal_pile[grep(paste0("/", sample), sif$Plasma_pile)], select = c(1,2))
  query <- GRanges(seqnames=p$chr, IRanges(start=(p$pos),end=(p$pos)))
  bed.ranges <- GRanges(seqnames=zzz$chr, IRanges(start=zzz$from,end=zzz$to))
  hits <- findOverlaps(query, bed.ranges)
  ZZZ <- cbind(p[queryHits(hits),], zzz[subjectHits(hits), c("log2", "type")])
  return(ZZZ)
}

plot_focality <- function(x, sample_name, gene_name){
  if (gene_name == "ERG_TMPRSS2") {
    z <- x
    plot(x = z$value, y = z$median.log2, ylim = c(-3,3), pch = 20, xlim = c(min(z$value), max(z$value)),
         xlab = "Genomic Position", ylab = "Uncorrected Median Log2R", las = 1, #main = paste0(sample_name, "\n", gene_name),
         sub = paste0("Log2R=", round(z$median.log2[nrow(z)], digits = 3))
    )
    abline(h = -0.5, lty = "dotted", col = "gray50")
    abline(h = 0.5, lty = "dotted", col = "gray50")
    segments(x0 = z$value[1], x1 = z$value[2], y0 = z$median.log2[1], y1 = z$median.log2[1], col = "black")
  }else{
    cols <- rep("black", nrow(x))
    cols[1] <- "red3"
    cols[nrow(x)] <- "turquoise4"
    z <- melt(x, id.vars = "median.log2", measure.vars = c("left_lim", "right_lim"))
    plot(x = z$value, y = z$median.log2, ylim = c(-3,3), pch = 20, xlim = c(min(z$value), max(z$value)),
         xlab = "Genomic Position", ylab = "Uncorrected Median Log2R", las = 1, #main = paste0(sample_name, "\n", gene_name),
         sub = paste0("Focal Log2R=", round(z$median.log2[1], digits = 3), "\tLog2R=", round(z$median.log2[nrow(z)], digits = 3))
    )
    #axis(1, at=z$value, labels = z$value)
    axis(1, at = c(x$left_lim[c(1,nrow(x))], x$right_lim[c(1,nrow(x))]), 
         labels = rep("", 4))
    #labels = c(x$left_lim[c(1,nrow(x))], x$right_lim[c(1,nrow(x))]))
    segments(x0 = c(x$left_lim[1], x$right_lim[1]), x1 = c(x$left_lim[1], x$right_lim[1]), 
             y0 = c(-3,-3), y1 = c(-2.95,-2.95), col = cols[1], lty = "dotted")
    segments(x0 = x$left_lim[1], x1 = x$right_lim[1], y0 = -2.95, y1 = -2.95, col = cols[1], lty = "dotted")
    text(x = mean(c(x$left_lim[1], x$right_lim[1])), y = -2.85, labels = x$right_lim[1] - x$left_lim[1],
         col = cols[1], cex = 0.7)
    
    segments(x0 = c(x$left_lim[nrow(x)], x$right_lim[nrow(x)]), x1 = c(x$left_lim[nrow(x)], x$right_lim[nrow(x)]), 
             y0 = c(-3,-3), y1 = c(-2.65,-2.65), col = cols[nrow(x)], lty = "dotted")
    segments(x0 = x$left_lim[nrow(x)], x1 = x$right_lim[nrow(x)], y0 = -2.65, y1 = -2.65, col = cols[nrow(x)], lty = "dotted")
    text(x = mean(c(x$left_lim[nrow(x)], x$right_lim[nrow(x)])), y = -2.55, labels = x$right_lim[nrow(x)] - x$left_lim[nrow(x)],
         col = cols[nrow(x)], cex = 0.7)
    abline(h = -0.5, lty = "dotted", col = "gray50")
    abline(h = 0.5, lty = "dotted", col = "gray50")
    # segments(x0 = x$left_lim[1], x1 = x$left_lim[1], y0 = x$median.log2[1], y1 = 3,
    #          lty = "dotted", col = cols[1])
    # segments(x0 = x$right_lim[1], x1 = x$right_lim[1], y0 = x$median.log2[1], y1 = 3,
    #          lty = "dotted", col = cols[1])
    # segments(x0 = x$left_lim[nrow(x)], x1 = x$left_lim[nrow(x)], y0 = x$median.log2[nrow(x)], y1 = 3,
    #          lty = "dotted", col = cols[nrow(x)])
    # segments(x0 = x$right_lim[nrow(x)], x1 = x$right_lim[nrow(x)], y0 = x$median.log2[nrow(x)], y1 = 3,
    #          lty = "dotted", col = cols[nrow(x)])
    segments(x0 = x$left_lim, x1 = x$right_lim, y0 = x$median.log2, y1 = x$median.log2, col = cols)
  }
}

c <- list.files(path = paste0("Analyses/", run_name, "/.focal_tables/"), pattern = "control", full.names = T)
dt <- rbindlist(lapply(c, function(x) {load(x); focal_control_genes}))

dir.create(paste0("Analyses/", run_name, "/Plots/MultiPanel"))

message("Control genes")
dir.create(paste0("Analyses/", run_name, "/Plots/MultiPanel/control_genes"))
mclapply(unique(dt$gene), function(g){
  message(g)
  pdf(file = paste0("Analyses/", run_name, "/Plots/MultiPanel/control_genes/", g, ".pdf"))
  for (s in unique(dt$sample)) {
    #layout(matrix(c(1,2,3,4,4,4), nrow = 6, ncol = 1, byrow = TRUE))
    #layout(matrix(c(1,2,3,3), nrow = 4, ncol = 1, byrow = TRUE))
    layout(matrix(c(1,2,3,4,4), nrow = 5, ncol = 1, byrow = TRUE))
    x <- dt[gene == g & sample == s]

    #plot AF germline
    p <- extract_snps(g, s, sif, coords, germ.distr, "Normal")
    z <- melt(x, id.vars = "median.log2", measure.vars = c("left_lim", "right_lim"))
    par(mar = c(0.1, 4, 0.1, 2) + 0.1)
    plot(x = p$pos, y = p$af, ylim = c(0,1), type = "h", xaxt='n', yaxt = 'n',
         xlim = c(min(z$value), max(z$value)), xlab = "", ylab = "AF Control", las = 1)
    text(min(z$value), 1, labels = paste0("AI ev.= ", beta.table[gene == g & sample == s]$evidence.n, "\n",
                                          "BetaN = ", round(beta.table[gene == g & sample == s]$beta.n, digits = 3)),
         adj = c(0,1), cex = 0.8)
    axis(2, at = seq(0,1,0.25), las = 1)
    abline(h = 0.5, col = "red", lty = "dashed")

    #plot AF tumor
    p <- extract_snps(g, s, sif, coords, germ.distr, "Tumor")
    z <- melt(x, id.vars = "median.log2", measure.vars = c("left_lim", "right_lim"))
    par(mar = c(0.1, 4, 0.1, 2) + 0.1)
    plot(x = p$pos, y = p$af, ylim = c(0,1), type = "h", xaxt='n', yaxt = 'n',
         xlim = c(min(z$value), max(z$value)), xlab = "", ylab = "AF ctDNA", las = 1)
    text(min(z$value), 1, labels = paste0("AI ev.= ", beta.table[gene == g & sample == s]$evidence, "\n",
                                          "BetaT = ", round(beta.table[gene == g & sample == s]$beta, digits = 3)),
         adj = c(0,1), cex = 0.8)
    axis(2, at = seq(0,1,0.25), las = 1)
    abline(h = 0.5, col = "red", lty = "dashed")

    #plot coverage
    par(mar = c(0.1, 4, 0.1, 2) + 0.1)
    pile <- extract_cov2(g, s, sif, bed.rc)
    pile$type <- ifelse(pile$type%in%c("introns","exons","intergenic"), pile$type, "mixed")
    cols_vector <- c(ggsci::pal_npg("nrc")(6)[c(1,5,6)], ggsci::pal_aaas()(3)[3])
    names(cols_vector) <- c("exons","intergenic","introns","mixed")
    pile$colors <- cols_vector[match(pile$type, names(cols_vector))]
    plot(x = pile$pos, y = pile$log2, ylim = c(-3,3), type = "h", xaxt='n', yaxt = 'n',
         xlim = c(min(z$value), max(z$value)), xlab = "", ylab = "Uncorrected Log2R", las = 1,
         col = pile$colors)
    axis(2, at = seq(-3,3,1), las = 1)
    for (h in seq(-3,3,.5)) {
      abline(h = h, col = "gray50", lty = "dotted")
    }
    legend(min(z$value), 3, legend = Hmisc::capitalize(names(cols_vector)),
           col = cols_vector, pch = 19, cex = 0.65, horiz = T, box.col = "black", bg = "white")

    par(mar = c(5, 4, 0.1, 2) + 0.1)
    plot_focality(x, s, g)

    mtext(paste0(s, "\n", Calls[sample == s & gene == g]$cn.call.corr),
          side = 1, line = 4, adj = 1, cex = 0.7)
    #mtext(paste0("TC = ", as.integer((1-tc[sample == s]$adm)*100), " %", ifelse(tc[sample == s]$is.manual == "y", "*", ""),
    mtext(paste0("TC = ", tc[sample == s]$tc, " %", ifelse(tc[sample == s]$is.manual == "y", "*", ""),
                 "\nPloidy = ", tc[sample == s]$ploidy),
          side = 1, line = 4, adj = 0, cex = 0.7)
  }
  dev.off()
}, mc.cores = cores.multipanel
)

t <- list.files(path = paste0("Analyses/", run_name, "/.focal_tables/"), pattern = "target", full.names = T)
dt <- rbindlist(lapply(t, function(x) {load(x); focal_target_genes}))

message("Target genes")
dir.create(paste0("Analyses/", run_name, "/Plots/MultiPanel/target_genes"))
mclapply(unique(dt$gene), function(g){
  message(g)
  pdf(file = paste0("Analyses/", run_name, "/Plots/MultiPanel/target_genes/", g, ".pdf"))
  for (s in unique(dt$sample)) {
    #layout(matrix(c(1,2,3,4,4,4), nrow = 6, ncol = 1, byrow = TRUE))
    #layout(matrix(c(1,2,3,3), nrow = 4, ncol = 1, byrow = TRUE))
    layout(matrix(c(1,2,3,4,4), nrow = 5, ncol = 1, byrow = TRUE))
    x <- dt[gene == g & sample == s]
    
    #plot AF germline
    p <- extract_snps(g, s, sif, coords, germ.distr, "Normal")
    z <- melt(x, id.vars = "median.log2", measure.vars = c("left_lim", "right_lim"))
    par(mar = c(0.1, 4, 0.1, 2) + 0.1)
    plot(x = p$pos, y = p$af, ylim = c(0,1), type = "h", xaxt='n', yaxt = 'n', 
         xlim = c(min(z$value), max(z$value)), xlab = "", ylab = "AF Control", las = 1)
    text(min(z$value), 1, labels = paste0("AI ev.= ", beta.table[gene == g & sample == s]$evidence.n, "\n",
                                          "BetaN = ", round(beta.table[gene == g & sample == s]$beta.n, digits = 3)), 
         adj = c(0,1), cex = 0.8)
    axis(2, at = seq(0,1,0.25), las = 1)
    abline(h = 0.5, col = "red", lty = "dashed")
    
    #plot AF tumor
    p <- extract_snps(g, s, sif, coords, germ.distr, "Tumor")
    z <- melt(x, id.vars = "median.log2", measure.vars = c("left_lim", "right_lim"))
    par(mar = c(0.1, 4, 0.1, 2) + 0.1)
    plot(x = p$pos, y = p$af, ylim = c(0,1), type = "h", xaxt='n', yaxt = 'n', 
         xlim = c(min(z$value), max(z$value)), xlab = "", ylab = "AF ctDNA", las = 1)
    text(min(z$value), 1, labels = paste0("AI ev.= ", beta.table[gene == g & sample == s]$evidence, "\n",
                                          "BetaT = ", round(beta.table[gene == g & sample == s]$beta, digits = 3)), 
         adj = c(0,1), cex = 0.8)
    axis(2, at = seq(0,1,0.25), las = 1)
    abline(h = 0.5, col = "red", lty = "dashed")
    
    #plot coverage
    par(mar = c(0.1, 4, 0.1, 2) + 0.1)
    pile <- extract_cov2(g, s, sif, bed.rc)
    pile$type <- ifelse(pile$type%in%c("introns","exons","intergenic"), pile$type, "mixed")
    cols_vector <- c(ggsci::pal_npg("nrc")(6)[c(1,5,6)], ggsci::pal_aaas()(3)[3])
    names(cols_vector) <- c("exons","intergenic","introns","mixed")
    pile$colors <- cols_vector[match(pile$type, names(cols_vector))]
    plot(x = pile$pos, y = pile$log2, ylim = c(-3,3), type = "h", xaxt='n', yaxt = 'n', 
         xlim = c(min(z$value), max(z$value)), xlab = "", ylab = "Uncorrected Log2R", las = 1,
         col = pile$colors)
    axis(2, at = seq(-3,3,1), las = 1)
    for (h in seq(-3,3,.5)) {
      abline(h = h, col = "gray50", lty = "dotted") 
    }
    legend(min(z$value), 3, legend = Hmisc::capitalize(names(cols_vector)), 
           col = cols_vector, pch = 19, cex = 0.65, horiz = T, box.col = "black", bg = "white")
    
    # par(mar = c(0.1, 4, 0.1, 2) + 0.1)
    # pile <- extract_cov(g, s, sif, coords, "Normal")
    # plot(x = pile$pos, y = pile$cov, ylim = c(-10,10), type = "h", xaxt='n', yaxt = 'n', 
    #      xlim = c(min(z$value), max(z$value)), xlab = "", ylab = "Cov. N", las = 1)
    # for (h in seq(500,3000,500)) {
    #   abline(h = h, col = "gray50", lty = "dotted") 
    # }
    # axis(2, at = seq(0,3000,1000), las = 1)
    # 
    # par(mar = c(0.1, 4, 0.1, 2) + 0.1)
    # pile <- extract_cov(g, s, sif, coords, "Tumor")
    # plot(x = pile$pos, y = pile$cov, ylim = c(0,3000), type = "h", xaxt='n', yaxt = 'n', 
    #      xlim = c(min(z$value), max(z$value)), xlab = "", ylab = "Cov. T", las = 1)
    # for (h in seq(500,3000,500)) {
    #   abline(h = h, col = "gray50", lty = "dotted") 
    # }
    # axis(2, at = seq(0,3000,1000), las = 1)
    
    par(mar = c(5, 4, 0.1, 2) + 0.1)
    plot_focality(x, s, g)
    
    #mtext(paste0(s), side=1, line=3, adj = 1, cex = 0.7)
    mtext(paste0(s, "\n", Calls[sample == s & gene == g]$cn.call.corr), 
          side = 1, line = 4, adj = 1, cex = 0.7)
    #mtext(paste0("TC = ", as.integer((1-tc[sample == s]$adm)*100), " %", ifelse(tc[sample == s]$is.manual == "y", "*", ""),
    mtext(paste0("TC = ", tc[sample == s]$tc, " %", ifelse(tc[sample == s]$is.manual == "y", "*", ""),
                 "\nPloidy = ", tc[sample == s]$ploidy),
          side = 1, line = 4, adj = 0, cex = 0.7)
  }
  dev.off()
}, mc.cores = cores.multipanel
)

# For ERG_TMPRSS2
message("Target genes (only ERG-TMPRSS2)")

g <- "ERG_TMPRSS2"
z <- as.data.table(rbind(coords[gene_symbol == g]$start, coords[gene_symbol == g]$end))
colnames(z) <- "value"

mclapply("ERG_TMPRSS2", function(g){
  message(g)
  pdf(file = paste0("Analyses/", run_name, "/Plots/MultiPanel/target_genes/", g, ".pdf"))
  for (s in unique(dt$sample)) {
    seg <- fread(paste0("Analyses/", run_name, "/Segmentation_focal.seg"))
    z$median.log2 <- seg[sample == s & gene == g]$log2ratio
    #layout(matrix(c(1,2,3,4,4,4), nrow = 6, ncol = 1, byrow = TRUE))
    #layout(matrix(c(1,2,3,3), nrow = 4, ncol = 1, byrow = TRUE))
    layout(matrix(c(1,2,3,4,4), nrow = 5, ncol = 1, byrow = TRUE))
    #x <- dt[gene == g & sample == s]
    #plot AF germline
    p <- extract_snps(g, s, sif, coords, germ.distr, "Normal")
    #z <- melt(x, id.vars = "median.log2", measure.vars = c("left_lim", "right_lim"))
    par(mar = c(0.1, 4, 0.1, 2) + 0.1)
    plot(x = p$pos, y = p$af, ylim = c(0,1), type = "h", xaxt='n', yaxt = 'n',
         xlim = c(min(z$value), max(z$value)), xlab = "", ylab = "AF Control", las = 1)
    text(min(z$value), 1, labels = paste0("AI ev.= ", beta.table[gene == g & sample == s]$evidence.n, "\n",
                                          "BetaN = ", round(beta.table[gene == g & sample == s]$beta.n, digits = 3)),
         adj = c(0,1), cex = 0.8)
    axis(2, at = seq(0,1,0.25), las = 1)
    abline(h = 0.5, col = "red", lty = "dashed")

    #plot AF tumor
    p <- extract_snps(g, s, sif, coords, germ.distr, "Tumor")
    #z <- melt(x, id.vars = "median.log2", measure.vars = c("left_lim", "right_lim"))
    par(mar = c(0.1, 4, 0.1, 2) + 0.1)
    plot(x = p$pos, y = p$af, ylim = c(0,1), type = "h", xaxt='n', yaxt = 'n',
         xlim = c(min(z$value), max(z$value)), xlab = "", ylab = "AF ctDNA", las = 1)
    text(min(z$value), 1, labels = paste0("AI ev.= ", beta.table[gene == g & sample == s]$evidence, "\n",
                                          "BetaT = ", round(beta.table[gene == g & sample == s]$beta, digits = 3)),
         adj = c(0,1), cex = 0.8)
    axis(2, at = seq(0,1,0.25), las = 1)
    abline(h = 0.5, col = "red", lty = "dashed")

    #plot coverage
    par(mar = c(0.1, 4, 0.1, 2) + 0.1)
    pile <- extract_cov2(g, s, sif, bed.rc)
    pile$type <- ifelse(pile$type%in%c("introns","exons","intergenic"), pile$type, "mixed")
    cols_vector <- c(ggsci::pal_npg("nrc")(6)[c(1,5,6)], ggsci::pal_aaas()(3)[3])
    names(cols_vector) <- c("exons","intergenic","introns","mixed")
    pile$colors <- cols_vector[match(pile$type, names(cols_vector))]
    plot(x = pile$pos, y = pile$log2, ylim = c(-3,3), type = "h", xaxt='n', yaxt = 'n',
         xlim = c(min(z$value), max(z$value)), xlab = "", ylab = "Uncorrected Log2R", las = 1,
         col = pile$colors)
    axis(2, at = seq(-3,3,1), las = 1)
    for (h in seq(-3,3,.5)) {
      abline(h = h, col = "gray50", lty = "dotted")
    }
    legend(min(z$value), 3, legend = Hmisc::capitalize(names(cols_vector)),
           col = cols_vector, pch = 19, cex = 0.65, horiz = T, box.col = "black", bg = "white")

    # par(mar = c(0.1, 4, 0.1, 2) + 0.1)
    # pile <- extract_cov(g, s, sif, coords, "Normal")
    # plot(x = pile$pos, y = pile$cov, ylim = c(-10,10), type = "h", xaxt='n', yaxt = 'n',
    #      xlim = c(min(z$value), max(z$value)), xlab = "", ylab = "Cov. N", las = 1)
    # for (h in seq(500,3000,500)) {
    #   abline(h = h, col = "gray50", lty = "dotted")
    # }
    # axis(2, at = seq(0,3000,1000), las = 1)
    #
    # par(mar = c(0.1, 4, 0.1, 2) + 0.1)
    # pile <- extract_cov(g, s, sif, coords, "Tumor")
    # plot(x = pile$pos, y = pile$cov, ylim = c(0,3000), type = "h", xaxt='n', yaxt = 'n',
    #      xlim = c(min(z$value), max(z$value)), xlab = "", ylab = "Cov. T", las = 1)
    # for (h in seq(500,3000,500)) {
    #   abline(h = h, col = "gray50", lty = "dotted")
    # }
    # axis(2, at = seq(0,3000,1000), las = 1)

    par(mar = c(5, 4, 0.1, 2) + 0.1)
    plot_focality(z, s, g)

    #mtext(paste0(s), side=1, line=3, adj = 1, cex = 0.7)
    mtext(paste0(s, "\n", Calls[sample == s & gene == g]$cn.call.corr), 
          side = 1, line = 4, adj = 1, cex = 0.7)
    #mtext(paste0("TC = ", as.integer((1-tc[sample == s]$adm)*100), " %", ifelse(tc[sample == s]$is.manual == "y", "*", ""),
    mtext(paste0("TC = ", tc[sample == s]$tc, " %", ifelse(tc[sample == s]$is.manual == "y", "*", ""),
                 "\nPloidy = ", tc[sample == s]$ploidy),
          side = 1, line = 4, adj = 0, cex = 0.7)
  }
  dev.off()
}, mc.cores = cores.multipanel
)
