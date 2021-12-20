library(CLONETv2)
library(TPES)
library(parallel)
library(GenomicRanges)
library(data.table)
library(stringr)

ConvertBetaTableToClonetv2Format <- function(beta.table, seg, bed){
  # Annotate segfile with gene names
  query <- GRanges(seqnames=seg$chr, IRanges(start=(seg$start),end=(seg$end)))
  bed.ranges <- GRanges(seqnames=bed$chr, IRanges(start=bed$start,end=bed$end))
  hits <- findOverlaps(query, bed.ranges)
  dt <- cbind(seg[queryHits(hits),], bed[subjectHits(hits), c("gene")])
  dt <- dt[!duplicated.data.frame(dt),]
  
  betas <- as.data.table(beta.table[,c("sample", "gene","beta","snps","snps.n","cov","cov.n","beta.n","evidence","evidence.n")])
  betas$sample <- gsub(x = betas$sample, pattern = ".snps", replacement = "")
  colnames(betas)[1] <- c("ID")
  
  dt$ID <- dt$sample
  data <- merge(dt, betas, by = c("ID", "gene"))
  data$ID <- NULL
  colnames(data) <- c("gene", "sample", "chr", "start", "end", "n.markers", "log2", "beta", "nsnp","snps.n", "cov","cov.n", "n_beta", "evidence","evidence.n")
  data$beta <- as.numeric(data$beta)
  data$nsnp <- as.numeric(data$nsnp)
  data$snps.n <- as.numeric(data$snps.n)
  data$cov <- as.numeric(data$cov)
  data$cov.n <- as.numeric(data$cov.n)
  data$n_beta <- as.numeric(data$n_beta)
  data$evidence <- as.numeric(data$evidence)
  data$evidence.n <- as.numeric(data$evidence.n)
  # Standard CLONETv2 beta.table
  #data <- data[,c("sample","chr","start","end","n.markers","log2","beta","nsnp","cov","n_beta")]
  # PCF_SELECT beta.table
  data <- data[,c("sample","gene","chr","start","end","n.markers","log2","beta","evidence","nsnp","cov","n_beta","evidence.n")]
  return(data)
}
ConvertAbemustoMAFforTPES <- function(sampleID, dir_abemus, optimalR=F, filter_snps=F, common_snps = NULL, af_control_thr = 0.01, cov_case_thr = 10){
  source(file.path(dir_abemus, "abemus_configure.R"))
  sif <- fread(filesif)
  file <- sif[plasma == sampleID]
  path_to_tab <- file.path(dir_abemus, "Results", file$patient)
  if (T %in% grepl(file$plasma, list.files(path_to_tab, pattern = "F3"))) {
    if (optimalR) {
      tb <- fread(paste0(path_to_tab, "/pmtab_F3_optimalR_", file$plasma, ".tsv"), select = c(2,3,12,19,5,27))
      if (filter_snps) {
        rsids <- readLines(common_snps)
        tb <- tb[!(dbsnp %in% rsids) & af_control < af_control_thr & cov_case >= cov_case_thr]
        tb <- tb[,1:4]
      }
      tb <- tb[,1:4]
    }else{
      tb <- fread(paste0(path_to_tab, "/pmtab_F3_", file$plasma, ".tsv"), select = c(2,3,12,19,5,27))
      if (filter_snp) {
        rsids <- readLines(common_snps)
        tb <- tb[!(dbsnp %in% rsids)& cov_case >= cov_case_thr] 
        tb <- tb[,1:4]
      }
      tb <- tb[,1:4]
    }
    if (!grepl("chr", tb$chr[1])) {
      tb$chr <- paste0("chr", tb$chr)
    }
    tb$start <- tb$pos
    tb$end <- tb$pos
    tb$ref.count <- tb$cov_case - tb$cov.alt
    tb$alt.count <- tb$cov.alt
    tb$pos <- NULL
    tb$cov_case <- NULL
    tb$cov.alt <- NULL
    tb$sample <- sampleID
    return(as.data.frame(tb))
  }else{
    message("No SNVs called by ABEMUS")
  }
}

ConvertAbemustoMAFLite <- function(tb){
  if (!grepl("chr", tb$chr[1])) {
    tb$chr <- paste0("chr", tb$chr)
  }
  tb$start <- tb$pos
  tb$end <- tb$pos
  tb$ref_allele <- tb$ref
  tb$alt_allele <- tb$alt
  j <- match(setdiff(colnames(tb), c("chr","start", "end", "ref_allele","alt_allele")), colnames(tb))
  tb[,j] <- NULL
  return(tb)
}

ConvertBetaTable <- function(beta.table){
  beta.table <- as.data.table(beta.table)
  beta.table$beta <- as.numeric(beta.table$beta)
  beta.table$beta.max <- as.numeric(beta.table$beta.max)
  beta.table$beta.min <- as.numeric(beta.table$beta.min)
  beta.table$snps <- as.numeric(beta.table$snps)
  beta.table$snps.n <- as.numeric(beta.table$snps.n)
  beta.table$cov <- as.numeric(beta.table$cov)
  beta.table$cov.n <- as.numeric(beta.table$cov.n)
  beta.table$beta.n <- as.numeric(beta.table$beta.n)
  beta.table$beta.n.max <- as.numeric(beta.table$beta.n.max)
  beta.table$beta.n.min <- as.numeric(beta.table$beta.n.min)
  beta.table$evidence <- as.numeric(beta.table$evidence)
  beta.table$evidence.n <- as.numeric(beta.table$evidence.n)
  return(beta.table)
}

match_TPs <- function(x, sample_match){
  sample_match$date <- paste0(stringr::str_split(sample_match$sample, "-", simplify = T)[,1], "-", sample_match$date)
  x$date <- sample_match$date[match(x$sample, sample_match$sample)]
  x$date[is.na(x$date)] <- x$sample[is.na(x$date)]
  x$sample <- x$date
  x$date <- NULL
  return(x)
}

rename_PRIMEsamples <- function(tab, matching){
  tab$orig_name <- stringr::str_split(tab$sample, "\\.", simplify = T)[,1]
  tab <- merge(tab, matching, all.x = T, by = "orig_name")
  tab$sample <- sapply(1:nrow(tab), function(i) gsub(tab$orig_name[i], tab$name[i], tab$sample[i]), simplify = T)
  tab$orig_name <- NULL
  tab$name <- NULL
  tab
}

annotate_GenomicRegions <- function(data, annotation, col_names = c()){
  suppressMessages(library(GenomicRanges))
  if (!length(intersect(c("chr", "start", "end"), colnames(data))) == 3) {
    stop("Data to be annotated needs chr, start and end columns")
  }else if (!length(intersect(c("chr", "start", "end"), colnames(annotation))) == 3){
    stop("Annotation data needs chr, start and end columns")
  }
  if (!(unique(grepl("chr", data$chr)) == unique(grepl("chr", annotation$chr)))) {
    stop("chr names need to be concordant")
  }
  if (length(col_names) == 0) {
    stop("Please specify names of the cols to be added to data")
  }
  query <- GRanges(seqnames=data$chr, IRanges(start=(data$start),end=(data$end)))
  bed.ranges <- GRanges(seqnames=annotation$chr, IRanges(start=annotation$start,end=annotation$end))
  hits <- findOverlaps(query, bed.ranges)
  dt <- cbind(data[queryHits(hits),], annotation[subjectHits(hits), ..col_names])
  dt <- dt[!duplicated.data.frame(dt),]
  return(dt)
}

Onco.AnnotateLabelsCalls <- function(aiTable, genes, genes_on_chrX){
  aiTable[cn.call.corr == "Gain", "cn.call.corr"] <- "DelGain"
  aiTable[cn.call.corr == "Amplification", "cn.call.corr"] <- "Gain"
  aiTable[cn.call.corr == "Uncertain HomoDel", "cn.call.corr"] <- "Uncertain_Homo"
  aiTable[cn.call.corr == "Uncertain Amplification", "cn.call.corr"] <- "Uncertain_Gain"
  aiTable[cn.call.corr == "GrayZone_Imb", "cn.call.corr"] <- "GrayZoneImb"
  aiTable[cn.call.corr == "GrayZone_NoImb", "cn.call.corr"] <- "GrayZoneNoImb"
  aiTable[cn.call.corr == "" | cn.call.corr == "Not enough SNPs", "cn.call.corr"] <- ""
  aiTable[cn.call.corr == "wt", "cn.call.corr"] <- ""
  
  if ("AR" %in% genes) {
    aiTable[gene %in% genes_on_chrX & cn.call.corr == "DelGain", "cn.call.corr"] <- "Gain"
  }
  aiTable
}
Onco.AnnotateLabelsCalls_lowTC <- function(aiTable, genes, genes_on_chrX){
  aiTable[cn.call.corr == "Gain", "cn.call.corr"] <- "DelGain"
  aiTable[cn.call.corr == "Amplification", "cn.call.corr"] <- "Gain"
  aiTable[cn.call.corr == "HomoDel", "cn.call.corr"] <- "Deletion"
  aiTable[cn.call.corr == "Uncertain HomoDel", "cn.call.corr"] <- "Deletion"
  aiTable[cn.call.corr == "Uncertain Amplification", "cn.call.corr"] <- "Gain"
  aiTable[cn.call.corr == "GrayZone_Imb", "cn.call.corr"] <- "GrayZoneImb"
  aiTable[cn.call.corr == "GrayZone_NoImb", "cn.call.corr"] <- "GrayZoneNoImb"
  aiTable[cn.call.corr == "" | cn.call.corr == "Not enough SNPs", "cn.call.corr"] <- ""
  aiTable[cn.call.corr == "wt", "cn.call.corr"] <- ""
  
  if ("AR" %in% genes) {
    aiTable[gene %in% genes_on_chrX & cn.call.corr == "DelGain", "cn.call.corr"] <- "Gain"
  }
  aiTable
}

correctDoubleGenes <- function(aiTable, genes, samples){
  double_genes <- c("ERCC1", "ERCC2", "MSH2", "MSH6", "RNF43", "RAD51C")
  rbindlist(lapply(samples, function(s){
    rbindlist(lapply(genes, function(g) {
      if (g %in% double_genes) {
        if (grep(g, double_genes) %% 2 == 1) {
          tmp <- aiTable[sample == s & gene == g][1] 
        }else{
          tmp <- aiTable[sample == s & gene == g][2] 
        }
      }else{
        tmp <- aiTable[sample == s & gene == g] 
      }
      tmp
    }
    ))
  }))
}

get_lesion_type <- function(i, AllelicImb.G, AllelicImb.T, log2R, AllelicImb.thr, log2R.thr, snps.no){
  if(length(AllelicImb.G) != length(AllelicImb.T)){
    stop("Lengths of input vectors differ")
  }
  if(is.na(snps.no[i])){
    return("Not enough SNPs")
  }else if(AllelicImb.G[i] >= AllelicImb.thr){
    return("Germline")
  }else if(is.na(log2R[i])){
    return(NA)
  }else{
    if (AllelicImb.T[i] >= AllelicImb.thr) {
      if(log2R[i] <= log2R.thr[1]){
        return("HemiDel")
      }else if(log2R[i] > log2R.thr[2]){
        return("Gain")
      }else if(log2R[i] <= 0.1 & log2R[i] >= -0.1){
        return("CNNL")
      }else{
        return("GrayZone_Imb")
      }
    }else{
      if(log2R[i] <= log2R.thr[1] & log2R[i] > -1){
        return("Uncertain HomoDel")
      }else if(log2R[i] <= -1){
        return("HomoDel")
        #}else if(log2R[i] > log2R.thr[2]){
      }else if(log2R[i] >= 0.5){
        return("Amplification")
      }else if(log2R[i] >= log2R.thr[2] & log2R[i] < 0.5){
        return("Uncertain Amplification")
      }else if(log2R[i] < 0.1 & log2R[i] > -0.1){
        return("wt")
      }else{
        return("GrayZone_NoImb")
      }
    }
  }
} 

get_lesion_type_v2 <- function(i, AllelicImb.G, AllelicImb.T, log2R.unc, log2R.corr, TC,
                               thr.imbalanceG = 0.5, 
                               thr.imbalanceT= 0.2, 
                               thr.TC = 0.15, 
                               thr.log2R_lowTC = 0.1, 
                               snps.no, use.corr.log2R = F){
  if(is.na(snps.no[i])){
    #return("Not enough SNPs")
    AllelicImb.G[i] <- 0
    AllelicImb.T[i] <- 0
  }
  if (use.corr.log2R) {
    log2R = log2R.corr[i]
  }else{
    log2R = log2R.unc[i]
  }
  if(length(AllelicImb.G) != length(AllelicImb.T)){
    stop("Lengths of input vectors differ")
  }else if(AllelicImb.G[i] >= thr.imbalanceG){
    return("Germline")
  }else if(is.na(log2R)){
    return(NA)
  }else{
    if (AllelicImb.T[i] >= thr.imbalanceT) {
      if(log2R <= -0.5){
        return("HemiDel")
      }else if(log2R > 0.3){
        return("Unb.Gain")
      }else if(log2R <= 0.1 & log2R >= -0.1){
        return("CNNL")
      }else{
        return("Imb_WTGrayZoneCN")
      }
    }else{
      if(TC[i] <= thr.TC | is.na(TC[i])){
        if (use.corr.log2R) {
          if (abs(log2R.unc[i]) < thr.log2R_lowTC) {
            if(log2R <= -0.5){
              return("Likely Loss")
            }else if(log2R > 0.3){
              return("Likely Gain")
            }else if(log2R <= 0.1 & log2R >= -0.1){
              return("wt")
            }else{
              return("NoImb_WTGrayZoneCN")
            }
          }else{
            log2R <- log2R.unc[i]
            if(log2R <= -0.5){
              return("Likely Loss")
            }else if(log2R > 0.3){
              return("Likely Gain")
            }else if(log2R <= 0.1 & log2R >= -0.1){
              return("wt")
            }else{
              return("NoImb_WTGrayZoneCN")
            }
          }
        }else{
          if(log2R <= -0.5){
            return("Likely Loss")
          }else if(log2R > 0.3){
            return("Likely Gain")
          }else if(log2R <= 0.1 & log2R >= -0.1){
            return("wt")
          }else{
            return("NoImb_WTGrayZoneCN")
          }
        }
      }else{
        if(log2R <= -0.5 & log2R > -1){
          return("Uncertain HomoDel")
        }else if(log2R <= -1){
          return("HomoDel")
        }else if(log2R >= 0.5){
          return("Bal.Gain")
        }else if(log2R > 0.3 & log2R < 0.5){
          return("Uncertain Bal.Gain")
        }else if(log2R <= 0.1 & log2R >= -0.1){
          return("wt")
        }else{
          return("NoImb_WTGrayZoneCN")
        }
      }
    }
  }
}

get_lesion_type_v3 <- function(i, AllelicImb.G, AllelicImb.T, log2R.unc, log2R.corr, TC, Gene, Chromosome,
                               thr.imbalanceG = 0.5, 
                               thr.imbalanceT = 0.2, 
                               thr.TC = 0.15, 
                               cutoff.table = NULL,
                               cutoff.significance = 0.005,
                               snps.no, use.corr.log2R = F){
  suppressMessages(library(data.table))
  
  if (is.null(cutoff.table)) {
    cutoff.table <- fread("/shares/CIBIO-Storage/CO/SPICE/downloads/PCF_SELECT/analyses/2020.02/Segmentation/Table_cutoffs.txt")
  }else{
    colnames(cutoff.table) <- c("gene", "type", "p_val", "value")
  }
  
  if(length(AllelicImb.G) != length(AllelicImb.T)){
    stop("Lengths of input vectors differ")
  }
  # if (log2R.unc > abs(cutoff.table[p_val == cutoff.significance & gene == Gene & type == "Gain"]$value)) {
  #   
  # }
  #### Define dynamic thresholds
  # Gains
  thr.gz <- abs(cutoff.table[p_val == cutoff.significance & gene == Gene[i] & type == "Gain"]$value)
  thr.uncbalgain <- thr.likelygain <- 0.3
  thr.balgain <- thr.unbgain <- 0.5
  if(thr.gz >= thr.uncbalgain & thr.gz < thr.balgain){
    thr.uncbalgain <- thr.likelygain <- thr.gz
  }else if (thr.gz >= thr.balgain){
    thr.uncbalgain <- Inf
    thr.balgain <- thr.unbgain <- thr.likelygain <- thr.gz
  }
  
  # Deletions
  thr.hemi <- thr.unchomo <- thr.likelyloss <- -0.5
  thr.homo <- -1
  if(-thr.gz <= thr.unchomo & thr.gz > thr.homo){
    thr.unchomo <- thr.likelyloss <- -thr.gz
  }else if (thr.gz <= thr.homo){
    thr.unchomo <- -Inf
    thr.homo <- thr.likelyloss <- -thr.gz
  }
  ####
  
  if(is.na(snps.no[i])){
    #return("Not enough SNPs")
    AllelicImb.G[i] <- 0
    AllelicImb.T[i] <- 0
  }
  
  ### Apply cutoff on uncorrected log2R
  if (Chromosome[i] == "X") {
    if(log2R.unc[i] >= -thr.gz & log2R.unc[i] <= thr.gz){
      return("wt")
    }
  }else{
    if(AllelicImb.G[i] >= thr.imbalanceG){
      return("Germline")
    }
    if(AllelicImb.T[i] < thr.imbalanceT & log2R.unc[i] >= -thr.gz & log2R.unc[i] <= thr.gz){
      return("wt")
    }
    if(AllelicImb.T[i] >= thr.imbalanceT & log2R.unc[i] >= -thr.gz & log2R.unc[i] <= thr.gz){
      return("CNNL")
    } 
  }
  ###
  
  if (use.corr.log2R) {
    log2R = log2R.corr[i]
  }else{
    log2R = log2R.unc[i]
  }
  
  ### Calls for chrX
  if (Chromosome[i] == "X") {
    if(is.na(log2R)){
      return(NA)
    }else{
      if (TC[i] <= thr.TC | is.na(TC[i])){
        if (log2R >= 0.5) {
          return("Likely Gain")
        }else if (log2R <= -1){
          return("Likely Loss")
        }
      }else{
        if (log2R >= 0.5) {
          return("Bal.Gain")
        }else if (log2R <= -1){
          return("Deletion on chrX")
        }
      }
    }
  }
  ###
  
  
  ### Calls on other chromosomes
  if(is.na(log2R)){
    return(NA)
  }else{
    if (AllelicImb.T[i] >= thr.imbalanceT) {
      if(log2R <= thr.hemi){
        return("HemiDel")
      }else if(log2R > thr.unbgain){
        return("Unb.Gain")
        # }else if(log2R <= thr.gz & log2R >= -thr.gz){
        #   return("CNNL")
      }else{
        return("Imb_WTGrayZoneCN")
      }
    }else{
      if(TC[i] <= thr.TC | is.na(TC[i])){
        # if (use.corr.log2R) {
        #   if (abs(log2R.unc[i]) < thr.log2R_lowTC) {
        #     if(log2R <= -0.5){
        #       return("Likely Loss")
        #     }else if(log2R > 0.3){
        #       return("Likely Gain")
        #     }else if(log2R <= 0.1 & log2R >= -0.1){
        #       return("wt")
        #     }else{
        #       return("NoImb_WTGrayZoneCN")
        #     }
        #   }else{
        #     log2R <- log2R.unc[i]
        if(log2R <= thr.likelyloss){
          return("Likely Loss")
        }else if(log2R > thr.likelygain){
          return("Likely Gain")
          # }else if(log2R <= thr.gz & log2R >= -thr.gz){
          #   return("wt")
        }else{
          return("NoImb_WTGrayZoneCN")
        }
        # }else{
        #   if(log2R <= -0.5){
        #     return("Likely Loss")
        #   }else if(log2R > 0.3){
        #     return("Likely Gain")
        #   }else if(log2R <= 0.1 & log2R >= -0.1){
        #     return("wt")
        #   }else{
        #     return("NoImb_WTGrayZoneCN")
        #   }
        # }
      }else{
        if(log2R <= thr.unchomo & log2R > thr.homo){
          return("Uncertain HomoDel")
        }else if(log2R <= thr.homo){
          return("HomoDel")
        }else if(log2R > thr.uncbalgain & log2R <= thr.balgain){
          return("Uncertain Bal.Gain")
        }else if(log2R > thr.balgain){
          return("Bal.Gain")
          # }else if(log2R <= thr.gz & log2R >= -thr.gz){
          #   return("wt")
        }else{
          return("NoImb_WTGrayZoneCN")
        }
      }
    }
  }
}

get_lesion_type_v4 <- function(i, AllelicImb.G, AllelicImb.T, log2R.unc, log2R.corr, TC, Gene, Chromosome,
                               thr.imbalanceG = 0.5, 
                               thr.imbalanceT = 0.2, 
                               thr.TC = 0.15, 
                               cutoff.table = NULL,
                               cutoff.significance = 0.005,
                               snps.no, use.corr.log2R = F){
  suppressMessages(library(data.table))
  
  if (is.null(cutoff.table)) {
    cutoff.table <- fread("/shares/CIBIO-Storage/CO/SPICE/downloads/PCF_SELECT/analyses/2020.02/Segmentation/Table_cutoffs.txt")
  }else{
    colnames(cutoff.table) <- c("gene", "type", "p_val", "value")
  }
  
  if(length(AllelicImb.G) != length(AllelicImb.T)){
    stop("Lengths of input vectors differ")
  }
  # if (log2R.unc > abs(cutoff.table[p_val == cutoff.significance & gene == Gene & type == "Gain"]$value)) {
  #   
  # }
  #### Define dynamic thresholds
  # Gains
  thr.gz <- abs(cutoff.table[p_val == cutoff.significance & gene == Gene[i] & type == "Gain"]$value)
  thr.uncbalgain <- thr.likelygain <- 0.3
  thr.balgain <- thr.unbgain <- 0.5
  if(thr.gz >= thr.uncbalgain & thr.gz < thr.balgain){
    thr.uncbalgain <- thr.likelygain <- thr.gz
  }else if (thr.gz >= thr.balgain){
    thr.uncbalgain <- Inf
    thr.balgain <- thr.unbgain <- thr.likelygain <- thr.gz
  }
  
  # Deletions
  thr.hemi <- thr.unchomo <- thr.likelyloss <- -0.5
  thr.homo <- -1
  if(-thr.gz <= thr.unchomo & thr.gz > thr.homo){
    thr.unchomo <- thr.likelyloss <- -thr.gz
  }else if (thr.gz <= thr.homo){
    thr.unchomo <- -Inf
    thr.homo <- thr.likelyloss <- -thr.gz
  }
  ####
  
  if(is.na(snps.no[i])){
    #return("Not enough SNPs")
    AllelicImb.G[i] <- 0
    AllelicImb.T[i] <- 0
  }
  
  ### Apply cutoff on uncorrected log2R
  if (Chromosome[i] == "X") {
    if(log2R.unc[i] >= -thr.gz & log2R.unc[i] <= thr.gz){
      return("wt")
    }
  }else{
    if(AllelicImb.G[i] >= thr.imbalanceG){
      return("Germline")
    }
    if(AllelicImb.T[i] < thr.imbalanceT & log2R.unc[i] >= -thr.gz & log2R.unc[i] <= thr.gz){
      return("wt")
    }
    # if(AllelicImb.T[i] >= thr.imbalanceT & log2R.unc[i] >= -thr.gz & log2R.unc[i] <= thr.gz){ ### in v3
    #   return("CNNL")                                                                          ### in v3
    # }                                                                                         ### in v3
  }
  ###
  
  if (use.corr.log2R) {
    log2R = log2R.corr[i]
  }else{
    log2R = log2R.unc[i]
  }
  
  ### Calls for chrX
  if (Chromosome[i] == "X") {
    if(is.na(log2R)){
      return(NA)
    }else{
      if (TC[i] <= thr.TC | is.na(TC[i])){
        if (log2R >= 0.5) {
          return("Likely Gain")
        }else if (log2R <= -1){
          return("Likely Loss")
        }
      }else{
        if (log2R >= 0.5) {
          return("Bal.Gain")
        }else if (log2R <= -1){
          return("Deletion on chrX")
        }
      }
    }
  }
  ###
  
  
  ### Calls on other chromosomes
  if(is.na(log2R)){
    return(NA)
  }else{
    if (AllelicImb.T[i] >= thr.imbalanceT) {
      if(log2R <= thr.hemi){
        return("HemiDel")
      }else if(log2R > thr.unbgain){
        return("Unb.Gain")
      }else if(log2R <= thr.gz & log2R >= -thr.gz){ ### not in v3
        return("CNNL")                              ### not in v3
      }else{
        return("Imb_WTGrayZoneCN")
      }
    }else{
      if(TC[i] <= thr.TC | is.na(TC[i])){
        # if (use.corr.log2R) {
        #   if (abs(log2R.unc[i]) < thr.log2R_lowTC) {
        #     if(log2R <= -0.5){
        #       return("Likely Loss")
        #     }else if(log2R > 0.3){
        #       return("Likely Gain")
        #     }else if(log2R <= 0.1 & log2R >= -0.1){
        #       return("wt")
        #     }else{
        #       return("NoImb_WTGrayZoneCN")
        #     }
        #   }else{
        #     log2R <- log2R.unc[i]
        if(log2R <= thr.likelyloss){
          return("Likely Loss")
        }else if(log2R > thr.likelygain){
          return("Likely Gain")
          # }else if(log2R <= thr.gz & log2R >= -thr.gz){
          #   return("wt")
        }else{
          return("NoImb_WTGrayZoneCN")
        }
        # }else{
        #   if(log2R <= -0.5){
        #     return("Likely Loss")
        #   }else if(log2R > 0.3){
        #     return("Likely Gain")
        #   }else if(log2R <= 0.1 & log2R >= -0.1){
        #     return("wt")
        #   }else{
        #     return("NoImb_WTGrayZoneCN")
        #   }
        # }
      }else{
        if(log2R <= thr.unchomo & log2R > thr.homo){
          return("Uncertain HomoDel")
        }else if(log2R <= thr.homo){
          return("HomoDel")
        }else if(log2R > thr.uncbalgain & log2R <= thr.balgain){
          return("Uncertain Bal.Gain")
        }else if(log2R > thr.balgain){
          return("Bal.Gain")
          # }else if(log2R <= thr.gz & log2R >= -thr.gz){
          #   return("wt")
        }else{
          return("NoImb_WTGrayZoneCN")
        }
      }
    }
  }
}

correct_log2 <- function(log2R, ploidy, tc){
  log2_shift <- -log2(ploidy/2)
  log2_pl_corr <- log2R - log2_shift
  return(log2(pmax((2^(log2_pl_corr) -
                      (1-tc)), 0)/tc)
  )
}

is.log2.cutoff.pass <- function(i, AllelicImbalance.table, cutoff.table = NULL, cutoff.significance = 0.005) {
  if (is.null(cutoff.table)) {
    cutoff.table <- fread("/shares/CIBIO-Storage/CO/SPICE/downloads/PCF_SELECT/analyses/2020.02/Segmentation/Table_cutoffs.txt")
  }else{
    colnames(cutoff.table) <- c("gene", "type", "p_val", "value")
  }
  
  if ((AllelicImbalance.table$all_log2[i] >= cutoff.table[gene == AllelicImbalance.table$gene[i] & p_val == cutoff.significance & type == "Gain"]$value)|(AllelicImbalance.table$all_log2[i] <= cutoff.table[gene == AllelicImbalance.table$gene[i] & p_val == cutoff.significance & type == "Loss"]$value)) {
    T
  }else{
    if (!is.na(AllelicImbalance.table$evidence[i]) & AllelicImbalance.table$evidence[i] >= thr.imbT) {
      T
    }else{
      F
    }
  }
}

compute_asP <- function(aiTable, panel.version = 2){
  size <- fread(paste0("AnnotationFiles/Panel_genes_size_v", panel.version, ".csv"))
  tb <- merge(aiTable, size)
  
  out <- data.frame(sample = unique(tb$sample))
  
  tb <- tb[!(is.na(cnA) | is.na(cnB))]
  to_merge <- rbindlist(lapply(unique(tb$sample), function(s){
    data.frame(sample = s, asP = round(sum((tb[sample == s]$cnA+tb[sample == s]$cnB)*tb[sample == s]$size)/sum(tb[sample == s]$size), digits = 2)) 
  }))
  
  if (nrow(to_merge) == 0) {
    out$asP <- NA
    data.table(out)
  }else{
    data.table(merge(out, to_merge, by = "sample") )
  }
}

adjustVAF <- function (SNVid, SNVtable.ext)
{
  thisSNV <- SNVtable.ext[SNVid, ]
  thisSNV$t_cov <- NA
  thisSNV$t_af <- NA
  thisSNV$n_admReads <- NA
  thisSNV$t_ref_count_corr <- NA
  thisSNV$t_af_corr <- NA
  thisSNV$cn.int <- NA
  thisSNV$CN_SNVmut <- NA
  thisSNV$VAFexp <- NA
  thisSNV$SNV.clonality <- NA
  thisSNV$SNV.clonality.int <- NA
  thisSNV$t_cov <- thisSNV$rc_alt_tumor + thisSNV$rc_ref_tumor
  if (thisSNV$t_cov == 0) {
    return(thisSNV)
  }
  thisSNV$t_af <- thisSNV$rc_alt_tumor/thisSNV$t_cov
  if (is.na(thisSNV$cnA) || is.na(thisSNV$rc_alt_tumor)) {
    return(thisSNV)
  }
  thisSNV$n_admReads <- round((2 * thisSNV$t_cov * thisSNV$adm)/(2 * 
                                                                   thisSNV$adm + (thisSNV$cnA + thisSNV$cnB) * (1 - thisSNV$adm)))
  thisSNV$t_ref_count_corr <- thisSNV$rc_ref_tumor - thisSNV$n_admReads
  thisSNV$t_af_corr <- thisSNV$rc_alt_tumor/(thisSNV$rc_alt_tumor + 
                                               thisSNV$t_ref_count_corr)
  thisSNV$cn.int <- thisSNV$cnA.int + thisSNV$cnB.int
  if (thisSNV$cn.int <= 0) {
    return(thisSNV)
  }
  possibleAF <- seq(1, thisSNV$cn.int, 1)/thisSNV$cn.int
  thisSNV$CN_SNVmut <- seq(1, thisSNV$cn.int, 1)[which.min(abs(possibleAF - 
                                                                 thisSNV$t_af_corr))]
  thisSNV$VAFexp <- (thisSNV$CN_SNVmut * (1 - thisSNV$adm))/(2 * 
                                                               thisSNV$adm + (thisSNV$cnA + thisSNV$cnB) * (1 - thisSNV$adm))
  thisSNV$SNV.clonality = ((thisSNV$t_af - thisSNV$VAFexp) * 
                             (2 * thisSNV$adm + (thisSNV$cnA + thisSNV$cnB) * (1 - 
                                                                                 thisSNV$adm)))/(thisSNV$CN_SNVmut * (1 - thisSNV$adm)) + 
    1
  thisSNV$SNV.clonality.int = ((thisSNV$t_af - thisSNV$VAFexp) * 
                                 (2 * thisSNV$adm + (thisSNV$cnA.int + thisSNV$cnB.int) * 
                                    (1 - thisSNV$adm)))/(thisSNV$CN_SNVmut * (1 - thisSNV$adm)) + 
    1
  return(thisSNV)
}

convertLog2toCN <- function(log2, chr = ""){
  if (!grepl("X", chr)) {
    2*2^log2
  }else{
    2^log2
  }
}
