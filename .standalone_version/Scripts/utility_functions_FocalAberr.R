library(data.table)
library(fitdistrplus)

# summary gene: median log2 computation -> grep instead of which == ..
annotate_bed_with_rc <- function(gene, window.length = 10000, bed, t.rc, n.rc){
  # if (gene == "ERG_TMPRSS2") {
  #   return(NULL)
  # }
  # if (gene == "RUNX1"){ # No flanking regions for RUNX1
  #   return(NULL)
  # }
  #cat(gene, "\n")
  ## Select subset of annotated bed covering the gene considered
  if(gene == "MYC") { # Also MYCN exists
    tmp_bed <- bed[which(bed$V4 == gene),]
  }else{
    tmp_bed <- bed[grep(gene, bed$V4),]
  }
  
  ## Add log2 values to subset of bed covering the gene considered
  ind <- which(names(t.rc) == gene)
  tmp_bed$log2 = log2(t.rc[[ind]]/n.rc[[ind]])
  tmp_bed$t.rc = t.rc[[ind]]
  tmp_bed$n.rc = n.rc[[ind]]
  tmp_bed$is.focal = ifelse(grepl(gene, tmp_bed$V6), "gene", "flanking")
  if (gene == "ERG_TMPRSS2") {
    tmp_bed[V6 == "ERG" | V6 == "TMPRSS2"]$is.focal = "gene"
  }
  tmp_bed$amplicon = as.character(1:nrow(tmp_bed))
  return(tmp_bed)
}

compute_log2_windows <- function(gene, window.length = 10000, bed, t.rc, n.rc){
  if (gene == "ERG_TMPRSS2") {
    return(NULL)
  }
  if (gene == "RUNX1"){ # No flanking regions for RUNX1
    return(NULL)
  }
  #cat(gene, "\n")
  ## Select subset of annotated bed covering the gene considered
  if(gene == "MYC") { # Also MYCN exists
    tmp_bed <- bed[which(bed$V4 == gene),]
  }else{
    tmp_bed <- bed[grep(gene, bed$V4),]
  }

  ## Add log2 values to subset of bed covering the gene considered
  ind <- which(names(t.rc) == gene)
  tmp_bed$log2 = log2(t.rc[[ind]]/n.rc[[ind]])
  tmp_bed$t.rc = t.rc[[ind]]
  tmp_bed$n.rc = n.rc[[ind]]

  ## Extract limits of the region covered exactly by the gene considered
  right_lim_gene <- tail(tmp_bed[grep(gene, V6), get("V3")], n = 1)
  left_lim_gene <- head(tmp_bed[grep(gene, V6), get("V2")], n = 1)

  size_amp <- sum(tmp_bed[grep(gene, V6), get("V3")]-tmp_bed[grep(gene, V6), get("V2")])

  summary_gene <- data.frame(left_lim = left_lim_gene, right_lim = right_lim_gene, size_region = right_lim_gene-left_lim_gene,
                             size_amp = size_amp,
                             n.amplicons = nrow(tmp_bed[which(tmp_bed$V3 >= left_lim_gene)[1]:tail(which(tmp_bed$V2 <= right_lim_gene), n = 1)]),
                             #median.log2 = median(tmp_bed[which(tmp_bed$V6 == gene),]$log2, na.rm = T),
                             median.log2 = median(tmp_bed[grep(gene, V6),]$log2, na.rm = T),
                             median.log2_left = NaN, median.log2_right = NaN)

  ## Iteratively increase window size for amplicon selection and log2 computation
  this_right_lim <- right_lim_gene
  this_left_lim <- left_lim_gene
  while (this_right_lim <= tail(tmp_bed$V3, n = 1) & this_left_lim >= tmp_bed$V2[1]) {
    # Update limits
    this_right_lim <- this_right_lim + window.length
    this_left_lim <- this_left_lim - window.length

    # Compute log2 and statistics adding flanking regions to gene
    ## Statistics
    l_lim <- tmp_bed[which(tmp_bed$V3 >= this_left_lim)[1], get("V2")]
    r_lim <- tmp_bed[tail(which(tmp_bed$V2 <= this_right_lim), n = 1), get("V3")]
    size_region <- r_lim - l_lim

    amp_bed_tmp <- tmp_bed[which(tmp_bed$V3 >= this_left_lim)[1]:tail(which(tmp_bed$V2 <= this_right_lim), n = 1)]
    size_amp <- sum(amp_bed_tmp$V3-amp_bed_tmp$V2)

    n.amplicons <- nrow(tmp_bed[which(tmp_bed$V3 >= this_left_lim)[1]:tail(which(tmp_bed$V2 <= this_right_lim), n = 1)])
    median.log2_left <- median(tmp_bed[which(tmp_bed$V3 >= this_left_lim)[1]:tail(which(tmp_bed$V2 <= left_lim_gene), n = 1), get("log2")], na.rm = T)
    median.log2_right <- median(tmp_bed[which(tmp_bed$V3 >= right_lim_gene)[1]:tail(which(tmp_bed$V2 <= this_right_lim), n = 1), get("log2")], na.rm = T)
    median.log2 <- median(tmp_bed[which(tmp_bed$V3 >= this_left_lim)[1]:tail(which(tmp_bed$V2 <= this_right_lim), n = 1), get("log2")], na.rm = T)

    summary_gene[nrow(summary_gene) + 1,] <- list(l_lim, r_lim, size_region, size_amp, n.amplicons, median.log2, median.log2_left, median.log2_right)
  }

  summary_gene$gene <- gene
  return(summary_gene)
}

is.focal.aberration <- function(f_control, f_target, n.rep = 10000){
  #diff_c <- unlist(sapply(unique(f_control$gene), function (x) abs(f_control[gene == x, get("median.log2")][1])-abs(tail(f_control[gene == x, get("median.log2")], n = 1)), simplify = T))
  diff_c <- unlist(sapply(unique(f_control$gene), function (x) abs(f_control[gene == x, get("median.log2")][1]-tail(f_control[gene == x, get("median.log2")], n = 1)), simplify = T))
  diff_c <- diff_c[!is.na(diff_c)]
  #diff_t <- unlist(sapply(unique(f_target$gene), function (x) abs(f_target[gene == x, get("median.log2")][1])-abs(tail(f_target[gene == x, get("median.log2")], n = 1)), simplify = T))
  diff_t <- unlist(sapply(unique(f_target$gene), function (x) abs(f_target[gene == x, get("median.log2")][1]-tail(f_target[gene == x, get("median.log2")], n = 1)), simplify = T))
  f <- fitdist(diff_c, distr = "norm")
  c.mean <- f$estimate[1]
  c.sd <- f$estimate[2]
  gene <- unlist(sapply(1:n.rep, function(i) {sim_c_diff <- rnorm(10000, mean = c.mean, sd = c.sd); names(diff_t[which(diff_t > max(sim_c_diff) | diff_t < min(sim_c_diff))])}))
  if(length(gene) == 0){
    return(data.table(gene = as.character(), N = as.numeric()))
  }else{
    out <- table(gene)/n.rep
    return(data.table(out))
  }
}

is.focal.aberration_strata <- function(f_control, f_target, n.rep = 10000, 
                                       background_control, background_target, stratify_by = 0.1){ # new parameters for stratification
  
  upper_limit_stratification <- ceiling(max(c(background_control$Median, background_target$Median), na.rm = T)*10)/10
  strata <- c(seq(0, upper_limit_stratification, by = stratify_by), Inf)
  
  diff_c <- unlist(sapply(unique(f_control$gene), function (x) abs(f_control[gene == x, get("median.log2")][1]-tail(f_control[gene == x, get("median.log2")], n = 1)), simplify = T))
  diff_c <- diff_c[!is.na(diff_c)]
  diff_t <- unlist(sapply(unique(f_target$gene), function (x) abs(f_target[gene == x, get("median.log2")][1]-tail(f_target[gene == x, get("median.log2")], n = 1)), simplify = T))
  gene <- c()
  for (i in 1:(length(strata)-1)) {
    chosen_controls <- background_control[abs(Median) > strata[i] & abs(Median) <= strata[i+1]]$gene
    chosen_targets <- background_target[abs(Median) > strata[i] & abs(Median) <= strata[i+1]]$gene
    x_diff_c <- diff_c[names(diff_c) %in% chosen_controls]
    x_diff_t <- diff_t[names(diff_t) %in% chosen_targets]
    if(length(x_diff_c) != 0 & length(x_diff_t) != 0){
      f <- fitdist(x_diff_c, distr = "norm")
      c.mean <- f$estimate[1]
      c.sd <- f$estimate[2]
      gene <- c(gene, unlist(sapply(1:n.rep, function(i) {sim_c_diff <- rnorm(10000, mean = c.mean, sd = c.sd); names(x_diff_t[which(x_diff_t > max(sim_c_diff) | x_diff_t < min(sim_c_diff))])}))) 
    }
  }
  if(length(gene) == 0){
    return(data.table(gene = as.character(), N = as.numeric()))
  }else{
    out <- table(gene)/n.rep
    return(data.table(out))
  }
}

check_continuity <- function(log2_vector, focal_type, flank_length = 2){
  if (focal_type == "gain") {
    # extract absolute maxima
    ind <- which.max(log2_vector)
    # find all local maxima
    loc <- which(diff(sign(diff(log2_vector)))==-2)+1
    # check for other local peaks inside the flank length
    if (ind == 1) {
      length(intersect(loc, ind:(ind+flank_length))) == 0
    }else{
      length(intersect(loc, (ind-flank_length):(ind+flank_length))) == 1
    }
  }else{
    # extract absolute minima
    ind <- which.min(log2_vector)
    # find all local minima
    loc <- which(diff(sign(diff(log2_vector)))==2)+1
    # check for other local peaks inside the flank length
    if (ind == 1) {
      length(intersect(loc, ind:(ind+flank_length))) == 0
    }else{
      length(intersect(loc, (ind-flank_length):(ind+flank_length))) == 1
    }
  }
}

compute_focal_region <- function(g, focal_target_genes){
  tmp <- focal_target_genes[gene == g]
  if(length(tail(tmp[gene == g, get("median.log2")], n = 1)) == 0){
    return(data.table(start.focal = NA, end.focal = NA, log2ratio.focal = NA))
  }
  if(tail(tmp[gene == g, get("median.log2")], n = 1) > tmp[gene == g, get("median.log2")][1]){ # for Focal Loss
    focal_log2 <- tmp$median.log2[1]
    m_left <- tmp$median.log2_left[which.min(tmp$median.log2_left)]
    m_right <- tmp$median.log2_right[which.min(tmp$median.log2_right)]
    
    ### add check for continuity
    if (!check_continuity(tmp$median.log2, "loss", flank_length = 2)) {
      focal_log2 <- 0
    }
    if (!check_continuity(tmp$median.log2_left, "loss", flank_length = 2)) {
      m_left <- 0
    }
    if (!check_continuity(tmp$median.log2_right, "loss", flank_length = 2)) {
      m_right <- 0
    }
    # if there is no continuity take the exact focal region
    if (focal_log2 == 0 & m_left == 0 & m_right == 0) {
      focal_log2 <- tmp$median.log2[1]
    }
    ###
    
    if (which.min(c(focal_log2, m_left, m_right)) == 1) {
      tmp <- tmp[1, c("left_lim", "right_lim", "median.log2")]
      colnames(tmp) <- c("start.focal", "end.focal", "log2ratio.focal")
      return(tmp)
    }else if(which.min(c(focal_log2, m_left, m_right)) == 2){ # Log2 lower considering left-hand flanking regions
      l_lim <- tmp[which.min(tmp$median.log2_left), get("left_lim")]
      r_lim <- tmp[1, get("right_lim")]
      median.log2 <- tmp[which.min(tmp$median.log2_left), get("median.log2_left")]
      return(data.table(start.focal = l_lim, end.focal = r_lim, log2ratio.focal = median.log2))
    }else{ # Log2 lower considering right-hand flanking regions
      l_lim <- tmp[1, get("left_lim")]
      r_lim <- tmp[which.min(tmp$median.log2_right), get("right_lim")]
      median.log2 <- tmp[which.min(tmp$median.log2_right), get("median.log2_right")]
      return(data.table(start.focal = l_lim, end.focal = r_lim, log2ratio.focal = median.log2))
    }
  }else{ # for Focal Gain
    focal_log2 <- tmp$median.log2[1]
    m_left <- tmp$median.log2_left[which.max(tmp$median.log2_left)]
    m_right <- tmp$median.log2_right[which.max(tmp$median.log2_right)]
    
    ### add check for continuity
    if (!check_continuity(tmp$median.log2, "gain", flank_length = 2)) {
      focal_log2 <- 0
    }
    if (!check_continuity(tmp$median.log2_left, "gain", flank_length = 2)) {
      m_left <- 0
    }
    if (!check_continuity(tmp$median.log2_right, "gain", flank_length = 2)) {
      m_right <- 0
    }
    # if there is no continuity take the exact focal region
    if (focal_log2 == 0 & m_left == 0 & m_right == 0) {
      focal_log2 <- tmp$median.log2[1]
    }
    ###
    
    if (which.max(c(focal_log2, m_left, m_right)) == 1) {
      tmp <- tmp[1, c("left_lim", "right_lim", "median.log2")]
      colnames(tmp) <- c("start.focal", "end.focal", "log2ratio.focal")
      return(tmp)
    }else if(which.max(c(focal_log2, m_left, m_right)) == 2){ # Log2 higher considering left-hand flanking regions
      l_lim <- tmp[which.max(tmp$median.log2_left), get("left_lim")]
      r_lim <- tmp[1, get("right_lim")]
      median.log2 <- tmp[which.max(tmp$median.log2_left), get("median.log2_left")]
      return(data.table(start.focal = l_lim, end.focal = r_lim, log2ratio.focal = median.log2))
    }else{ # Log2 higher considering left-hand flanking regions
      l_lim <- tmp[1, get("left_lim")]
      r_lim <- tmp[which.max(tmp$median.log2_right), get("right_lim")]
      median.log2 <- tmp[which.max(tmp$median.log2_right), get("median.log2_right")]
      return(data.table(start.focal = l_lim, end.focal = r_lim, log2ratio.focal = median.log2))
    }
  }
}

# load("Analyses/FocalLog2_Cornell_blacklist_test/focal_tables/55_cfDNA3_target.RData")
# s <- "110_cfDNA3"
# load(paste0("Analyses/FocalLog2_Cornell_blacklist_test/focal_tables/", s, "_target.RData"))
# seg <- fread("Analyses/allele_specific_call_Cornell/segmentation.seg")
# seg[sample == s & prob.focal > 0]
# focal_prob_table <- seg[sample == s & prob.focal > 0, c("gene", "prob.focal")]
# colnames(focal_prob_table) <- c("gene", "N")
# tmp <- focal_target_genes[gene == "CDKN2A"]
# tmp <- focal_target_genes[gene == "TP53"]

annotate_seg_focal <- function(seg, focal_prob_table, focal_target_genes){
  if (F %in% is.na(unique(focal_prob_table$N))) {
    out <- rbindlist(lapply(focal_prob_table$gene, FUN = compute_focal_region, focal_target_genes = focal_target_genes))
    out$gene <- focal_prob_table$gene
    out$prob.focal <- focal_prob_table$N
    seg <- merge(seg, out, all.x = T, by = "gene")
    seg <- seg[, c("sample", "chr", "start", "end", "n.markers", "log2ratio", "gene", 
                   "start.focal", "end.focal", "log2ratio.focal", "prob.focal")]
    return(seg)
  }else{
    seg$start.focal <- NA
    seg$end.focal <- NA
    seg$log2ratio.focal <- NA
    seg$prob.focal <- NA
    seg <- seg[, c("sample", "chr", "start", "end", "n.markers", "log2ratio", "gene", 
                   "start.focal", "end.focal", "log2ratio.focal", "prob.focal")]
    return(seg)
  }
}
