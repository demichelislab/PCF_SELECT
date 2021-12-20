library(data.table)
library(parallel)
library(dplyr)

source("PATHS.R")

sif <- fread(path_sif)
# sif$Normal <- gsub("/shares/CIBIO-Storage/CO", "", sif$Normal)
# sif$Plasma <- gsub("/shares/CIBIO-Storage/CO", "", sif$Plasma)

localMaxima <- function(x) {
  y <- diff(c(-Inf, x)) > 0L
  y <- cumsum(rle(y)$lengths)
  y <- y[seq(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

extract_peaks <- function(dens){
  ydens <- dens$y
  xdens <- dens$x
  peaks <- xdens[localMaxima(ydens)]
  return(peaks)
}

findPeak <- function (AF, bw, RMB=0.47) {
  dens_cha <- stats::density(AF, bw = bw, na.rm = T)
  p_cha <- extract_peaks(dens_cha)
  return(p_cha[which(abs(p_cha-RMB) == min(abs(p_cha-RMB)))])
}

out <- mclapply(1:nrow(sif), function(i) {
  germ <- fread(sif$Normal[i])
  germ <- germ[af >= 0.2 & af <= 0.8]
  germ_shift <- findPeak(germ[,get("af")], bw = "SJ", RMB = 0.47)
  
  plas <- fread(sif$Plasma[i])
  plas <- plas[rsid %in% germ$rsid]
  plas_shift <- findPeak(plas[,get("af")], bw = "SJ", RMB = 0.47)
  
  #s_name <- strsplit(sif$Plasma[i], split = "/", fixed = 7)[[1]][7]
  s_name <- basename(sif$Plasma[i])
  s_name <- gsub(".snps", "", s_name)
  data.frame("sample" = s_name, "germ_shift" = germ_shift, "plas_shift" = plas_shift, "germ_median" = median(germ$af), "plas_median" = median(plas$af))
}, mc.cores = cores.peakcorrection)

peaks_shifts <- rbindlist(out)
write.table(peaks_shifts, file = path_shift_table, sep = "\t", quote = F, col.names = T,
            row.names = F)

if(compute_GermlineModel){
  shift <- fread(path_shift_table)
  
  computeGermlineDistributions <- function(germline,label="germline",snps.list=c(),shift)
  {
    snps = c()
    for(i in 1:length(germline))
    {
      cat(germline[i],"\n")
      geno = fread(germline[i], data.table = F)
      geno = geno[which(geno$af > 0.2 & geno$af < 0.8),]
      # Apply peak correction
      geno$af = geno$af+(0.5-shift$germ_shift[i])
      #geno$af = geno$af+(0.5-median(geno$af))
      geno = geno[which(!geno$rsid%in%snps.list),]
      snps = union(snps,paste(geno$chr,geno$pos,geno$rsid,sep="-"))
      stats = geno[,c("rsid","af","cov")]
      if(i==1)
      {
        stats.tot = stats
      } else
      {
        stats.tot = full_join(stats.tot,stats,by="rsid")
      }
      colnames(stats.tot)[2*i] = paste("af",germline[i],sep="-")
      colnames(stats.tot)[2*i+1] = paste("cov",germline[i],sep="-")
    }
    
    snps = data.frame(chr=sapply(snps,function(x) strsplit(x,"\\-")[[1]][1]),
                      pos=sapply(snps,function(x) strsplit(x,"\\-")[[1]][2]),
                      rsid=sapply(snps,function(x) strsplit(x,"\\-")[[1]][3]),stringsAsFactors = FALSE)
    isort = match(stats.tot$rsid,snps$rsid)
    snps = snps[isort,]
    cat(all(snps$rsid==stats.tot$rsid),"\n")
    
    af = unlist(stats.tot[,seq(2,ncol(stats.tot),by=2)])
    cov = unlist(stats.tot[,seq(3,ncol(stats.tot),by=2)])
    
    cov.inter = quantile(cov,prob=seq(0,1,0.1),na.rm = T)
    cov.inter[1] = 0
    af.sd.cov = sapply(1:(length(cov.inter)-1),function(i) sd(af[which(cov>=cov.inter[i]&cov<cov.inter[i+1])],na.rm = T))
    rsid.af.stats = data.frame(rsid = stats.tot$rsid, chr = snps$chr, pos = as.numeric(snps$pos),
                               af.mean = apply(stats.tot[,seq(2,ncol(stats.tot),by=2)],1,function(x) mean(x,na.rm = T)),
                               af.cv = apply(stats.tot[,seq(2,ncol(stats.tot),by=2)],1,function(x) sd(x,na.rm = T))/
                                 apply(stats.tot[,seq(2,ncol(stats.tot),by=2)],1,function(x) mean(x,na.rm = T)),
                               af.freq = apply(stats.tot[,seq(2,ncol(stats.tot),by=2)],1,function(x) length(which(!is.na(x)))/length(x)),
                               af.alt = apply(stats.tot[,seq(2,ncol(stats.tot),by=2)],1,function(x) length(which(!is.na(x)))),
                               cov.mean = apply(stats.tot[,seq(3,ncol(stats.tot),by=2)],1,function(x) mean(x,na.rm = T)),
                               stringsAsFactors = FALSE)
    
    return(list(rsid.af.stats,af.sd.cov,cov.inter))
  }
  
  germ.distr <- computeGermlineDistributions(germline = unique(sif$Normal), label = "germline", shift = shift)
  save(germ.distr, file = path_germ_distr, version = 2)
}
