library(data.table)
library(parallel)

source("PATHS.R")

# SEGMENTATION
sif = fread(path_sif, data.table = F)

pileup.dir = path_dirRC

res = mclapply(1:nrow(sif), function(i)
{
  cat(i,"\n")
  # normal.name = gsub(".bam","",basename(sif$Normal.Bam.Name[i]))
  # tumor.name = gsub(".bam","",basename(sif$Tumor.Bam.Name[i]))
  normal.name = gsub(".snps","",basename(sif$Normal[i]))
  tumor.name = gsub(".snps","",basename(sif$Plasma[i]))
  bed = fread(path_bed,data.table=F)
  # control.genes = fread(path_controlgenes,data.table=F)
  # control.genes = unique(gsub("CONTROL_GENE\\|","",control.genes[grep("CONTROL_GENE",control.genes[,4]),4]))
  control.genes = fread(path_annotations)
  control.genes = control.genes[PCF_SELECT == "control_gene" & chromosome_name != "X"]$hgnc_symbol
  normal.rc = fread(paste(pileup.dir,normal.name,".rc",sep=""),data.table=F)
  tumor.rc = fread(paste(pileup.dir,tumor.name,".rc",sep=""),data.table=F)

  ## GC content normalization
  normal.rc$gc.rounded = (round((normal.rc$gc*10)*2)/2)/10
  tumor.rc$gc.rounded = (round((tumor.rc$gc*10)*2)/2)/10
  normal.rc$rc = sapply(1:nrow(normal.rc),function(x) (normal.rc$rc[x]*mean(normal.rc$rc))/mean(normal.rc$rc[which(normal.rc$gc.rounded==normal.rc$gc.rounded[x])]))
  tumor.rc$rc = sapply(1:nrow(tumor.rc),function(x) (tumor.rc$rc[x]*mean(tumor.rc$rc))/mean(tumor.rc$rc[which(tumor.rc$gc.rounded==tumor.rc$gc.rounded[x])]))

  normal.rc.X = normal.rc[which(bed[,1]=="chrX"),]
  normal.rc = normal.rc[which(bed[,1]!="chrX"),]
  tumor.rc.X = tumor.rc[which(bed[,1]=="chrX"),]
  tumor.rc = tumor.rc[which(bed[,1]!="chrX"),]
  bed.X = bed[which(bed[,1]=="chrX"),]
  bed = bed[which(bed[,1]!="chrX"),]
  genes.X = unique(bed.X[,4])
  genes = unique(bed[,4])
  genes.X = unique(strsplit(paste(genes.X,collapse="|"),"\\|")[[1]])
  genes = unique(strsplit(paste(genes,collapse="|"),"\\|")[[1]])
  all(normal.rc.X$from==bed.X[,2]+1)
  all(normal.rc.X$to==bed.X[,3])
  all(normal.rc$from==bed[,2]+1)
  all(normal.rc$to==bed[,3])

  # control.rc.normal = median(sapply(control.genes,function(x) median(normal.rc$rc[which(sapply(bed[,4],function(z) x%in%strsplit(z,"\\|")[[1]]))])),na.rm = T)
  # control.rc.tumor = median(sapply(control.genes,function(x) median(tumor.rc$rc[which(sapply(bed[,4],function(z) x%in%strsplit(z,"\\|")[[1]]))])),na.rm = T)
  # control.rc.normal.X = median(sapply(c("TEX11","HDAC8"),function(x) median(normal.rc.X$rc[which(bed.X[,4]==x)])),na.rm = T)
  # control.rc.tumor.X = median(sapply(c("TEX11","HDAC8"),function(x) median(tumor.rc.X$rc[which(bed.X[,4]==x)])),na.rm = T)
  control.rc.normal = mean(normal.rc$rc, na.rm = T)
  control.rc.tumor = mean(tumor.rc$rc, na.rm = T)

  normal.genes.rc = lapply(sapply(genes,function(x) normal.rc$rc[which(sapply(bed[,4],function(z) x%in%strsplit(z,"\\|")[[1]]))]),function(y) y/control.rc.normal)
  tumor.genes.rc = lapply(sapply(genes,function(x) tumor.rc$rc[which(sapply(bed[,4],function(z) x%in%strsplit(z,"\\|")[[1]]))]),function(y) y/control.rc.tumor)
  genes.log2ratio = log2(sapply(1:length(normal.genes.rc),function(x) median(tumor.genes.rc[[x]]/normal.genes.rc[[x]],na.rm = T)))
  # hist(genes.log2ratio,breaks=39,xlim=c(-2,2))

  normal.genes.rc.X = lapply(sapply(genes.X,function(x) normal.rc.X$rc[which(bed.X[,4]==x)]),function(y) y/control.rc.normal)
  tumor.genes.rc.X = lapply(sapply(genes.X,function(x) tumor.rc.X$rc[which(bed.X[,4]==x)]),function(y) y/control.rc.tumor)
  genes.log2ratio.X = log2(sapply(1:length(normal.genes.rc.X),function(x) median(tumor.genes.rc.X[[x]]/normal.genes.rc.X[[x]],na.rm = T)))
  # hist(genes.log2ratio.X,breaks=39,xlim=c(-2,2))

  if (store_amplicons_log2R) {
    dir.create(path_amplicons_log2r)

    tmp = sapply(1:length(normal.genes.rc),function(x) {tumor.genes.rc[[x]]/normal.genes.rc[[x]]})
    names(tmp) = genes
    out <- rbindlist(lapply(1:length(tmp), function(x){
      data.frame(gene_symbol = names(tmp)[x], log2R = log2(tmp[[x]]), sample = tumor.name)
    }))
    out$chr <- normal.rc$chr
    out$from <- normal.rc$from
    out$to <- normal.rc$to
    out <- out[,c(4:6, 1:3)]

    tmp.X = sapply(1:length(normal.genes.rc.X),function(x) {tumor.genes.rc.X[[x]]/normal.genes.rc.X[[x]]})
    names(tmp.X) = genes.X
    out.X <- rbindlist(lapply(1:length(tmp.X), function(x){
      data.frame(gene_symbol = names(tmp.X)[x], log2R = log2(tmp.X[[x]]), sample = tumor.name)
    }))
    out.X$chr <- normal.rc.X$chr
    out.X$from <- normal.rc.X$from
    out.X$to <- normal.rc.X$to
    out.X <- out.X[,c(4:6, 1:3)]

    out_amp_log2r <- as.data.table(rbind(out, out.X))
    save(out_amp_log2r, file = paste0(path_amplicons_log2r, tumor.name, ".RData"), version = 2)
  }

  seg.file = t(sapply(1:length(genes),function(x)
    c(tumor.name,bed[which(sapply(bed[,4],function(z) genes[x]%in%strsplit(z,"\\|")[[1]]))[1],1],
      min(bed[which(sapply(bed[,4],function(z) genes[x]%in%strsplit(z,"\\|")[[1]])),2]),
      max(bed[which(sapply(bed[,4],function(z) genes[x]%in%strsplit(z,"\\|")[[1]])),2]),NA,genes.log2ratio[x])))
  seg.file.X = t(sapply(1:length(genes.X),function(x)
    c(tumor.name,bed.X[which(bed.X[,4]==genes.X[x])[1],1],min(bed.X[which(bed.X[,4]==genes.X[x]),2]),max(bed.X[which(bed.X[,4]==genes.X[x]),2]),NA,genes.log2ratio.X[x])))

  seg.file = rbind(seg.file,seg.file.X)
  colnames(seg.file) = c("sample","chr","start","end","n.markers","log2ratio")
  seg.file[,2] = gsub("chr","",seg.file[,2])

  seg.file = as.data.table(seg.file)

  return(seg.file)
}, mc.cores = cores.readdepth)

seg.file = rbindlist(res)
fwrite(seg.file, file = path_segfile, sep="\t", row.names=F, quote=F)

if (store_amplicons_log2R) {
  log2r_amplicons <- rbindlist(lapply(list.files(path_amplicons_log2r, full.names = T), function(x){
    load(x)
    out_amp_log2r
  }))
  fwrite(log2r_amplicons, file = paste0(path_amplicons_log2r, "Segmentation_amplicons_log2r.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
}

## Beta Computation
source(path_betafunctions)

#load("/elaborazioni/sharedCO/PCF_SELECT/Analysis/BetaComputation/Germ_Distr_UCL.RData")
load(path_germ_distr)
#load("/elaborazioni/sharedCO/PCF_SELECT/Analysis/BetaComputation/Germ_Distr_Cornell.RData")
tmp = germ.distr[[1]]

tmp = tmp[which(tmp$af.mean>0.35&tmp$af.mean<0.65),]
tmp = tmp[which((tmp$af.cv<=0.1&tmp$af.alt>1)),]
tmp = tmp[which(tmp$af.freq<0.8),]
germ.distr[[1]] = tmp

coords = fread(path_panelgenescoords,data.table=F)
coords = coords[which(coords$chr!="X"),]

pshift <- fread(path_shift_table, header = T)

res = mclapply(1:nrow(sif),function(i, shift = pshift)
{
  cat(sif$Plasma[i],"\n")
  # Normal
  p = fread(sif$Normal[i])
  p = p[which(p$af<0.8 & p$af>0.2),]
  
  ##### Apply shifting and record rsids to discard
  p$af <- p$af+(0.5-shift$germ_shift[i])
  
  plas <- fread(sif$Plasma[i])
  plas$af <- plas$af+(0.5-shift$plas_shift[i])
  rsids_to_exclude <- plas[af < 0 | af > 1, get("rsid")]
  if (length(rsids_to_exclude) != 0) {
    p <- p[!rsid %in% rsids_to_exclude]
  }
  #####
  
  het.snps = p$rsid
  list.beta = list()
  for(j in 1:nrow(coords))
  {
    afs = p$af[which(p$chr==coords$chr[j]&p$pos>=coords$start[j]&p$pos<=coords$end[j])]
    rsid = p$rsid[which(p$chr==coords$chr[j]&p$pos>=coords$start[j]&p$pos<=coords$end[j])]
    cov = p$cov[which(p$chr==coords$chr[j]&p$pos>=coords$start[j]&p$pos<=coords$end[j])]
    afs[which(afs<0.5)] = 1-afs[which(afs<0.5)]
    
    list.beta[[length(list.beta)+1]] <- computeBeta(cov = cov, afs = afs, rsid = rsid, germ.distr = germ.distr, p.thr = 0.01, paired = T, times = ai.nrep)
  }
  
  normal.tab = matrix(unlist(list.beta),ncol=length(list.beta[[1]]),byrow = T)
  
  # Plasma
  p = fread(sif$Plasma[i])
  p$af = p$af+(0.5-shift$plas_shift[i])
  p = p[which(p$rsid%in%het.snps),]
  list.beta = list()
  for(j in 1:nrow(coords))
  {
    afs = p$af[which(p$chr==coords$chr[j]&p$pos>=coords$start[j]&p$pos<=coords$end[j])]
    rsid = p$rsid[which(p$chr==coords$chr[j]&p$pos>=coords$start[j]&p$pos<=coords$end[j])]
    cov = p$cov[which(p$chr==coords$chr[j]&p$pos>=coords$start[j]&p$pos<=coords$end[j])]
    afs[which(afs<0.5)] = 1-afs[which(afs<0.5)]
    
    list.beta[[length(list.beta)+1]] <- computeBeta(cov = cov, afs = afs, rsid = rsid, germ.distr = germ.distr, p.thr = 0.01, paired = T, times = ai.nrep)
  }
  
  plasma.tab = matrix(unlist(list.beta),ncol=length(list.beta[[1]]),byrow = T)
  
  tab = cbind(gsub(".snps","",basename(sif$Plasma[i])),
              coords[,1],coords[,2],coords[,3],coords[,4],
              plasma.tab,normal.tab)
  colnames(tab) = 
    c("sample","chr","start","end","gene",
      "beta","beta.max","beta.min","snps","evidence","cov",
      "beta.n","beta.n.max","beta.n.min","snps.n","evidence.n","cov.n")
  
  return(tab)
  
},mc.cores=cores.betacomputation)

beta.table = do.call(rbind,res)
save(beta.table, file=path_beta_table,compress=T, version = 2)
#save(beta.table, file="/SPICE/downloads/PCF_SELECT/analyses/2019.12/in_vitro_dilutions-ucl/betaTable_modelCornell.RData",compress=T)
