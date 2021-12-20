# Adapted from /SPICE/downloads/PCF_SELECT/analyses/2020.02/Segmentation/Cornell/20200221_Segmentation.R

library(data.table)
library(parallel)
library(fitdistrplus)

source("PATHS.R")

# SEGMENTATION
sif = fread(path_sif, data.table=F)
pileup.dir <- path_dirRC
#
source("Scripts/utility_functions_FocalAberr.R")

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
  gl = fread(path_annotations)
  control.genes = gl[PCF_SELECT == "control_gene" & chromosome_name != "X"]$hgnc_symbol
  normal.rc = fread(paste(pileup.dir,normal.name,".rc",sep=""),data.table=F)
  other.genes = gl[PCF_SELECT == "other" & chromosome_name != "X"]$hgnc_symbol
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
  
  focal = T
  bed.focal <- fread(path_bed_focal)
  
  # x <- rbindlist(lapply(genes, annotate_bed_with_rc, bed = bed.focal,
  #                       window.length = 10000, t.rc = tumor.genes.rc, n.rc = normal.genes.rc))
  # x$sample <- tumor.name
  
  if (focal == T) {
    dir.create(paste0("Analyses/", run_name, "/.focal_tables"))
    
    seg.file$gene <- c(genes, genes.X)
    focal_target_genes <- rbindlist(lapply(setdiff(genes, control.genes), compute_log2_windows, bed = bed.focal,
                                           window.length = 10000, t.rc = tumor.genes.rc, n.rc = normal.genes.rc))
    focal_target_genes[sapply(focal_target_genes, is.infinite)] <- 0
    focal_control_genes <- rbindlist(lapply(control.genes, compute_log2_windows, bed = bed.focal,
                                            window.length = 10000, t.rc = tumor.genes.rc, n.rc = normal.genes.rc))
    focal_control_genes[sapply(focal_control_genes, is.infinite)] <- 0
    focal_prob_table <- is.focal.aberration(f_control = focal_control_genes, f_target = focal_target_genes, n.rep = 10000)
    focal_prob_table <- as.data.table(rbind(focal_prob_table, data.table(gene = setdiff(genes, c(control.genes, other.genes, focal_prob_table$gene)),
                                                                         N = NA)))
    seg.file <- annotate_seg_focal(seg.file, focal_prob_table, focal_target_genes)

    focal_control_genes$sample <- tumor.name
    focal_target_genes$sample <- tumor.name
    save(focal_target_genes, file = paste0("Analyses/", run_name, "/.focal_tables/",
                                           tumor.name, "_target.RData"), version = 2)
    save(focal_control_genes, file = paste0("Analyses/", run_name, "/.focal_tables/",
                                           tumor.name, "_control.RData"), version = 2)
  }
  
  return(as.data.table(seg.file))
}, mc.cores = cores.readdepth)

seg.file = rbindlist(res)
fwrite(seg.file, file = paste0("Analyses/", run_name, "/Segmentation_focal.seg"), sep="\t", row.names=F, quote=F)

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
  gl = fread(path_annotations)
  control.genes = gl[PCF_SELECT == "control_gene" & chromosome_name != "X"]$hgnc_symbol
  normal.rc = fread(paste(pileup.dir,normal.name,".rc",sep=""),data.table=F)
  other.genes = gl[PCF_SELECT == "other" & chromosome_name != "X"]$hgnc_symbol
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

  bed.focal <- fread(path_bed_focal)

  x <- rbindlist(lapply(genes, annotate_bed_with_rc, bed = bed.focal,
                        window.length = 10000, t.rc = tumor.genes.rc, n.rc = normal.genes.rc))
  x$sample <- tumor.name
  
  return(x)
}, mc.cores = cores.readdepth)

bed.rc = rbindlist(res)
save(bed.rc,  file = paste0("Analyses/", run_name, "/.focal_tables/bed_with_rc.RData"))
