library(data.table)
library(parallel)

source("PATHS.R")

load(path_germ_distr)
load(path_beta_table)

beta.table <- ConvertBetaTable(beta.table)
seg <- fread(path_segfile, select = 1:6)
seg$chr <- paste0("chr", seg$chr)
bed <- fread(path_bed, select = 1:4)
colnames(bed) <- c("chr", "start", "end", "gene")
bt <- ConvertBetaTableToClonetv2Format(beta.table, seg, bed)

sif = fread(path_sif,data.table=F)
sif$Normal_pile <- gsub(".snps", ".pileup", sif$Normal)
sif$Plasma_pile <- gsub(".snps", ".pileup", sif$Plasma)
pileup.dir = path_dirRC
sif$Normal_rc <- paste0(pileup.dir, basename(gsub("snps", "rc", sif$Normal)))
sif$Plasma_rc <- paste0(pileup.dir, basename(gsub("snps", "rc", sif$Plasma)))
coords <- fread(path_panelgenescoords)

tc <- fread(paste0("Analyses/", run_name, "/tc_estimations_CLONETv2.tsv"), header = T)
pl_table <- tc[, c("sample", "ploidy")]
tc$adm <- 1-tc$tc
tc$adm.min <- 1-tc$tc
tc$adm.max <- 1-tc$tc
admixture_table <- tc[, c("sample", "adm", "adm.min", "adm.max")]

samples <- unique(beta.table$sample)

scna_table <- rbindlist(lapply(1:length(samples), function(i){
  ploidy_table <- as.data.frame(pl_table[sample == samples[i]])
  adm_table <- as.data.frame(admixture_table[sample == samples[i]])
  # 
  # #### Compute Admixture when CLONETv2 return NA
  # if (is.na(adm_table$adm)) {
  #   bt_sample <- bt[which(bt$sample == samples[i]), ]
  #   log2_shift <- -log2(ploidy_table$ploidy/2)
  #   bt_sample$log2.plcorr <- bt_sample$log2 - log2_shift
  #   beta <- mean(bt_sample[evidence != 0 & log2.plcorr < -0.05]$beta)
  #   if (is.nan(beta)) {
  #     beta <- NA
  #   }
  #   adm_table$adm <- beta/(2-beta)
  # }
  # ####
  
  scna_table <- compute_allele_specific_scna_table(as.data.frame(bt[which(bt$sample == samples[i]), ]), 
                                                   ploidy_table, adm_table)
  scna_table$tc <- round(1-adm_table$adm, digits = 2)
  scna_table$ploidy <- ploidy_table$ploidy
  if (!is.na(adm_table$adm)) {
    log2_shift <- round(-log2(ploidy_table$ploidy/2), digits = 3)
    log2_pl_corr <- scna_table$log2 - log2_shift
    scna_table$log2.corr <- log2(pmax((2^(log2_pl_corr) -
                                         (1-scna_table$tc)), 0)/scna_table$tc)
  }
  scna_table
}))

beta.table <- scna_table
beta.table$tc <- NULL
beta.table$ploidy <- NULL
beta.table <- merge(beta.table, tc[,c("sample", "tc", "ploidy")])

# info for chrX
seg <- fread(paste0("Analyses/", run_name, "/Segmentation_focal.seg"), select = c(1,7,2,6))
seg.X <- seg[chr == "X"]
seg.X$chr <- NULL

####### ADD SNVs calls to output ######
abemus_folder <- paste0("ABEMUS/", run_name, "/")
onco <- fread(paste0(abemus_folder, "Results/table_mutations_nocommonSNPs.oncotator.txt"))
res <- fread(paste0(abemus_folder, "Results/table_mutations_nocommonSNPs.tsv"), select = c(1:3,7:11))
colnames(res)[2:3] <- c("Chromosome", "Start_position")
onco$sample <- res$sample
dt <- onco[,c(287,1,2,5,6,7,9,14,40:42)]
dt <- merge(dt, res, by = c("sample", "Chromosome", "Start_position"))
####### ADD SNVs calls to output ######

extract_snps <- function(gene, sample, sif, coords, germ.distr, sample_type){
  ### extract het snps limited to gene
  j <- which(coords$gene_symbol == gene)
  p = fread(sif$Normal[grep(paste0("/", sample, ".snps"), sif$Plasma)])
  p = p[which(p$chr==coords$chr[j]&p$pos>=coords$start[j]&p$pos<=coords$end[j])]
  p = p[which(p$af<0.8 & p$af>0.2),]
  het.snps = p$rsid
  
  if (sample_type == "Normal") {
    ### extract AF of hetsnps in germ sample
    p = p[which(p$rsid%in%germ.distr[[1]][,1]),] 
    return(p)
  }else{
    ### extract AF of hetsnps in tumor sample
    p = fread(sif$Plasma[grep(paste0("/", sample, ".snps"), sif$Plasma)])
    p = p[which(p$rsid%in%het.snps),]
    p = p[which(p$rsid%in%germ.distr[[1]][,1]),]
    return(p)
  }
}

#### Target ####
t <- list.files(path = paste0("Analyses/", run_name, "/.focal_tables/"), pattern = "target", full.names = T)
c <- list.files(path = paste0("Analyses/", run_name, "/.focal_tables/"), pattern = "control", full.names = T)
dt <- as.data.table(rbind(rbindlist(lapply(t, function(x) {load(x); focal_target_genes})),
                          rbindlist(lapply(c, function(x) {load(x); focal_control_genes}))))

message("Compute table")
ai_log_list <- mclapply(unique(beta.table$gene), function(g){
  message(g)
  out <- rbindlist(lapply(unique(beta.table$sample), function(s) {
    if(g == "ERG_TMPRSS2"){
      seg <- fread(paste0("Analyses/", run_name, "/Segmentation_focal.seg"))
      x <- data.frame("median.log2" = seg[sample == s & gene == g]$log2ratio,
                      "median.log2_right" = NA,
                      "median.log2_left" = NA
      )
    }else{
      x <- dt[gene == g & sample == s]
    }
    nrows <- nrow(beta.table[gene == g & sample == s])
    data.frame("sample" = s,
               "gene" = g,
               "evidence_n" = beta.table[gene == g & sample == s]$evidence.n,
               "beta_n" = round(beta.table[gene == g & sample == s]$n_beta, digits = 3),
               "evidence" = beta.table[gene == g & sample == s]$evidence,
               "beta" = round(beta.table[gene == g & sample == s]$beta, digits = 3),
               "focal_log2" = rep(round(dplyr::first(x$median.log2), digits = 3), nrows),
               #"all_log2" = round(dplyr::last(x$median.log2), digits = 3),
               "all_log2" = round(beta.table[gene == g & sample == s]$log2, digits = 3),
               "all_log2_right" = rep(round(dplyr::last(x$median.log2_right), digits = 3), nrows),
               "all_log2_left" = rep(round(dplyr::last(x$median.log2_left), digits = 3), nrows),
               "informative_snps" = rep(nrow(extract_snps(g, s, sif, coords, germ.distr, "Normal")), nrows),
               "snps" = beta.table[gene == g & sample == s]$nsnp,
               "tc" = beta.table[gene == g & sample == s]$tc,
               "ploidy" = beta.table[gene == g & sample == s]$ploidy,
               "cnA" = beta.table[gene == g & sample == s]$cnA,
               "cnB" = beta.table[gene == g & sample == s]$cnB,
               "cnA.int" = beta.table[gene == g & sample == s]$cnA.int,
               "cnB.int" = beta.table[gene == g & sample == s]$cnB.int
    )
  }
  ))
  # fwrite(out, file = paste0("Analyses/allele_specific_call_Cornell_ploidy2/target_genes/",
  #                           g, ".csv"), sep = ",", quote = F, row.names = F, col.names = T)
  out
}, mc.cores = cores.aitablegen
)

ai_log_table <- rbindlist(ai_log_list)
fwrite(ai_log_table, file = paste0("Analyses/", run_name, "/ai_log2_table.csv"), 
       sep = ",", quote = F, row.names = F, col.names = T)

# for genes on chrX
ai_log_list.X <- mclapply(unique(seg.X$gene), function(g){
  message(g)
  out <- rbindlist(lapply(unique(seg.X$sample), function(s) {
    x <- seg.X[gene == g & sample == s]
    data.frame("sample" = s,
               "gene" = g,
               "evidence_n" = NA,
               "beta_n" = NA,
               "evidence" = NA,
               "beta" = NA,
               "focal_log2" = NA,
               "all_log2" = round(x$log2ratio, digits = 3),
               "all_log2_right" = NA,
               "all_log2_left" = NA,
               "informative_snps" = NA,
               "snps" = NA,
               "tc" = unique(beta.table[sample == s]$tc),
               "ploidy" = unique(beta.table[sample == s]$ploidy),
               "cnA" = NA,
               "cnB" = NA,
               "cnA.int" = NA,
               "cnB.int" = NA
    )
  }
  ))
  # fwrite(out, file = paste0("Analyses/FocalLog2_Cornell_blacklist_test/allelic_imbalance_log2_tables/target_genes/",
  #                           g, ".csv"), sep = ",", quote = F, row.names = F, col.names = T)
  out
}, mc.cores = cores.aitablegen 
)

ai_log_table.X <- rbindlist(ai_log_list.X)
fwrite(ai_log_table.X, file = paste0("Analyses/", run_name, "/ai_log2_table.chrX.csv"), 
       sep = ",", quote = F, row.names = F, col.names = T)
