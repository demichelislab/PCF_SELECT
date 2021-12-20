library(data.table)
library(parallel)

source("PATHS.R")

#FOLDER_PCF_ANALYSES <- "/shares/CIBIO-Storage/CO/SPICE/downloads/PCF_SELECT/analyses/2020.02/"
#FOLDER_PCF_DATA <- "/shares/CIBIO-Storage/CO/SPICE/downloads/PCF_SELECT/"
#setwd(FOLDER_PCF_ANALYSES)
#source("Scripts/utility_functions.R")

#run_name <- "Vancouver"

#seg_path <- paste0("Analyses/", run_name, "/Segmentation.seg")
#bed_path <- "BedFiles/Nimblegen_Regions.bed"
#beta.table_path <- paste0("Analyses/", run_name, "/betaTable.RData")
output_path <- paste0("Analyses/", run_name, "/tc_estimations_CLONETv2.tsv")

seg_tb <- fread(path_segfile, header = T, stringsAsFactors = F)
seg_tb$chr <- paste0("chr", seg_tb$chr)
samples <- as.character(unique(seg_tb$sample))
bed <- fread(path_bed, select = c(1:4))
colnames(bed) <- c("chr", "start", "end", "gene")
load(path_beta_table)
gl <- fread(path_annotations)

#### Automatic TC estimation using CLONETv2

beta.table <- ConvertBetaTableToClonetv2Format(beta.table = beta.table, seg = seg_tb[,1:6], bed = bed)
beta.table <- beta.table[gene %in% gl[PCF_SELECT == "yes"]$hgnc_symbol]

adm <- matrix(nrow = 0, ncol = 4)
out <- mclapply(1:length(samples), function(i) {
  ploidy_table <- compute_ploidy(beta_table = beta.table[which(beta.table$sample == samples[i]), ])
  adm_table <- compute_dna_admixture(beta_table = beta.table[which(beta.table$sample == samples[i]), ], ploidy_table = ploidy_table)
  cbind(adm_table, "ploidy" = ploidy_table$ploidy)
}, mc.cores = cores.TCestimation)

CLONETv2_table <- rbindlist(out, fill = T)
CLONETv2_table$tc <- 1-CLONETv2_table$adm

#### Manual TC estimation
# polyploid <- c("4_cfDNA2", "PM63-SP158_ctDNA")
# CLONETv2_table[ploidy > 2.15 & !(sample %in% polyploid)]$ploidy <- 2
tc <- CLONETv2_table
bt <- beta.table
pl_table <- tc[, c("sample", "ploidy")]

admixture_table <- tc[, c("sample", "adm", "adm.min", "adm.max")]

samples <- unique(beta.table$sample)

tc_manual <- rbindlist(lapply(1:length(samples), function(i){
  ploidy_table <- as.data.frame(pl_table[sample == samples[i]])
  adm_table <- as.data.frame(admixture_table[sample == samples[i]])
  
  #### Compute Admixture when CLONETv2 return NA
  bt_sample <- bt[which(bt$sample == samples[i]), ]
  log2_shift <- -log2(ploidy_table$ploidy/2)
  bt_sample$log2.plcorr <- bt_sample$log2 - log2_shift
  beta <- mean(bt_sample[evidence != 0 & log2.plcorr < -0.05]$beta)
  if (is.nan(beta)) {
    beta <- NA
  }
  adm_table$tc_manual <- 1-round(beta/(2-beta), digits = 2)
  adm_table[,c("sample", "tc_manual")]
}))

CLONETv2_table$tc_CLONETv2 <- CLONETv2_table$tc
TC <- merge(CLONETv2_table, tc_manual)
TC$tc <- ifelse(is.na(TC$tc), TC$tc_manual, TC$tc)
TC$adm <- 1-TC$tc

fwrite(TC, file = output_path, row.names = F, col.names = T, quote = F)
