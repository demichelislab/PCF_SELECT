# Apply scaling factor and filters to ABEMUS results and output one annotated mutations table (filtering common SNPs or not)

library(data.table)
library(parallel)

FOLDER_ANALYSES <- PATH TO PCF_SELECT MAIN FOLDER
cores <- 40 

source(paste0(FOLDER_ANALYSES, "ABEMUS/utility/ABEMUS_functions.R"))

abemus_folder <- "/shares/CIBIO-Storage/CO/SPICE/data/PSMA/PCF_SELECT/ABEMUS/PSMA/"
pacbamfolder_bychrom <- "/shares/CIBIO-Storage/CO/SPICE/data/PSMA/PaCBAM_byChrom/"
bed_path <- paste0(FOLDER_ANALYSES, "BedFiles/Nimblegen_Regions_v2.bed")
vcf_path <- gsub(".bed", ".vcf", bed_path)
use_optimal_R <- T
min_cov_case <- 50
min_cov_control <- 50
max_af_control <- 0.01
thr.MAF <- 0.1

### remove when package is updated
tab_optimal_R <- read.table(paste0(FOLDER_ANALYSES, "ABEMUS/utility/table_optimal_R.tsv"), header = T)
###

setwd(abemus_folder)
if(!file.exists("samples_info_file_rpa.tsv")){
  sif <- fread("samples_info_file.tsv", data.table = F)
  sif <- get_case_mean_coverage(tabindex = sif, pacbamfolder_bychrom = pacbamfolder_bychrom)
  sif$tabcalls_f1 <- paste0(abemus_folder, "Results/", sif$patient,"/pmtab_F1_", sif$plasma, ".tsv")
  sif$tabcalls_f2 <- paste0(abemus_folder, "Results/", sif$patient,"/pmtab_F2_", sif$plasma, ".tsv")
  sif$tabcalls_f3 <- paste0(abemus_folder, "Results/", sif$patient,"/pmtab_F3_", sif$plasma, ".tsv")
  write.table(sif, paste0(abemus_folder, "samples_info_file_rpa.tsv"), row.names = F, col.names = T, quote = F, sep = "\t")
}else{
  sif <- fread(paste0(abemus_folder, "samples_info_file_rpa.tsv"))
  sif$tabcalls_f1 <- paste0(abemus_folder, "Results/", sif$patient,"/pmtab_F1_", sif$plasma, ".tsv")
  sif$tabcalls_f2 <- paste0(abemus_folder, "Results/", sif$patient,"/pmtab_F2_", sif$plasma, ".tsv")
  sif$tabcalls_f3 <- paste0(abemus_folder, "Results/", sif$patient,"/pmtab_F3_", sif$plasma, ".tsv")
  }

tsize = get_target_size(bed_path)
su2c <- fread(paste0(FOLDER_ANALYSES, "ABEMUS/utility/snv_list_SU2C2_PCF_SELECT.tsv"))
su2c$sign <- paste(su2c$chr, su2c$start, sep = ":")

tabindex <- as.data.frame(sif)
if (use_optimal_R) {
  tabindex <- apply_scaling_factor(tabindex, use.optimal.R = T, target_size = tsize)
  write.table(tabindex, file = "tabindex_optimalR.tsv", sep = "\t", col.names = T, row.names = F, quote = F)
}else{
  tabindex <- apply_scaling_factor(tabindex, R = 1, use.optimal.R = F, target_size = tsize)
  write.table(tabindex, file = "tabindex.tsv", sep = "\t", col.names = T, row.names = F, quote = F)
}

res <- mclapply(1:nrow(tabindex), function(i) {
  tab <- fread(tabindex$tabcalls_f3_optimalR[i])
  tab <- tab[cov_case >= min_cov_case & cov_control >= min_cov_control]
  tab$sign <- paste(tab$chr, tab$pos, sep = ":")
  tab$in.su2c <- ifelse(tab$sign %in% su2c$sign, 1, 0)
  tab$sample <- tabindex$plasma[i]
  tab
  }, mc.cores = cores)
res_raw <- rbindlist(res)
### Re-annotate SNPs
vcf <- fread(vcf_path, select = 1:3)
colnames(vcf) <- c("chr", "pos", "dbsnp.n")
res_raw <- merge(res_raw, vcf, all.x = T, by = c("chr", "pos"))
res_raw$dbsnp <- res_raw$dbsnp.n
res_raw$dbsnp.n <- NULL
res <- res_raw[,c("sample", "chr", "pos", "ref", "alt", "dbsnp", "af_case", "cov_case", "af_control", "cov_control","bperr", "in.su2c")]

### filter out common SNPs ###
rsids <- readLines(paste0(FOLDER_ANALYSES, "ABEMUS/utility/rsids_common_", thr.MAF, ".txt"))
res_raw <- res_raw[!(dbsnp %in% rsids) & af_control < max_af_control]
res <- res[!(dbsnp %in% rsids) & af_control < max_af_control]

fwrite(res_raw, file = paste0(abemus_folder, "Results/table_mutations.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
fwrite(res, file = paste0(abemus_folder, "Results/table_mutations_nocommonSNPs.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)
