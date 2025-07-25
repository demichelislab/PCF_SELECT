library(data.table)
library(parallel)

source("PATHS.R")

tostring <- function(x){paste0("'", x, "'")}

cores = 50
out.spia <- paste0("Analyses/", run_name, "/SPIA/") 
out.spia.vcf <- paste0("Analyses/", run_name, "/SPIA/vcfs/")
out.spia.sif <- paste0("Analyses/", run_name, "/SPIA/sifs/")
out.spia.res <- paste0("Analyses/", run_name, "/SPIA/results/")
out.spia.config <- paste0("Analyses/", run_name, "/SPIA/config/")

out.spia.sif_comb <- paste0("Analyses/", run_name, "/SPIA/sifs_combinations/")
out.spia.res_comb <- paste0("Analyses/", run_name, "/SPIA/results_combinations/")
out.spia.config_comb <- paste0("Analyses/", run_name, "/SPIA/config_combinations/")

dir.create(out.spia)
dir.create(out.spia.vcf)
dir.create(out.spia.sif)
dir.create(out.spia.res)
dir.create(out.spia.config)
dir.create(out.spia.sif_comb)
dir.create(out.spia.res_comb)
dir.create(out.spia.config_comb)


path_SPIA = "Scripts/SPIA/SPIA.R"
path_SPIA_vcf_generation = "Scripts/SPIA/SPIA_generate_genotype_vcf.sh"
path_SPIA_functions = "Scripts/SPIA/SPIAfunctions.R"
path_SPIA_config = "Scripts/SPIA/config_file.R"

source(path_SPIA_functions)

min_depth_of_coverage=50
heterozygous_perc=0.2
snps_list_file="Scripts/SPIA/snps_spia_PCFv2-hg19_nochr.vcf"

sif <- fread(path_sif)

mclapply(1:nrow(sif), function(i){
  normal <- basename(sif$Normal)
  plasma <- basename(sif$Plasma[i])
  
  #### Retrieve VCF
  ### Normal
  output_prefix = paste0(out.spia.vcf, gsub(".snps", "", normal))
  # input_file = sif$Normal[i]
  # system(paste("sh", path_SPIA_vcf_generation, min_depth_of_coverage, heterozygous_perc, output_prefix, snps_list_file, input_file, basename(output_prefix), sep = " "))
  vcf.normal = paste0(output_prefix, ".vcf")
  
  ### Plasma
  output_prefix = paste0(out.spia.vcf, gsub(".snps", "", plasma))
  # input_file = sif$Plasma[i]
  # system(paste("sh", path_SPIA_vcf_generation, min_depth_of_coverage, heterozygous_perc, output_prefix, snps_list_file, input_file, basename(output_prefix), sep = " "))
  vcf.plasma = paste0(output_prefix, ".vcf")
  ####
  
  #### create SIF for SPIA
  vcfFileList <- normalizePath(paste0(out.spia.sif_comb, gsub(".snps", ".txt", plasma)))
  fileConn <- file(vcfFileList)
  writeLines(c(normalizePath(vcf.normal), normalizePath(vcf.plasma)), fileConn)
  close(fileConn)
  ####
  
  ## SPIA table
  outSPIAtable_file = paste0(normalizePath(out.spia.res_comb), "/", gsub(".snps", "", plasma), "_spia.tsv")
  ## SPIA plot file name
  SPIAplot_file = paste0(normalizePath(out.spia.res_comb), "/", gsub(".snps", "", plasma), "_spia_plot.pdf")
  
  ### Prepare config file
  ## create copy of config file
  path_config <- paste0(normalizePath(out.spia.config_comb), "/", gsub(".snps", "", plasma), "_config_file.R")
  system(paste("cp", path_SPIA_config, path_config, sep = " "))
  write(paste0("SPIAplot_file = ", tostring(SPIAplot_file)), file = path_config, append=TRUE)
  write(paste0("SPIAfunctions_location = ", tostring(normalizePath(path_SPIA_functions))), file = path_config, append=TRUE)
  write(paste0("vcfFileList = ", tostring(vcfFileList)), file = path_config, append=TRUE)
  write(paste0("outSPIAtable_file = ", tostring(outSPIAtable_file)), file = path_config, append=TRUE)
  
  ### Run SPIA
  system(paste("Rscript", normalizePath(path_SPIA), normalizePath(path_config), "2>&1 | tee", paste0(normalizePath(out.spia.res_comb), "/", gsub(".snps", "", plasma), "_spia.log"), sep = " "))
}, mc.cores = cores)


spia.report <- rbindlist(lapply(list.files(paste0("Analyses/", run_name, "/SPIA/results/"), pattern = ".tsv", full.names = T), function(i) fread(i)))
fwrite(spia.report, paste0(out.spia, "SPIA_summary.csv"), sep = ";", quote = F, col.names = T, row.names = F)