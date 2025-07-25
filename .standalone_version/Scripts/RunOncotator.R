library(data.table)
library(R.utils)
library(parallel)

FOLDER_PCF_ANALYSES <- "/shares/CIBIO-Storage/CO/SPICE/downloads/PCF_SELECT/analyses/2020.02/"
FOLDER_PCF_DATA <- "/shares/CIBIO-Storage/CO/SPICE/downloads/PCF_SELECT/"
setwd(FOLDER_PCF_ANALYSES)
source("Scripts/utility_functions.R")
source("Scripts/parameters.R")

ANALYSIS_FOLDER <- "Analyses/mutect_cornell/"
cores <- 6

mutect_files <- list.files(paste0(ANALYSIS_FOLDER, "Mutect/"), pattern = "vcf$", full.names = T)

mclapply(1:length(mutect_files), function(i){
  x <- fread(mutect_files[i], skip = as.integer(system(paste0("less ", mutect_files[i], " | grep '##' | wc -l"), intern = T)), header = T)
  x <- x[FILTER == "PASS" & `#CHROM` %in% c(1:22,"X","Y")]
  colnames(x)[c(1,2,4,5)] <- c("chr", "pos", "ref", "alt")
  maflite <- ConvertAbemustoMAFLite(x)
  
  # Run oncotator on ABEMUS mutation list
  # NOTE: run after ABEMUS_ApplyScalingFactor.R
  input_file <- gsub("vcf", "maflite.tsv", mutect_files[i])
  write.table(maflite, file = input_file,
              sep = "\t", quote = F, row.names = F, col.names = T)
  input_file <- paste0(FOLDER_PCF_ANALYSES, input_file)
  
  
  output_file = gsub(input_file, pattern = 'maflite.tsv$',replacement ='oncotator.txt')
  res_dir='/shares/CIBIO-Storage/CO/SPICE/spice_marco/pipeline/pipeline_extra/resources/'
  oncotator_db=file.path(res_dir,"oncotator/oncotator_v1_ds_Jan262015/")
  canonical_file=file.path(res_dir,"oncotator/oncotator_v1_ds_Jan262015/tx_exact_uniprot_matches.txt")
  genome_build='hg19'
  input_format='MAFLITE'
  output_format='TCGAMAF'
  oncotator_singularity = "/shares/CIBIO-Storage/CO/SPICE/spice_marco/pipeline/pipeline_extra/tools/oncotator/versions/1.9.6.1/oncotator"
  
  cmd = paste(oncotator_singularity,'-v',
              '--input_format',input_format,
              '--output_format',output_format,
              '--db-dir',oncotator_db,
              '-c',canonical_file,
              input_file,
              output_file,
              genome_build
  )
  
  system("cd $(readlink -f $PWD)")
  cat(cmd,'\n')
  system(cmd)
}, mc.cores = cores)
