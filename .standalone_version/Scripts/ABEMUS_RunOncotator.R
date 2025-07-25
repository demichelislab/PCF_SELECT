# Run oncotator on ABEMUS mutation list
# NOTE: run after ABEMUS_ApplyScalingFactor.R
library(data.table)

FOLDER_ANALYSES <- "PATH TO PCF_SELECT MAIN FOLDER"

abemus_folder <- "PATH TO ABEMUS FOLDER/run_name/"
exclude_common_snps <- T

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

if(exclude_common_snps){
  m <- fread(paste0(abemus_folder, "Results/table_mutations_nocommonSNPs.tsv"))
  maflite <- ConvertAbemustoMAFLite(m)
  write.table(maflite, file = paste0(abemus_folder, "Results/table_mutations_nocommonSNPs_MAFLITE.tsv"), 
              sep = "\t", quote = F, row.names = F, col.names = T)
  input_file <- paste0(abemus_folder, "Results/table_mutations_nocommonSNPs_MAFLITE.tsv")
}else{
  m <- fread(paste0(abemus_folder, "Results/table_mutations.tsv"))
  maflite <- ConvertAbemustoMAFLite(m)
  write.table(maflite, file = paste0(abemus_folder, "Results/table_mutations_MAFLITE.tsv"),
              sep = "\t", quote = F, row.names = F, col.names = T)
  input_file <- paste0(abemus_folder, "Results/table_mutations_MAFLITE.tsv")
}

output_file = gsub(input_file,pattern = '_MAFLITE.tsv$',replacement = '\\.oncotator.txt')
# res_dir='/shares/CIBIO-Storage/CO/SPICE/spice_marco/pipeline/pipeline_extra/resources/'
res_dir = "PATH TO ONCOTATOR RESOURCES"
oncotator_db=file.path(res_dir,"oncotator/oncotator_v1_ds_Jan262015/")
canonical_file=file.path(res_dir,"oncotator/oncotator_v1_ds_Jan262015/tx_exact_uniprot_matches.txt")
genome_build='hg19'
input_format='MAFLITE'
output_format='TCGAMAF'
# oncotator_singularity = "/shares/CIBIO-Storage/CO/SPICE/spice_marco/pipeline/pipeline_extra/tools/oncotator/versions/1.9.6.1/oncotator"
oncotator_singularity = PATH TO ONCOTATOR SINGULARITY

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
