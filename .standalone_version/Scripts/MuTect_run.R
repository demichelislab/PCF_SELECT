library(data.table)
library(parallel)

FOLDER_PCF_ANALYSES <- "/shares/CIBIO-Storage/CO/SPICE/downloads/PCF_SELECT/analyses/2020.02/"
FOLDER_PCF_DATA <- "/shares/CIBIO-Storage/CO/SPICE/downloads/PCF_SELECT/"
setwd(FOLDER_PCF_ANALYSES)

tumor_only_mode <- F
path_analysis <- "Analyses/CellLines_Paola/"
cores <- 6

sif <- fread("SampleInfoFiles/CellLines_Paola.tsv")
sif$Plasma <- gsub("PaCBAM_ctDNA/", "", sif$Plasma)
sif$Plasma <- gsub("snps", "bam", sif$Plasma)
sif$Normal <- gsub("PaCBAM/", "", sif$Normal)
sif$Normal <- gsub("snps", "bam", sif$Normal)
dir.create(paste0(path_analysis, "Mutect"))

reference_fasta = "/shares/CIBIO-Storage/CO/elaborazioni/sharedCO/PCF_Project/Alignment_Pipeline/human_g1k_v37.fasta"
if (tumor_only_mode) {
  tumor_bams = sif$Plasma 
}else{
  tumor_bams = sif$Plasma
  normal_bams = sif$Normal
}

#Export tools
system("source /shares/CIBIO-Storage/CO/scratch/sharedCO/tools/export_tools")

mclapply(1:nrow(sif), function(i) {
  tumor_bam <- tumor_bams[i]
  if (!tumor_only_mode) {
    normal_bam <- normal_bams[i] 
  }
  
  output_mutect = paste0(FOLDER_PCF_ANALYSES, path_analysis, "Mutect/", gsub("bam", "vcf.gz", basename(tumor_bam)))
  
  if(tumor_only_mode){
    cmd = paste0("gatk4 Mutect2",
                 ' -R ', reference_fasta,
                 ' -I ', tumor_bam,
                 ' -O ', output_mutect,
                 ' 2>&1 | tee ', gsub("vcf.gz", "log", output_mutect)
    )
    cmd_filter = paste0("gatk4 FilterMutectCalls",
                        ' -R', reference_fasta,
                        ' -O', gsub(".vcf.gz", ".filtered.vcf", output_mutect),
                        ' -V', output_mutect)
  }else{
    cmd = paste0("gatk4 Mutect2",
                 ' -R ', reference_fasta,
                 ' -I ', tumor_bam,
                 ' -I ', normal_bam,
                 ' --tumor-sample ', gsub(".bam", "", basename(tumor_bam)),
                 ' --normal-sample ', gsub(".bam", "", basename(normal_bam)),
                 ' -O ', output_mutect,
                 ' 2>&1 | tee ', gsub("vcf.gz", "log", output_mutect)
    )
    cmd_filter = paste0("gatk4 FilterMutectCalls",
                        ' -R', reference_fasta,
                        ' -O', gsub(".vcf.gz", ".filtered.vcf", output_mutect),
                        ' -V', output_mutect)
  }
  
  #system("cd $(readlink -f $PWD)")
  cat(cmd,'\n')
  system(cmd)
  system(cmd_filter)
}, mc.cores = cores)
