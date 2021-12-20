## Convert gene-based approach sif to ABEMUS sif format
library(data.table)
library(stringr)

cat(paste("\n[",Sys.time(),"]\tConverting SIF file to ABEMUS SIF format","\n"))

args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=3){
  message("\nERROR:\t3 arguments required")
  message("\nUSAGE:\tRscript ConvertSIFtoABEMUSformat.R sif_path bam_path out_path\n")
  cat(paste("\t[1] sif_path\tPath to SIF file to convert",
            "\t[2] bam_path\tPath to folder containing BAM files",
            "\t[3] out_path\tOutput folder",sep="\n"),"\n\n")  
  quit()
}

sif_path <- args[1]
bam_path <- args[2]
out_path <- args[3]
out_path <- paste0(file.path(out_path), "samples_info_file.tsv")

s <- fread(sif_path)
tmp <- str_split(s$Plasma, "/", simplify = T)[,ncol(str_split(s$Plasma, "/", simplify = T))]
n_plas_bam <- gsub(".snps", ".bam", tmp)
n_plas <- gsub(".bam", "", n_plas_bam)
tmp <- str_split(s$Normal, "/", simplify = T)[,ncol(str_split(s$Normal, "/", simplify = T))]
n_germ_bam <- gsub(".snps", ".bam", tmp)
n_germ <- gsub(".bam", "", n_germ_bam)
pat <- gsub("_germline", "", n_germ)
n_germ_bam <- paste0(file.path(bam_path), n_germ_bam)
n_plas_bam <- paste0(file.path(bam_path), n_plas_bam)
out <- data.frame("patient" = pat, 
                  "plasma" = n_plas, "plasma.bam" = n_plas_bam,
                  "germline" = n_germ, "germline.bam" = n_germ_bam)
fwrite(out, file = out_path, sep = "\t", quote = F, row.names = F, col.names = T)