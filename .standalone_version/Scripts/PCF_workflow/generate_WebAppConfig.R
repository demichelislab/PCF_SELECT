library(data.table)

source("PATHS.R")

path_local <- gsub("/shares/CIBIO-Storage/CO", "", FOLDER_PCF_ANALYSES)

sif <- fread(path_sif)

out <- data.frame(sample = gsub(".snps", "", basename(sif$Plasma)),
                  betaTable = paste0(path_local, path_beta_table),
                  segTable = paste0(path_local, "Analyses/", run_name, "/Segmentation_focal.seg"),
                  tcTable = paste0(path_local, "Analyses/", run_name, "/tc_estimations_CLONETv2_automatic.tsv"),
                  bedFile = paste0(path_local, "BedFiles/", basename(path_bed)),
                  annotationTable = paste0(path_local, path_annotations),
                  outTable = paste0(path_local, "Analyses/", run_name, "/tc_estimations_CLONETv2.tsv"))

### Update TC tables
tcTable = paste0(path_local, "Analyses/", run_name, "/tc_estimations_CLONETv2_automatic.tsv")
outTable = paste0(path_local, "Analyses/", run_name, "/tc_estimations_CLONETv2.tsv")

tmp <- fread(outTable)
tc <- tmp
tmp$notes <- ""
tmp$inspected <- "N"
fwrite(tc, tcTable, sep = ",", quote = F, col.names = T, row.names = F)
fwrite(tmp, outTable, sep = ",", quote = F, col.names = T, row.names = F)

fwrite(out, file = paste0(FOLDER_PCF_ANALYSES, "Analyses/", run_name, "/PCF_SELECT_WebApp.config"), sep = "\t", quote = F, col.names = T, row.names = F)