FOLDER_PCF_ANALYSES <- "/shares/CIBIO-Storage/CO/SPICE/downloads/PCF_SELECT/analyses/2020.02/"
FOLDER_PCF_DATA <- "/shares/CIBIO-Storage/CO/SPICE/downloads/PCF_SELECT/"
setwd(FOLDER_PCF_ANALYSES)
source("Scripts/utility_functions.R")
source("Scripts/parameters.R")

site <- "UCLv2"
run_name <- "TR064"

Normal <- list.files(paste0(FOLDER_PCF_DATA, "PaCBAM/", site), pattern = run_name, full.names = T)
Normal <- Normal[grepl("snps", Normal)]

Plasma <- list.files(paste0(FOLDER_PCF_DATA, "PaCBAM_ctDNA/", site), pattern = run_name, full.names = T)
Plasma <- Plasma[grepl("snps", Plasma)]

sif <- data.frame(Plasma = Plasma, Normal = Normal)

write.table(sif, file = paste0("SampleInfoFiles/", run_name, ".tsv"), sep = "\t", quote = F, row.names = F, col.names = T)