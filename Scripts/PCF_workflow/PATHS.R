run_name <- "run_name"

### Absolute path to folder where the repository is cloned
FOLDER_PCF_ANALYSES <- "/home/francesco.orlando/PCF_SELECT_github/"
### Absolute path to folder containing data folders (e.g. BAMs, PaCBAM, PaCBAM_RC)
FOLDER_PCF_DATA <- "/home/francesco.orlando/PCF_SELECT_github/Data/"

setwd(FOLDER_PCF_ANALYSES)
source("Scripts/utility_functions.R")
source("Scripts/parameters.R")
path_betafunctions <- "Scripts/betaFunctions.R"
path_dirRC <- paste0(FOLDER_PCF_DATA, "/PaCBAM_RC/")
path_bed <- "BedFiles/Nimblegen_Regions_v2.bed" ### depends on panel version
path_bed_focal <- "BedFiles/Nimblegen_Regions_v2_FocalLog2.bed" ### depends on panel version
path_controlgenes <- "BedFiles/Panel_v2_20190503.bed" ### depends on panel version
path_sif <- paste0("SampleInfoFiles/", run_name, ".tsv")
path_segfile <- paste0("Analyses/", run_name, "/Segmentation.seg")
path_amplicons_log2r <- paste0("Analyses/", run_name, "/.amplicons_log2R/")
path_germ_distr <- paste0("GermlineModels/", run_name, ".RData") ### pre-computed reference models can be used
path_beta_table <- paste0("Analyses/", run_name, "/betaTable.RData")
path_panelgenescoords <- "BedFiles/Panelv2_Genes_Coordinates_hg19.bed" ### depends on panel version
path_shift_table <- paste0("Analyses/", run_name, "/peaks_shift.tsv")
path_annotations <- "AnnotationFiles/GeneList_PCF-SELECT_Annotation.tsv"
path_abemus <- paste0("ABEMUS/", run_name, "/")

compute_GermlineModel <- F
store_amplicons_log2R <- F

cores.peakcorrection = 40
cores.readdepth = 40
cores.betacomputation = 40
cores.TCestimation = 40
cores.aitablegen = 40
cores.multipanel = 40
cores.SPIA = 40
