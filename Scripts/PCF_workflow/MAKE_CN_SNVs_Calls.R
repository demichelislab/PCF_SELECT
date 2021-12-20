library(data.table)
library(parallel)

source("PATHS.R")

sif = fread(path_sif, data.table=F)

gl <- fread(path_annotations)
control.genes <- gl[PCF_SELECT == "control_gene"]$hgnc_symbol
# control.genes = fread("BedFiles/Panel_20170516.bed",data.table=F)
# control.genes = unique(gsub("CONTROL_GENE\\|","",control.genes[grep("CONTROL_GENE",control.genes[,4]),4]))
# control.genes = c(control.genes, c("TEX11", "HDAC8"))

# Load ABEMUS calls
abemus_folder <- path_abemus
oncotator <- fread(paste0(abemus_folder, "Results/table_mutations_nocommonSNPs.oncotator.txt"))
tab_muts <- fread(paste0(abemus_folder, "Results/table_mutations_nocommonSNPs.tsv"), select = c(1:3,7:11))
colnames(tab_muts)[2:3] <- c("Chromosome", "Start_position")
oncotator$sample <- tab_muts$sample
dt <- oncotator[,c(287,1,2,5,6,7,9,14,40:42)]
dt <- merge(dt, tab_muts, by = c("sample", "Chromosome", "Start_position"))
dt$sign <- paste(dt$Chromosome, dt$Start_position, sep = ":")

# Consider only Non-Synonymous mutations
# apply stringent filters on minumum coverage in case sample
# NB: applied filters are ABEMUS with optimal R, retaining only non-common SNPs according to dbSNP
dt <- dt[Variant_Classification %in% c("Nonsense_Mutation", "Missense_Mutation") & cov_case >= 50]

# annotate SU2C2 SNVs
su2c2 <- fread("ABEMUS/utility/snv_list_SU2C2_PCF_SELECT.tsv")
su2c2$sign <- paste(su2c2$chr, su2c2$start, sep = ":")

dt$in.su2c2 <- ifelse(dt$sign %in% su2c2$sign, T, F)
fwrite(dt, paste0("Analyses/", run_name, "/SNVs_calls.csv"), sep = ",",
       col.names = T, row.names = F, quote = F)

# Load AI tables
ai <- fread(paste0("Analyses/", run_name, "/ai_log2_table.csv"))
ai$gene_type <- ifelse(ai$gene %in% control.genes, "control", "target")
ai$gene_type <- ifelse(ai$gene %in% gl[PCF_SELECT == "other"]$hgnc_symbol, "other", ai$gene_type)
ai.X <- fread(paste0("Analyses/", run_name, "/ai_log2_table.chrX.csv"))
ai.X$gene_type <- ifelse(ai.X$gene %in% control.genes, "control", "target")

ai <- as.data.table(rbind(ai, ai.X))

ai$log2_cutoff_pass <- sapply(1:nrow(ai), is.log2.cutoff.pass,
                              AllelicImbalance.table = ai
)

ai$chr <- gl[match(ai$gene, gl$hgnc_symbol)]$chromosome_name
ai <- ai[!is.na(chr)]
ai$all_log2.corr <- ifelse(ai$log2_cutoff_pass, correct_log2(ai$all_log2, ai$ploidy, ai$tc), ai$all_log2)
ai$all_log2.corr <- ifelse(is.na(ai$tc), NA, ai$all_log2.corr)

### Dynamic thresholds for Log2R
#ai$cn.call.corr.old <- ai$cn.call.corr
ai$cn.call.corr <- sapply(1:nrow(ai), get_lesion_type_v4,
                          AllelicImb.G = ai$evidence_n, 
                          AllelicImb.T = ai$evidence, 
                          log2R.unc = ai$all_log2,
                          log2R.corr = ai$all_log2.corr,
                          Gene = ai$gene,
                          Chromosome = ai$chr,
                          TC = ai$tc,
                          snps.no = ai$snps,
                          use.corr.log2R = T
)

### add SNVs info to table
ai$non_syn_snv <- 0
ai$snv_in_su2c <- 0
for(i in 1:nrow(dt)) {
  s <- dt[i,]$sample
  g <- dt[i,]$Hugo_Symbol
  in.su2c <- dt[i,]$in.su2c2
  ai[sample == s & gene == g]$non_syn_snv <- ai[sample == s & gene == g]$non_syn_snv + 1
  if (in.su2c) {
    ai[sample == s & gene == g]$snv_in_su2c <- ai[sample == s & gene == g]$snv_in_su2c + 1
  }
}

# fwrite(ai, "Analyses/allele_specific_call/Cornell_CN_SNVs_calls.csv", sep = ",", quote = F,
#        row.names = F, col.names = T)
fwrite(ai, paste0("Analyses/", run_name, "/CN_SNVs_calls.csv"), sep = ",", quote = F,
       row.names = F, col.names = T)

