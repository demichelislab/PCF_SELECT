library(data.table)
library(parallel)

source("PATHS.R")

m <- fread(paste0(path_abemus, "Results/table_mutations_nocommonSNPs.tsv"))
onco <- fread(paste0(path_abemus, "Results/table_mutations_nocommonSNPs.oncotator.txt"))
onco$sample <- m$sample
dt <- onco[,c(287,1,2,5,6,7,9,14,40:42)]

res <- fread(paste0(path_abemus, "Results/table_mutations_nocommonSNPs.tsv"), select = c(1:3,7:11))
colnames(res)[2:3] <- c("Chromosome", "Start_position")
dt <- merge(dt, res, by = c("sample", "Chromosome", "Start_position"))
dt$sign <- paste(dt$Chromosome, dt$Start_position, sep = ":")
dt$patient <- gsub("_ctDNA-.*$", "", dt$sample)
dt$rc_alt_tumor <- round(dt$af_case*dt$cov_case, digits = 0)
dt$rc_ref_tumor <- round((1-dt$af_case)*dt$cov_case, digits = 0)

SNVTable <- dt
SNVTable$gene <- SNVTable$Hugo_Symbol

panel <- fread(path_panelgenescoords)
aiTable <- as.data.table(rbind(fread(paste0("Analyses/", run_name, "/ai_log2_table.csv")), fread(paste0("Analyses/", run_name, "/ai_log2_table.chrX.csv"))))
aiTable$adm <- 1- aiTable$tc
aiTable$Hugo_Symbol <- aiTable$gene
aiTable$chr <- sapply(1:nrow(aiTable), function (i) panel[match(aiTable$gene[i], panel$gene_symbol)]$chr, simplify = T)
aiTable$start <- sapply(1:nrow(aiTable), function (i) panel[match(aiTable$gene[i], panel$gene_symbol)]$start, simplify = T)
aiTable$end <- sapply(1:nrow(aiTable), function (i) panel[match(aiTable$gene[i], panel$gene_symbol)]$end, simplify = T)

d <- merge(SNVTable, aiTable, by = c("sample", "gene"), all.x = T)

x <- rbindlist(lapply(1:nrow(d), function(i) adjustVAF(i, SNVtable.ext = d)))
SNVTable <- x[Variant_Classification %in% c("Nonsense_Mutation", "Missense_Mutation")]

cols_to_exclude <- c("Hugo_Symbol.x", "Hugo_Symbol.y", "t_af","chr", "start", "end", "adm", "sign")
out_table <- as.data.frame(SNVTable)[,which(!colnames(SNVTable) %in% cols_to_exclude)]

### Test full impairment ###
out_table$is.fullimpaired <- ifelse((out_table$cn.int == out_table$CN_SNVmut) | (out_table$cnA.int == 0 & out_table$cnB.int == 0), 1, 0)
###

fwrite(out_table, file = paste0("Analyses/", run_name, "/SNVs_calls_corrected.csv"), sep = ";", quote = F, row.names = F, col.names = T)