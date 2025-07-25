#ABEMUS functions

get_case_mean_coverage <- function (tabindex, pacbamfolder_bychrom) 
{
  pacbam_list <- list.files(pacbamfolder_bychrom, full.names = TRUE)
  chrom_mean_cov <- function(cp) {
    chrom_cov <- fread(input = cp, sep = "\t", skip = 1, 
                       stringsAsFactors = FALSE, select = 9, data.table = FALSE)
    return(mean(as.numeric(chrom_cov[, 1])))
  }
  get_pacbam_case_pileups <- function(pat, pacbam_list) {
    pacbam_case_folder <- pacbam_list[grep(pattern = paste("^", pat, "$", sep = ""), x = basename(pacbam_list))]
    pacbam_case_pileups <- list.files(file.path(pacbam_case_folder, 
                                                "pileup"), full.names = TRUE)
    mean_covs <- as.numeric(sapply(pacbam_case_pileups, FUN = chrom_mean_cov))
    return(mean(mean_covs, na.rm = T))
  }
  tabindex$case_mean_coverage <- as.numeric(sapply(tabindex[, 2], get_pacbam_case_pileups, pacbam_list))
  return(tabindex)
}

apply_scaling_factor <- function (tabindex, R = 1, use.optimal.R = FALSE, target_size = NA)
{
  if (use.optimal.R) {
    if (!"case_mean_coverage" %in% colnames(tabindex)) {
      message(paste("[", Sys.time(), "]\tError. column 'case_mean_coverage' must be present in tabindex when use.optimal.R = TRUE"))
      stop()
    }
    if (is.na(target_size)) {
      message(paste("[", Sys.time(), "]\tError. 'target_size' needed when use.optimal.R = TRUE"))
      stop()
    }
    tabindex[, paste0("tabcalls_f3", "_optimalR")] <- NA
    tabindex[, paste0("tabcalls_f3", "_optimalR_used")] <- NA
    for (i in 1:nrow(tabindex)) {
      a <- fread(input = tabindex$tabcalls_f3[i], stringsAsFactors = FALSE, data.table = FALSE)
      sizes <- sort(unique(tab_optimal_R$target_Mbp))
      closest_target_size <- sizes[which(abs(sizes - target_size) == 
                                           min(abs(sizes - target_size)))]
      bombanel <- tab_optimal_R[which(tab_optimal_R$target_Mbp == closest_target_size), ]
      covs <- bombanel$mean_coverage
      closest_coverage <- covs[which(abs(covs - tabindex$case_mean_coverage[i]) == 
                                       min(abs(covs - tabindex$case_mean_coverage[i])))]
      optR <- bombanel$scalingfactor[which(bombanel$mean_coverage == closest_coverage)]
      a$filter.pbem_coverage <- a$filter.pbem_coverage * optR
      a$pass.filter.pbem_coverage <- 0
      a$pass.filter.pbem_coverage[which(a$af_case >= a$filter.pbem_coverage)] <- 1
      a <- a[which(a$pass.filter.pbem_coverage == 1), ]
      out.name <- gsub(basename(tabindex$tabcalls_f3[i]), 
                       pattern = "pmtab_F3_", replacement = paste0("pmtab_F3_optimalR_"))
      out.path <- gsub(tabindex$tabcalls_f3[i], pattern = basename(tabindex$tabcalls_f3[i]), 
                       replacement = out.name)
      write.table(x = a, file = out.path, quote = FALSE, 
                  sep = "\t", row.names = FALSE, col.names = TRUE)
      tabindex[i, paste0("tabcalls_f3", "_optimalR")] <- out.path
      tabindex[i, paste0("tabcalls_f3", "_optimalR_used")] <- optR
    }
    message(paste("[", Sys.time(), "]\talright."))
    return(tabindex)
  }
  else {
    tabindex[, paste0("tabcalls_f3", "_R", R)] <- NA
    for (i in 1:nrow(tabindex)) {
      a <- fread(input = tabindex$tabcalls_f3[i], stringsAsFactors = FALSE, data.table = FALSE)
      a$filter.pbem_coverage <- a$filter.pbem_coverage * R
      a$pass.filter.pbem_coverage <- 0
      a$pass.filter.pbem_coverage[which(a$af_case >= a$filter.pbem_coverage)] <- 1
      a <- a[which(a$pass.filter.pbem_coverage == 1), ]
      out.name <- gsub(basename(tabindex$tabcalls_f3[i]), 
                       pattern = "pmtab_F3_", replacement = paste0("pmtab_F3_R", R, "_"))
      out.path <- gsub(tabindex$tabcalls_f3[i], pattern = basename(tabindex$tabcalls_f3[i]), 
                       replacement = out.name)
      write.table(x = a, file = out.path, quote = FALSE, 
                  sep = "\t", row.names = FALSE, col.names = TRUE)
      tabindex[i, paste0("tabcalls_f3", "_R", R)] <- out.path
    }
    message(paste("[", Sys.time(), "]\talright."))
    return(tabindex)
  }
}

get_target_size <- function (targetbed, Mbp = TRUE) 
{
  bed <- fread(input = targetbed, colClasses = list(character = 1), 
               data.table = FALSE, stringsAsFactors = FALSE, header = FALSE)
  if (Mbp) {
    return(sum(bed$V3 - bed$V2)/1e+06)
  }
  else {
    return(sum(bed$V3 - bed$V2))
  }
}

