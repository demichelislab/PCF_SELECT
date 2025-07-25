library(data.table)
library(parallel)

source("PATHS.R")

load(path_beta_table)
beta.table <- ConvertBetaTable(beta.table)

pdf(paste0("Analyses/", run_name, "/Plots/check_control_samples.pdf"))
for (s in unique(beta.table$sample)) {
  plot(beta.table[sample == s]$evidence, beta.table[sample == s]$evidence.n, xlim = c(0,1), ylim = c(0,1), pch = 16,
       xlab = "AI ev. in tumor sample", ylab = "AI ev. in matched control", main = s)
}
dev.off()
