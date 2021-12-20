suppressMessages(library(colorRamps))

thr.tc = 0.15

thr.imbG = 0.5

thr.imbT = 0.2

thr.log2R_lowtc = 0.1

ai.nrep = 100 ### Default used in PCF_SELECT
ai.nrep = 1 ### For testing

### Evidence Color Mapping
ColorSchema <- colorRampPalette(c("yellow", "red"))(100)
# colfunc(100)
# plot(rep(1,100),col=colfunc(100),pch=19,cex=3)

#ColorSchema <- grDevices::heat.colors(n = 100, rev = T)
colorMapping <- c("white", ColorSchema)
names(colorMapping) <- as.character(seq(0,1,0.01))
rm(ColorSchema)
###

### Chromosomes colors
# chr.cols <- terrain.colors(24, alpha = .8)
chr.cols <- primary.colors(24, steps = 3)
names(chr.cols) <- c(1:22, "X", "Y")
chr.cols.withchr <- chr.cols
names(chr.cols.withchr) <- paste0("chr", names(chr.cols.withchr))
