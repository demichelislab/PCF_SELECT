#!/usr/bin/env Rscript
######################################
##Set current directory
#setwd(".")

######################################
##Load parameters: 
cmd_args <- commandArgs(trailingOnly = T)

if (length(cmd_args) != 1){
  cat(paste("[",Sys.time() ,"] Usage: SPIA.R <config_file>\n"))
  stop("Usage: SPIA.R <config_file>\n")
}

######################################
##Load configuration file
configFile <- cmd_args[1]
cat(paste("[",Sys.time() ,"] Configuration file: ",configFile,"\n",sep=""))
source(configFile)


######################################
##Open log file if needed
if( ! print_on_screen ){
  today <- format(Sys.time(), "%a_%b_%d")
  logFile = paste(today,".SPIA.log",sep="")
  sink(logFile)
}


######################################
##Plot intro message
SPIA.version<-1.1
cat(paste("[",Sys.time() ,"] Running SPIA analysis v",SPIA.version,"\n",sep=""))
cat(paste("[",Sys.time() ,"] Working directory: ",getwd(),"\n",sep=""))
cat(paste("[",Sys.time() ,"] Starting analysis\n",sep=""))


source(SPIAfunctions_location)


SPIAtable <- createSPIAtableFromVCFlist(vcfFileList,vcfSuffix)

if(saveGenotype){
  write.table(SPIAtable,genotypeTable_file,sep="\t",row.names=T,col.names=T,quote=F)
}

SPIAtest.param = list(Pmm = Pmm, nsigma = nsigma, Pmm_nonM = Pmm_nonM, nsigma_nonM = nsigma_nonM, PercValidCall=PercValidCall)
SPIAtest.result = SPIATest(SPIAtable,test.param=SPIAtest.param,verbose=verbose)

cat("Saving analsis results\n")
outTable <- SPIAtest.result$SPIAresult
write.table(outTable,outSPIAtable_file,sep="\t",row.names=F,col.names=T,quote=F)

if (saveSPIAplot){
  pdf(file=SPIAplot_file,width=11,height=8.50,useDingbats = FALSE )
  SPIAPlot(SPIAanalysis=SPIAtest.result)  
  #dev.off()
}

#end
cat(paste("[",Sys.time() ,"] SPIA analysis complete.\n",sep=""))
if ( ! print_on_screen ){
  sink()
}
