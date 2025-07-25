#################################
#
# Function to manage VCF file and 
# create the table for SPIA analysis
#
#################################


##Utility function
paste0 <- function(...){
  return(paste(...,sep=""))
}



#################################
# Take a field from a VCF file and search
# for GT. Return a code compatible with SPIAssay, i.e., 
# SNPs are coded as 0,1,2,NA (AA,BB,AB,NoCall) 
getGT <- function(field){
  
  val <- regexpr("0/0",field)
  if ( val == 1 ){return(0) }
  
  val <- regexpr("0\\|0",field)
  if ( val == 1 ){return(0) }
    
  val <- regexpr("0/1",field)
  if ( val == 1 ){return(2) }
  
  val <- regexpr("0\\|1",field)
  if ( val == 1 ){return(2) }
  
  val <- regexpr("1/0",field)
  if ( val == 1 ){return(2) }
  
  val <- regexpr("1\\|0",field)
  if ( val == 1 ){return(2) }
    
  val <- regexpr("1/1",field)
  if ( val == 1 ){return(1) }
  
  val <- regexpr("1\\|1",field)
  if ( val == 1 ){return(1) }
  
  return(NA)
  
}


#################################
# Take a vcf where GT is encoded in the INFO field
# and return a DF whith one column with the genotype
# and with rownames equal to dbSNP id
loadVCFgenotype <- function(file_name){
  
  
  #Load VCF columns name
  con  <- file(file_name, open = "r")
  
  #say if colnames are found
  colnamesFound <- F
  

  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0  && !colnamesFound ) {
    #cat(oneLine)
    
    if ( regexpr("##",oneLine) == 1 ){ next }
    
    if ( regexpr("#",oneLine) == 1 ){ 
      colnames <- strsplit(oneLine,'\t')[[1]]
      colnamesFound <- T  
    }
    
    if ( regexpr("#",oneLine) != 1 ){
      break
    }
       
  } 
  
  close(con)
  
  if (!colnamesFound || length(colnames) < 10){
    err<-paste0("VCF ",file_name," with less than 10 columns\n")
    cat(paste0("[",Sys.time(),"] Error: " ,err))
    stop(err)
  }
  
  #load the vcf data
  vcfFile <- read.csv(file=file_name,sep="\t",header=F,comment.char="#",as.is=T)
  
  #set vcf data colnames
  colnames(vcfFile) <- colnames
  
  #dbSNP id will be used as key
  #check if it is unique
  if (length(unique(vcfFile$ID)) < length(vcfFile$ID)){
    err<-paste0("VCF file ",file_name," has not unique elements in the ID field. They will not be considered.\n")
    cat(paste0("[",Sys.time(),"] Warning: " ,err))
    warning(err)
    
    #remove duplicates
    vcfFile <- vcfFile[which(!(duplicated(vcfFile$ID) | duplicated(vcfFile$ID,fromLast=T))),]
    
  }
  
  #retrive the genotype info
  genotypeDF <- data.frame(row.names=vcfFile$ID)

  
  #parse columns
  for(i in c(10:length(colnames))){
    genotypeDF[,colnames[i]] <- sapply(X=vcfFile[,i],FUN=getGT)  
  }  
  return(genotypeDF)
}

#################################
# Take a file where the first column contain the list 
# of VCF files with the genotype and create a table
# with a row for each SNP and a column for each sample.
# File names are used as sample names. 
# Optionally a suffix can be removed from vcf file names
createSPIAtableFromVCFlist <- function(vcfFileList,vcfSuffix=""){
  
  #check existence
  if (  file.access(vcfFileList,4) != 0){
    err <- paste0("VCF list ",vcfFileList," does not exist or it is not readable\n")
    cat(paste0("[",Sys.time(),"] Error: " ,err))
    stop(err)
  }
  
  fileList <- read.csv(vcfFileList,header=F,sep="\t",as.is=T)
  
  #check emptyness
  if (length(fileList) == 0){
    err <- paste0("VCF list ",vcfFileList," is empty")
    cat(paste0("[",Sys.time(),"] Error: " ,err))
    stop(err)  
  }
  
  #use the first file to fix the number of rows of the DF
  firstVcfFile <- fileList[,1][1]
  SPIAtable <- loadVCFgenotype(firstVcfFile)
  
  #if other vcf file are available
  for( vcfFile in tail(fileList[,1],n=-1 )  ){
    
    vcf <-loadVCFgenotype(vcfFile)
    
    #check if SNPs are comaptible
    if (! all(rownames(vcf) %in% rownames(SPIAtable))){
      err <- paste0("The SNPs in the VCF file ",vcfFile," do not match those of VCF file ",firstVcfFile,". It will be ignored.\n")
      cat(paste0("[",Sys.time(),"] Warning: " ,err))
      warning(err)
      next
    }
    
    SPIAtable[rownames(vcf),colnames(vcf)] <- vcf

  }
  
  return(SPIAtable)
  
}



#################################
#
# Functions to compute the SPIA distance
#
#################################


####################################################################
#
# ComputingSPIADistance
#  Input: - two vectors of SNPs, one for each cell line
#           SNPs are coded as 0,1,2,NA (AA,BB,AB,NoCall)
#         - significative digits (default 5)
#         - verbose mode
#  Output: a row with informations of about the distance
#           1. SPIA distance  
#           2. number of valid calls
#           3. number of total calls
#           4. number of calls where one of the two SPNs are not available
#           5. number of calls where both SNPs are not available
#           6. number of calls where SNP change from {AA,BB} to AB or from AB to {AA,BB}
#           7. number of calls where SNP change from AA to BB or from BB to AA
#
####################################################################

ComputingSPIADistance<-function(vector1, vector2, namesVect1, nameVect2, defaultDigists = 5, verbose = FALSE){
  
  ## Control that the two vectors has the same length 
  if(length(vector1)!=length(vector2)){
    cat(paste0("Warning: in function ComputingSPIADistance input vectors",nameVect1," and ",nameVect2," must have the same length\n"))
    return(-1)
  }  
  
  vector1 <- factor(vector1)
  vector2 <- factor(vector2)
  levels(vector1) <- c(0,1,2)
  levels(vector2) <- c(0,1,2)
  
  ##Build contingency table of vector1 and vector2
  SummaryTable<-table(vector1,vector2,exclude=NULL);
  
  lenVet<-length(vector1);
  
  ##COMPUTATION START
  if (verbose) {cat("SPIA: computing distance:\n")}
  
  ##Both SNPs are not available
  counterBothNA <- SummaryTable[4,4]
  if (verbose) {cat(paste0("SPIA: number of calls where both SNPs are not available:",counterBothNA,"\n")) }
  
  ##One of the two SPNs are not available
  counterOneNA <- sum(SummaryTable[4,c(1:3)], SummaryTable[c(1:3),4])
  if (verbose) {cat(paste0("SPIA: number of calls where one of the two SNPs are not available:",counterOneNA,"\n")) }
  
  ##From AA to BB or from BB to AA
  counterDiffAvsBorBvsA <- SummaryTable["1","0"]+SummaryTable["0","1"]
  if (verbose) {cat(paste0("SPIA: number of calls where AA becomes BB or BB becomes AA:",counterDiffAvsBorBvsA,"\n")) }
  
  ##From {AA,BB} to AB or from AB to {AA,BB}
  counterDiffAorBvsABorviceversa <- SummaryTable["2","0"]+
    SummaryTable["2","1"]+
    SummaryTable["0","2"]+
    SummaryTable["1","2"]
  if (verbose) {cat(paste0("SPIA: number of calls where AA or BB become AB or vice versa:",counterDiffAorBvsABorviceversa,"\n")) }                                    
  
  ##From {AA,BB} to AB
  counterDiffABvsAorB <- SummaryTable["0","2"] + SummaryTable["1","2"]
  if (verbose) {cat(paste0("SPIA: number of calls where AA or BB become AB:",counterDiffABvsAorB,"\n")) }
  
  ##Both SNPs are Homozygous
  counterBothHomoz <- SummaryTable["0","0"] + SummaryTable["1","1"]
  if (verbose) {cat(paste0("SPIA: number of calls where both SNPs are homozygous:",counterBothHomoz,"\n")) }
  
  ##Both SNPs are Hetherozygous
  counterBothHeter<-SummaryTable["2","2"]
  if (verbose) {cat(paste0("SPIA: number of calls where both SNPs are heterozygous:",counterBothHeter,"\n")) }
  
  ##Computing SPIA Distance
  somma<-counterDiffAorBvsABorviceversa+counterDiffAvsBorBvsA
  counter<-counterDiffAorBvsABorviceversa+counterDiffAvsBorBvsA+counterBothHomoz+counterBothHeter;
  distance<-somma/counter;
  if (verbose) {cat(paste0("SPIA: distance is ",distance,"\n")) }
  
  ##Return    
  return(c(signif(distance,digits=defaultDigists),counter,lenVet,counterOneNA,counterBothNA,counterDiffAvsBorBvsA,counterDiffAorBvsABorviceversa,counterDiffABvsAorB,counterBothHomoz,counterBothHeter))
  
}


######################################################################
#
# getSPIALimits: computes the thresholds that define the limits for the SPIA test
# 
# Input: 1. N: number of valid calls
#        2. Pmm: probability of mismatch in a matching population (dafault 0.1)
#        3. nsigma: parameter that characterize the limit for mm (dafault 2)
#        4. Pmm_nonM: probability of mismach in a non matching population (dafault 0.6)
#        5. nsigma_nonM: parameter that characterize the limit for mm_nonM (dafault 3)
#        6. verbose: print verbose information on error
#
######################################################################

getSPIALimits<-function(N,Pmm = 0.1, nsigma = 2, Pmm_nonM = 0.6, nsigma_nonM = 3, verbose = FALSE){
  
  #Pmm has to be less than Pmm_nonM 
  if (Pmm >= Pmm_nonM){    
    cat("Warning: in function getSPIALimits the probability of mismatch in a matching population has to be less than the 
        mismatch probability in a mismatching population")    
    return(-1)
  }
  
  ##Gaussian approximationfor the binomial distribution B(N,Pmm)
  sigmasq_mm<-(N*Pmm*(1-Pmm)) 
  media_mm<-N*Pmm
  
  ##Gaussian approximationfor the binomial distribution B(N,Pmm_nonM)
  sigmasq_mm_nonM<-(N*Pmm_nonM*(1-Pmm_nonM))
  media_mm_nonM<-N*Pmm_nonM
  
  ##Compute limits
  k.inf <- media_mm+(nsigma*sqrt(sigmasq_mm))  
  k.sup <- max(media_mm_nonM-(nsigma_nonM*sqrt(sigmasq_mm_nonM)),0)
  
  
  if (k.inf > k.sup){
    cat("Warning: in function getSPIALimits the inf limit is greater than the sup limit.\n")
    cat("Warning: verify your input parameters:\n")
    cat("\t P mismatch in a matching population:", Pmm,"\n")
    cat("\t P mismatch in a mismatching population:", Pmm_nonM,"\n")
    cat("\t the limit for the mismatch in a matching population is ", k.inf,"\n")
    cat("\t the limit for the mismatch in a mismatching population is ", k.sup,"\n")
    return(list(liminf=k.inf/N,limsup=k.sup/N,err=TRUE ))
  }
  
  return(list(liminf=k.inf/N,limsup=k.sup/N,err=FALSE ))
}





############################################################################
#
# Compute SPIA test
#  Input: 1. x: a matrix with a column for each cell line and a row for each SNP
#            in the SPIA format (use toSPIAData before)
#         2. row.names: specify if the fisrt column contains the name of the SNPs
#         3. test.prob: specify if probabilistic test has to be performed
#         4. test.param: specify the parameters of the test
#            - test.param$Pmm: SNP probability of mismatch in a matching population (dafault 0.1)
#            - test.param$nsigma: area limit for Pmm
#            - test.param$Pmm_nonM: SNP probability of mismach in a non matching population (dafault 0.6)
#            - test.param$nsigma_nonM: area limit for Pmm_nonM
#            - test.param$PercValidCall: percentage of valid call to consider the test valid
#         5. verbose: print verbose information
#  Output: a matric with a line for each cell line and with columns with the
#          informationss about distances 
#
#############################################################################

SPIATest<-function (x, row.names = FALSE, test.prob = TRUE,  
                    test.param = list(Pmm = 0.1, nsigma = 2, Pmm_nonM = 0.6, nsigma_nonM = 3, PercValidCall=0.9), verbose = FALSE) 
{  
  
  #Print welcome message  
  cat("Input data:\n")
  
  
  #If row.names is set, the first column is removed from x
  if (row.names){
    rowNames <- x[,1]
    x <- x[,c(2:dim(x)[2])]
  }
  
  ##Nr is the number of SNP
  Nr <- nrow(x)
  ##Nc is the number of cell lines
  Nc <- ncol(x)  
  ##cuples gives the number of pairs: if there is 4 cell lines  
  ##       a total of (4 * 3 )/2 = 6 combinations are possible
  couples<-(Nc * (Nc - 1)/2)
  
  
  #if (verbose){
  #Print the statistics about the data
  cat(paste0("Number of samples: ",Nc,"\n"))
  cat(paste0("Number of pairs: ",couples,"\n"))
  cat(paste0("Number of SNPs: ",Nr,"\n"))
  #}
  
  if (test.prob){
    #Create the array that will store the results of the distance computing
    result<-array(,c(couples,13));  
    colnames(result)<-c("Sample_1","Sample_2","Distance","SPIA_Score","SNP_available","Total_SNP","One_SNP_NA","Bot_SNP_NA","Diff_AvsB_or_BvsA","Diff_AorBvsAB_or_vic","DiffABvsAorB","counterBothHomoz","counterBothHeter");
  }  else  {
    #Create the that store the results of the distance computing with the result of the statistical test
    result<-array(,c(couples,12));  
    colnames(result)<-c("Sample_1","Sample_2","Distance","SNP_available","Total_SNP","One_SNP_NA","Bot_SNP_NA","Diff_AvsB_or_BvsA","Diff_AorBvsAB_or_vic","DiffABvsAorB","counterBothHomoz","counterBothHeter");
  }
  
  #Verify if cell line names are available
  if (!is.null(colnames(x))){
    CLnames <- colnames(x)
  } else if (!is.null(names(x)) ) {
    CLnames <- names(x)
  } else {
    CLnames <- c(1:dim(x)[2])
  }      
  #Print dots when computing (each 20 pairs)
  if (couples < 20){
    steps <- 1
  } else {
    steps <- as.integer(couples/20) + 1
  }
  #The counter is used to trace the current pair of cell line processed by the for cycle
  counter<-0;
  for(i in (1:(Nc-1))){  
    for(j in ((i+1):Nc)){                 
      #Compute distance and other information for the cell i vs cell j
      dist<-ComputingSPIADistance(x[,i],x[,j],CLnames[i],CLnames[j],verbose=verbose)
      
      #Save the result in result and increment counter        
      result[(counter+1),1] <- CLnames[i]
      result[(counter+1),2] <- CLnames[j]
      
      #If the probabilistic test is not eables
      if (!test.prob){
        #save the current value of the distance
        result[(counter+1),c(3:12)]<-dist  
      }  else  {
        #save the current value of the distance and verify the probabilist test
        result[(counter+1),3]<-dist[1]
        #call the probabilistic test
        limits <- getSPIALimits(dist[2],test.param$Pmm,test.param$nsigma,test.param$Pmm_nonM,test.param$nsigma_nonM)
        #verify if there are problems with the limits defined by the binomial distributions
        if (limits$err)
        {
        #  message("SPIA error: the parameters of the probabilistic probabilistic test define ")
        #  message("  a negative 'similar' region.")
        #  message(paste("  The limit for the mismatch in a matching population is ",limits$liminf,sep=""))
        #  message(paste("  The limit for the mismatch in a mismatching population is ",limits$limsup,sep=""))
          return(-1)
        }
        
        #save the result of the probabilistic test
        if ( dist[2] / dist[3] < test.param$PercValidCall){
          result[(counter+1),4] <-"Not Valid"
        } else 
          if ( dist[1] < limits$liminf ){
            result[(counter+1),4] <-"Similar"
          } else 
            if ( dist[1] > limits$limsup ){
              result[(counter+1),4] <-"Different"
            } else {
              result[(counter+1),4] <-"Uncertain"
            }                    
        result[(counter+1),c(5:13)]<-dist[c(2:10)]
      }
      
      #print one dot each n steps
      if ((counter %% steps) == 0) { cat(".") }
      counter<-counter+1
    }
  }
  cat("\n")  
  return(list(SPIAresult = result, parameters = test.param, input.param = list(N_samples = Nc, N_SNPs = Nr, testDone = test.prob)))
}





#########################
#
# Plot functions
#
########################






############################################################################
#
# SPIAPlot
#  Input: 1. SPIAanalysis: the result of SPIATest function
#  Output: a plot of the cell line pairs versus the SPIA distance with 
#          information about the probabilistic test
#
#############################################################################

SPIAPlot <- function(SPIAanalysis, titleplot = NULL, colors = NULL, xlab=NULL, ...){  
  
  #if title is not defined
  if (is.null(titleplot)){
    titleplot<-paste("Pair-wise comparison of ", SPIAanalysis$input.param$N_samples," samples on ", SPIAanalysis$input.param$N_SNPs, " SNPs",sep="")  
  }
  
  if (is.null(colors)){
    colors <- c(rgb(54,144,192,maxColorValue=255),  #uncertain
                rgb(65,174,118,maxColorValue=255),  #similar
                rgb(215,48,31,maxColorValue=255),  #different
                rgb(37,37,37,maxColorValue=255)  #no computed
      )    
  }
  
  if(is.null(xlab)){
    xlab="Index of Sample Pairs"
  }
  
  baseCol <- rgb(37,37,37,maxColorValue=255)
  
  plot(NA,
       xlim=c(0,dim(SPIAanalysis$SPIAresult)[1]),
       ylim=c(0,1),
       axes=F,xlab="",ylab="",
       ...)
  
  axis(1,col=baseCol,col.axis=baseCol)
  axis(2,col=baseCol,col.axis=baseCol)
  title(main=titleplot,ylab="Distance D (% of discordant genotype calls)",xlab=xlab,col.main=baseCol,col.lab=baseCol)
         
  #check if the SPIAanalysis contain the probabilistic test
  if (SPIAanalysis$input.param$testDone)
    #compute the color depending on the result of the statistical test
  {      
    for(i in c(1:dim(SPIAanalysis$SPIAresult)[1])){
      pointConf <- switch(SPIAanalysis$SPIAresult[i,4], Uncertain = c(colors[1],19), Similar = c(colors[2],19), Different = c(colors[3],20), c(colors[4],19))    
      points(i,SPIAanalysis$SPIAresult[i,3],col=pointConf[1],pch=as.integer(pointConf[2]),cex=0.7);     
    }
    #plot the legend that describe the color code
    legendText <- c("SPIA TEST: uncertain", "SPIA TEST: match", "SPIA TEST: different", paste("< ",SPIAanalysis$parameters$PercValidCall,"% of available calls",sep=""))
    legend(0,1,legendText, col = colors, pch=c(19,19,20,19),cex=0.6)                            
  } else {
    #use the same color for each point
    for(i in c(1:dim(SPIAanalysis$SPIAresult)[1])){      
      points(i,SPIAanalysis$SPIAresult[i,3],col="black",pch=20,cex=0.7);     
    }
  }
      
}
