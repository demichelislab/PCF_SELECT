###########################
#
# Configuration file for SPIA
#
###########################

#### Parameters for the SPIA statistical test

## Mean and standard deviation of the probability
## that two sample that match have a SNPs with different genotype
Pmm = 0.1 
nsigma = 2 

## Mean and standard deviation of the probability
## that two sample that do not match have a SNPs with different genotype
Pmm_nonM = 0.6
nsigma_nonM = 5

## Minimum percentage of valid calls to perform the SPIA statistical test
PercValidCall=0.7

########
#### Output files

## Save SPIA plot
saveSPIAplot = T

## Save SPIA genotype (for debugging)
saveGenotype = T

## SPIA plot file name
genotypeTable_file = ""

########
#### Other parameters

#print verbose information
verbose = F

#print output on screen (if False it create a log file)
print_on_screen = T
