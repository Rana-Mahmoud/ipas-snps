# This file will contain script that apply iPAS on SNPs data to quantify the activity of a given pathway
setwd("~/Desktop/Rs-Work/ipas-snps")

# First : lets simulate SNPs data
#---------------------------------
snps.sourcePath = "generate-snps-data.R"
source(snps.sourcePath)

# number of target snps 
# "NB : Must be divisable by 4"
n.snps <- 120

# get pathogenecity matrix
d.snps.pathogenicity <- get.Snps.matrix(n.snps) # diseased
n.snps.pathogenicity <- get.Snps.matrix(n.snps) # normal
#colnames(d.snps.pathogenicity)

# get snps names
snps.info <- as.data.frame( get.Snps.list(n.snps))
#head(snps.names)

# update matric snps names
colnames(d.snps.pathogenicity) = snps.info$snpname
colnames(n.snps.pathogenicity) = snps.info$snpname

# Second : Lets map SNPs to genes 
#----------------------------------
snps.genes.hits <- map.snps.to.genes(snps.info)
snps.genes.hits 
# !!! # issue > number of hits is very small so it seems that
#  SNPs and gene annotations are not for the same reference build 
# (ie. both shoud have coordinates with respect to hg19/GRCh37)
# but SNPs generated are from :
# https://web.stanford.edu/class/gene210/files/data/culprit.txt
# which i cant find its annotation gtf file !

## then column names of the SNPs matricies should be renamed with corresponding gene to apply iPAS and rum correctly 

# Third : Lets Set iPAS inputs 
#--------------------------------
# diseased snps matrix
#tumorFile = "Beer_10samples.txt"
diseased.snps = t(d.snps.pathogenicity)
View(d.snps.pathogenicity)

# normal/ref snps matrix
#refFile = "LUAD_nREF.txt"
ref.snps = t(n.snps.pathogenicity)
View(n.snps.pathogenicity)

# have 2 options : 1st try the same reactoms file
pathwayListFile= "reactome2.txt"
# option 2 : generat our own pathways mapped to current genes 
# check this link 
# https://bioconductor.org/packages/devel/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html
# !!! BUT $ The output will not be similar to given pathway file 


# set modes
whoami = "snpstestRun1"
runMode = "nRef"

# keep results in the same repo 
bgDistributionPath = ""

# Source iPAS library
sourcePath = "iPAS_library.r"
source(sourcePath )

# Apply Method
#iPAS( tumorFile, refFile, pathwayListFile, whoami, runMode, bgDistributionPath  )
iPAS( diseased.snps, ref.snps, pathwayListFile, whoami, runMode, bgDistributionPath  )
