# This file will contain script that apply iPAS on SNPs data to quantify the activity of a given pathway

# Set inputs 
#--------------
# diseased snps matrix
tumorFile = "Beer_10samples.txt"
diseased.snps = t()

# normal/ref snps matrix
refFile = "LUAD_nREF.txt"
ref.snps = t()

# have 2 options : 1st try the same reactoms file
pathwayListFile= "reactome2.txt"
# option 2 : generat our own pathways mapped to current genes 
# check this link 
# https://bioconductor.org/packages/devel/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html

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
