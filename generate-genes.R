# this file generates genes related to snps list
# ``````````````````````````````````````````````````
# First : lets simulate SNPs data
#---------------------------------
snps.sourcePath = "generate-snps-data.R"
source(snps.sourcePath)

# number of target snps
n.snps <- 150

# get pathogenecity matrix
snps.pathogenicity <- get.Snps.matrix(n.snps)
#colnames(snps.pathogenicity)
  
# get snps names
snps.names <- get.Snps.list(n.snps)
#head(snps.names)
  
# update matric snps names
colnames(snps.pathogenicity) = snps.names