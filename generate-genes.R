# Target : build a mapping of SNPs to pathways 
#  - by first looking for closest genes in the vicinity of a SNP, 
#  - or via eQTL mapping (expression quantitative trait loci)
# ``````````````````````````````````````````````````
# First : lets simulate SNPs data
#---------------------------------
snps.sourcePath = "generate-snps-data.R"
source(snps.sourcePath)

# number of target snps 
# "NB : Must be divisable by 4"
n.snps <- 120

# get pathogenecity matrix
snps.pathogenicity <- get.Snps.matrix(n.snps)
#colnames(snps.pathogenicity)
  
# get snps names
snps.info <- as.data.frame( get.Snps.list(n.snps))
#head(snps.names)
  
# update matric snps names
colnames(snps.pathogenicity) = snps.info$snpname
# ````````````````````````````````````````````````````
# Secound : generate genes from current snps
#--------------------------------------------
# map snps to genes 
# http://jef.works/blog/2016/12/06/mapping-snps-and-peaks-to-genes-in-R/

# Try documentation example 





















