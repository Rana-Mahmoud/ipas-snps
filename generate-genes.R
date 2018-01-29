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

# Try documentation example I want to map a set of SNPs to genes.

# Sample snps
snps <- c("1:11873", "1:69100", "1:752761")

# Sample genes from GTF
gtfFile <- '~/Desktop/Rs-Work/Homo_sapiens.GRCh37.75.gtf'
gtf <- read.table(gtfFile, header=F, stringsAsFactors=F, sep='\t', nrows=1000) # limit number of rows for testing
gtf.gene <- gtf[gtf[,3]=="gene", c(1,4,5)]
gene.names <- unlist(lapply(gtf[gtf[,3]=="gene", 9], function(x) {
  y <- strsplit(x, ';')[[1]][2]
  gsub(' gene_name ', '', y)
}))
rownames(gtf.gene) <- gene.names
# contains (gene.name , chromosome , start.pos.range, end.pos.range)
head(gtf.gene)

# ---------------------------------------

# 



















