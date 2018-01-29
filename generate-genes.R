# Target : build a mapping of SNPs to pathways 
#  - by first looking for closest genes in the vicinity of a SNP, 
#  - or via eQTL mapping (expression quantitative trait loci)
# ``````````````````````````````````````````````````
setwd("~/Desktop/Rs-Work/ipas-snps")
library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
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

# Define a few helper functions

#' Convert from string to range
#' 
#' @param pos A vector of strings ex. chr1 2938302 2938329
#' @param delim Delimiter for string splitting
#' @param region Boolean of whether region or just one position
#'
#' @returns Dataframe of ranges
#' 
string2range <- function(pos, delim=' ', region=TRUE) {
  posp <- as.data.frame(do.call(rbind, strsplit(pos, delim)))
  posp[,1] <- posp[,1]
  posp[,2] <- as.numeric(as.character(posp[,2]))
  if(region) {
    posp[,3] <- as.numeric(as.character(posp[,3]))
  } else {
    posp[,3] <- posp[,2]
  }
  return(posp)
}

#' Convert from ranges to GRanges
#' 
#' @param df Dataframe with columns as sequence name, start, and end
#' 
#' @returns GRanges version 
#' 
range2GRanges <- function(df) {
  require(GenomicRanges)
  require(IRanges)
  gr <- GenomicRanges::GRanges(
    seqnames = df[,1],
    ranges=IRanges(start = df[,2], end = df[,3])
  )
  return(gr)
}

# ---------------------------------------

# convert SNPs to GRanges
snps.ranges <- string2range(snps, delim=":", region=FALSE)
head(snps.ranges)

snps.granges <- range2GRanges(snps.ranges)
names(snps.granges) <- snps
head(snps.granges)

# convert genes to GRanges
gtf.granges <- range2GRanges(gtf.gene)
names(gtf.granges) <- gene.names
head(gtf.granges)

# ---------------------------------------

# Now that we have our two GRanges objects,
# we can easily overlap them using GenomicRanges::findOverlaps.

r1 <- snps.granges
r2 <- gtf.granges
overlap <- GenomicRanges::findOverlaps(r1, r2)

# explor the output
slotNames(overlap)
#[1] "from"            "to"              "nLnode"          "nRnode"         
#[5] "elementMetadata" "metadata"    

# where ,
# from == "queryHits"
# to == "subjectHits"

# make vector of SNPs to genes
#hits <- names(r2)[slot(overlap, "subjectHits")]
hits <- names(r2)[slot(overlap, "to")]

#names(hits) <- names(r1)[slot(overlap, "queryHits")]
names(hits) <- names(r1)[slot(overlap, "from")]

hits <- as.matrix(hits)
View(hits)








