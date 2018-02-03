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

# ---------------------------------------

# Define a few helper function

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
# function that maps snps to nearest genes
map.snps.to.genes <- function(snps.info)
{
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
  
  # clean snps data 
  snps.pos <- as.numeric( as.matrix(snps.info$position))
  snps.chr <- as.numeric( as.matrix(snps.info$chr))
  
  # create range with epositions of current snps
  snps.pos.range <- data.frame("chr"= snps.chr,"startRange"=snps.pos,"endRange"=snps.pos)
  head(snps.pos.range)
  
  # convert SNPs to GRanges
  snps.granges <- range2GRanges(snps.pos.range)
  names(snps.granges) <- snps.info$snpname
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
  overlap
  #Hits object with 2 hits and 0 metadata columns:
  #  queryHits subjectHits
  # <integer>   <integer>
  #  [1]         1          51
  # [2]         1          52
  #-------
  #  queryLength: 120 / subjectLength: 69
  
  slotNames(overlap)
  #[1] "from"            "to"              "nLnode"          "nRnode"         
  #[5] "elementMetadata" "metadata"    
  
  # where ,
  # from == "queryHits"
  # to == "subjectHits"
  
  # make vector of SNPs to genes
  target.genes <- names(r2)[slot(overlap, "to")]
  target.snps <-names(r1)[slot(overlap, "from")]
  
  # check results
  mapped.hits <- as.matrix(cbind(target.snps,target.genes),colnames("SNPs","genes"))
  View(mapped.hits)
  
  return(mapped.hits)
}

hits <- map.snps.to.genes(snps.info)
# Output
#     target.snps target.genes   
#[1,] "rs3094315" "RP11-206L10.9"
#[2,] "rs3094315" "RP11-206L10.8"


# To map genes to pathways using Reactom , Check this package
# https://bioconductor.org/packages/devel/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html

