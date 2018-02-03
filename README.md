# iPAS-SNPs: Adaption of iPAS method to SNPs pathoginicity scores
This Repo created to adapt iPAS statistical method to SNPs pathoginicity scores to quantify the activity of a given pathway.

# Files
- generate-snps-data.R 
  - Simulates SNPs pathogenecity matrix data using `Scrim` package
  - Simulates SNPs using `TriadSim` package 

- generate-genes.R
  - Target : build a mapped list of SNPs to genes
  - by looking for closest genes in the vicinity of a SNP using this [Blog](http://jef.works/blog/2016/12/06/mapping-snps-and-peaks-to-genes-in-R/)

- snps-iPAD.R
  - Apply iPAS on SNPs data to quantify the activity of a given pathway
  
- iPAS_library.r
  - Its the iPAS library that applies the statistical analysis

- reactome2.txt
  - Its the sample file of pathways mapped to genes from REACTOME
  
# Current issues

> hits between genes and SNPs is very low `snps.genes.hits` so it seems that SNPs and gene annotations are not for the same reference build (ie. both shoud have coordinates with respect to hg19/GRCh37) , but SNPs generated are from this [ref](https://web.stanford.edu/class/gene210/files/data/culprit.txt)
  - In which I cant find its annotation gtf file !

> generation of our own pathways mapped to current genes, found this [ReactomePA package](https://bioconductor.org/packages/devel/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html)
  - BUT The output will NOT be similar to given pathway file  