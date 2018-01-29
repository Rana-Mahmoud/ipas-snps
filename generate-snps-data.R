# Here we will do 2 steps :
# 1- generate SNPs pathogenecity matrix data using "Scrim" package
# 2- generate SNPs using "TriadSim" packages 
# ----------------------------
# Load needed libraries
library(scrime)
library(TriadSim)

# Step 1 : SNPs matrix
get.Snps.matrix <- function(n.snps)
{
  # Simulate a data set containing 2000 observations (1000 cases
  # and 1000 controls) and 100 SNPs, where one three-way and two 
  # two-way interactions are chosen randomly to be explanatory 
  # for the case-control status.
  obs <- 2000 # number of observations
  #n.snps <- 100 # number of snps
  intr.vec <- c(3, 2, 2) # interactions 
  # call genrate function
  sim.snps <- simulateSNPs(obs , n.snps , intr.vec )
  # matrix of data 
  #View(sim.snps$data)
  return(sim.snps$data)
}

# Step 2 : generate SNPs 
get.Snps.list <- function(n.snp)
{
  curr.dir <- getwd()
  # access results and take subset of snps names = n.snp
  target.path <- paste0(curr.dir,"/snps-results/triad1.bim")
  # check if file exist
  if (!file.exists(target.path))
  {
    m.file <- paste(system.file(package = "TriadSim"),'/extdata/pop1_4chr_mom',sep='')
    f.file <- paste(system.file(package = "TriadSim"),'/extdata/pop1_4chr_dad',sep='')
    k.file <- paste(system.file(package = "TriadSim"),'/extdata/pop1_4chr_kid',sep='')
    input.plink.file <- c(m.file, f.file, k.file)
    
    
    #, out.put.file=paste(curr.dir,"/snps-results/",'triad',sep='')
    
    TriadSim(input.plink.file
             , out.put.file=paste(curr.dir,"/snps-results/",'triad',sep='')
             , fr.desire=0.05
             , pathways=list(1:4,5:8)
             , n.ped=1000
             , N.brk=3
             , target.snp=NA
             , P0=0.001
             , is.OR=FALSE
             , risk.exposure= 1
             , risk.pathway.unexposed=c(1.5, 2)
             , risk.pathway.exposed=c(1.5, 2)
             , is.case=TRUE
             , e.fr=NA
             , pop1.frac=NA
             , P0.ratio=1
             , rcmb.rate
             , no_cores=NA)
    
  }
  # list of current snps chromosomes
  chromosomes.list <- c("1","8","17","20")
  # quarter the n.snps
  q.n.snps <- n.snps/4
  # create empty matrix of (snps.name, chromosome.number , position)
  snps.mat <- matrix(0,nrow = 1,ncol = 3)
  colnames(snps.mat) <-  c("snpname","chr","position")
  head(snps.mat)
  for (i in chromosomes.list)
  {
    curr.chr.path <- paste0(curr.dir,"/snps-results/triad",i,".bim")
    snps.file <- read.table(curr.chr.path)
    head(snps.file)
    # subset quarter the n.snps
    q.snps.info <- snps.file[1:q.n.snps,c(2,1,4)]
    colnames(q.snps.info) <- c("snpname","chr","position")
    head(q.snps.info)
    # append to main snpa list
    snps.mat <- rbind(snps.mat,q.snps.info)
    head(snps.mat)
  }
  #dim(snps.mat)
  #head(snps.mat)
  snps.info <- as.matrix(snps.mat[1:n.snps+1,])
  #head(snps.info)
  #tail(snps.info)
  #dim(snps.info)
  return(snps.info)
}

##### Generate SNPS list and SNPs pathogenicity scores matrix
# Needed number of snps
#n.snps <- 100 # number of snps


#simulate.snps.data <- function(n.snps = 100)
#{
  # get pathogenecity matrix
#  snps.pathogenicity <- get.Snps.matrix(n.snps)
  #colnames(snps.pathogenicity)
  
  # get snps names
#  snps.names <- get.Snps.list(n.snps)
#  head(snps.names)
  
  # update matric snps names
#  colnames(snps.pathogenicity) = snps.names
  
#  return(c(snps.names,snps.pathogenicity))
#}

