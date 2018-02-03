######################################################################################
# Individualized Pathway Aberrance Score
# Author : TaeJin Ahn, Eunjin Lee
# Date : 2013.11.25
# Contact : taejin.ahn@gmail.com
######################################################################################

library(MASS)
library(preprocessCore)

doNormalization <- function(data, ref){

	all.data = c()
	
	if(identical(ref,data)){
		print("iPAS preprocessing: quantile normalization of all reference samples" )

		tmp <- data
		q.tmp <- normalize.quantiles(as.matrix(tmp)); rownames(q.tmp) <- rownames(tmp)
		z.tmp <- t(scale(t(q.tmp)))
		all.data <- as.matrix(z.tmp) ; 
		colnames(all.data) <- colnames(data)
	}else {
		print("iPAS preprocessing: proceeding 1 (one tumor) + N (normal references) quantile normalization" )
	
		all.data <- c()
		for( i in 1:ncol(data)){
			tmp <- cbind(data[,i],ref); 	
			id <- colnames(data)[i];
			q.tmp <- normalize.quantiles(as.matrix(tmp)); 
			rownames(q.tmp) <- rownames(tmp)
				
			normalQ = q.tmp[,-1] 
			#the first colum is cancer, second to the last are normal
			u = apply( normalQ, 1, mean)
			sd = apply( normalQ, 1, sd )

			z.tmp = (q.tmp - u ) / sd
			all.data <- cbind(all.data, as.matrix(z.tmp[,1])); 
			colnames(all.data)[i] <- id
			
			if( i%%10 ==0 ) {
				print( paste( "iPAS preprocessing: " , i , " out of ", ncol(data), " samples have finished", sep="" ) ) 
			}
		}
	}	
	return(all.data)
}

iPASpreprocessing <- function(tumors, refs,runMode="nRef"){
	z.data = c();
	print("#iPAS preprocessing#")
	if(runMode =="cohort") {
		print("non reference based normalization, using entire test cohort data set")
		
		rownames(tumors) <- tumors[,1] ; 
		tumors <- tumors[,-1]
		
		#using entire cohort for normalization
		tumors_ <- doNormalization(tumors, tumors) 
		refs_ <- tumors_
		
		z.data <- list(refs_, tumors_); 
		names(z.data) <- c("reference", "tumors")

		return( z.data )
	} else{
		
		rownames(refs) <- refs[,1] ; refs <- refs[,-1]
		rownames(tumors) <- tumors[,1] ; tumors <- tumors[,-1]

		refs_ <- doNormalization(refs, refs)
		tumors_ <- doNormalization(tumors, refs)
	
		z.data <- list(refs_, tumors_); 
		names(z.data) <- c("reference", "tumors")

		return( z.data )
	}
	return ( z.data)
}

iPAScalculation <- function( tumors, reference, gene.list , pathwayName, bgDistributionPath){

	GS.ref <- reference[na.omit(match(gene.list, rownames(reference))),]
	GS.tumors <- tumors[na.omit(match(gene.list, rownames(tumors))),]
	
	centroid <- rowMeans(GS.ref) ; 
	normal.cov <-  cov(t(GS.ref))

	geneset.size = dim(GS.ref)[1]

	if( nchar (bgDistributionPath ) == 0  ){			
		AveZ_comp <- doAveZ(GS.ref)
		ED_comp <- doEuclidean(GS.ref, centroid, geneset.size) 
		MD_comp <- doMahalanobis(GS.ref, centroid, normal.cov, geneset.size) 
		
		normal.train.fisher.ref <- c(); normal.train.bin.ref <- c(); normal.train.GSEA.ref <- c()
		for( i in 1:ncol(reference)){
			normal.train.fisher.ref <- c(normal.train.fisher.ref, doFisherScore(reference[,i],gene.list))
			normal.train.GSEA.ref <- c(normal.train.GSEA.ref, getESfastest(reference,gene.list,i))
		}
		write.table(AveZ_comp, file=paste("AVGZ_bg_",pathwayName, ".txt" ), sep="\t", row.names=F, quote=F)
		write.table(ED_comp, file=paste("ED_bg_",pathwayName, ".txt") , sep="\t", row.names=F, quote=F)
		write.table(MD_comp, file=paste("MD_bg_",pathwayName, ".txt") , sep="\t", row.names=F, quote=F)
		write.table(normal.train.fisher.ref, file=paste("Fisher_bg_",pathwayName, ".txt") , sep="\t", row.names=F, quote=F)
		write.table(normal.train.GSEA.ref, file=paste("GSEA_bg_",pathwayName, ".txt") , sep="\t", row.names=F, quote=F)
	}else{
		#print("loading precalculated background distribution")	
		currentWD = getwd()
		setwd(bgDistributionPath)
		AveZ_comp = as.numeric( unlist( read.table(paste("AVGZ_bg_",pathwayName, ".txt" ), sep="\t", header=T) ) )
		ED_comp = as.numeric( unlist( read.table(paste("ED_bg_",pathwayName, ".txt" ), sep="\t", header=T) ) )
		MD_comp = as.numeric( unlist(read.table(paste("MD_bg_",pathwayName, ".txt" ), sep="\t", header=T) ) )
		normal.train.fisher.ref = as.numeric( unlist( read.table(paste("Fisher_bg_",pathwayName, ".txt" ), sep="\t", header=T) ) )
		normal.train.GSEA.ref = as.numeric( unlist( read.table(paste("GSEA_bg_",pathwayName, ".txt" ), sep="\t", header=T) ) )
		setwd(currentWD)
	}

	all.tumor <- c()
	for( i in 1:ncol(tumors)){
		test.tumor.all <- as.matrix(tumors[,i])
		test.tumor.data <- as.matrix(GS.tumors[,i])
	
		test.tumor.aveZ.p <- doAveZ.P(test.tumor.data,AveZ_comp)
		test.tumor.fisher.p <- doFisher.P(test.tumor.all, gene.list, normal.train.fisher.ref)
		test.tumor.ED.p <- doEuclidean.P(test.tumor.data, centroid, ED_comp, geneset.size)
		test.tumor.MD.p <- doMahalanobis.P(test.tumor.data, centroid, normal.cov, MD_comp, geneset.size)
		test.tumor.GSEA.p <- doGSEA.P(test.tumor.all,gene.list,1,normal.train.GSEA.ref)
		
		result <- t(as.matrix(c(test.tumor.aveZ.p,test.tumor.fisher.p,  test.tumor.ED.p, test.tumor.MD.p, test.tumor.GSEA.p)))
		rownames(result) <- colnames(tumors)[i]

		all.tumor <- rbind(all.tumor, result)
	}
	
	all.result <- list(all.tumor, normal.cov);
	names(all.result) <- c("tumor","normal.cov")
	return(all.result)
}	

doAveZ <- function(all.data){ 
	aveZscore <- colMeans(all.data); 
	return(aveZscore)
	}
	
doAveZ.P <- function(data, comp){
	aveZ <- doAveZ(data)
	aveZ.p.n <- pnorm(aveZ, mean=mean(comp),sd=sd(comp),lower.tail=FALSE) ; if(aveZ.p.n > 0.5) aveZ.p.n <- 1-aveZ.p.n
	tmp <- c(aveZ, comp);
	aveZ.p.e <- 1-(rank(tmp)[1]/length(tmp)) ; if(aveZ.p.e > 0.5) aveZ.p.e <- 1-aveZ.p.e
	aveZ.z <- (aveZ - mean(comp)) / sd(comp)
	aveZ.result <- c(aveZ, aveZ.p.n, aveZ.p.e); 
	names(aveZ.result) <- c("AveZ.score", "AveZ.normalDist.P", "AveZ.empirical.P")
	return(aveZ.result)
}

	
doFisherScore <- function(all.data,gene.list){
	all.data <- na.omit(all.data)
	n.ALL <- length(all.data)
	p.tmp <- pnorm(all.data) ; for( i in 1:length(p.tmp)){if(p.tmp[i]>0.5) { p.tmp[i] <- 1-p.tmp[i]}} 
	n.FG <- length(p.tmp[p.tmp<0.05])
	all.data <- as.matrix(all.data)
	data <- as.matrix(all.data[na.omit(match(gene.list,rownames(all.data))),])
	n.inNetG <- nrow(data)
	p.data <- pnorm(data) ; for(t in 1:nrow(p.data)){ if(p.data[t,] > 0.5) p.data[t,] <- as.numeric(1-p.data[t,]) }
	n.FG_inNetG <- length(p.data[which(p.data<0.05),])
	test.table <- t(matrix(c(n.FG_inNetG, n.FG-n.FG_inNetG, n.inNetG-n.FG_inNetG, n.ALL-n.FG-(n.inNetG-n.FG_inNetG)), nrow=2))
	dimnames(test.table) <- list(c("Focus","~Focus"),c("InNet","~InNet"))
	fisher.p <- fisher.test(test.table,alternative="greater")$p.value
	return(fisher.p)
}

doFisher.P <- function(all.data, gene.list, fisher.ref){
	fisher.p <- doFisherScore(all.data, gene.list)

	tmp <-  c(fisher.p,fisher.ref)
	fisher.p.e <- rank(tmp)[1]/length(tmp)

	fisher.score <- -log10(fisher.p) ;if(fisher.score==Inf) {fisher.score <- 30}
	fisher.ref1 <- -log10(fisher.ref); 
	fisher.ref2 <- fisher.ref1[fisher.ref1!=Inf]; fisher.ref1[fisher.ref1==Inf] <- 30; 
	fisher.p.n <- pnorm(fisher.score, mean=mean(fisher.ref2),sd=sd(fisher.ref1),lower.tail=FALSE)
	fisher.z = ((fisher.score) - mean(fisher.ref2)) / sd(fisher.ref1)

	fisher.result <- c(fisher.p, fisher.p.n, fisher.p.e); names(fisher.result) <- c("fisher.score","fisher.normalDist.P","fisher.empirical.P")
	return(fisher.result)
}

doEuclidean <- function(data, centroid, size){
	ED <- c()
	for( t in 1:ncol(data)){ ED[t] <- dist(rbind(data[,t], centroid),method="euclidean")}
	names(ED) <- colnames(data)
	ED = ED / size
	return(ED)
	}

doEuclidean.P <- function(data, centroid, comp, size){
	ED <- doEuclidean (data, centroid, size)
	ED.p.n <- pnorm(ED, mean=mean(comp),sd=sd(comp),lower.tail=FALSE) 
	tmp <-  c(ED,comp); 
	ED.p.e <- 1-(rank(tmp)[1]/length(tmp))
	ED.z = (ED - mean(comp)) / sd(comp)
	
	ED.result <- c(ED, ED.p.n, ED.p.e); names(ED.result) <- c("ED.score", "ED.normalDist.P","ED.empirical.P")
	return(ED.result)
}

doMahalanobis <- function(data, centroid, covar, size){
	#MD <- mahalanobis(t(data), centroid, covar)
	MD <- mahalanobis(t(data), centroid, ginv(covar), inverted=TRUE)
	MD = sqrt(MD) / size
	return(MD)
}

doMahalanobis.P <- function(data, centroid, covar, comp, size){
	MD <- doMahalanobis ( data, centroid, covar, size)
	MD.p.n <- pnorm(MD, mean=mean(comp),sd=sd(comp),lower.tail=FALSE) 
	tmp <-  c(MD,comp); 
	MD.p.e <- 1-(rank(tmp)[1]/length(tmp))
	MD.z = (MD - mean(comp) )/ sd(comp)
	
	MD.result <- c(MD, MD.p.n, MD.p.e); names(MD.result) <- c("MD.score","MD.normalDist.P","MD.empirical.P")
	return(MD.result)
}

getESfastest = function(data, gene.set, test_col_idx) {
	sorted <- as.data.frame(data[ order(data[,test_col_idx], decreasing=TRUE), ])
	data1 <- cbind(rownames(sorted),sorted[,test_col_idx]);
	data1 <- as.data.frame(data1); 
	data1[,2] <- as.numeric(as.character(data1[,2])); 
	colnames( data1 ) = c("gene","stat") 

	Nr_tmp <- as.numeric(as.character(data1[data1$gene%in%gene.set,]$stat))
	Nr = sum( abs(Nr_tmp) )
	N = dim(data1)[1]
	Nh = length(gene.set)

	maxES = -10000
	hitIdxs = which( data1$gene%in%gene.set ==TRUE )
	for( k in 1:length(hitIdxs) ){
		hit_cursor = hitIdxs[1:k]
		total_cnt_of_miss = (hitIdxs[k] - 1) - ( length(hitIdxs[1:k]) -1 )
		Phit = sum( abs(data1[hit_cursor,2]) ) / Nr
		Pmiss =  total_cnt_of_miss * ( 1 / (N -Nh) )
		ES = Phit - Pmiss
		if( ES > maxES ) {maxES = ES}
	}
	return (maxES)
}

doGSEA.P <- function(data, gene.set, test_col_idx, normal.GSEA.ref){
	test.score <- getESfastest(data, gene.set, test_col_idx)
	GSEA.p.n <- pnorm(test.score, mean=mean(normal.GSEA.ref),sd=sd(normal.GSEA.ref),lower.tail=FALSE)
	tmp <-  c(test.score,normal.GSEA.ref)
	GSEA.p.e <- 1-(rank(tmp)[1]/length(tmp))
	GSEA.z = ( test.score - mean( normal.GSEA.ref) )/ sd(normal.GSEA.ref)
	
	GSEA.result <- c(test.score, GSEA.p.n, GSEA.p.e); 
	names(GSEA.result) <- c("GSEA.score","GSEA.normalDist.P","GSEA.empirical.P")
	return(GSEA.result)
}
	


# replaced with direct input matrix of snps data
#iPAS = function( tumorFile, refFile, pathwayListFile, whoami, runMode , bgDistributionPath )
iPAS = function( all.tumor, all.normal, pathwayListFile, whoami, runMode , bgDistributionPath ){

	outputName = paste( whoami, runMode, Sys.Date(), sep="_" )

	pathway.data <- read.delim( pathwayListFile, header=F,sep="\t")
	
	n.pathway= nrow(pathway.data)

	# replaced with direct input matrix
	#all.normal <- read.delim( refFile , header=T)
	#all.tumor <- read.delim( tumorFile , header=T, sep="\t")

	colnames(all.tumor)[1] <- "Gene.Symbol"
	gene.list.tmp <- intersect(all.normal$Gene.Symbol, all.tumor$Gene.Symbol)
	gene.list <- as.data.frame(gene.list.tmp) ; colnames(gene.list) <- "Gene.Symbol"
	all.tumor <- merge(gene.list, all.tumor)
	ref.normal <- merge(gene.list, all.normal)

	#set output path
	write.path <- paste(getwd(), "/", outputName, "_iPAS",sep=""); 
	dir.create(write.path)

	#preprocessing individual tumor samples
	result <- iPASpreprocessing(all.tumor, ref.normal , runMode)
	reference <- result$reference
	tumors <- result$tumors


	for( i in 1:n.pathway){

		tmp <- t(pathway.data[i,])
		pathway.id <- tmp[1]
		
		print( paste("iPAS ", i, "/", n.pathway," : " , pathway.id , sep="") )
		p.pathway.id = pathway.id
		#handling bug : if pathway id is too long, some machines can not write the result
		#cutting the pathway name which is longer than 70 charcters
		if( nchar( pathway.id) > 70 ){
			p.pathway.id = substring(pathway.id,1,70)
		}
		
		gene.list.tmp <- tmp[-c(1,2)]	
		gene.list <- gene.list.tmp[gene.list.tmp!=""]

		result <- iPAScalculation (tumors, reference, gene.list , pathway.id, bgDistributionPath  )
				
		t.result <- as.matrix(result$tumor); 
		t.result <- cbind(gsub(".rma.Signal","",rownames(t.result)), t.result);
		rownames(t.result) <- c(); 
		t.result = cbind( p.pathway.id, t.result)
		colnames(t.result)[1] <- "Pathway.ID"
		colnames(t.result)[2] <- "Sample.ID"
		if ( i == 1 ){
			write.table(t.result, file=paste(write.path,"/iPASresult.txt",sep=""),row.names=FALSE, sep="\t",quote=FALSE, append=F)
		}else{
			write.table(t.result, file=paste(write.path,"/iPASresult.txt",sep=""),row.names=FALSE, col.names=F, sep="\t",quote=FALSE, append=T)
		}
	}
}
