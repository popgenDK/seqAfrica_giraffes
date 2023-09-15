# /usr/bin/env R
require(Relate)

#### geno
##geno is the genotypes
## genotypes coded as NA,0,1,2
## N x M
## N number of individuals
## M number of sites

#### maf
## minimum maf

#### mis
## maximum missingness
LDdecay <- function(geno,pos,chr,maf=0.05,mis=0.05,depth=50){

    if(missing(pos))
        pos <- 1:ncol(geno)
    if(missing(chr))
        chr <- rep(-1,length(pos))
    
    
    freq <- colMeans(geno,na.rm=T)/2
    missingness <-  colMeans(is.na(geno))

    keepSNP <- freq>maf & freq < 1-maf & missingness < mis

    if(sum(keepSNP)<2)
        stop(paste("only",sum(keepSNP)," SNPs passed filter\n"))

    keepPos <- pos[keepSNP]
    keepChr <- chr[keepSNP]
    ld <- ld.snp3(geno[,keepSNP]+1,back=depth)

    snpR2 <-rowMeans(ld$rmisc,na.rm=T)
    res <- reform(ld$rmisc,keepPos,keepChr)

    
    list(r2=res,ldResults=ld,positions=keepPos,keepSNP=keepSNP,SNPr2=snpR2)
}

reform<-function(x,pos,chr){

    depth <-nrow(x)
    M <- length(pos)
    if(missing(chr))
        chr<- rep(-1,length(pos))
   res <- matrix(NA,ncol=4,nrow=length(x))

    colnames(res)<-c("r2","pos1","pos2","chr")

    res[,1] <- as.vector(x)
    res[,2] <- rep(pos[-M],each=depth)
    ord2<- rep(2:M-1,each=depth)+1:depth
    res[,3] <- pos[ord2]
    res[,4] <- rep(chr[-M],each=depth)
    sameChr <- res[,4] == chr[ord2]
    keep <- sameChr & !is.na(res[,"pos2"])
    
    res[keep,]

}


makeBin <- function(x,max=1,n=100){

    res <- x$r2
    seq <- seq(0,max,by=max/n)
    bin<-cut(res[,3]-res[,2],seq)
    r2bin<-tapply(res[,1],bin,mean,na.rm=T)
    list(r2bin=r2bin,seq=seq[-1])
}
