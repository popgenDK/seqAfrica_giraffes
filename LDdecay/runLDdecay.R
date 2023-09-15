# LDdecay curves for the giraffe populations

# description https://github.com/aalbrechtsen/LDdecay
# https://htmlpreview.github.io/?https://github.com/aalbrechtsen/LDdecay/blob/master/LDdecay-exported.html

DATA=data

# PLINK v1.90b6

# thin data to reduce comp time
plink --bfile $DATA --thin 0.1 --chr 1 --out chr1 --make-bed

library(Relate)
#BiocManager::install("snpStats")
library(snpStats)
source("/home/albrecht/Rfun/LDdecay.R")

if(FALSE){
    #to install
    library(devtools)
    install_github("aalbrechtsen/relate")
}

plink<-
function(plinkFile){
    pl <- snpStats::read.plink(plinkFile)
    ind<-pl$fam[,1]
    snp<-pl$map[,2]
    geno<-as.integer(as.integer(pl$genotypes)-1)
    dim(geno)<-c(length(ind),length(snp))
    geno[geno==-1]<-NA
    rownames(geno)<-ind
    colnames(geno)<-colnames(pl)
    bim<-read.table(paste0(plinkFile,".bim"),as.is=T,header=F)
    fam<-read.table(paste0(plinkFile,".fam"),as.is=T,header=F)
    list(geno=geno,bim=bim,fam=fam,pl=pl)
}

plinkFile <- "/path/to/"
pl <- plink(plinkFile)
popInfo <- pl$fam[,2]
table(popInfo)

# smallest sample size Kordofan=4
# run each pop with n=4

#### Nubian n=8
NUBgenoAll <- pl$geno[popInfo=="Nubian",]
NUBgeno <- head(NUBgenoAll,4)
NUBRes <- LDdecay(NUBgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)
NUBBin <- makeBin(NUBRes,max=15,n=1000)

#### Kordofan n=4
KORgenoAll <- pl$geno[popInfo=="Kordofan",]
KORgeno  <- head(KORgenoAll,4)
KORRes <- LDdecay(KORgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)
KORBin <- makeBin(KORRes,max=15,n=1000)

#### Masai n=6
MASgenoAll <- pl$geno[popInfo=="Masai",]
MASgeno <- head(MASgenoAll,4)
MASRes <- LDdecay(MASgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)
MASBin <- makeBin(MASRes,max=15,n=1000)

#### Southern African Central n=8
SACgenoAll <- pl$geno[popInfo=="Southern_African_Central",]
SACgeno <- head(SACgenoAll,4)
SACRes <- LDdecay(SACgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)
SACBin <- makeBin(SACRes,max=15,n=1000)

#### Angolan n=12
ANGgenoAll <- pl$geno[popInfo=="Angolan",]
ANGgeno <- head(ANGgenoAll,4)
ANGRes <- LDdecay(ANGgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)
ANGBin <- makeBin(ANGRes,max=15,n=1000)

#### Masai_Selous n=5
MSEgenoAll <- pl$geno[popInfo=="Masai_Selous",]
MSEgeno <- head(MSEgenoAll,4)
MSERes <- LDdecay(MSEgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)
MSEBin <- makeBin(MSERes,max=15,n=1000)

#### Reticulated n=12
RECgenoAll <- pl$geno[popInfo=="Reticulated",]
RECgeno <- head(RECgenoAll,4)
RECRes <- LDdecay(RECgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)
RECBin <- makeBin(RECRes,max=15,n=1000)

#### Masai_Thornicrofts n=5
MTHgenoAll <- pl$geno[popInfo=="Masai_Thornicrofts",]
MTHgeno <- head(MTHgenoAll,4)
MTHRes <- LDdecay(MTHgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)
MTHBin <- makeBin(MTHRes,max=15,n=1000)

#### Southern_African n=12
SOAgenoAll <- pl$geno[popInfo=="Southern_African",]
SOAgeno <- head(SOAgenoAll,4)
SOARes <- LDdecay(SOAgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)
SOABin <- makeBin(SOARes,max=15,n=1000)

#### West_African n=5
WEAgenoAll <- pl$geno[popInfo=="West_African",]
WEAgeno <- head(WEAgenoAll,4)
WEARes <- LDdecay(WEAgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)
WEABin <- makeBin(WEARes,max=15,n=1000)

dist <- WEARes$r2[,"pos2"] - WEARes$r2[,"pos1"]
maxDist <- tapply(dist,WEARes$r2[,"pos1"],max)
target<- 10.0
maxDist[maxDist>target] <- target
pdf("Depth5000_West_African_LDdecay_hist.pdf")
hist(maxDist,br=100,col="red",
xlab="Distance (truncated) in Mb", main="West_African")
dev.off()

#### PLOT ALL TOGETHER ####
# colors
Nubian 184, 106, 10 #b86a0a
Kordofan 225, 146, 20 #e19214
Masai 118, 29, 70 #761d46
Southern_African_Central 46, 61, 138 #2e3d8a
Angolan 23, 30, 62 #171e3e
Masai_Selous 158, 37, 108 #9e256c
Reticulated 190, 26, 14 #be1a0e
Masai_Thornicrofts 174, 108, 168 #ae6ca8
Southern_African 61, 123, 191 #3d7bbf
West_African 250, 205, 56 #facd38

# order of labels
# West- Kordofan - Nubian - Reticulated - Masai - Masai Selous - Masai Thornicroft - 
# Southern - Southern Central - Angolan

pdf("Depth5000_LDdecay_n4.pdf")
plot(NUBBin$seq,NUBBin$r2bin,type="l",lwd=3,col='#b86a0a',ylab='mean r2',
	xlab='distance (Mb)',ylim=c(0,1.0), main='Depth5000_LDdecay_n4')
 	lines(KORBin$seq,KORBin$r2bin,lwd=3,col='#e19214')
 	lines(MASBin$seq,MASBin$r2bin,lwd=3,col='#761d46')
	lines(SACBin$seq,SACBin$r2bin,lwd=3,col='#2e3d8a')
	lines(ANGBin$seq,ANGBin$r2bin,lwd=3,col='#171e3e')
 	lines(MSEBin$seq,MSEBin$r2bin,lwd=3,col='#9e256c')
	lines(RECBin$seq,RECBin$r2bin,lwd=3,col='#be1a0e')
 	lines(MTHBin$seq,MTHBin$r2bin,lwd=3,col='#ae6ca8')
	lines(SOABin$seq,SOABin$r2bin,lwd=3,col='#3d7bbf')
 	lines(WEABin$seq,WEABin$r2bin,lwd=3,col='#facd38')
 legend("topright",fill=c('#facd38','#e19214','#b86a0a','#be1a0e','#761d46',
 	'#9e256c','#ae6ca8','#3d7bbf','#2e3d8a','#171e3e'),
 	c("West African","Kordofan","Nubian","Reticulated","Masai","Masai Selous",
 	"Masai Thornicroft","Southern African","Southern African Central","Angolan"))
dev.off()

# n=8
#### Nubian n=8
NUBgenoAll <- pl$geno[popInfo=="Nubian",]
NUBgeno <- head(NUBgenoAll,8)
NUBRes <- LDdecay(NUBgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)
NUBBin <- makeBin(NUBRes,max=15,n=1000)

#### Southern African Central n=8
SACgenoAll <- pl$geno[popInfo=="Southern_African_Central",]
SACgeno <- head(SACgenoAll,8)
SACRes <- LDdecay(SACgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)
SACBin <- makeBin(SACRes,max=15,n=1000)

#### Angolan n=12
ANGgenoAll <- pl$geno[popInfo=="Angolan",]
ANGgeno <- head(ANGgenoAll,8)
ANGRes <- LDdecay(ANGgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)
ANGBin <- makeBin(ANGRes,max=15,n=1000)

#### Reticulated n=12
RECgenoAll <- pl$geno[popInfo=="Reticulated",]
RECgeno <- head(RECgenoAll,8)
RECRes <- LDdecay(RECgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)
RECBin <- makeBin(RECRes,max=15,n=1000)

#### Southern_African n=12
SOAgenoAll <- pl$geno[popInfo=="Southern_African",]
SOAgeno <- head(SOAgenoAll,8)
SOARes <- LDdecay(SOAgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)
SOABin <- makeBin(SOARes,max=15,n=1000)

pdf("Depth5000_LDdecay_n8.pdf")
plot(NUBBin$seq,NUBBin$r2bin,type="l",lwd=3,col='#b86a0a',ylab='mean r2',
	xlab='distance (Mb)',ylim=c(0,1.0), main='Depth5000_LDdecay_n8')
	lines(SACBin$seq,SACBin$r2bin,lwd=3,col='#2e3d8a')
	lines(ANGBin$seq,ANGBin$r2bin,lwd=3,col='#171e3e')
	lines(RECBin$seq,RECBin$r2bin,lwd=3,col='#be1a0e')
	lines(SOABin$seq,SOABin$r2bin,lwd=3,col='#3d7bbf')
 legend("topright",fill=c('#b86a0a','#be1a0e','#3d7bbf','#2e3d8a','#171e3e'),
 	c("Nubian","Reticulated","Southern African","Southern African Central","Angolan"))
dev.off()


# n=10
#### Angolan n=12
#ANGgenoAll <- pl$geno[popInfo=="Angolan",]
ANGgeno <- head(ANGgenoAll,10)
ANGRes <- LDdecay(ANGgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)
ANGBin <- makeBin(ANGRes,max=15,n=1000)

#### Reticulated n=12
#RECgenoAll <- pl$geno[popInfo=="Reticulated",]
RECgeno <- head(RECgenoAll,10)
RECRes <- LDdecay(RECgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)
RECBin <- makeBin(RECRes,max=15,n=1000)

#### Southern_African n=12
#SOAgenoAll <- pl$geno[popInfo=="Southern_African",]
SOAgeno <- head(SOAgenoAll,10)
SOARes <- LDdecay(SOAgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)
SOABin <- makeBin(SOARes,max=15,n=1000)

pdf("Depth5000_LDdecay_n10.pdf")
plot(ANGBin$seq,ANGBin$r2bin,type="l",lwd=3,col='#171e3e',ylab='mean r2',
	xlab='distance (Mb)',ylim=c(0,1.0), main='Depth5000_LDdecay_n10')
	lines(RECBin$seq,RECBin$r2bin,lwd=3,col='#be1a0e')
	lines(SOABin$seq,SOABin$r2bin,lwd=3,col='#3d7bbf')
 legend("topright",fill=c('#be1a0e','#3d7bbf','#171e3e'),
 	c("Reticulated","Southern African","Angolan"))
dev.off()

# only angolan. masai selous and reticulated

#### Angolan n=12
ANGgenoAll <- pl$geno[popInfo=="Angolan",]
ANGgeno <- head(ANGgenoAll,5)
ANGRes <- LDdecay(ANGgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)
ANGBin <- makeBin(ANGRes,max=15,n=1000)

#### Masai_Selous n=5
MSEgenoAll <- pl$geno[popInfo=="Masai_Selous",]
MSEgeno <- head(MSEgenoAll,5)
MSERes <- LDdecay(MSEgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)
MSEBin <- makeBin(MSERes,max=15,n=1000)

#### Reticulated n=12
RECgenoAll <- pl$geno[popInfo=="Reticulated",]
RECgeno <- head(RECgenoAll,5)
RECRes <- LDdecay(RECgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)
RECBin <- makeBin(RECRes,max=15,n=1000)


#### PLOT ALL TOGETHER ####
# colors
Nubian 184, 106, 10 #b86a0a
Kordofan 225, 146, 20 #e19214
Masai 118, 29, 70 #761d46
Southern_African_Central 46, 61, 138 #2e3d8a
Angolan 23, 30, 62 #171e3e
Masai_Selous 158, 37, 108 #9e256c
Reticulated 190, 26, 14 #be1a0e
Masai_Thornicrofts 174, 108, 168 #ae6ca8
Southern_African 61, 123, 191 #3d7bbf
West_African 250, 205, 56 #facd38

# order of labels
# West- Kordofan - Nubian - Reticulated - Masai - Masai Selous - Masai Thornicroft - 
# Southern - Southern Central - Angolan

pdf("Depth5000_LDdecay_n5_ANG_MSE_REC.pdf")
plot(ANGBin$seq,ANGBin$r2bin,type="l",lwd=3,col='#b86a0a',ylab='mean r2',
	xlab='distance (Mb)',ylim=c(0,1.0), main='Depth5000_LDdecay_n5')
 #	lines(KORBin$seq,KORBin$r2bin,lwd=3,col='#e19214')
# 	lines(MASBin$seq,MASBin$r2bin,lwd=3,col='#761d46')
#	lines(SACBin$seq,SACBin$r2bin,lwd=3,col='#2e3d8a')
	lines(ANGBin$seq,ANGBin$r2bin,lwd=3,col='#171e3e')
 	lines(MSEBin$seq,MSEBin$r2bin,lwd=3,col='#9e256c')
	lines(RECBin$seq,RECBin$r2bin,lwd=3,col='#be1a0e')
 #	lines(MTHBin$seq,MTHBin$r2bin,lwd=3,col='#ae6ca8')
#	lines(SOABin$seq,SOABin$r2bin,lwd=3,col='#3d7bbf')
 #	lines(WEABin$seq,WEABin$r2bin,lwd=3,col='#facd38')
 legend("topright",fill=c('#be1a0e',
 	'#9e256c','#171e3e'),
 	c("Reticulated","Masai Selous",
"Angolan"))
dev.off()

# END
