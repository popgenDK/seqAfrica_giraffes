# /usr/bin/env R
source("LDdecay_functions.R")

#Load libraries
library(Relate)
library(snpStats)
if(FALSE){
    library(devtools)
    install_github("aalbrechtsen/relate")
}

#pl <- plink

#Plotting colours
cols <- c("Ghana"="#4E052D","Togo"="#7d064a","Nigeria"="#cb1a7e","Cameroon"="#F359A1","Eq Guinea"="#F89AA6","Gabon"="#FCB85F","DR Congo"="#D0E25A","Uganda"="#09C189","Ethiopia"="#0F8964","Tanzania"="#5DD9E4","Zimbabwe"="#40A0D8","South Africa"="#1347A2","Madagascar"="#0a1c6c")

#1.Cameroon
CamgenoAll <- pl$geno[popInfo=="Cameroon",]
Camgeno <- head(CamgenoAll,5)
CamRes <- LDdecay(Camgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=2000)

#Calculate LD and plot physical distance in Mb (Max 5 Mb and 75 bins)
CamBin <- makeBin(CamRes,max=5,n=75)

#2. Eq_Guinea
EqGgenoAll <- pl$geno[popInfo=="Eq_Guinea",]
EqGgeno <- head(EqGgenoAll,8)
EqGRes <- LDdecay(EqGgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)

#Calculate LD and plot physical distance in Mb (Max 5 Mb and 75 bins)
EqGBin <- makeBin(EqGRes,max=5,n=75)

#3. Madagascar
MadgenoAll <- pl$geno[popInfo=="Madagascar",]
Madgeno <- head(MadgenoAll,32)
MadRes <- LDdecay(Madgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)

#Calculate LD and plot physical distance in Mb (Max 5 Mb and 75 bins)
MadBin <- makeBin(MadRes,max=5,n=75)

#4. Tanzania
TzgenoAll <- pl$geno[popInfo=="Tanzania",]
Tzgeno <- head(TzgenoAll,5)
TzRes <- LDdecay(Tzgeno,pos=pl$bim[,4]/1e6,chr=pl$bim[,1],maf=0.05,mis=0.01,depth=5000)

#Calculate LD and plot physical distance in Mb (Max 5 Mb and 75 bins)
TzBin <- makeBin(TzRes,max=5,n=75)

#First bin
CamBin$r2bin[1] <- 0.5
EqGBin$r2bin[1] <- 0.5
MadBin$r2bin[1] <- 0.5
TzBin$r2bin[1] <- 0.5

#Plot
pdf("LDdecay_n5_5mb_5000depth.pdf", height=6, width=8)
par(mar = c(5, 5, 5, 5))
plot(CamBin$seq,CamBin$r2bin,type="l",lwd=5,col="#F359A1",ylab=expression("Mean r"^2),
    xlab='Distance (Mb)',ylim=c(0.17,0.4),xlim=c(0.01,5),
    cex.axis = 1.3, cex.lab = 1.3)
    lines(EqGBin$seq,EqGBin$r2bin,lwd=5,col="#F89AA6")
    lines(MadBin$seq,MadBin$r2bin,lwd=5,col="#0a1c6c")
    lines(TzBin$seq,TzBin$r2bin,lwd=5,col="#5DD9E4")
 legend("topright",fill=c("#F359A1","#F89AA6","#0a1c6c","#5DD9E4"),
     c("Cameroon (n=5)","Eq Guinea (n=5)","Madagascar (n=5)","Tanzania (n=5)"))
dev.off()
