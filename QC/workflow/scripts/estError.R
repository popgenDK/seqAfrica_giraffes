## copy from ANGSD/angsd repository
bases<-c("A","C","G","T")
b <- c(bases,"N")
########### do not change ################3
l<-commandArgs(TRUE)
getArgs<-function(x,l)
  unlist(strsplit(grep(paste("^",x,"=",sep=""),l,val=T),"="))[2]
Args<-function(l,args){
 if(! all(sapply(strsplit(l,"="),function(x)x[1])%in%names(args))){
   cat("Error -> ",l[!sapply(strsplit(l,"="),function(x)x[1])%in%names(args)]," is not a valid argument")
   q("no")
 }
 arguments<-list()
 for(a in names(args))
   arguments[[a]]<-getArgs(a,l)
 
 if(any(!names(args)%in%names(arguments)&sapply(args,is.null))){
   cat("Error -> ",names(args)[!names(args)%in%names(arguments)&sapply(args,is.null)]," is not optional!\n")
   q("no")
 }
 for(a in names(args))
   if(is.null(arguments[[a]]))
     arguments[[a]]<-args[[match(a,names(args))]]
 
 
 arguments
}

print.args<-function(args,des){
  if(missing(des)){
    des<-as.list(rep("",length(args)))
    names(des)<-names(args)
  }
  cat("->  Needed arguments:\n")
  mapply(function(x)cat("\t",x,":",des[[x]],"\n"),cbind(names(args)[sapply(args,is.null)]))
  cat("->  Optional arguments (defaults):\n")
  mapply(function(x)cat("\t",x," (",args[[x]],")",":",des[[x]],"\n"),cbind(names(args)[!sapply(args,is.null)]))
  q("no")
}

## choose your parameters and defaults
## NULL is an non-optional argument, NA is an optional argument with no default, others are the default arguments
args<-list(file=NULL,
           out="errorEst",
           indNames="ind",
           nIter=100,
           subset=NA,
           main="Error rate using an outgroup and a high quality genome",
           height=7,
           width=12,
           srt=90,
           cex=1,
           doPng=FALSE,
           ylim1=NA,
           ylim2=NA
           )
#if no argument aree given prints the need arguments and the optional ones with default
des<-list(file="the ancError File",
          out="Name of the out files",
          indNames="postFix, file with names or comma seperated names of individuals",
          nIter="Numer of optimazation attemps",
          subset="comma seperated numbers of the individuals to include",
          width="width of the pdf",
          height="height of the pdf",
          srt="angle of ind names",
          cex="scale of names",
          doPng="Make png instead of pdf",
          ylim1="Set max ylim for barplot1",
          ylim2="Set max ylim for barplot2"
          )

#######get arguments and add to workspace
### do not change
if(length(l)==0) print.args(args,des)
attach(Args(l,args))
args <- commandArgs(TRUE)
if(length(args)==0){
  cat(" Arguments: output prefix\n")
  q("no")
}

###################################
tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")



width=as.numeric(width)
height=as.numeric(height)

if(!is.na(ylim1))
    ylim1=as.numeric(ylim1)
if(!is.na(ylim2))
    ylim2=as.numeric(ylim2)


palette(tol21rainbow[c(1:8,10:21)])
srt<-as.numeric(srt)


nIter<-as.integer(nIter)
cex=as.numeric(cex)
cat("----------\nfile: ",file," out: ",out," nIter: ",nIter," subset: ",subset,"\n width: ",width," height: ",height," srt: ",srt," cex: ",cex," doPng: ",doPng,"\n-----------\n")

r<-as.matrix(read.table(file))

if(!is.na(subset)){
  subset<- as.integer(unlist(strsplit(subset,",")))
  print(subset)
  r<-r[subset,]
}

nInd<-nrow(r)
cat("Number of individuals read:",nInd,"\n")
{
  if(length(grep(",",indNames))>0)
    indNames<-unlist(strsplit(indNames,","))
  else  {
    options("warn"=-1)
    try(indNames<-basename(scan(indNames,what="theFuck")),silent=TRUE)
    options("warn"=0)
  }
}

if(length(indNames)==1&nInd>1){
  indNames<-paste(indNames,1:nInd,sep="")
}
cat("Ind names:\n")
print(indNames)

if(length(indNames)!=nInd){
  cat("Error: Wrong number of ind Names\n")
  q("no")
}

getMat<-function(x){
  m<-array(0,dim=c(5,5,5),dimnames=list(b,b,b))
  for(s in 0:4)
    for(p in 0:4)
      for(a in 0:4)
        m[a+1,p+1,s+1]<-x[a*25+p*5+s+1]
  m
}

logLike<-function(x,Xch,Pch){
  eMat<-matrix(0,4,4)
  eMat[-c(1,6,11,16)]<-x
  diag(eMat)<-1-rowSums(eMat)
  P <- Pch %*% eMat
  ll <- -sum(log(P)*Xch)
  return(ll)
}

res<-NULL
for(j in 1:nInd){
  m<-getMat(r[j,])
  print(m)

  ##remove if missing
  m<-m[-5,-5,-5]
  
  Pch<-matrix(0,4,4)
  for(i in 1:4)
    Pch<-Pch+m[,,i]
  print(Pch)
  Pch<-Pch/rowSums(Pch)
  cat("Matrix of allele change probabilities.\nA diagonal close to 1 denotes a low probability of errors.\n\n")
  print(Pch)
  
  Xch<-matrix(0,4,4)
  for(i in 1:4)
      Xch<-Xch+m[,i,]
  cat("\n\nObserved allele changes:\n\n")
  print(Xch)
  maxErr=0.02;

  conv<-list(par=1,objective=Inf)
      
  for(i in 1:nIter){
      Tempconv <- nlminb(runif(12)/100,logLike,upper=rep(maxErr,12),lower=rep(1e-10,12),Xch=Xch,Pch=Pch)
      if( any(Tempconv$par > maxErr*0.99) )
          maxErr<-min(1,maxErr*10);
      if(Tempconv$objective<conv$objective){
          conv<-Tempconv
      
      }
  }
  
  res<-rbind(res,conv$par)
}

getover<-function(r,nInd){
  over<-NULL
  for(j in 1:nInd){
    m<-getMat(r[j,])
    
    ##remove if missing
    m<-m[-5,-5,-5]
    Xch<-matrix(0,4,4)
    for(i in 1:4)
      Xch<-Xch+m[,i,]
    
    Pch<-matrix(0,4,4)
    for(i in 1:4)
      Pch<-Pch+m[,,i]
    
    N1dot<-sum(Xch)-sum(diag(Xch))
    Ndot1<-sum(Pch)-sum(diag(Pch))
    Ndot0<-sum(diag(Pch))
    err<-(N1dot-Ndot1)/(Ndot0-Ndot1)
    over<-c(over,err)
  }
  over
}
over<-getover(r,nInd)
nam<-paste(rep(bases,4),"->",rep(bases,each=4))[-c(1,6,11,16)]

{
  if(doPng==FALSE){
    pdf(paste(out,".pdf",sep=""),width=width,height=height)
  }  else
    bitmap(paste(out,".png",sep=""),w=width,h=height,res=300)
}

{
if(is.na(ylim1))
    barplot(res,beside=T,col=1:nInd,names=nam,main=main,ylab="error rate")
else
    barplot(res,beside=T,col=1:nInd,names=nam,main=main,ylab="error rate",ylim=c(0,ylim1))
}
legend("top",paste(indNames,round(over*100,2),"%"),fill=1:nInd,bty="n")
dev.off()

colnames(res)<-nam
rownames(res)<-indNames

write.table(res,file=paste(out,".txt",sep=""),sep="\t")
{
if(doPng==FALSE)
    pdf(paste(out,"Overall.pdf",sep=""),w=width,h=height)
else
    bitmap(paste(out,"Overall.png",sep=""),w=width,h=height,res=300)
}

par(mar=par()$mar+c(3,0,0,0))
{
if(is.na(ylim2)){
    h<-barplot(over,col=1:nInd,names=NULL,main=main,ylab="error rate")
}else{
    h<-barplot(over,col=1:nInd,names=NULL,main=main,ylab="error rate",ylim=c(0,ylim2))
}
}
text(h,rep(-max(over)/50,length(h)),indNames,xpd=T,srt=srt,adj=1,cex=cex)
dev.off()
{
    if(doPng==FALSE){
        cat("figure:",paste(out,".pdf",sep=""),"\n")
        cat("figure2:",paste(out,"Overall.pdf",sep=""),"\n")
    }else{
        cat("figure:",paste(out,".png",sep=""),"\n")
        cat("figure2:",paste(out,"Overall.png",sep=""),"\n")
    }

}
cat("table:",paste(out,".txt",sep=""),"\n")

write.table(cbind(paste(indNames,round(over*100,4),"%")),file=paste(out,".txt",sep=""),sep="\t",append=T,col=F,row=F)

chrFile<-paste(file,"Chr",sep="")

if(file.info(chrFile)$size!=0 & !is.na(file.info(chrFile)$size)){

    overList<-list()
    con<-file(chrFile,"r")
    while(length(cchr<-scan(con,nlines=1,what="theFuck"))>0){
        r<-read.table(con,nrow=nInd)
        r<-as.matrix(r)
        overList[[cchr[2]]]<-getover(r,nInd)       
    }
    close(con)
    
    rbind.list<-function (l) 
    {
        n <- length(l)
        results <- NULL
        for (i in 1:n) {
            results <- rbind(results, l[[i]])
        }
        results
    }
    res<-rbind.list(overList)
    chrNames<-names(overList)
    rownames(res) <- chrNames
    colnames(res) <- indNames
    pdf(paste(out,"OverChr.pdf",sep=""),w=width,h=height)
    for(i in 1:nInd){
        dotchart(res[,i],chrNames,xlab="Error rates",col=3,main=indNames[i])
        # dotchart(res[res[,i]>0,i],chrNames,xlab="Error rates (positive)",col=3,main=indNames[i])
    }
    dev.off()
    write.table(res,file=paste(out,"OverChr.txt",sep=""),qu=F,row=T,col=T)
}

#warnings()
