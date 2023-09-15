#!/usr/bin/env Rscript
####### ROH version 0.1.1 ####### 
# Fill adjacent ROH region with different color.
################################
####### ROH version 0.1 ####### 
# This code was made for visualzing ROH
# It requries plink file as input #
################################
# Author: Xiaodong Liu (UCPH), Shixu He (UCPH)

#### Load libraries ####
source("http://popgen.dk/albrecht/open/online.R")
library(snpStats)
library(windowscanr)
library(tidyverse)
library("RColorBrewer")
suppressPackageStartupMessages(library("argparse"))
#### end of loading

#### functions ####
extract_xcoord <- function(x,y,v1,lower.bound, upper.bound) {
  output=list()
  xcoords = c()
  ycoords =c()
  ystart = y[1]
  dec.factor = lower.bound
  incr.factor = upper.bound
  for (i in v1) {
    if (i==1 ) {
      xcoords = append(xcoords,x[i])
      ycoords = append(ycoords,ystart-dec.factor)
    } else if (i==v1[length(v1)]) {
      xcoords = append(xcoords,x[i])
      ycoords = append(ycoords,ystart+incr.factor )
    } else {
      xend = x[i]
      j= i+1
      xstart = x[j]
      #cat(xstart)
      xcoords = append(xcoords,c(xend,xstart))
      ycoords = append(ycoords,c(ystart+incr.factor ,ystart-dec.factor))
    }
  }
  output[[1]] = xcoords
  output[[2]] = ycoords
  return(output)
}

prop_het <- function(x) {
  return(length(which(x==1)) / length(x[!is.na(x)]))
}

snp_count <- function(x) {
  return( length(x[!is.na(x)]) )
}

het_count <- function(x) {
  return( length(which(x==1)) )
}

hom_count <- function(x) {
  return( length(which(x==0|x==2) ) )
}
#### end of functions

#### parse argument ####
# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-p", "--plinkfile", type="character",
                    help="plink bed file",metavar="string")
parser$add_argument("-s", "--sample", type="character", 
                    help="sample to plot",
                    metavar="string")
parser$add_argument("-w","--windowsize",type="integer",default=100000,
                    help ="sliding window size for density calculation [100000 by default]",
                    metavar="interger")
parser$add_argument("--homfile", type="character",
                    metavar="string",
                    help = "hom output by plink")


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

# print some progress messages to stderr if "quietly" wasn't requested
if ( args$verbose ) { 
  write("writing some verbose information while processing..\n", stderr()) 
}

filename = args$plinkfile
individual = args$sample
homfile = args$homfile
windowsize = args$windowsize
#### end of argunemt parser

#### read plink file ####

plinkfile = gsub(".bed", "", filename )
dat = plink(plinkFile =plinkfile)

write("finish plink file reading\n", stderr())

if( args$verbose) {
  cat(plinkfile, "was read into R.\n")
}


i = which(dat$fam[,2] ==individual)
print(i)
genotype <-dat$geno[i,]

t2 <-cbind(dat$bim[,c(1,4)], genotype)
names(t2) <- c("chr","pos","gt")
x=dat$bim$V4
y0=as.factor((dat$bim$V1))
y1 = factor(y0,levels=unique(y0))
levels( y1 ) <- 1:length(unique(y1))
space.chr = 5
y=as.numeric(y1) *space.chr

t2$chr <- as.numeric(y1)
#print(head(t2))
#### end of reading plink file


#### calculate SNP/ hetero density ####
output <- winScan(t2,groups="chr",position="pos",values="gt",win_size=100000,funs=c("snp_count","het_count","hom_count"),cores =5)
output$gt_snp_count = output$gt_het_count + output$gt_hom_count
output$het_ratio = output$gt_het_count/output$gt_snp_count
if( args$verbose) {
  cat("Density calculation was done.\n")
}

#### end of calculating density

#### plot ####
if( args$verbose) {
  cat("start plotting.\n")
}

hom_gts=which(t2$gt==2 | t2$gt==0)
het_gts=which(t2$gt==1)
xhets<-x[het_gts]
yhets<-y[het_gts]
xhoms<-x[hom_gts]
yhoms<-y[hom_gts]

rect.height = max(yhoms)/2*0.02
upper.bound  = ifelse(rect.height > 0.5, 0.5,  rect.height)
lower.bound  = ifelse(rect.height > 0.5, 0.5,  rect.height)

## plot rectangles for 
triangle.d = data.frame(x1 =numeric(),
                        x2 = numeric(),
                        y1 = numeric(),
                        y2 = numeric())

for (chr in unique(yhoms)) {
  index = which(yhoms == chr)
  yhom= yhoms[index]
  xhom1 = xhoms[index]
  
  xhom2 = c(xhom1[2:length(xhom1)],NA)
  dist = xhom2 - xhom1
  gap.index =which(dist>=500000)
  rect.index = c(1,gap.index,length(xhom1))
  ( xy.coords = extract_xcoord(xhom1,yhom,rect.index,lower.bound,upper.bound) )
  
  for (rect in seq(1,length(xy.coords[[1]]),2) ) {
    #rect(xy.coords[[1]][rect],xy.coords[[2]][rect],xy.coords[[1]][rect+1],xy.coords[[2]][rect+1],col="blue",border = NA)
    triangle.d =  triangle.d %>% add_row(x1= as.numeric(xy.coords[[1]][rect]), x2 =as.numeric(xy.coords[[1]][rect+1]),y1 = xy.coords[[2]][rect], y2 = xy.coords[[2]][rect+1])
    # triangle.d[nrow(triangle.d+1),] = c(xy.coords[[1]][rect], xy.coords[[1]][rect+1], xy.coords[[2]][rect], xy.coords[[2]][rect+1])
  }
  
}

g1 = ggplot() + 
  scale_x_continuous(name="x") + 
  scale_y_continuous(name="y") +
  geom_rect(data=triangle.d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="blue", fill=alpha("blue", 1))

segment.d = data.frame(x = xhets, y =yhets-lower.bound, xend = xhets,yend = yhets+upper.bound )
g2 = g1 + geom_segment(data=segment.d, mapping= aes(x=x, y =y, xend =xend,yend=yend),col="red", linewidth = 0.03)

## plot density
start = 1.5
( cap.density = quantile(output$gt_snp_count,0.9,na.rm = T) )
#hist(output$het_ratio)
( cap.proportion = quantile(output$het_ratio,0.5,na.rm = T) )

g3 = g2
for (i in 1:length(unique(output$chr))) {
  output.chr <- output[output$chr==i,]
  #  lines(output.chr$win_mid,(output.chr$gt_snp_count)/110+start+6*(i-1) ) ;
  
  density = ifelse( (output.chr$gt_snp_count/cap.density) >1, 1, output.chr$gt_snp_count/cap.density)
  #lines(output.chr$win_mid,density+start+space.chr*(i-1),col="purple"   )
  baseline1= start+space.chr*(i-1)
  density.d = data.frame(x = output.chr$win_mid, func1 = rep(baseline1, length(output.chr$win_mid)),
                         func2= density+baseline1)
  #density.d = rbind(density.d, density.d.chr)
  g3 = g3 + geom_ribbon(data = density.d,mapping = aes(x=x,ymax = func2, ymin = func1), fill = "purple", alpha = 1)
  ###polygon( c(0,output.chr$win_mid,max(output.chr$win_end)), c(baseline1, density+baseline1,baseline1),col="purple",border=NA)
  
  
  
  baseline2 = start+1.25+space.chr*(i-1)
  proportion = ifelse( output.chr$het_ratio>cap.proportion, 1, output.chr$het_ratio/cap.proportion)
  proportion.d = data.frame(x = output.chr$win_mid,func1= rep(baseline2,length(output.chr$win_mid)),
                            func2 = proportion+baseline2)
  g3 = g3 + geom_ribbon(data = proportion.d,mapping = aes(x=x,ymax = func2, ymin = func1), fill = "grey55", alpha = 1)
  #lines(output.chr$win_mid,output.chr$het_ratio + start+1.25+space.chr*(i-1),col="black"   )
  ###lines(output.chr$win_mid, baseline2+proportion,col="green")
  ###lines(output.chr$win_mid,rep(baseline2,length(output.chr$win_mid)), col="green")
  ###polygon(c(output.chr$win_mid, rev(output.chr$win_mid)), c(baseline2+proportion, rev(rep(baseline2,length(output.chr$win_mid)) )), border=NA, col="grey55")
  #polygon(c(0,output.chr$win_mid,max(output.chr$win_end)),c(baseline2, baseline2+proportion,baseline2),col="grey55" ,border=NA)
}

### deal with inferred ROH by plink
hom.table <- read.table(homfile, header=T)
hom.table.ind <- hom.table[hom.table$FID==individual & hom.table$KB>=1000,]
hom.table.ind$CHR=as.character(hom.table.ind$CHR)

xref= unique(y) # reference for change chr name to chr number
names(xref) =as.character(unique(y0))
hom.table.ind$clean <- ifelse(is.na(xref[hom.table.ind$CHR]), hom.table.ind$CHR, xref[hom.table.ind$CHR])
hom.table.ind$clean = hom.table.ind$clean + upper.bound * 1.8

#g4 = g3 +  geom_segment(data=hom.table.ind, mapping= aes(x=POS1, y =clean, xend =POS2,yend=clean),col="black", linewidth = 0.6)

colist <- rep(c("black", "darkgray"),floor(length(hom.table.ind$POS1)/2) )

if (length(hom.table.ind$POS1) %% 2){
    colist <- c(colist, "black")
}

g4 = g3 +  geom_segment(data=hom.table.ind, mapping= aes(x=POS1, y =clean, xend =POS2,yend=clean),col=colist, linewidth = 0.6)
g4

### print the figure in a file
myblanktheme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "white"),
                     axis.title=element_text(size=14,face="bold"))


bitmap(paste0(individual,".ROH.Density.png"),h=14,w=20,res=300)
g4 + myblanktheme +  scale_y_continuous(name = NULL, # remove y title
                                        breaks = unique(yhoms), 
                                        labels = paste0("chr",unique(t2$chr))) +
  labs(y="", x="Positions",title=individual)
dev.off()
