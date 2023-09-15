library(tidyverse)
library(ggplot2)
library(data.table)
library(ggthemes)
library("RColorBrewer")

# bim_df
# chr, var, morgan, coord/pos (1-based), A1 (minor), A2(major)
#
# return
# data.frame(chr, chrLen)
get_chr_length <- function(bim_df){
        if (! all(c("chr", "pos") %in% colnames(bim_df))){
                stop("bim dataframe doesn't have required columns \"chr\" or \"pos\"")
                return(NULL)
        }
        as.data.frame(bim_df %>% group_by(chr) %>% summarise(len=max(pos)-min(pos)+1))
}

# chrLen_df
# chr, len
#
# return
# <integer>
est_geno_len <- function(chrLen_df){        
        if (! all(c("chr", "len") %in% colnames(chrLen_df))){
                stop("chromosome length dataframe doesn't have required columns \"chr\" or \"len\"")
                return(NULL)
        }
        sum(chrLen_df$len)
}

# input
# roh_len_list: from plink .hom file, column "KB"
# ranges_upper: length range upper inclusive boundary.
# ranges_name: name of each range after grouped
#
# return: 
# data.frame() with different length range per column
roh_len_genome_prop <-function(roh_lens, geno_len=1000, ranges_upper=c(1000, 2000, 5000, 10000), ranges_name=c("lt1m", "rg1m2m", "rg2m5m", "rg5m10m", "gt10m")){
        if (! (length(ranges_upper)+1)==length(ranges_name)){
                stop("number input ranges uppers doesn't match number of input ranges names. len(upper)+1==len(name)")
                return(NULL)
        }
	n=geno_len/1000 # unit: KB
        lt1m=sum(roh_lens[roh_lens<=1000])/n # <=1M
	rg1m2m=sum(roh_lens[roh_lens>1000 & roh_lens<=2000])/n # (1M, 2M]
        rg2m5m=sum(roh_lens[roh_lens>2000 & roh_lens<=5000])/n # (2M, 5M]
        rg5m10m=sum(roh_lens[roh_lens>5000 & roh_lens<=10000])/n # (5M, 10M]
        gt10m=sum(roh_lens[roh_lens > 10000])/n # (10M, +inf)
        return(setNames(data.frame(matrix(c(lt1m, rg1m2m, rg2m5m, rg5m10m, gt10m), nrow=1)), ranges_name))
}

args <- commandArgs(trailingOnly=TRUE)

homfile <- args[1] # "all.hom"
bimfile <- args[2] # "filtered.maf5.geno5.bim"

geno_len <- est_geno_len(get_chr_length(fread(bimfile, header=FALSE, col.names=c("chr", "var", "morgan", "pos", "A1", "A2"))))

hom_df <- fread(homfile, header=TRUE)

roh_prop <- as.data.frame(hom_df %>% group_by(IID) %>% summarize(roh_len_genome_prop(KB, geno_len = geno_len)))


# normal plot
rownames(roh_prop) <- roh_prop$IID
roh_prop$IID <- NULL
plotdata <- t(roh_prop)
#pdf("ROH_barplot.pdf", width=18, height=10)
png("ROH_barplot.png", width=1800, height=800)
par(mar=c(24, 12, 8, 10)+0.1, mgp=c(5,1,0))
#**!!!!change space to vector to have a "grouped" view!!!**
barplot(plotdata, col=brewer.pal(n = nrow(plotdata), name = "Blues"), border="blue", width=4, space=c(4, 2), las=2, ylab="Genome proportion", cex.axis=1.5, cex.names=2.2, cex.lab=3, ylim=c(0, 0.3), main="ROH Proportion on Genome", cex.main=3)
#abline(h=seq(0.05, 0.25, 0.05), lwd=0.5, lty=2)
legend('topright', legend = c("<=1 Mbp", "1~2 Mbp", "2~5 Mbp", "5~10 Mbp", ">10 Mbp") ,bty='o', pch = 22, col="blue", pt.cex=3, pt.bg=brewer.pal(n = nrow(plotdata), name = "Blues"), xpd = T, x.intersp=0.8, cex=2)
dev.off()


stop()

# ggplot version

plotdata <- pivot_longer(
  roh_prop,
  cols = c("lt1m","rg1m2m","rg2m5m","rg5m10m","gt10m"),
  names_to = "range",
  values_to = "geno_prop",
  values_transform = list(geno_prop = as.numeric)
)
plotdata$group <- sub("_filt", "", plotdata$IID)
plotdata$range <- factor(plotdata$range, levels=c("lt1m","rg1m2m","rg2m5m","rg5m10m","gt10m"))

ggplot(plotdata, aes(fill= range, y = geno_prop, x= IID, group=group)) +
  geom_bar(position="stack", stat="identity") +
  scale_y_continuous(labels = scales::percent_format(acurracy=1)) +
  scale_fill_manual(values = brewer.pal(n = length(unique(plotdata$range)), name = "Blues")) +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(), 
    axis.line.y = element_line(),
    
    ) + labs(y="Genome proportion (%)", x="", title="ROH Proportion on Genome")
  ylab("Genom") +
  xlab("")

# horizontal plot
par(mar=c(8,12,4,6))
barplot(t(roh_prop[, 2:dim(roh_prop)[2]]), col=brewer.pal(n = dim(roh_prop)[2]-1, name = "Blues"), border="blue", space=0.5, horiz=TRUE, las=1, xlab="Genom proportion (%)", cex.axis=1.5, cex.names=1.5, cex.lab=1.5, xlim=c(0, 0.3), main="ROH Proportion on Genome", cex.main=4)
abline(v=seq(0.05, 0.25, 0.05), lwd=0.5, lty=2)
