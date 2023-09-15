require(reshape2)
require(plotly)
require(admixtools)
whereDir <- function(){
    # function to get directory where scripts are, so accompanying functions can be sourced even if the script is run from outside
    # made by krishang
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    needle <- "--file"
    match <- grep(needle, cmdArgs)
    tf <- unlist(strsplit(cmdArgs[match], "="))[2]
    d <- dirname(tf)
    return(d)
}

fillmat <- function(df){
    newdf <- df
    colchange <- c(1,2,5,6,7,8,9,10,11,12)
    newcol <- c(2,1,6,5,8,7,10,9,12,11)
    temp <- newdf[,newcol]
    newdf[,colchange] <- temp
    newdf[,c(4,11,12)] <- newdf[,c(4,11,12)] * -1
    colnames(newdf) <- colnames(df)
    rbind(df, newdf)
}

d <- whereDir()
# d <- "../scripts/"
source('/home/krishang/projects/DNA/africa1kg/qpgraphs/scripts/qpgraphRfuns.R')

args <- commandArgs(trailing=TRUE)
# args <- c("allpop_tests.txt.Rdata", "allpop_tests.txt", "res.pdf")
rdataFile <- args[1]
compsFile <- args[2]
admixStart <- as.numeric(args[3])+1
admixEnd <- as.numeric(args[4])+1

load(rdataFile)

compsall_temp <- data.table::fread(compsFile, data.table=F)

print(table(keep <- !is.na(compsall_temp$p)))
compsall_temp <- compsall_temp[keep,]
compsall <- fillmat(compsall_temp)
compsall <- compsall[order(compsall$score_test1,compsall$score_test2),]

#filter out those with temporal conflicts

source('/home/albrecht/projects/admixer2022/scripts/timepolice.R')
idxs_caught_by_police = c()
for (i in 1:length(all_fits)){
        if(length( res<-timepolice(as.data.frame(all_fits[[i]]$edges)) )>0 ){
                idxs_caught_by_police<-append(idxs_caught_by_police,i)
        }
}

library(tidyverse) # cause conflicts with timepolice: Error: C stack usage  7969828 is too close to the limit
compsall=filter(compsall,!(idx1%in%idxs_caught_by_police | idx2%in%idxs_caught_by_police))


source('/home/users/xiaodong/Documents/Project/Wildebeest/findgraph/summary/NestedQpgraph.R')
mylist <- vector(mode="list", length=length(unique(compsall$adm1)))
for (nadmix in seq(0,max(unique(compsall$adm1)))){
        compsMax<-subset(compsall,adm1<=nadmix & adm2<=nadmix)
        comps <- subset(compsMax, p>0.05 & idx1==compsMax[1,1])
        #print(paste0('n: ',nadmix))
        #print(paste0('best: ', compsMax[1,1]))
        #print(head(comps,2))
        if(nrow(comps)==0){
                comps <- subset(compsall,  adm1==nadmix)
                mylist[[nadmix+1]] = list(all_fits[[comps[1,1]]]$edges)
                names(mylist[[nadmix+1]]) = paste0('GraphID',as.character(comps[1,1]))

        }
        else{
                compsidx1 <- comps[1,1] # ifelse(comps[1,"score_test1"] > comps[1,"score_test2"], comps[1,1], comps[2,2])
                IDs=c(compsidx1)
                mylist[[nadmix+1]] = list(all_fits[[compsidx1]]$edges)
                if(nrow(comps)>1){
                        for (idx in 1:nrow(comps)){
                                idx2 = comps[idx,"idx2"]
                                mylist[[nadmix+1]][[idx+1]] = all_fits[[idx2]]$edges
                                IDs=c(IDs,idx2)
                        }
                }
                names(mylist[[nadmix+1]]) = paste0('GraphID',as.character(IDs))
        }
}

mylist2=mylist[admixStart:admixEnd]

all = data.frame()

#for (OldGraphIndex in head(seq_along(mylist2),-1) ){
#	for	(newGraghIndex in tail(seq_along(mylist2),-1)){
#		if (newGraghIndex > OldGraphIndex){
#			for (newGraghIndex2 in 1:length(mylist2[[newGraghIndex]]) ){
#                		newGraghID = names(mylist2[[newGraghIndex]])[newGraghIndex2]
#				nestedResult = data.frame(nested=FALSE)
#				#print(newGraghID)
#				#print(newGraghIndex)
#				#print(newGraghIndex2)
#				for (OldGraphIndex2 in 1:length(mylist2[[OldGraphIndex]]) ){
#					oldGraphID = names(mylist2[[OldGraphIndex]])[OldGraphIndex2]
#					#print(oldGraphID)
#					#print(OldGraphIndex)
#					#print(OldGraphIndex2)
#					nestedResult = nested(as.data.frame(mylist2[[OldGraphIndex]][[OldGraphIndex2]]),as.data.frame(mylist2[[newGraghIndex]][[newGraghIndex2]]))
#					#print(nestedResult)
#					if (nestedResult$nested){
#						#bugTrace
#						nestedResult = add_column(nestedResult,admix1Graph=oldGraphID,.after=1)
#						nestedResult = add_column(nestedResult,admix2Graph=newGraghID,.after=3)
#						#nestedResult = add_column(nestedResult,admix2Graph=newGraghID,.after=2)
#						break 
#					}
#					else{
#						nestedResult = nested2(as.data.frame(mylist2[[OldGraphIndex]][[OldGraphIndex2]]),as.data.frame(mylist2[[newGraghIndex]][[newGraghIndex2]]))
#						nestedResult = add_column(nestedResult,admix1Graph=oldGraphID,.after=1)
#						nestedResult = add_column(nestedResult,admix2Graph=newGraghID,.after=3)
#						#nestedResult = add_column(nestedResult,admix2Graph=newGraghID,.after=2)
#					}
#				}
#				all = rbind(all,nestedResult)
#			}
#			
#
#		#	dd=lapply(mylist2[[newGraghIndex]],function(newGragh,oldGraph) nested(as.data.frame(oldGraph),as.data.frame(newGragh)),oldGraph=m)%>%reduce(rbind)
#		#	dd=add_column(dd,admix1Graph=m_name,.after=1)
#		#	dd=add_column(dd,admix2Graph=names(mylist2[[newGraghIndex]]),.after=3)
#			#print(dd)
#		#	all=rbind(all,dd)
#		}
#	}
#} 
for (OldGraphIndex in head(seq_along(mylist2),-1) ){
	m = mylist2[[OldGraphIndex]][[1]]
	m_name = names(mylist2[[OldGraphIndex]])[1]
	#print(paste0('Index1: ',OldGraphIndex))
	for	(newGraghIndex in tail(seq_along(mylist2),-1)){
		if (newGraghIndex > OldGraphIndex){
			#print(paste0('Migration1: ',OldGraphIndex-1, ' ,Migration2: ', newGraghIndex-1))
			#print(paste0('length: ', length(mylist2[[newGraghIndex]])))
			#all=rbind(all,lapply(mylist2[[newGraghIndex]],function(newGragh,oldGraph) nested(as.data.frame(oldGraph),as.data.frame(newGragh)),oldGraph=oldGraph) %>% reduce(rbind))
			dd=lapply(mylist2[[newGraghIndex]],function(newGragh,oldGraph) nested2(as.data.frame(oldGraph),as.data.frame(newGragh)),oldGraph=m)%>%reduce(rbind)
			dd=add_column(dd,admix1Graph=m_name,.after=1)
			dd=add_column(dd,admix2Graph=names(mylist2[[newGraghIndex]]),.after=3)
			#print(dd)
			all=rbind(all,dd)
		}
	}
}
all

#all %>%group_by(admix1,admix2)%>%summarize(n=length(range),supported=length(!is.na(range)),min=min(range),max=max(range))

separator = ";"
num.split.range = max(sapply(all$range, function(x) {lengths(regmatches(x, gregexpr(separator, x)))}))

mymin = function(x) {
    return(min(x,na.rm=T))
}
mymax = function(x) {
    return(max(x,na.rm=T))
}

if (num.split.range >= 1) { # if at least 2 admix in range column
    split.col.names = c( paste("ADMIX",0:num.split.range,sep="_"))
#    summary.min =  all %>% separate(range,sep=";",into=split.col.names,convert=F,remove=F) %>%
#        group_by(admix1,admix2) %>%
#        summarize(n=length(range),supported=length(!is.na(range)), across(c(split.col.names),function(x) {min(x,na.rm=T)})) %>% 
#        as.data.frame()
#    names(summary.min)[5:ncol(summary.min)] = c(paste0("minimum_admix",0:num.split.range))
                                
    summary.all =  all %>% separate(range,sep=";",into=split.col.names,convert=F,remove=F) %>%
        group_by(admix1,admix2) %>%
        summarize(n=length(range),supported=length(which(!is.na(range))), across(c(split.col.names),list(min=mymin,max=mymax ))) %>%
         as.data.frame()

                        
} else {
    summary.all =  all %>%  group_by(admix1,admix2) %>%
	summarize(n=length(range),supported=length(which(!is.na(range))), ADMIX_min=min(range),ADMIX_max=max(range))  %>%
        as.data.frame()
}

cat("\n\nHere is the summary of the results:\n")
print(summary.all)
write_tsv(summary.all,"nested/nestedAdmix.summary.txt",col_names=T)
write_excel_csv(summary.all,file="nested/nestedAdmix.summary.csv",col_names=T)
