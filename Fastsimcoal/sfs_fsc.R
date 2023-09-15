library(genio)
### Load data
df <- read_plink("/projects/seqafrica/data/giraffe/qpGraph/7pops/Okapi_GCam_variable_sites_mergeAAme-OJoh_nomultiallelics_noindels_10dp_2het_nodup.polyIngroup.nonpolyOutgroup")
poplabel <- read.table("/projects/seqafrica/people/rdc143/sfs/ind51.7pops.wOutgroups.txt")


fam <- as.data.frame(df$fam)
fam$pop <- poplabel[, 2][match(fam[,1], poplabel[, 1])]

X <- as.matrix(df$X[, !is.na(fam$pop)])

popinfo <- fam[!is.na(fam$pop), ]

print(paste0("No. of sites: ", nrow(X)))

missing <- apply(is.na(X), 2, mean)
out <- subset(popinfo, pop %in% c("Okapi", "Pronghorn"))[, 1]


###Polarize genos in relation to outgroup
keep <- !apply(is.na(X[, out]), 1, any)
X <- X[keep, ]


flip <- rowSums(X[, out] ) == 6
print("Sites needing to be flipped:")
table(flip)

X[flip, ] <- 2 - X[flip, ]


### Check for conflicts between outgroup inds
keep <- rowSums(X[, out] ) == 0
X <- X[keep, ]

print(paste0("No. of sites after removing NA's and conflicting sites in outgroup: ", nrow(X)))


### check conflicts between popinfo and fam
print("Inds in pop label file:")
table(k <- popinfo[, 1] %in% poplabel)

pops <- unique(popinfo$pop)


#5d
pops5d <- c("West_African", "Kordofan", "Southern_African", "Southern_Central", "Reticulated")
#pops5d <- c("West_African", "Kordofan", "Reticulated", "Masai", "Southern_Central") # alt_pops

tab <- table(rowSums(X[, which(popinfo$pop == pops5d[1])]),
             rowSums(X[, which(popinfo$pop == pops5d[2])]),
             rowSums(X[, which(popinfo$pop == pops5d[3])]),
             rowSums(X[, which(popinfo$pop == pops5d[4])]),
             rowSums(X[, which(popinfo$pop == pops5d[5])]))

tab <- aperm(tab, 5:1) # change order of dims to correspond with fsc/dadi,
                       # flattening array in r iterates over dim 1 first, should be iterating through dim 4 first, then 3 etc.

tab[1] <- 0
tab[length(tab)] <- 0
tab[1] <- sum(tab) / (850879 / 83314370) # Fraction of variable sites from another sfs estimate with more thorough filtering


file <- file("../fsc/giraffe_5pops_DSFS.obs", open = "wt")
writeLines("1 observations. No. of demes and sample sizes are on next line", file)
cat(paste0(length(pops5d), "\t"), file = file)
write(dim(tab), file) # maybe these should be in reverse order again, but fsc doesnt seem to read the dims from here anyway
write(as.vector(tab),ncol=prod(dim(tab)), file)
close(file)
