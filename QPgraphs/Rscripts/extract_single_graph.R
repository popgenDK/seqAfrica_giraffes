require(admixtools)
require(igraph)
require(dplyr)
require(readr)

args <- commandArgs(trailing=TRUE)

rdataFile <- args[1]           # "/home/users/xiaodong/Documents/Project/African1kg/buffalo/qpgraph/result_Zimbabwe/tests.txt.Rdata"
treeidx <- as.integer(args[2]) # 145
outFile <- args[3]             # "~/to_anders_timepolice_nasty.txt"

load(rdataFile)
write_delim(all_fits[[treeidx]]$edges, file=outFile)

