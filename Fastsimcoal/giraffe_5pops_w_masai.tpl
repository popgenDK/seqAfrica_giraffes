//Parameters for the coalescence simulation program : fsimcoal2.exe
5 samples to simulate :
//Population effective sizes (number of genes)
NW
NK
NR
NMA
NSC
//Samples sizes and samples age
10
8
24
10
16
//Growth rates  : negative growth implies population expansion
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
8 historical event
TM1 3 2 MIGR1 1 0 0
TM1 3 3 1 FACM1 0 0
TM2 2 1 MIGR2 1 0 0
TM2 2 2 1 FACM2 0 0
T1 1 0 1 FAC1 0 0
T2 4 3 1 FAC2 0 0
TM3 2 3 1 FACM3 0 0
T3 3 0 1 FAC3 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ  1   0   2.12E-08 OUTEXP
