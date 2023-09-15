//Parameters for the coalescence simulation program : fsimcoal2.exe
5 samples to simulate :
//Population effective sizes (number of genes)
NW
NK
NSA
NSC
NR
//Samples sizes and samples age
10
8
24
16
24
//Growth rates  : negative growth implies population expansion
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
6  historical event
T1 2 3 1 FAC1 0 0
TA$ 4 1 MIGR 1 0 0
TA$ 4 4 0 FACA 0 0
T2 1 0 1 FAC2 0 0
TB 4 3 1 FACB 0 0
T3 3 0 1 FAC3 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ  1   0   2.12E-08 OUTEXP
