//Parameters for the coalescence simulation program : fastsimcoal.exe
2 samples to simulate :
//Population effective sizes (number of genes)
NEU
NAF
//Samples sizes and samples age
SAMPLE_SIZE
SAMPLE_SIZE
//Growth rates	: negative growth implies population expansion
R1
0
//Number of migration matrices : 0 implies no migration between demes
3
//Migration matrix 0
0 MIG_AFEU
MIG_AFEU 0
//Migration matrix 1
0 MIG_AFB
MIG_AFB 0
//Migration matrix 2
0 0
0 0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
3 historical events
TEU 0 0 0 RES1 0 1
TOOA 0 1 1 1 0 2 // move all genes from deme 0 (EU) to 1 (AF)
TAF 1 1 0 RES2 0 2 
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 1e-8 OUTEXP
