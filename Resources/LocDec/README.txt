Here are the sources files of softwares used in the paper Local decoding:

compare (CompareDistance) returns the Pearson coefficient between two distances

align (ComputeDistanceAln) returns a distance matrix of proportion of mismatches 
from a multiple alignment

local (ComputeDistanceAln) returns a distance matrix of trivial or Pham dissimilarity 

Typical uses:

align -f n -s d -i 10 test.clustal testA.nex
//returns the matrix of proportions of mistmaches between the sequences of the alignment test.clustal

local -f n -s d -d l -m 0 -o 10 test.fasta testL.nex
//returns the distance matrix distL10 between the sequences of the set of sequences test.fasta

local -f n -s d -d s -m 0 -o 10 test.fasta testB.nex
//returns the distance matrix distB10 between the sequences of the set of sequences  test.fasta

local -f n -s d -d l -m 1 -o 10 test.fasta testP.nex
//returns the distance matrix distP between the sequences of the set of sequences  test.fasta

compare testA.nex testL.nex
//returns the Mantel score between testA.nex and testL.nex
