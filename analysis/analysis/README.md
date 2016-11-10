# Dotaligner HPC clustering pipeline

*Reads a set of .fasta files, folds them with the partition function algorithm, extracts the base pairing probabilities, then performs all-vs all pairwise comparions in a HPC job array*

Usage: 

`./launcher.sh file.fasta`

Dependencies (version tested): 
- DotAligner (0.3)
- ViennaRNA (2.1.3)
- LocaRNA (1.8.11)
- R (3.2.3)
