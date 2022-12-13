# Biopythontest
The code here was to test out biopython for SNP calling and a couple other applications.
## Functionality
This script takes in fasta sequences and returns:

	1. file with the SNPs found for each isolate and the number of SNPs that were found (.txt)
	2. fasta file containing the SNP sequences for each isolate
	3. Phylogenetic tree build using the SNP data (.png)
	4. Aligned sequences (.aln)
	5. Aligned snp seuences (.aln)
	6. two dnd files that were outputs from clustalw2

This can be used to find SNPs bewteen genes or genomic sequences
The first sequences in the sequences.fasta will be used at the reference genomic sequence
The sequences included here are all from the same fungal species
The tree was constructed using clustalw2
