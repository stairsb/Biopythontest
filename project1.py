from Bio.Seq import Seq
from Bio import SeqIO
import os
from itertools import combinations
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as pyplot
from Bio import Phylo


#aligns sequences in fasta file
cmd = "clustalw2 -infile=sequences.fasta -type=DNA -outfile=fungi.aln"
os.system(cmd)

#converts the contents of .aln to fasta format
from Bio import AlignIO
align = AlignIO.read("fungi.aln", "clustal")
f = format(align, "fasta")
with open("aligned_sequences.fasta", "w") as output_handle:
	output_handle.write(f)
	output_handle.close()

#adds the SNP and location to a list for each
SNPlist = []
locationlist = []
def addsnpstofile(SNP, location):
	SNPlist.append(SNP)
	locationlist.append(str(location))

#opens the sequences file and get the SNP and location
y=0
alignment = AlignIO.read("aligned_sequences.fasta", "fasta")
seq_names = list(SeqIO.parse("aligned_sequences.fasta", "fasta"))

#finds snps location and records number of snps found
def getsnploc(num_seq):	
	y=0
	for r in range(0,len(alignment[1].seq)):
		if alignment[0,r] != alignment[num_seq,r]:
			if alignment[0,r] != "-" and alignment[num_seq,r] != "-":
				locationlist.clear()
				y=y+1
#				print(r, alignment[0,r], alignment[2,r], y)
				addsnpstofile(alignment[num_seq,r],y)

#converts list to a string
def listToString(s): 
    str1 = ""   
    for ele in s: 
        str1 += ele   
    return str1 

#adds the sequence header, list of snps, and total snps found to a file for each sequence
num_s = 0
for seque in alignment:
	getsnploc(num_s)
	if num_s > 0:
		with open('SNPS_withtotal.txt', 'a') as f:
			f.write("%s\n" % seq_names[num_s].id)
			f.write("%s\n" % listToString(SNPlist))
			f.write("%s\n" % listToString(locationlist))
			SNPlist = []
			locationlist = []
			f.close()	
	num_s += 1

#print out a fasta file
num_s2 = 0
for seque in alignment:
	getsnploc(num_s2)
	if num_s2 > 0:
		with open('phylo.fasta', 'a') as f:
			f.write("%s\n" % seq_names[num_s2].id)
			f.write("%s\n" % listToString(SNPlist))
			SNPlist = []
			f.close()
	num_s2 += 1

cmd2 = "sed -i 's/R/>/g' phylo.fasta"
os.system(cmd2)


#aligns fasta sequences
cmd1 = "clustalw2 -infile=phylo.fasta -type=DNA -outfile=phylo.aln"
os.system(cmd1)	

#builds a phylogenetic tree
tree = Phylo.read('phylo.dnd','newick')
Phylo.draw(tree)
pyplot.savefig("SNPs_tree.png")

