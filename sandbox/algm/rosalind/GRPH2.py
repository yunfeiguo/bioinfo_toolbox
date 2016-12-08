from Bio import SeqIO
from collections import defaultdict
def matchK(i,j, k):
    '''take two SeqIO objects
    return true if last k char of 
    first obj match first k char of
    second obj
    i and j must not be equal
    '''
    if i.seq == j.seq:
	return False
    if i.seq[len(i.seq) - k:len(i.seq)] == j.seq[0:k]:
	return True
    return False
file = 'data/rosalind_grph.txt'
fa = SeqIO.index(file, 'fasta')
k = 3
kmerPrefixDict = defaultdict(list)
kmerSuffixDict = defaultdict(list)
for i in fa:
    kmerPrefixDict[str(fa[i].seq[0:k])].append(i)
    kmerSuffixDict[str(fa[i].seq[len(fa[i].seq) - k : len(fa[i].seq)])].append(i)
for kmer in kmerPrefixDict:
    if kmer in kmerSuffixDict:
	for n1 in kmerSuffixDict[kmer]:
	    for n2 in kmerPrefixDict[kmer]:
		if fa[n1].seq == fa[n2].seq:
		    continue
		else:
		    print(fa[n1].id + ' ' + fa[n2].id)
