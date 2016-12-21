from Bio import SeqIO
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
for i in fa:
    for j in fa:
	if matchK(fa[i],fa[j], k):
	    print(fa[i].id + ' ' + fa[j].id)
