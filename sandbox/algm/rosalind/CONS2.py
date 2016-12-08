from Bio import SeqIO

aCount = []
tCount = []
cCount = []
gCount = []
length = -1    
count = None
for seq in SeqIO.parse('data/rosalind_cons.txt','fasta'):
    if length == -1:
	length = len(seq.seq)
	count = [0]*length
	for i in range(length):
	    count[i] = {'A':0,'T':0,'C':0,'G':0}
    else:
	assert(length == len(seq.seq))
    #assume all sequences are of same lengths
    for i in range(length):
	c = str.upper(str(seq.seq[i]))
	count[i][c] += 1

consensus = ""	
for i in range(length):
    maxCount = max([count[i]['A'],count[i]['T'],count[i]['C'],count[i]['G']])
    aCount.append(str(count[i]['A']))
    tCount.append(str(count[i]['T']))
    cCount.append(str(count[i]['C']))
    gCount.append(str(count[i]['G']))
    if maxCount == count[i]['A']:
	consensus += 'A'
    elif maxCount == count[i]['T']:
	consensus += 'T'
    elif maxCount == count[i]['C']:
	consensus += 'C'
    else:
	consensus += 'G'
print(consensus)
print 'A: ' + ' '.join(aCount)
print 'C: ' + ' '.join(cCount)
print 'G: ' + ' '.join(gCount)
print 'T: ' + ' '.join(tCount)
