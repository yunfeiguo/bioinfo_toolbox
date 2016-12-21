from collections import Counter
from Bio import SeqIO
def arrayToString(a):
    s = []
    for i in a:
	s.append(str(i))
    return " ".join(s)
def transpose(a):
    '''transpose a 2D matrix
    '''
    b = []
    for i in range(len(a[0])):
	row = []
	for j in range(len(a)):
	    row.append(a[j][i])
	b.append(row)
    return b

aCount = []
tCount = []
cCount = []
gCount = []
sequences = []
letterMatrix = []
length = 0
for seq in SeqIO.parse('data/rosalind_cons.txt','fasta'):
    if len(sequences) == 0:
	length = len(str(seq.seq))
    sequences.append(str.upper(str(seq.seq)))
    #assume all sequences have equal length
    assert(length == len(str(seq.seq)))

for i in range(len(sequences)):
    row = []
    for j in range(length):
        row.append(sequences[i][j])
    letterMatrix.append(row)

letterMatrix = transpose(letterMatrix)    
for i in range(length):
    count = Counter(letterMatrix[i])
    aCount.append(count['A'])
    tCount.append(count['T'])
    cCount.append(count['C'])
    gCount.append(count['G'])
	
consensus = ""	
for i in range(length):
    maxCount = max([aCount[i],tCount[i],cCount[i],gCount[i]])
    if maxCount == aCount[i]:
	consensus += 'A'
    elif maxCount == tCount[i]:
	consensus += 'T'
    elif maxCount == cCount[i]:
	consensus += 'C'
    else:
	consensus += 'G'
print(consensus)
print 'A: ' + arrayToString(aCount)
print 'C: ' + arrayToString(cCount)
print 'G: ' + arrayToString(gCount)
print 'T: ' + arrayToString(tCount)
