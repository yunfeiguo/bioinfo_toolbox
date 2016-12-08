from Bio import SeqIO
def isCommonSubstr(s, c, fa):
    '''check all contigs except c
    see if s is their common substring
    '''
    for i in fa:
	if i == c:
	    continue
	else:
	    try:
		str(fa[i].seq).index(s)
	    except:
	        return False
    return True
file = 'data/rosalind_lcsm.txt'
fa = SeqIO.index(file,'fasta')
shortestContig = None
minLen = 0
#find shortest contig
for i in fa:
    if shortestContig is None:
	shortestContig = i
	minLen = len(fa[shortestContig].seq)
    elif len(fa[i].seq) < minLen:
	shortestContig = i
	minLen = len(fa[i].seq)
#iterate over all possible substrings of shortest contig
#output the longest common substring
shortestContigSeq = str(fa[shortestContig].seq)
maxCommonSubstr = 0
result = None
for i in range(len(shortestContigSeq)):
    for j in range(i,len(shortestContigSeq)):
	if isCommonSubstr(shortestContigSeq[i:j+1], shortestContig, fa):
	    if j+1-i > maxCommonSubstr:
		result = shortestContigSeq[i:j+1]
		maxCommonSubstr = j+1-i
print(result)		
