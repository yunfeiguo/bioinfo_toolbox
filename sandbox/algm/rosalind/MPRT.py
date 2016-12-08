import urllib
import regex as re

def getNextID(file):
    fh = open(file,'r')
    for i in fh:
	yield i.rstrip()
def getFasta(id):
    fa = urllib.urlopen('http://www.uniprot.org/uniprot/{0}.fasta'.format(id))
    seq = fa.read()
    seq = ''.join(seq.split('\n')[1:])
    seq = seq.upper()
    return seq
file = 'data/rosalind_mprt.txt'
regex = re.compile('N[GAVLIMCFYWHKRQNEDST][ST][GAVLIMCFYWHKRQNEDST]')
for id in getNextID(file):
    seq = getFasta(id)
    pos = []
    #why overlapped = True?? Because motifs may overlap
    for match in re.finditer(regex, seq, overlapped = True):
	pos.append(str(match.start() + 1))
    if len(pos) > 0:
	print(id)
	print(' '.join(pos))
