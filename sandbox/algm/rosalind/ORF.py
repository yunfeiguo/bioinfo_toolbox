from Bio import SeqIO
from Bio.Seq import Seq
import regex as re
from sets import Set
def readGeneticCode(file):
    '''read genetic code from a file
    codons are upper case letters
    '''
    fh = open(file,"r")
    codon = {}
    for i in fh:
	s = i.rstrip()
	fields = s.split(" ")
	key = ""	
	for j in fields:
	    if len(j) == 3:
		key = str.upper(j)
	    elif len(j) == 1 or j == 'Stop':
		codon[key] = j	    
    return codon

def getAllProtein(s, codon, result):
    '''take an RNA string
    translate it to all possible
    proteins based on codon, and 
    add them to a set named result
    '''
    startRegex = re.compile('AUG', overlapped = True)
    allStart = re.finditer(startRegex, s)
    for i in allStart:
	protein = translate(s, i.start(), len(s), codon)
	if protein != "":
	    result.add(protein)
    return

def translate(s, start, end, codon):
    '''take an RNA string, beginning at start codon,
    ending at stop codon, translate it into protein sequence
    '''
    protein = []
    for i in range(start, end, 3):
	if codon[s[i:i+3]] != 'Stop':
	    #if there is no more codon, then just return ""
	    if i + 6 > end:
		return ""
	    protein.append(codon[s[i:i+3]])
	else:
	    break
    return ''.join(protein)

codon = readGeneticCode("data/genetic_code.txt")
file = 'data/rosalind_orf.txt'
fa = SeqIO.read(file,'fasta')
fa = str(fa.seq)
fa = fa.replace('T','U') #to RNA
fa = Seq(fa)
rcfa = fa.reverse_complement()
result = Set()
getAllProtein(str(fa), codon, result)
getAllProtein(str(rcfa), codon, result)
for i in result:
    print(i)
