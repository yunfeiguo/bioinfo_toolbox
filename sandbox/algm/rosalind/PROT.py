#!/usr/bin/env python
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

geneticCode = readGeneticCode("data/genetic_code.txt")
fh = open("data/rosalind_prot.txt","r")
rna = str.upper(next(fh))
rna = rna.rstrip()
aa = ""
for i in range(0,len(rna) - 2, 3):
    newAA = geneticCode[rna[i:(i+3)]]
    if newAA == 'Stop':
	aa += " "
    else:
	aa += newAA
print(aa)    
