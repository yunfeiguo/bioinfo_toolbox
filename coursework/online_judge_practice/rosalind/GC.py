#!/usr/bin/env python
from collections import Counter

def calcGC(seq):
    count = Counter(seq)
    #print(count)
    currentGC = 1.0 * (count['G'] + count['g'] + count['c'] + count['C'])
    currentGC /= count['G'] + count['g'] + count['c'] + count['C'] + count['A'] + count['a'] + count['T'] + count['t']
    return currentGC*100

fh = open('data/rosalind_gc.txt','r')
id = ""
gc = 0
seq = ""
for i in fh:
    line = i.rstrip()
    #print(line)
    if line.startswith('>'):
	if seq != "":
	    currentGC = calcGC(seq)
	    if currentGC > gc:
	        id = currentID
	        gc = currentGC
	currentID = line[1:]
	seq = ""
    else:
	seq += line
fh.close()	
if seq != "":
    currentGC = calcGC(seq)
    if currentGC > gc:
        id = currentID
        gc = currentGC
print(id)
print(gc)
