#!/usr/bin/env python
fh = open("data/rosalind_dna.txt",'r')
countA = 0
countC = 0
countG = 0
countT = 0
for i in fh:
    for c in i:
	if c == 'a' or c == 'A':
	    countA += 1
	elif c == 'c' or c == 'C':
	    countC += 1
	elif c == 'g' or c == 'G':
	    countG += 1
	elif c == 't' or c == 'T':
	    countT += 1
	else:
	    pass
print("{0} {1} {2} {3}".format(countA, countC, countG, countT))
