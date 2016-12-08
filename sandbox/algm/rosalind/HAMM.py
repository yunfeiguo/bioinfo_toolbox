#!/usr/bin/env python
def hammingDistance(x, y):
    '''given two strings
    calculate their hamming
    distance
    '''
    assert(len(x) == len(y))
    distance = 0
    for i in range(len(x)):
	if x[i] != y[i]:
	    distance += 1 
    return distance
fh = open('data/rosalind_hamm.txt','r')
line1 = next(fh)
line2 = next(fh)
line1 = line1.rstrip()
line2 = line2.rstrip()
print(hammingDistance(line1, line2))
