#!/usr/bin/env python
fh = open('data/rosalind_iprb.txt')
line = next(fh)
fields = line.split(' ')
k = int(fields[0])
m = int(fields[1])
n = int(fields[2])
N = k + m + n
#1 - complimentary probability
pA = 1 - (0.25*m*(m - 1) + m * n + n*(n-1))*1.0/N/(N-1)
print(pA)
