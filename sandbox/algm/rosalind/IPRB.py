#!/usr/bin/env python
from math import factorial
def choose(n, k):
    assert(n >= k)
    assert(k >= 0)
    return factorial(n) / factorial(n - k) /factorial(k)
fh = open('data/rosalind_iprb.txt')
line = next(fh)
fields = line.split(' ')
k = int(fields[0])
m = int(fields[1])
n = int(fields[2])
N = k + m + n
pAA = 1 - choose(N - k, 2)*1.0/choose(N, 2)
pAaAa = choose(m, 2)*1.0/choose(N,2)
pAaaa = choose(m, 1)*choose(n, 1)*1.0/choose(N,2)
pA = pAA*1 + 0.75*pAaAa + 0.5*pAaaa
print(pA)
