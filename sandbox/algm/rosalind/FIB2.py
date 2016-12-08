#!/usr/bin/env python
def fib(n, k = 1):
    assert(n >= 1)
    previous, current = 0, 1
    for i in range(n - 1):
	previous, current = current, k*previous + current
    return current
fh = open('data/rosalind_fib.txt','r')
args = fh.next()
args = args.split(" ")
n = int(args[0])
k = int(args[1])
print(fib(n,k))    
