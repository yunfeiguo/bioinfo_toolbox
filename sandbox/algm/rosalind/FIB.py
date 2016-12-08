#!/usr/bin/env python
fh = open('data/rosalind_fib.txt','r')
args = fh.next()
args = args.split(" ")
n = int(args[0])
k = int(args[1])
previous = 1
current = 1
for i in range(3,n+1):
    tmp = current
    current = k*previous + current
    previous = tmp
print(current)    
