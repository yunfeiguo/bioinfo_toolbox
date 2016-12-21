#!/usr/bin/env python
fh = open('data/rosalind_subs.txt')
s = next(fh)
s = s.rstrip()
t = next(fh)
t = t.rstrip()
positions = ""
assert(len(t) <= len(s))
for i in range(0, len(s) - len(t)):
    if t == s[i:(i+len(t))]:
	positions += str(i + 1)
	positions += " "
print(positions)	
