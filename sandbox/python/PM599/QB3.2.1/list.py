#!/usr/bin/python
li = ['a','b','c']
print li
print li[:]
print li[0],li[1],li[2]
li.append('d')
print li
li = li + ['e']
print li
print li[-1]
print li[-1:]
print li[-2:]
print li[0:3:2]
li.extend(['e'])
print li
li.append(['e'])
print li

del li[-1]
print 'del li[-1]'
print li

print 'li.pop'
print li.pop()

print li
li = li[:-1]
print 'li=li[:-1]'
print li
