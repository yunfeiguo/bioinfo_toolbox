#!/usr/bin/env python
def reverse(s):
    '''reverse string'''
    return(s[::-1])
def isPalind(s):
    return(s==reverse(s))
s = raw_input('Please enter:')
s = s.replace(' ','')
s = s.replace('.','')
s = s.replace(',','')
s = s.lower()
print(s)
if isPalind(s):
    print('Palindromic')
else:
    print('not palindromic')
