#!/usr/bin/env python

def longestCommonPrefix(strs):
    if len(strs) > 0:
	#sample will the string we extract prefix
    	sample = strs[0]
    	longestCommonPrefix = []
	for pos in range(len(sample)):
	    targetLetter = sample[pos]
	    end = False
	    for i in strs:
	        if pos > len(i)-1 or i[pos] != targetLetter:
	    	    end = True
	    if end:
	        break
	    else:
	        longestCommonPrefix.append(targetLetter)
        return ''.join(longestCommonPrefix)	    
    else:
        return ''
print "longestCommonPrefix([])"
print longestCommonPrefix([])    
print "longestCommonPrefix([''])"    
print longestCommonPrefix([''])    
print "longestCommonPrefix(['',''])"
print longestCommonPrefix(['',''])    
print "longestCommonPrefix(['abc','abcde','efg','efd','lkm'])"
print longestCommonPrefix(['abc','abcde','efg','efd','lkm'])    
print "longestCommonPrefix(['abc','abcde','abcefg','abcefd','abclkm'])"
print longestCommonPrefix(['abc','abcde','abcefg','abcefd','abclkm'])    
