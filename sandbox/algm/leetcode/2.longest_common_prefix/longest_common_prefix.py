#!/usr/bin/env python

def longestCommonPrefix(strs):
    prefixDict = {}
    for i in strs:
	for j in range(len(i)):
	    prefix = i[:(j+1)]
	    if prefix in prefixDict:
		prefixDict[prefix] += 1
	    else:
		prefixDict[prefix] = 1
    longestCommonPrefix = ''
    for i in prefixDict:
	#when count of a prefix occurence equal to # of elements
	#in strs, this is a common prefix
	if prefixDict[i] == len(strs):
	    if len(i) > len(longestCommonPrefix):
		longestCommonPrefix = i
    return longestCommonPrefix
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
