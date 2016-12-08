#!/usr/bin/env python

class hashSet(object):
    def __init__(self,numBuckets):
	if not type(numBuckets) == int:
	    raise ValueError
	if numBuckets <= 0:
	    raise ValueError
	self.hash={}
	for i in range(numBuckets):
	    self.hash[i]=[]
	self.numBuckets=numBuckets
    def hashValue(self,e):
	if not type(e) == int:
	    raise ValueError
	return e % self.numBuckets
    def member(self,e):
	if not type(e) == int:
	    raise ValueError
	if e in self.hash[self.hashValue(e)]:
	    return True
	else:
	    return False
    def insert(self,e):
	if not type(e)==int:
	    raise ValueError
	else:
	    self.hash[self.hashValue(e)].append(e)
    def remove(self,e):
	if not type(e) == int:
	    raise ValueError
	if not e in self.hash[self.hashValue(e)]:
	    raise ValueError
	self.hash[self.hashValue(e)].remove(e)
    def getNumBuckets(self):
	return self.numBuckets
    def __str__(self):
	sum=0
	for i in self.hash.keys():
	    sum+=len(self.hash[i])
	print "There are "+str(self.numBuckets)+" buckets and "+str(sum)+" elements."
