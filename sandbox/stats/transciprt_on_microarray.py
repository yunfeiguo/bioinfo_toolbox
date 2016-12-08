#!/usr/bin/env python
import math
from sets import Set
from random import randint
def nchoosek(n, k):
    combinatorial = 1
    for i in range(k):
	combinatorial *= (n - i) / (i + 1)
'''
we want to calcualte the number of holes required to 
make sure >99% probability that all n kinds of 
transcripts bind to m holes. Each kind of transcript
has an infinite number of copies.

how to formulate this?
P(exactly n kinds of transcripts bind) = 1 - P(exactly n - 1 kinds bind) - ... - P(exactly 1 kind binds)
P(x kinds of transcripts) = choose(n,x)*((n-x)/n)^m

the probabilities are difficult to calculate, an alternative method is Monte Carlo. We sample from a discrete distribution {1,2,...,n}, probability of each value is equal to 1/n. We repeat m times, simulating binding of m transcripts. At then end of this experiment, we know how many different transcripts bind. We carry out 10,000 experiments, and can calculate frequency of having all n transcripts.
'''
nExperiment = 100
nTranscript = 20000
nHole = nTranscript*15
countForAllTranscript = 0 #number of times seeing all n transcripts
for i in range(1,nExperiment+1):
    seen = Set()
    for j in range(nHole):
	seen.add(randint(1,nTranscript))
    if len(seen) == nTranscript:
	countForAllTranscript += 1
    if i % 100 == 0:
	print("{0} experiments done.".format(i))
print("{0} holes, {1} transcripts".format(nHole, nTranscript))
print("observed {0} times ({1}%) of all different types of transcripts binding".format(countForAllTranscript,100.0*countForAllTranscript/nExperiment))
print("in {0} experiments.".format(nExperiment))
