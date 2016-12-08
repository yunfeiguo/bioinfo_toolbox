#!/usr/bin/env python

#initialize list to store sequence
protSeq = []
#open pdb file
f1 = open('2Q6H.pdb', 'r')
#loop over lines in file
for next in f1:
#identify lines that contain sequences
 if next[:6] == 'SEQRES':
 #strip away white space and
 #convert line into list
  line = next.strip().split()
 #delete descriptor information
 #at beginning of each line
  del line[:4]
 #loop over amino acids in line
  for aa in line:
 #add to sequence list
    protSeq.append(aa)
    #close file
f1.close()

print len(protSeq)
dict = {}
for i in protSeq:
    if i in dict:
	dict[i] = dict[i] + 1
    else:
	dict[i] = 1
for i in dict.keys():
    print i,str(dict[i])
