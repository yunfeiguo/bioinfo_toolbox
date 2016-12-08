import collections
myList = ['a','b','c','a','d']
c = collections.Counter(myList)
print c
print c.most_common() #sort by frequency

myTuple = [(1,'abc'),(1,'a'),(2,'b'),(3,'c')]
d = collections.defaultdict(list)

for n,letter in myTuple:
    d[n].append(letter)
print d    

