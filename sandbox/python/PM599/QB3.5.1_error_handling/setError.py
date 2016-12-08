mylist = [1,2,3,4,5,6]
print mylist
myset = set(mylist)
print myset
myset.discard(1)
print myset
try:
    myset.remove(1)
except KeyError:    
    print "oops no such key"
finally:    
    print myset
