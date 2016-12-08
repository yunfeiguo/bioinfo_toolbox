try:
     #fh = open("xxyxco.txt","r")
     pass
except IOError:
     print "oops, an IO error"
else:
     print "no error"
finally:
     print "final clause"

