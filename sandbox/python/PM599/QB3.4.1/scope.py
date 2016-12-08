def myFunc(a):
    a += 1
    print 'in myFunc: a = %d' % a
b = 1
print b
myFunc(b)
print b

def myFunc2(a):
    a[0] += 1
    print 'in myFunc: a[0] = %d' % a[0]
b = [1,2]
print b
myFunc2(b) #here b is passed as reference, so change to its elements is persitent
print b

def myFunc3(a):
    a[0] += 1
    print 'in myFunc: a[0] = %d' % a[0]
b = {0:1,1:2}
print b
myFunc3(b) #here b is passed as reference, so change to its elements is persitent
print b
