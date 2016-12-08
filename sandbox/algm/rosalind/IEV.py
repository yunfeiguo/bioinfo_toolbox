file = 'data/rosalind_iev.txt'
fh = open(file,'r')
input = next(fh)
input = input.rstrip()
n1,n2,n3,n4,n5,n6 = input.split(' ')
n1 = int(n1)
n2 = int(n2)
n3 = int(n3)
n4 = int(n4)
n5 = int(n5)
n6 = int(n6)
#AA-AA AA
#AA-Aa AA,Aa
#AA-aa Aa
#Aa-Aa AA,Aa,aa
#Aa-aa Aa,aa
#aa-aa aa
countAA = n1*2 + n2*2*0.5 + n4*2*0.25
countAa = n2*2*0.5 + n3*2 + n4*2*0.5 + n5*2*0.5
print(countAA + countAa)
