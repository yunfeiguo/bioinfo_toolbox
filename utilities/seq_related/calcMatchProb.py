import math
def ncombo(n,m):
    """
    calculate combinatorial number
    """
    return(math.factorial(n)/math.factorial(n-m)/math.factorial(m))
def probMatch(m,l):
    """
    given length of sequence l,
    number of matches m,
    calculate probability seeing m or more matches
    assume, a,t,c,g are independent, equal prob.
    """
    p = 0
    for i in range(m,l+1,1):
        #all cases for i matches in l
        numerator = float(ncombo(l,i)*(4**i)*((16-4)**(l-i)))
        #all cases in l regardless of matches
        denominator = float(16**l)
	p += numerator/denominator
    return(p)
print(probMatch(11,26))
