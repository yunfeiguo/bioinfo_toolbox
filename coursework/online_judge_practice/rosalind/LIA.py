import math
def nchoosek(n, k):
    '''return combination number
    for choosing k elements from n elements
    '''
    combinatorial = 1
    for i in range(k):
	combinatorial *= (n - i)*1.0 / (i + 1)
    return combinatorial
file = 'data/rosalind_lia.txt'
fh = open(file,'r')
input = next(fh)
input = input.rstrip()
k,n = input.split(' ')
k = int(k)
n = int(n)
result = 0
total = math.pow(2,k)
'''
for each generation, we don't know how many offsprings for each genotype
but we know for sure, the probability of A allele is equal to prob(a), and 
prob(B) = prob(b). Therefore prob(Aa) = prob(Bb) = 1/2, prob(AaBb) = 1/4.
Then the rest of the problem becomes a binomial distribution.
'''
for i in range(n,int(total)+1):
    result += nchoosek(total,i)*math.pow(0.25,i)*math.pow(0.75,total-i)
print(result)
