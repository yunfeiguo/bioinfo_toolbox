import sys
def mortalFibonacci(n, m):
    assert(n > 0)
    assert(m > 1)
    result = [0] * (n + 1)
    result[0] = 0
    result[1] = 1
    #before we reach m+1 months
    #conventional Fibonacci series
    for i in range(2, m + 1):
	result[i] = result[i - 1] + result[i - 2]
    if n >= m + 1:	
	result[m + 1] = sum(result[1:m])
    #o(n) = o(n - 2) + o(n-3)+...+o(n-m)
    #o(n-1) = o(n-3) + ... + o(n-1-m)
    #o(n) = o(n-1)+o(n-2) - o(n-1-m) for n >= m+2
    for i in range(m + 2, n + 1):
	#rabbits from last month
	result[i] = result[i - 1] + result[i - 2] - result[i - m - 1]
    return result[n]

if len(sys.argv) != 3:
    print('Usage: <n> <m>')
    sys.exit(1)
n = int(sys.argv[1])
m = int(sys.argv[2])
print(mortalFibonacci(n, m))
