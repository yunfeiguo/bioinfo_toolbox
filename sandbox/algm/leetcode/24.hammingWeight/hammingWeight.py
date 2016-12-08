class Solution(object):
    def hammingWeight(self, n):
        """
        :type n: int
        :rtype: int
        """
        if n is None or n <= 0 or (n % 1) != 0:
            #must be unsigned integer
            return(0)
        else:
            s = 0
            for i in str(bin(n))[2:]:
                s += int(i)
            return(s)
x = Solution()
print(x.hammingWeight(0))
print(x.hammingWeight(1))
print(x.hammingWeight(11))
print(x.hammingWeight(-1))
