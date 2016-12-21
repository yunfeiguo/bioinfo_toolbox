class Solution(object):
    def trailingZeroes(self, n):
        """
        :type n: int
        :rtype: int
        """
        if n is None or n <= 0:
            return(0)
        else:
            assert(n % 1 == 0)
            #the number of trailing 0 is 
            #equal to number of 5s in multiplication
            maxFive = int(math.floor(math.log(n)/math.log(5)))
            countZero = 0
            for i in range(1,maxFive+1,1):
                countZero += math.floor(n/(5**i))
            countZero = int(countZero)
            return(countZero)
