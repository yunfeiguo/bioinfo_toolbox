class Solution:
    """    
    @param n: Given the range of numbers
    @param k: Given the numbers of combinations
    @return: All the combinations of k numbers out of 1..n   
    """
    def combine(self, n, k):      
        # write your code here  
        if n < 1 or k < 1 or n < k:
            return([])
        if n == k:
            return([range(1,n+1)])
        elif k == 1:
            return([[i] for i in xrange(1,n+1)])
        #k must not be 1 from now on
        ret = self.combine(n-1,k)
        containN = self.combine(n-1,k-1) #cannot be empty as k>=2, n>k
        for i in containN:
            i.append(n)
        ret.extend(containN)
        return(ret)
