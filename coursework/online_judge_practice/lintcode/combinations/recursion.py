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
        ret = []
        for i in xrange(n,k-1,-1):
            cur = self.combine(i-1,k-1)
            for j in cur:
                j.append(i)
            ret.extend(cur)
        return(ret)
