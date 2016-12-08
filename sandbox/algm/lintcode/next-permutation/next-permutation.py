class Solution:
    # @param num :  a list of integer
    # @return : a list of integer
    def nextPermutation(self, num):
        # write your code here
        if len(num) <= 1:
            return(num)
        else:
            switchIdx = None
            n = len(num)
            for i in xrange(n-1):
                if num[i] < num[i+1]:
                    switchIdx = i
            if switchIdx is None:
                #num is sorted in descending order
                num.reverse()
            else:
                num[switchIdx+1:] = sorted(num[switchIdx+1:])
                i = switchIdx+1
                while num[i] <= num[switchIdx] and i < n - 1:
                    #make sure i is NOT out of range
                    i += 1
                tmp = num[switchIdx]
                num[switchIdx] = num[i]
                num[i] = tmp
            return(num)
