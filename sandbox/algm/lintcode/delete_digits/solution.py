class Solution:
    """
    @param A: A positive integer which has N digits, A is a string.
    @param k: Remove k digits.
    @return: A string
    """
    def DeleteDigits(self, A, k):
        # write you code here
        if A is None:
            return('0')
        else:
            i = 0
            n = len(A)
            while k > 0:
                if i == n - 1 or int(A[i]) > int(A[i+1]):
                    #we are at the end or encounter reverse sorted numbers
                    A = A[:i] + A[i+1:]
                    k -= 1
                    n -= 1
                    if i > 0:
                        i -= 1#go back for comparison of i-1's new neighbor
                else:
                    i += 1
            i = 0
            #remove leading zeros
            while i < n and A[i] == '0':
                i += 1
            A = A[i:]
            if A == '':
                return('0')
            else:
                return(A)
x = Solution()
print(x.DeleteDigits('178542',4))
print(x.DeleteDigits('1782',4))
print(x.DeleteDigits('',0))
print(x.DeleteDigits('1',0))
print(x.DeleteDigits('100',1))
print(x.DeleteDigits('90249',2))
