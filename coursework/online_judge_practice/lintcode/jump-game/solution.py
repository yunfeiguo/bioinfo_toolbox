class Solution:
    # @param A, a list of integers
    # @return a boolean
    def canJump(self, A):
        # write your code here
        if A is None or len(A) <= 1:
            return(True)
        else:
            i = 0
            remain = A[0]
            n = len(A)
            while remain > 0 and i < n - 1:
                remain -= 1
                i += 1
                remain = max(A[i],remain)
            return(i == n - 1)
