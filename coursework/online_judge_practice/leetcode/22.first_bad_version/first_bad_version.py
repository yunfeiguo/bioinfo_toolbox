# The isBadVersion API is already defined for you.
# @param version, an integer
# @return a bool
# def isBadVersion(version):

class Solution(object):
    def firstBadVersion(self, n):
        """
        :type n: int
        :rtype: int
        """
        if n is None or n <= 0:
            return(None)
        else:#n>0
            assert(n % 1 == 0)
            low = 1
            high = n
            middle = (1+n)/2
            #use bisection search
            while high - low > 1:
                print(str(low)+","+str(high))
                if bool(isBadVersion(low)) ^ bool(isBadVersion(high)):
                    middle = (low+high)/2
                    if isBadVersion(middle):
                        high = middle
                    else:
                        low = middle
                else:#low and high have same boolean
                    if isBadVersion(low) is True:
                        return(low)
                    else:
                        return(high)
            #at this point, low == high or low == high - 1
            if isBadVersion(low):
                return(low)
            elif isBadVersion(high):
                return(high)
            else:
                return(None)
def isBadVersion(self,n):
    if n >= 2:
        return(True)
    else:
        return(False)
