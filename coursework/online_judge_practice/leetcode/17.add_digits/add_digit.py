class Solution(object):
    def addDigits(self, num):
        """
        :type num: int
        :rtype: int
        """
        assert(num >= 0)
        if num == 0:
            return(0)
        else:
            mod = num % 9
            if mod == 0:
                return(9)
            else:
                return(mod)
