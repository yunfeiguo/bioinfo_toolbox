class Solution(object):
    def titleToNumber(self, s):
        """
        :type s: str
        :rtype: int
        """
        if s is None or len(s) == 0:
            return(0)
        else:
            return(self.alphaStr2deci(s))
    def alphaStr2deci(self, s):
        """take str/unicode s of length>1
        convert it to decimal
        """
        pos = 0
        r = 0
        for i in range(len(s)-1,-1,-1):
            r += (26**pos)*self.alpha2deci(s[i])
            pos += 1
        return(r)
    def alpha2deci(self, c):
        """take single str/unicode
        and return decimal number
        assume only uppercase letter
        """
        if(len(c) == 0):
            return(0)
        else:
            return(ord(c)-ord('A')+1)
