class Solution(object):
    def addBinary(self, a, b):
        """
        :type a: str
        :type b: str
        :rtype: str
        """
        if a is None:
            return(b)
        elif b is None:
            return(a)
        elif len(a) == 0:
            return(b)
        elif len(b) == 0:
            return(a)
        else:
            aDeci = self.binStr2Deci(a)
            bDeci = self.binStr2Deci(b)
            return(self.deci2binStr(aDeci+bDeci))
    def binStr2Deci(self,a):
        """here you can safely assume
        a has at least one number and convert
        it to decimal
        """
        assert(type(a) is str or type(a) is unicode)
        aDeci = 0
        pos = 0
        for i in range(len(a)-1,-1,-1):
            aDeci += (2**pos)*int(a[i])
            pos += 1
        return(aDeci)
    def deci2binStr(self,aDeci):
        """take a decimal number and 
        convert it to binary string
        """
        assert(type(aDeci) is int or type(aDeci) is long)
        r = ""
        if aDeci <= 0:
            return("0")
        else:
            while aDeci>0:
                i = aDeci % 2
                aDeci = aDeci / 2
                r = str(i) + r
            return(r)
