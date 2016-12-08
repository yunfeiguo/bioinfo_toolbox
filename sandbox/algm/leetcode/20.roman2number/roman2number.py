class Solution(object):
    def romanToInt(self, s):
        """
        :type s: str
        :rtype: int
        """
        if s is None or len(s) == 0:
            return(0)
        else:
            Roman = {'I':1, 'V':5, 'X':10, 'L':50, 'C':100, 'D':500, 'M':1000}
            n = 0
            last = 0
            for i in s:
                current = Roman[i]
                if last != 0 and current > last: #we got something like IV, XL
                    n -= last
                    n += current - last
                    last = current
                    continue
                #normally, we would just add whatever we got
                n += current
                last = current
            return(n)
