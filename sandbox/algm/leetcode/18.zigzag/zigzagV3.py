import math
class Solution(object):
    def convert(self, s, numRows):
        """
        :type s: str
        :type numRows: int
        :rtype: str
        """
        if s is None or len(s) == 0 or numRows <= 0:
            return("")
        elif numRows == 1:
            return(s)
        elif numRows == 2:
            r = ""
            #1 3 5 7
            #2 4 6 8
            for i in range(0,len(s),2):
                r += s[i]
            for i in range(1,len(s),2):
                r += s[i]
            return(r)
        else:
            r = self.fillMatrix(s,numRows)
            return(self.outputMatrix(r))
    def fillMatrix(self,s,n):            
        #pattern
        #0       8
        #1     7 9
        #2   6   10     ...
        #3 5     11  13
        #4       12
        l = len(s)
        totalVLine = int(math.floor((l-1)/(2*n-2))) + 1
        r = ['']*n*2*totalVLine #row major layout
        goingVertical = True
        moveRight = False
        currentRow = 0
        currentCol = 0
        for i in range(0,len(s),1):
            mod = i % (2*n - 2)
            if mod == 0:
                #on top
                goingVertical = True
                if currentRow != 0:
                    currentCol += 1
                currentRow = 0
            elif mod == n-1:
                #on bottom
                currentRow = n-1
                goingVertical = False
                moveRight = True
            else:
                if goingVertical:
                    currentRow += 1
                else:
                    currentRow -= 1
                    if moveRight:
                        currentCol += 1
                        moveRight = False
            r[currentRow*2*totalVLine + currentCol] = s[i]
        return(r)
    def outputMatrix(self,s):
        '''output 1D array one by one
        skipping empty elements
        '''
        r = ""
        for i in s:
            if i is not '':
                r += i
        return(r)
