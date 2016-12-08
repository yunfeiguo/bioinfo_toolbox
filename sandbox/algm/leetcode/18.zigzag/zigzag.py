class Solution(object):
    def convert(self, s, numRows):
        """
        :type s: str
        :type numRows: int
        :rtype: str
        """
        if s is None or len(s) == 0 or numRows == 0:
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
            #pattern
            #0       8
            #1     7 9
            #2   6   10     ...
            #3 5     11  13
            #4       12
            #so here we just need a function
            #that converts index in s to index in r
            r = ['']*len(s)
            for i in range(0,len(s),1):
                r[self.convertIdx(i,numRows,len(s))] = s[i]
            return(r)
    def convertIdx(self,i,n,l):
        '''take an index in input str
        convert it to index in output str
        in a zigzag fashion
        '''
        r = 0
        #determine if the index appears in vertical lines
        mod = i % (2*n-2) #every cycle consists of 2n-2 elements
        if mod >= n: #not in vertical lines
            #calculate its index based on the element on vertical line
            #before it, at the same horizontal line
            r = self.convertIdx(i-2*(mod-(n-1)),n,l)+1
        else:
            #vertical line is column number
            #mod is row number
            #distance between adjacent elements on the same col can be easily calculated
            currentVLine = math.floor(i/(2*n-2)) #line number begins with 0
            totalVLine = math.floor((l-1)/(2*n-2))
            if mod == 0:
                #top of vertical line
                r = currentVLine
            else:
                if currentVLine == totalVLine:
                    return(self.convertIdx(i-2*mod,n,l)+1)
                if mod == 1:
                    r = self.convertIdx(i-1,n,l) + (totalVLine-currentVLine) + 2*currentVLine + 1
                else:
                    #it is necessary to adjust for last col
                    if (l-1) % (2*n-2) >= mod - 1: #l-1 because start is 0
                        adjust = 0
                    else:
                        adjust = -1
                    if mod == n-1:
                        r = self.convertIdx(i-1,n,l) + 2*(totalVLine-currentVLine) + currentVLine + 1
                    else:
                        r = self.convertIdx(i-1,n,l) + 2*(totalVLine-currentVLine) + 2*currentVLine + 1
                    r += adjust
        r = int(r)
        return(r)

#test
convertIdx(2,3,4)
