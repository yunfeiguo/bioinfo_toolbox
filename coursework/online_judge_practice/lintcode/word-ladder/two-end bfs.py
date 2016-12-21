import string
class Solution:
    # @param start, a string
    # @param end, a string
    # @param dict, a set of string
    # @return an integer
    def ladderLength(self, start, end, dict):
        # write your code here
        #two-end BFS
        if start == end:
            return(2)
        sHead = [start]
        sEnd = [end]
        dist = 2
        while sHead and sEnd:
            if len(sHead) > len(sEnd):
                search = sEnd
                backup = sHead
            else:
                search = sHead
                backup = sEnd
            nextSearch = []
            for i in search:
                wordList = list(i)                
                for j in range(len(i)):
                    oldLetter = wordList[j]
                    for k in string.lowercase:
                        wordList[j] = k
                        newWord = ''.join(wordList)
                        if newWord in backup:
                            return(dist)
                        elif newWord in dict:
                            dict.discard(newWord)
                            nextSearch.append(newWord)
                    wordList[j] = oldLetter
            if len(sHead) > len(sEnd):
                sEnd = nextSearch
            else:
                sHead = nextSearch
            dist += 1
        return(0)
                        
                


