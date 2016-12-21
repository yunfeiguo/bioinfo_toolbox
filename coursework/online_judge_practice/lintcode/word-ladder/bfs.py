import string
class Solution:
    # @param start, a string
    # @param end, a string
    # @param dict, a set of string
    # @return an integer
    def ladderLength(self, start, end, dict):
        # write your code here
        #use BFS
        dist = 0
        if start == end:
            return(2)
        from Queue import Queue
        toVisit = Queue()
        toVisit.put(start)
        dict.add(end)
        while not toVisit.empty():
            n = toVisit.qsize()
            dist += 1
            for i in range(n):
                cur = toVisit.get()
                if cur == end:
                    return(dist)
                else:
                    addNeighbor(cur,dict,toVisit)
        return(0)
def addNeighbor(cur,dict,toVisit):
    dict.discard(cur) #circles may exist, remove it so we won't visit it again
    curInList = list(cur)
    for i in range(len(cur)):
        oldLetter = curInList[i]
        for j in string.lowercase:
            curInList[i] = j
            word = ''.join(curInList)
            if word in dict:
                toVisit.put(word)
                dict.discard(word) #do NOT visit a word twice!
        curInList[i] = oldLetter
    return(None)


