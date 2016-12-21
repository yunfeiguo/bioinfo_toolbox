import collections
class Solution:
    # @param board, a list of lists of 1 length string
    # @param words: A list of string
    # @return: A list of string
    def wordSearchII(self, board, words):
        # write your code here
        #first construct trie to store the words
        dict = trieNode()
        for i in words:
            dict.insert(i)
        #then examine the board
        n = len(board)
        if n == 0:
            return([])
        m = len(board[0])
        if m == 0:
            return([])
        self.ret = []
        self.n = n
        self.m = m
        self.board = board
        self.visited = [False]*n
        for i in xrange(n):
            self.visited[i] = [False]*m
        for i in xrange(n):
            for j in xrange(m):
                self.searchNeighbor('',i,j,dict)
        from sets import Set
        return(list(Set(self.ret)))
    def searchNeighbor(self,s,i,j,oldHit):
        #examine all possible strings
        #starting at board[i][j]
        #extend by one letter every time
        if i < 0 or j < 0 or i >= self.n or j >= self.m or self.visited[i][j]:
            return
        t = self.board[i][j]
        s += t
        #guiding the search, not for final output
        if t in oldHit.child:
            hit = oldHit.child[t]
            self.visited[i][j] = True
            #print('partial hit'+s)
            self.searchNeighbor(s,i-1,j,hit)
            self.searchNeighbor(s,i+1,j,hit)
            self.searchNeighbor(s,i,j-1,hit)
            self.searchNeighbor(s,i,j+1,hit)
            self.visited[i][j] = False
            if hit.val:
                self.ret.append(s)
class trieNode:
    def __init__(self):
        self.val = False
        self.child = collections.defaultdict(trieNode)
    def insert(self,s):
        '''insert a string s
        into trie
        '''
        n = len(s)
        if n > 0:
            cur = self
            for i in s:
                cur = cur.child[i]
            cur.val = True
    def search(self,s):
        '''take a string s
        return the node where s ends
        or False if s not in trie
        '''
        hit = self.partialHit(s)
        if hit is None:
            return(False)
        else:
            return(hit.val)

    def partialHit(self,s):
        '''
        return None if not partial hit found
        otherwise return the node where partial hit ends
        '''
        if not s:
            return(None)
        else:
            cur = self
            for i in s:
                if i in cur.child:
                    cur = cur.child[i]
                else:
                    return(None)
            return(cur) #could partial or full hit

def demo():	
	x = Solution()
	print(x.wordSearchII(["aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa","aaaaaaaaaaaaaaaaaaaaaaaaaaaaab"], ["baaaaaaaaaaaaa","a","aa","aaaa","aaaax","abaaabbaz"]))
'''            
x = trieNode()
print(x.search('s'))
print(x.insert('s'))
print(x.insert('shell'))
print(x.search('s'))
print(x.search('she'))
'''
import cProfile
cProfile.run('demo()')
