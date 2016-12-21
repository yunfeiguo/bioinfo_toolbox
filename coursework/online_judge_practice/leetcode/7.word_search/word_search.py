from sets import Set
debug = False
visited = []
class Solution:
    # @param {character[][]} board
    # @param {string} word
    # @return {boolean}
    def exist(self, board, word):
        if board == None or word == None:
            return False
        #create a 2D array to mark if visited
        for i in range(len(board)):
            visited.append([])
            for j in range(len(board[i])):
                visited[i].append(False)
        for i in range(len(board)):
            for j in range(len(board[i])):
        	    if found(board,word,i,j):
        	        return True
        return False
def found(board,word,i,j):
    '''takes board and word,look at i,j
    if board[i][j] != word[pos], then return False
    otherwise look at i+-1 and j+-1
    '''
    if len(word) == 1:
    	if i<0 or j < 0 or i < len(board) or j < len(board[i]) or visited[i][j] or word[pos] != board[i][j]:
            return False
	else:
	    return True
    else len(word) == 0:
	return False
    else:
    	visited[i][j] = True
    	if found(board,word[1:len(word)-1],pos+1,i+1,j) or found(board,word,pos+1,i-1,j) or found(board,word,pos+1,i,j+1) or found(board,word,pos+1,i,j-1):
        	return True
    	visited[i][j] = False #uncheck it because when we visit this element again, the path will be different
    	return False
def main():
    if __name__ == "__main__": 	
	board = [
		'ABCD',
		'EFGH',
		'IJKL',
		]
	board2 = ["CAA","AAA","BCD"]
	word4 = "AAB"
	word1 = 'BFG'
	word2 = 'AEL'
	word3 = 'BFGF'
	board3 = ["ABCE","SFCS","ADEE"]
	word5 = "SEE"
	board4 = ["ABCE","SFES","ADEE"]
	word6 = "ABCESEEEFS"
	x = Solution()
	print word1
	print x.exist(board,word1)
	print word2
	print x.exist(board,word2)
	print word3
	print x.exist(board,word3)
	print word4
	print x.exist(board2,word4)
	print word5
	print board3
	print x.exist(board3,word5)
	print word6
	print board4
	print x.exist(board4,word6)
main()		
