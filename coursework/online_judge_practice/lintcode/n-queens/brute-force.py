import copy
class Solution:
    """
    Get all distinct N-Queen solutions
    @param n: The number of queens
    @return: All distinct solutions
    """
    def solveNQueens(self, n):
        # write your code here
        # use brute force to check all
        # possible solutions
        ret = []
        board = []
        for i in range(n):
            board.append(['.']*n)
        recursiveSolvNQueen(ret,board,0,n)
        joinAll(ret)
        return(ret)
def recursiveSolvNQueen(ret, board, r, n):
    if r == n:
        ret.append(copy.deepcopy(board))
        return(None)
    for c in range(n):
        board[r][c] = 'Q'
        if isValid(board,r,c,n):
            recursiveSolvNQueen(ret, board, r+1,n)
        board[r][c] = '.'
    return(None)
def isValid(board,r,c,n):
    #no need to check horizontal
    #check vertical
    for i in range(0,r):
        if board[i][c] == 'Q':
            return(False)
    #check 45deg diagonal
    i = r - 1
    j = c + 1
    while i >= 0 and j < n:
        if board[i][j] == 'Q':
            return(False)
        i -= 1
        j += 1
    #check 135deg diagonal
    i = r - 1
    j = c - 1
    while i >= 0 and j >= 0:
        if board[i][j] == 'Q':
            return(False)
        i -= 1
        j -= 1
    return(True)
def joinAll(container):
    for i in range(len(container)):
        for j in range(len(container[i])):
            container[i][j] = ''.join(container[i][j])
    return(None)
