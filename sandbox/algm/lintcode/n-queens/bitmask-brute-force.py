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
        # use binary flag to determine validity
        # use 1 to denote no Queen, 0 otherwise
        # there are 2n-1 45degree diagonals
        colFlag = [True]*n
        deg45Flag = [True]*(2*n-1)
        deg135Flag = [True]*(2*n-1)
        board = []
        ret = []
        for i in range(n):
            board.append(['.']*n)
        recurSolvNQueen(ret, board, colFlag, deg45Flag, deg135Flag, 0, n)
        joinAll(ret)
        return(ret)
def recurSolvNQueen(ret, board, colFlag, deg45Flag, deg135Flag, r, n):
    if r == n:
        ret.append(copy.deepcopy(board))
        return(None)
    for c in range(n):
        if colFlag[c] and deg45Flag[c+r] and deg135Flag[n-1-c+r]: #valid placement
            colFlag[c] = deg45Flag[c+r] = deg135Flag[n-1-c+r] = False
            board[r][c] = 'Q'
            recurSolvNQueen(ret, board, colFlag, deg45Flag, deg135Flag, r+1, n)
            board[r][c] = '.'
            colFlag[c] = deg45Flag[c+r] = deg135Flag[n-1-c+r] = True
    return(None)
def joinAll(container):
    for i in range(len(container)):
        for j in range(len(container[i])):
            container[i][j] = ''.join(container[i][j])
    return(None)        
