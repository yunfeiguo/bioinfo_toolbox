class Solution:
    """
    @param S: A set of numbers.
    @return: A list of lists. All valid subsets.
    """
    def subsetsWithDup(self, S):
        # write your code here
        #iterative solution
        ret = [[]]
        S.sort()
        prevSize = 0
        for i in range(len(S)):
            if i > 0 and S[i] == S[i-1]:
                #duplicate
                copyFrom = prevSize
            else:
                copyFrom = 0
            prevSize = len(ret)
            for j in range(copyFrom,prevSize):
                ret.append(ret[j][:]+[S[i]])
        return(ret)
