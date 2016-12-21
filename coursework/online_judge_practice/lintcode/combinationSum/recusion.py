class Solution:
    # @param candidates, a list of integers
    # @param target, integer
    # @return a list of lists of integers
    def combinationSum(self, candidates, target):
        # write your code here
        if not candidates:
            return([])
        elif target < min(candidates):
            return([])
        elif len(candidates) == 1 and target % candidates[0] == 0:
            return([candidates*(target/candidates[0])])
        candidates.sort(reverse = True)
        ret = []
        for i in range(len(candidates)):
            n = candidates[i]
            if target == n:
                containN = [[n]]
            elif target < n:
                containN = []
            else:
                containN = self.combinationSum(candidates[i:],target - n)
                for j in containN:
                    j.append(n)
                #if containN is empty, then we cannot find any combination
                #that ends with n
            ret.extend(containN)
        return(ret)
