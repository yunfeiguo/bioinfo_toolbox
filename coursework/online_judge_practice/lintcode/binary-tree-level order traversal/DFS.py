"""
Definition of TreeNode:
class TreeNode:
    def __init__(self, val):
        self.val = val
        self.left, self.right = None, None
"""


class Solution:
    """
    @param root: The root of binary tree.
    @return: Level order in a list of lists of integers
    """
    def levelOrder(self, root):
        # write your code here
        #try BFS first
        if root is None:
            return([])
        l = [[root,0]] 
        ret = []
        while l:
            cur = l.pop()
            if len(ret) < cur[1]+1:
                ret.append([cur[0].val])
            else:
                ret[cur[1]].append(cur[0].val)
            if cur[0].right:
                l.append([cur[0].right,cur[1]+1])
            if cur[0].left:
                l.append([cur[0].left,cur[1]+1])
        return(ret)
                
            
