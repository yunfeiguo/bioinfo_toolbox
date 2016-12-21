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
    @return: Preorder in ArrayList which contains node values.
    """
    def preorderTraversal(self, root):
        # write your code here
        queue = [root]
        ret = []
        while queue:
            cur = queue.pop()
            if cur:
                ret.append(cur.val)
                queue.append(cur.right)
                queue.append(cur.left)
        return(ret)            
            
        

