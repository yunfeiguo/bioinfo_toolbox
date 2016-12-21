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
    @return: An integer
    """ 
    def maxDepth(self, root):
        # write your code here
        if root is None:
            return(0)
        elif root.left is None and root.right is None:
            #leaf node
            return(1)
        else:
            return(1+max(self.maxDepth(root.left),self.maxDepth(root.right)))

