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
    @return: True if this Binary tree is Balanced, or false.
    """
    def isBalanced(self, root):
        # write your code here
        if root is None:
            return(True)
        leftDepth = getDepth(root.left)
        rightDepth = getDepth(root.right)
        return(abs(leftDepth - rightDepth) <= 1 and self.isBalanced(root.left) and self.isBalanced(root.right))
def getDepth(root):
    if root is None:
        return(0)
    else:
        return(1+max(getDepth(root.left),getDepth(root.right)))
