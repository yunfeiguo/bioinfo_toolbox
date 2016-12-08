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
        #use DFS search
        #each node is accessed once, O(N)
        return(dfsHeight(root) != -1)
def dfsHeight(root):
    if root is None:
        return(0)
    left = dfsHeight(root.left)
    if left == -1:
        return(-1)
    right = dfsHeight(root.right)
    if right == -1:
        return(-1)
    if abs(left - right) > 1:
        return(-1)
    else:
        return(max(left,right) + 1)
        

