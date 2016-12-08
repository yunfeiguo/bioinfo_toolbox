"""
Definition of TreeNode:
class TreeNode:
    def __init__(self, val):
        self.val = val
        self.left, self.right = None, None
"""
class Solution:
    """
    @param root: The root of the binary search tree.
    @param value: Remove the node with given value.
    @return: The root of the binary search tree after removal.
    """    
    def removeNode(self, root, value):
        # write your code here
        if root is None:
            return(root)
        oldRoot = root
        if root.val == value:
            return(rm(root))
        pre = None
        while root:
            if root.val < value:
                pre = root
                root = root.right
            elif root.val > value:
                pre = root
                root = root.left
            else:
                newSubRoot = rm(root)
                if value > pre.val:
                    pre.right = newSubRoot
                else:
                    pre.left = newSubRoot
                break
        return(oldRoot)
def rm (root):
    if root.left is None and root.right is None:
        return(None)
    elif root.left is None:
        return(root.right)
    elif root.right is None:
        return(root.left)
    else:
        if root.right.left is None:
            root.right.left = root.left
            return(root.right)
        else:
            pre = root.right
            dive = pre.left
            while dive.left:
                pre = dive
                dive = dive.left
            pre.left = None
            root.val = dive.val
            return(root)
