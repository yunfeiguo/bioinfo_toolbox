"""
Definition of TreeNode:
class TreeNode:
    def __init__(self, val):
        self.val = val
        self.left, self.right = None, None
"""


class Solution:
    """
    @param preorder : A list of integers that preorder traversal of a tree
    @param inorder : A list of integers that inorder traversal of a tree
    @return : Root of a tree
    """
    def buildTree(self, preorder, inorder):
        # write your code here
        assert(len(preorder) == len(inorder))
        if len(preorder) == 0:
            return(None)
        elif len(preorder) == 1:
            return(TreeNode(preorder[0]))
        else:
            root = TreeNode(preorder[0])
            i = 0
            while inorder[i] != root.val:
                i += 1
            root.left = self.buildTree(preorder[1:i+1],inorder[0:i])
            root.right = self.buildTree(preorder[i+1:],inorder[i+1:])
            return(root)
