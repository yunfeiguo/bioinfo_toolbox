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
    @param node: insert this node into the binary search tree.
    @return: The root of the new binary search tree.
    """
    def insertNode(self, root, node):
        # write your code here
        if root is None:
            return(node)
        else:
            anchor = root
            parent = root
            while anchor != None:
                #loop will stop at None
                parent = anchor
                assert(anchor.val != node.val)
                if anchor.val < node.val:
                    anchor = anchor.right
                else:
                    anchor = anchor.left
            if parent.val < node.val:
                parent.right = node
            elif parent.val >= node.val:
                parent.left = node
            return(root)

