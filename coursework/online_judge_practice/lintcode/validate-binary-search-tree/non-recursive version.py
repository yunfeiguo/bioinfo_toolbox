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
    @return: True if the binary tree is BST, or false
    """  
    def isValidBST(self, root):
        # write your code here
        #use in-order traversal
        #gives serialized BST in ascending order
        toVisit = [root]
        s = [] #s stands for serialization
        seen = {}
        while toVisit:
            cur = toVisit.pop()
            if cur is None:
                continue
            if cur in seen and seen[cur] == True:
                s.append(cur.val)
            else:
                toVisit.append(cur.right)
                toVisit.append(cur)
                toVisit.append(cur.left)
                seen[cur] = True
            if len(s) > 1 and s[-1] <= s[-2]:
                return(False)
        return(True)
