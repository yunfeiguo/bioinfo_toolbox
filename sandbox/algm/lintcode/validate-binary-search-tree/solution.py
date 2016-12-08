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
        if root is None:
            return(True)
        else:
            ret = True
            if root.left != None:
                ret &= (root.val > root.left.val) & self.isValidBST(root.left)
                #check if all nodes on the right > root
                ret &= self.allCompare(root.left,root.val,'<')
            if root.right != None:
                ret &= (root.val < root.right.val) & self.isValidBST(root.right)
                #check if all nodes on the left < root
                ret &= self.allCompare(root.right,root.val,'>')
            return(ret)
    def allCompare(self,root,val,op):
        queue = [root]
        ret = True
        while queue:
            cur = queue.pop()
            if op == '>':
                ret &= cur.val > val
            elif op == '<':
                ret &= cur.val < val
            if cur.left:
                queue.append(cur.left)
            if cur.right:
                queue.append(cur.right)
            if not ret:
                break
        return(ret)
            
        
