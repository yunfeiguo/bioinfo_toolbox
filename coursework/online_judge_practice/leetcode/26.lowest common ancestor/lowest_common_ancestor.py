# Definition for a binary tree node.
class TreeNode(object):
    def __init__(self, x):
        self.val = x
        self.left = None
        self.right = None
class Solution(object):
    def lowestCommonAncestor(self, root, p, q):
        """
        :type root: TreeNode
        :type p: TreeNode
        :type q: TreeNode
        :rtype: TreeNode
        """
        #assume p and q are in the tree
	pLeft = False
        pRight = False
        qLeft = False
        qRight = False
        if root is None or p is None or q is None:
            return(None)
        elif root == p or root == q:
            return(root)
        else:
            pLeft = self.isDescendent(root.left,p)
            qRight = self.isDescendent(root.right,q)
            if pLeft and qRight:
                return(root)
            pRight = self.isDescendent(root.right,p)
            qLeft = self.isDescendent(root.left,q)
            if pRight and qLeft:
                return(root)
            if pRight and qRight:
                return(self.lowestCommonAncestor(root.right,p,q))
            if pLeft and qLeft:
                return(self.lowestCommonAncestor(root.left,p,q))
            #p or q not in the tree
            return(None)
    def isDescendent(self,root,x):
        """
        decide whether x is descendent of root (including itself)
        """
        if root is None or x is None:
            return(False)
        elif root == x:
            return(True)
        else:
            #it's impossible to have left and right both as ancestors
            return(self.isDescendent(root.left,x) ^ self.isDescendent(root.right,x))
x = TreeNode(1)
x.left = TreeNode(2)
x.right = TreeNode(3)
sol = Solution()
print(sol.isDescendent(x,x))
print(sol.isDescendent(x,x.left))
print(sol.isDescendent(x,x.right))
print(sol.lowestCommonAncestor(x,x,x).val)
print(sol.lowestCommonAncestor(x,x,x.left).val)
print(sol.lowestCommonAncestor(x,x,x.right).val)   
print(sol.lowestCommonAncestor(x,None,x))
print(sol.lowestCommonAncestor(x,TreeNode(10),x.left))
print(sol.lowestCommonAncestor(x,x.left,x.right).val)
print(sol.lowestCommonAncestor(x,x.left,x.left).val)
