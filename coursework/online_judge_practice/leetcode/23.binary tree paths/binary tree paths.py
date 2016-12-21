# Definition for a binary tree node.
# class TreeNode(object):
#     def __init__(self, x):
#         self.val = x
#         self.left = None
#         self.right = None

class Solution(object):
    def binaryTreePaths(self, root):
        """
        :type root: TreeNode
        :rtype: List[str]
        """
        if root is None:
            return([])
        else:
            allPaths = self.binaryTreePaths(root.left)
            allPaths.extend(self.binaryTreePaths(root.right))
            for i in range(0,len(allPaths),1):
                allPaths[i] = str(root.val) + '->' + str(allPaths[i])
            if root.left is None and root.right is None:# this is a leaf node
                allPaths.append(str(root.val))
            return(allPaths)
class TreeNode(object):
    def __init__(self,x):
        self.val = x
        self.left = None
        self.right = None
x = TreeNode(1)
y = Solution()
print(y.binaryTreePaths(x))
x.left = TreeNode(2)
print(y.binaryTreePaths(x))
x.right = TreeNode(3)
print(y.binaryTreePaths(x))
x.left.left = TreeNode(4)
x.left.right = TreeNode(5)
print(y.binaryTreePaths(x))
x.right.left = TreeNode(9)
print(y.binaryTreePaths(x))
   
    
    
