#!/usr/bin/env python
# Definition for a binary tree node.
from collections import deque
class TreeNode:
    def __init__(self, x):
        self.val = x
        self.left = None
        self.right = None
    def printTree(self):
	nextQ = deque([self])
	while nextQ:
	    q = deque(nextQ)
	    nextQ.clear()
	    while q:
		current = q.popleft()
	        if current == None:
		    continue
	        else:
		    if current.left != None:
		        nextQ.append(current.left)
		    if current.right != None:
		        nextQ.append(current.right)
		    print current.val," ",
            print "\n"
class Solution:
    # @param {TreeNode} root
    # @return {TreeNode}
    def invertTree(self, root):
        toVisit = [root]
        #for each node that has at least one child,
        #switch left and right
        while toVisit:
            currentNode = toVisit.pop()
            if currentNode == None:
                continue
            else:
                if currentNode.left == None and currentNode.right == None :
                    continue #do nothing
                else:
                    tmp = currentNode.left
                    currentNode.left = currentNode.right
                    currentNode.right = tmp
                    toVisit.append(currentNode.left)
                    toVisit.append(currentNode.right)
        return(root)            
root = TreeNode(1)
root.left = TreeNode(2)
root.right = TreeNode(3)
root.left.left = TreeNode(4)
root.left.right = TreeNode(5)
root.right.left = TreeNode(6)
root.right.right = TreeNode(7)
sol = Solution()
print "empty tree inverted"
print sol.invertTree(None)
print "rootonly tree inverted"
sol.invertTree(TreeNode(1)).printTree()
root.printTree()
print "inverted"
sol.invertTree(root).printTree()
