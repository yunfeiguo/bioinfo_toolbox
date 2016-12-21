#!/usr/bin/env python
# Definition for a binary tree node.
# class TreeNode:
#     def __init__(self, x):
#         self.val = x
#         self.left = None
#         self.right = None
class Solution:
    # @param {TreeNode} root
    # @return {integer}
	def maxDepth(self, root):
	    toVisit = []
	    current = []
	    depth = 0
	    if root == None:
		#empty tree
		return depth
	    else:
		toVisit.append(root)
	    while toVisit:
		current = toVisit
		toVisit = []
		while current:
			target = current.pop() #get last element
			if target.left != None:
			    toVisit.append(target.left)
			if target.right != None:
			    toVisit.append(target.right)
		depth += 1
	    return depth
