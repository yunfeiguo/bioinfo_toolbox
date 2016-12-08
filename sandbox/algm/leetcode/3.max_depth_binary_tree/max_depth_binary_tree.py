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
	    if root == None:
		return 0
	    elif root.left == None and root.right == None:
		return 1
	    else:
		return 1 + max(self.maxDepth(root.left),self.maxDepth(root.right))
