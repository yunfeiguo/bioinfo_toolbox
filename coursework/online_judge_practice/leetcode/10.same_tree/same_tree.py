# Definition for a binary tree node.
# class TreeNode:
#     def __init__(self, x):
#         self.val = x
#         self.left = None
#         self.right = None

class Solution:
    # @param {TreeNode} p
    # @param {TreeNode} q
    # @return {boolean}
    def isSameTree(self, p, q):
	#use iterative rather than recursive method
        if p == None and q == None:
            return True
        elif (p != None and q == None) or (p == None and q != None):
            return False
        else:
	    q1 = [p]
	    q2 = [q]
	    result = True
	    while len(q1) >= 1 and len(q2) >= 1:
		#check if values are equal
		current1 = q1.pop()
		current2 = q2.pop()
		if current1 != None and current2 != None:
		    pass
		elif current1 == None and current2 == None:
		    continue
		else:
		    #only one of them is None, unequal depth
		    return False
		if current1.val != current2.val:
		    return False
		q1.append(current1.left)
		q1.append(current1.right)
		q2.append(current2.left)
		q2.append(current2.right)
	    if len(q1) == len(q2):
		return True
	    else:
		return False
