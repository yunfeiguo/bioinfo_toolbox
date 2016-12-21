# Definition for singly-linked list.
# class ListNode:
#     def __init__(self, x):
#         self.val = x
#         self.next = None

class Solution:
    # @param {ListNode} node
    # @return {void} Do not return anything, modify node in-place instead.
    def deleteNode(self, node):
	if node is None:
	    pass #do nothing if given nothing
	elif node.next is None:
	    pass #this is the tail, do nothing
	else:
	    #since this is singly linked list, we can't get the previous node
	    #but we can modify the current node because its reference is
	    #stored in previous node's next
	    #node = node.next #this is wrong!!, the original reference is gone
	    node.next = node.next.next
	    node.val = node.next.val
	    #python should now consider the old node.next as garbage, because no reference to it exists.
