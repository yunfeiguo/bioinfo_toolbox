 #Definition for singly-linked list.
class ListNode:
     def __init__(self, x):
         self.val = x
         self.next = None
     def printList(self):
	 nxt = self
	 while nxt != None:
	     print nxt.val
	     nxt = nxt.next

class Solution:
    # @param {ListNode} head
    # @return {ListNode}
    def reverseList(self, head):
        if head == None:
            #empty list
            return head
        elif head.next == None:
            #one element list, return head
            return head
        else:
            nxt = head.next
            prev = head
            #traverse entire list, reverse link one by one
            while nxt != None:
                tmp = nxt.next
                nxt.next = prev
                prev = nxt
                nxt = tmp
	    head.next = None
            return prev
one = ListNode(1)
two = ListNode(2)
one.next = two
one.printList()
sol = Solution()
revOne = sol.reverseList(one)
revOne.printList()
