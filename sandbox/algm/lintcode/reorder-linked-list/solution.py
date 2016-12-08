"""
Definition of ListNode
class ListNode(object):

    def __init__(self, val, next=None):
        self.val = val
        self.next = next
"""
class Solution:
    """
    @param head: The first node of the linked list.
    @return: nothing
    """
    def reorderList(self, head):
        # write your code here
        #break the list into two halves
        #reverse the second half
        #weaven two lists
        n = getLen(head)
        if n <= 2:
            return(head)
        newHead = head
        for i in range((n-1)/2):
            newHead = newHead.next
        oldTail = newHead
        newHead = newHead.next        
        oldTail.next = None
        newHead = rev(newHead)
        head = join(head,newHead)
        return(head)
def getLen(head):
    n = 0
    while head:
        n += 1
        head = head.next
    return(n)
def rev(head):
    #reverse a singly linked list
    dummy = ListNode(0)
    cur = head
    dummy.next = cur
    if head is None or head.next is None:
        return(head)
    cur = cur.next
    while cur:
        tmp = dummy.next
        dummy.next = cur
	#this must be done first!
        cur = cur.next
        dummy.next.next = tmp
    head.next = None
    return(dummy.next)
def join(h1,h2):
    dummy = ListNode(0)
    dummyHead = dummy
    while h1 or h2:
        if h1:
            dummyHead.next = h1
            h1 = h1.next
            dummyHead = dummyHead.next
        if h2:
            dummyHead.next = h2
            h2 = h2.next
            dummyHead = dummyHead.next
    return(dummy.next)
