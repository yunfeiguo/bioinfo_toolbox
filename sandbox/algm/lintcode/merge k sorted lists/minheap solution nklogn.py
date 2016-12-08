"""
Definition of ListNode
class ListNode(object):

    def __init__(self, val, next=None):
        self.val = val
        self.next = next
"""
class Solution:
    """
    @param lists: a list of ListNode
    @return: The head of one sorted list.
    """
    def mergeKLists(self, lists):
        # write your code here
        #use minheap
        import heapq
        dummy = ListNode(0)
        head = dummy
        h = []
        for i in lists:
            if i:
                heapq.heappush(h,(i.val,i))
        while h:
            cur = heapq.heappop(h)
            head.next = cur[1]
            if cur[1].next:
                heapq.heappush(h,(cur[1].next.val,cur[1].next))
            head = head.next
        return(dummy.next)


    
    
        


