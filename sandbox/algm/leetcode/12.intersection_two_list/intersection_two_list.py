# Definition for singly-linked list.
# class ListNode:
#     def __init__(self, x):
#         self.val = x
#         self.next = None
class Solution:
    # @param two ListNodes
    # @return the intersected ListNode
    def getIntersectionNode(self, headA, headB):
        if headA is None or headB is None:
            return None
        else:
            #two lists intersect at the end
            #if they don't have equal length, then, the dif = abs(len1-len2)
            #intersection doesn't appear in first dif nodes
            dif = 0
            if self.count(headA) > self.count(headB):
                dif = self.count(headA)-self.count(headB)
                while dif > 0:
                    headA = headA.next
                    dif -= 1
            elif self.count(headA) < self.count(headB):
                dif = self.count(headB) - self.count(headA)
                while dif > 0:
                    headB = headB.next
                    dif -= 1
            nodeA = headA
            nodeB = headB
            while nodeA != nodeB:
                nodeA = nodeA.next
                nodeB = nodeB.next
            return nodeA
    def count(head):
        '''count number of nodes'''
        current = head
        n = 0
        while current is not None:
            n += 1
            current = current.next
        return n
