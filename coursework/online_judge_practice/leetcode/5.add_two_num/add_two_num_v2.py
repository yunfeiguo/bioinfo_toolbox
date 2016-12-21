# Definition for singly-linked list.
# class ListNode:
#     def __init__(self, x):
#         self.val = x
#         self.next = None

class Solution:
    # @param {ListNode} l1
    # @param {ListNode} l2
    # @return {ListNode}
    def addTwoNumbers(self, l1, l2):
        addOne = False
        result = None
        return_obj = None
        while l1 != None or l2 != None:
            increase = 0
            if l1 != None:
                increase += l1.val
                l1 = l1.next
            if l2 != None:
                increase += l2.val
                l2 = l2.next
            if addOne:
                increase += 1
            if increase >= 10:
                addOne = True
                increase -= 10
            else:
                addOne = False
            if result == None:
                result = ListNode(increase)
                returnObj = result
            else:
                result.next = ListNode(increase)
                result = result.next
        if addOne:
            if result == None:
                result = ListNode(1)
                returnObj = result
            else:
                result.next = ListNode(1)
        return returnObj
