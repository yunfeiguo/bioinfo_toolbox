#!/usr/bin/env python
class Solution:
    # @param {integer[]} nums
    # @return {integer}
    def singleNumber(self, nums):
        if len(nums) == 0:
            return None
        else:
	    '''
            some numbers appear twice, the one number appears once
            bitwise, 1s can appear either 0 time or 2 times for the paired numbers
            so we can use one bit to track the number of times 1 appears
            0(0 time)->1(1 time)->0(2 times)
            the the truth table is
            count   input   output(new count)
            0       0       0
            0       1       1
            1       1       0
            1       0       1
            to get the above truth table, we only need one operation, the XOR!
            The magic of XOR is, when one operand is 0, the result is same as the other operand
            so when 1 either appears 0 time or 2 times, the corresponding bit is 0, and 0 XOR the
            single number, will give us the single number.
            '''
	    ones = 0
	    for i in nums:
		ones = ones ^ i
	    return ones
test = Solution()
print "[1]"
print test.singleNumber([1])
print "[1,2,1]"
print test.singleNumber([1,2,1])
print "[1,2,1,3,5,5,3]"
print test.singleNumber([1,2,1,3,5,5,3])
