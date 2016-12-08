#!/usr/bin/env python
class Solution:
    # @param {integer[]} nums
    # @return {integer}
    def singleNumber(self, nums):
        if len(nums) == 0:
            return None
        else:
            #easiest solution
            #just use hash
            count={}
            for i in nums:
                if i in count:
                    count[i] += 1
                else:
                    count[i] = 1
            for i in count:
		#all elements appear twice except for one
		#so we return the element doesn't appear twice
                if count[i] != 2:
                    return i
