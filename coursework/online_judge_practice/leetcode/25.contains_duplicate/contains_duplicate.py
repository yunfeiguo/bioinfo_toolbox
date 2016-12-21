class Solution(object):
    def containsDuplicate(self, nums):
        """
        :type nums: List[int]
        :rtype: bool
        """
        if nums is None or len(nums) == 0:
            return(False)
        else:
            seen = {}
            for i in nums:
                if i in seen:
                    return(True)
                else:
                    seen[i] = True
            return(False)
x = Solution()
print(x.containsDuplicate([]))
print(x.containsDuplicate([1]))
print(x.containsDuplicate([1,3]))
print(x.containsDuplicate([1,3,1]))
