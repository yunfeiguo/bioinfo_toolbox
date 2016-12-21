class Solution:
    """
    @param nums: A list of integers
    @return: An integer denote the sum of maximum subarray
    """
    def maxSubArray(self, nums):
        # write your code here
        maxSoFar = nums[0]
        maxEndingHere = nums[0]
        for i in range(1,len(nums)):
            maxEndingHere = max(nums[i],maxEndingHere+nums[i])
            maxSoFar = max(maxEndingHere,maxSoFar)
        return(maxSoFar)
x = Solution()
print(x.maxSubArray(range(1000)))    
