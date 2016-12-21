class Solution:
	def maxSubArray(self,nums):
		maxSoFar = 0
		n = len(nums)
		for i in range(n):
		    for j in range(n):
			subSum = sum(nums[i:j+1])
			maxSoFar = max(subSum,maxSoFar)
		return(maxSoFar)
x = Solution()
print(x.maxSubArray(range(1000)))
