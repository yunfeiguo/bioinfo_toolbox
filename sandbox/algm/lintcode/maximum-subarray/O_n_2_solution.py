class Solution:
	def maxSubArray(self,nums):
		maxSoFar = 0
		n = len(nums)
		for i in range(n):
		    subSum = 0
		    for j in range(i,n):
			subSum += nums[j]
			maxSoFar = max(subSum,maxSoFar)
		return(maxSoFar)
x = Solution()
print(x.maxSubArray(range(1000)))
