class Solution:
	def maxSubArray(self,nums):
		maxSoFar = 0
		n = len(nums)
		cumSum = [nums[0]]*n
		for i in range(1,n):
		    cumSum[i] = cumSum[i-1]+nums[i]
		for i in range(n):
		    for j in range(i,n):
			subSum = 0
			if i == 0:
			    subSum = cumSum[j]
			else:
	  		    subSum = cumSum[j] - cumSum[i-1]
			maxSoFar = max(subSum,maxSoFar)
		return(maxSoFar)
x = Solution()
print(x.maxSubArray(range(1000)))
