class Solution:
	def maxSubArray(self,nums):
		n = len(nums)
		if n == 0:
		    return(0)
		elif n == 1:
		    if nums[0] > 0:
			return(nums[0])
		    else:
			return(0)
		else:
		    max1 = self.maxSubArray(nums[:n/2])
		    max2 = self.maxSubArray(nums[n/2:])
		    maxCenter = 0 #maximum of subarray that crosses center
		    cumSumLeft = [nums[n/2-1]]*(n/2)
		    cumSumRight = [nums[n/2]]*(n-n/2)
		    for i in range(n/2-2,-1,-1):
			cumSumLeft[i] = cumSumLeft[i+1] + nums[i]
		    maxLeft = max(cumSumLeft)
		    for j in range(n/2+1,n):
			cumSumRight[j-n/2] = cumSumRight[j-1-n/2] + nums[j]
		    maxRight = max(cumSumRight)
		    maxCenter = maxLeft+maxRight
		    return(max(max1,max2,maxCenter))
x = Solution()
print(x.maxSubArray(range(1000)))
