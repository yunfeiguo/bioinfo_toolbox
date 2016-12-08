class Solution(object):
    def moveZeroes(self, nums):
        """
        :type nums: List[int]
        :rtype: void Do not return anything, modify nums in-place instead.
        """
        #so the idea is that scan the entire array
        #every time we see a zero
        #we swap it with next non-zero
        #until the end of array
        if nums is None or len(nums) == 0:
            pass
        else:
            for i in range(0,len(nums)):
                if nums[i] == 0:
                    if i == len(nums)-1:
                        break #stop if we are the end
                    j = 0
                    for j in range(i+1,len(nums)):
                        if nums[j] != 0:
                            break
                    if nums[j] == 0:
                        break #we have reached the end
                              #yet there is no more non-zero
                    else:
                        self.swap(i,j,nums)
                    
    def swap(self,i,j,nums):
        tmp = nums[i]
        nums[i] = nums[j]
        nums[j] = tmp
x = Solution()
testCase = [[],None,[0],[1],[0,1],[1,0],[0,1,0],[0,0,1],[1,0,0]]
for y in testCase:
    x.moveZeroes(y)
print(testCase)    
