#!/usr/bin/env python

def twoSum( nums, target):
    return_list = []
    dict = {}
    for i in range(len(nums)):
	n = nums[i]
	if n in dict:
	    dict[n].append(i)
	else:
	    dict[n] = [i]
    for i in range(len(nums)):
	n = nums[i]
	remain = target - n
	if remain in dict:
	    if n == remain:
		if len(dict[n]) < 2:
		    continue
	    return_list.append(dict[n].pop())
	    return_list.append(dict[remain].pop())
	    break
    return_list.sort();
    for i in range(len(return_list)):
	return_list[i] += 1
    return return_list	
print 'twoSum([2,7,11,20],31)'
print twoSum([2,7,11,20],31)		
print 'twoSum([2,7,11,20],51)'
print twoSum([2,7,11,20],51)
print 'twoSum([],31)'
print twoSum([],31)		
