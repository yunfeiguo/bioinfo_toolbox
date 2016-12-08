"""
Definition of TreeNode:
class TreeNode:
    def __init__(self, val):
        self.val = val
        self.left, self.right = None, None
"""


class Solution:
    """
    @param root: The root of binary tree.
    @return: Level order in a list of lists of integers
    """
    def levelOrder(self, root):
        # write your code here
        #try BFS first
        if root is None:
            return([])
        from Queue import Queue
        q = Queue()
        q.put(root)
        ret = []
        while not q.empty():
            curLvl = []
            n = q.qsize()
            for i in range(n):
                e = q.get()
                curLvl.append(e.val)
                if e.left:
                    q.put(e.left)
                if e.right:
                    q.put(e.right)
            ret.append(curLvl)
        return(ret)
                
            

