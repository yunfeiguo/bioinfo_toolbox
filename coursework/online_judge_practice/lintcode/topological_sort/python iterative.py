# Definition for a Directed graph node
# class DirectedGraphNode:
#     def __init__(self, x):
#         self.label = x
#         self.neighbors = []

class Solution:
    """
    @param graph: A list of Directed graph node
    @return: A list of integer
    """
    def topSort(self, graph):
        # write your code here
        from collections import deque
        incomingEdge = countIncomingEdge(graph)
        toVisit = deque(graph)
        ret = []
        while toVisit:
            cur = toVisit.popleft()
            if incomingEdge[cur] == 0:
                #ret.append(cur.label)
                ret.append(cur)
                for i in cur.neighbors:
                    incomingEdge[i] -= 1
            else:
                toVisit.append(cur)
        return(ret)
def countIncomingEdge(graph):
    ret = {}
    for cur in graph:
        if not cur in ret:
            ret[cur] = 0
        for i in cur.neighbors:
            if i in ret:
                ret[i] += 1
            else:
                ret[i] = 1
    return(ret)
        
