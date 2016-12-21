class MinStack:
    # initialize your data structure here.
    def __init__(self):
        self.stack = []
        self.min = None

    # @param x, an integer
    # @return an integer
    def push(self, x):
        x = int(x)
        self.stack.append(x)
        if self.min == None:
            self.min = x
        elif x < self.min:
            self.min = x
        return self.stack[len(self.stack)-1]
    # @return nothing
    def pop(self):
        if len(self.stack) == 0:
            pass
        else:
            x = self.stack.pop()
            if len(self.stack) == 0:
                self.min = None
            elif x == self.min:
                self.min = min(self.stack)

    # @return an integer
    def top(self):
        if len(self.stack) == 0:
            pass
        else:
            return self.stack[len(self.stack)-1]

    # @return an integer
    def getMin(self):
        return self.min
x = MinStack()    
print "push",x.push(2147483646)
print "push",x.push(2147483646)
print "push",x.push(2147483647)
print "top",x.top()
print "pop",x.pop()
print "min",x.getMin()
print "pop",x.pop()
print "push",x.push(2147483647)
print "top",x.top()
print "min",x.getMin()
print "push",x.push(-2147483648)
print x.top()
print x.getMin()
print x.pop()
print x.getMin()
