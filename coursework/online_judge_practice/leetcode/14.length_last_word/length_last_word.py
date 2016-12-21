class Solution:
    # @param {string} s
    # @return {integer}
    def lengthOfLastWord(self, s):
        if s is None or len(s) == 0:
            return(0)
        elif len(s) == 1:
            if s.isspace():
                return(0)
            else:
                return(1)
        else:
            idx = range(0,len(s))
            idx.reverse()
            count = 0
            trailing = True
            for i in idx:
                c = str(s[i])
                if c.isspace():
                    print(str(i)+"th char is space")
                    if trailing: #trailing space does not influence our counting
                        continue
                    else:
                        break
                else:
                    if trailing:
                        trailing = False
                    count += 1
            return(count)
