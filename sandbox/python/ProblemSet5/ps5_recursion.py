#!/usr/bin/env python
# 6.00x Problem Set 5
#
# Part 2 - RECURSION

#
# Problem 3: Recursive String Reversal
#
def reverseString(aStr):
    """
    Given a string, recursively returns a reversed copy of the string.
    For example, if the string is 'abc', the function returns 'cba'.
    The only string operations you are allowed to use are indexing,
    slicing, and concatenation.
    
    aStr: a string
    returns: a reversed string
    """
    if len(aStr) == 1:
	return aStr
    else:
	return aStr[-1]+reverseString(aStr[:-1])

#
# Problem 4: X-ian
#
def x_ian(x, word):
    """
    Given a string x, returns True if all the letters in x are
    contained in word in the same order as they appear in x.

    >>> x_ian('eric', 'meritocracy')
    True
    >>> x_ian('eric', 'cerium')
    False
    >>> x_ian('john', 'mahjong')
    False
    
    x: a string
    word: a string
    returns: True if word is x_ian, False otherwise
    """
    if len(x)==0:
	return True
    if len(word) < len(x):
	return False
    else:
	if word[0] == x[0]:
	    return x_ian(x[1:],word[1:])
	else:
	    return x_ian(x,word[1:])

#
# Problem 5: Typewriter
#
def insertNewlines(text, lineLength):
    """
    Given text and a desired line length, wrap the text as a typewriter would.
    Insert a newline character ("\n") after each word that reaches or exceeds
    the desired line length.

    text: a string containing the text to wrap.
    line_length: the number of characters to include on a line before wrapping
        the next word.
    returns: a string, with newline characters inserted appropriately. 
    """
    if len(text)<=lineLength:
	return text
    else:
	i=findSpace(text,lineLength)
	i+=1
	result=text[:i]+"\n"
	return result+insertNewlines(text[i:],lineLength)

def findSpace(text,n):
    """
    Given text and an integer n, find a space with index not smaller than n-1

    text: a string
    n: an int
    return: index for the space, not smaller than n-1
    """
    i=n-1
    if text[i]== ' ':
	return i
    else:
	return findSpace(text,n+1)
if __name__ == '__main__':
    print insertNewlines('While I expect new intellectual adventures ahead, nothing will compare to the exhilaration of the world-changing accomplishments that we produced together.', 15)
    print insertNewlines('Nuh-uh! We let users vote on comments and display them by number of votes. Everyone knows that makes it impossible for a few persistent voices to dominate the discussion.', 20)
