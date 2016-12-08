class ShortInputException(Exception):
    def __init__(self,length,atleast):
        Exception.__init__(self)
        self.length = length
        self.atleast = atleast
try:
    text = raw_input('enter:')
    if len(text)<3:
        raise ShortInputException(len(text),3)
except EOFError:
    print('EOF')
except KeyboardInterrupt:
    print('keyboard')
except ShortInputException as se:
    print('Expect at least {}'.format(se.atleast))
else:
    print('ok')
