class Stack
{
    public:
	void Push(int value); //Push an integer, checking for overflow
	int top;		//index of the top of the stack
	int stack[10];		//the elements the stack
};

void
Stack::Push(int value)
{
    ASSERT(top<10);	//stack should never overflow
    stack[top++]=value;
}
