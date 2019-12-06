void
Stack::Push(int value)
{
    ASSERT(top<10);	//stack should never overflow
    stack[top++]=value;
}
