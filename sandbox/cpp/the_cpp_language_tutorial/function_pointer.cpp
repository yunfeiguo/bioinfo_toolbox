#include<iostream>
using namespace std;
int addition(int x, int y) {
    return x + y;
}
int subtract(int x, int y) {
    return x - y;
}
int operation(int x, int y, int (*function)(int, int)) {
	return (*function)(x, y);
}
int functionPointer() {
    cout << "1 + 2 = " << operation(1, 2, addition) << endl;
    cout << "1 - 2 = " << operation(1, 2, subtract) << endl;
}
