#include<iostream>
using namespace std;
void doubleNumber(int& a) {
    a *= 2;
}
int passByRef() {
    int x (3);
    cout << "Before " << x << endl;
    doubleNumber(x);
    cout << "After " << x << endl;
    return 0;
}
