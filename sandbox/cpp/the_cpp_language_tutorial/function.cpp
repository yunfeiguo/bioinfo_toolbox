#include<iostream>
using namespace std;
int addition(int a, int b) {
    return a + b;
}
int function() {
    int z (addition(3,2));
    cout << "The result is " << z << endl;
    return 0;
}
