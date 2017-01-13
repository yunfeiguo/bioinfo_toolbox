#include<iostream>
using namespace std;
int myPlus(int a, int b) {
    return a + b;
}
int myPlus(int a, int b, int c) {
    return a + b + c;
}
int overloading() {
    cout << "two numbers " << myPlus(1,2) << endl;
    cout << "three numbers " << myPlus(1,2,3) << endl;

}
