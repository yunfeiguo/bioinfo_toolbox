#include<iostream>
using namespace std;
int array() {
    int num1 [5];
    int num2 [] = {1,2,3};
    //uninitialized values
    for (int i = 0; i < sizeof(num1)/sizeof(num1[0]); i++) {
	cout << "num1 " << i << ":" << num1[i] << endl;
    }
    for (int i = 0; i < sizeof(num2)/sizeof(num2[0]); i++) {
	cout << "num2 " << i << ":" << num2[i] << endl;
    }
    //segmentaiton fault
    //cout << "num1 " << 10000 << "(out of bound):" << num1[10000] << endl;
    return 0;
}
