#include<iostream>
using namespace std;
int dynamicMemory() {
    int **p;
    p = new int*[2000000000];
    p = new (nothrow) int*[1];
    if (p == 0) {
	cout << "p is 0" << endl;
    }
    delete [] p;
}
