#include<iostream>
#include<new>
using namespace std;
int dynamicMemory() {
    int **p;
    p = new int*[2000000000];
    p = new (nothrow) int*[1];
    if (p == 0) {
	cout << "p is 0" << endl;
    }
    delete [] p;

    int n,i;
    int *numbers;
    cout << "how many numbers you want?" << endl;
    cin >> n;
    numbers = new (nothrow) int[n];
    if (numbers == 0) {
        cout << "allocation fail" << endl;
        return 1;
    }
    for (i = 0; i < n; i++) {
	cout << "enter the " << i + 1 << "th number:" << endl;
	cin >> numbers[i];
    }
    cout << "you entered: " << endl;
    for (i = 0; i < n; i++) {
	cout << numbers[i] << " ";
    }
    return 0;
}
