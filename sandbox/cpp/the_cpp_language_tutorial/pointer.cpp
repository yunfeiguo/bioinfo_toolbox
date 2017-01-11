#include<iostream>
using namespace std;
void increase(void *p, int size) {
    if (size == sizeof(char)) {
	(*(char *)p)++;
    } else if (size == sizeof(int)){
	(*(int *)p)++;
    }
}
int main() {
    char c ('x');
    int n (5);
    cout << "before\n";
    cout << "c " << c << " n " << n << endl;
    increase(&c, sizeof(c));
    increase(&n, sizeof(n));
    cout << "after\n";
    cout << "c " << c << " n " << n << endl;
}
