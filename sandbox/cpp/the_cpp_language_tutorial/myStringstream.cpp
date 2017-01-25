#include <iostream>
#include <string>
#include <sstream>
using namespace std;
int myStringstream() {
    int n;
    float x;
    string s;
    cout << "integer:\n";
    getline(cin, s);
    stringstream(s) >> n;
    cout << "float: \n";
    getline(cin, s);
    stringstream(s) >> x;

    cout << "result " << x << " x " << n << endl;
    return 0;
}
