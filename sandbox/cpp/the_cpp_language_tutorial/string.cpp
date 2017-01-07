#include <iostream>
#include <string>
using namespace std;
int main() {
    string s ("x");
    cout << "initial:";
    cout << s << endl;
    cout << "now enter yours:\n";
    cin >> s;
    cout << "you enetered:\n";
    cout << s << endl;
    cout << "another input:\n";
    getline (cin, s);
    cout << "you enetered:\n";
    cout << s << endl; //you can only get newline here, why? because previously cin stops at whitespace
    return 0;
}
