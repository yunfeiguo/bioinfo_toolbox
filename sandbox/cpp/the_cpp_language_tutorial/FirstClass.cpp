//
// Created by guoy28 on 1/18/17.
//

#include "FirstClass.h"
#include<iostream>
using namespace std;

SecondClass::SecondClass(int a, int b) {
    x = a;
    y = b;
}

int SecondClass::area() {
    return x*y;
}

SecondClass::~SecondClass() {
    cout << "deleting..." << endl;
}