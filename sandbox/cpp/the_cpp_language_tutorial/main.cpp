//
// Created by guoy28 on 1/13/17.
//
#include "pointer.h"
#include "hello.h"
#include "array.h"
#include "dynamic_memory.h"
#include "myStringstream.h"
#include "myString.h"
#include "overloading.h"
#include "FirstClass.h"
#include<iostream>
using namespace std;

int main() {
    if (0) {
        hello();
        myArray();
        dynamicMemory();
        pointer();
        myStringstream();
        myString();
        overloading();
        FirstClass x;
        x.setValues(2,3);
        cout << "area: " << x.area() << endl;
    }
    SecondClass *y = new SecondClass(3,4);
    cout << "area: " << y->area() << endl;
    delete y;
    SecondClass z(3,4);
    cout << "area: " << z.area() << endl;
//desctructor will be automatically called once z is out of scope.
}

