//
// Created by guoy28 on 1/13/17.
//
#include "OperatorOverloading.h"
#include "pointer.h"
#include "hello.h"
#include "array.h"
#include "dynamic_memory.h"
#include "myStringstream.h"
#include "myString.h"
#include "overloading.h"
#include "FirstClass.h"
#include "ThisUsage.h"
#include "ClassStaticVariable.h"
#include "FriendExample.h"
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

        SecondClass *y = new SecondClass(3,4);
        cout << "area: " << y->area() << endl;
        delete y;
        SecondClass z(3,4);
        cout << "area: " << z.area() << endl;
//desctructor will be automatically called once z is out of scope.
        OperatorOverloading one(1,2);
        OperatorOverloading two(2,3);
        OperatorOverloading three;
        three = one + two;
        cout << "add one and two using overloaded +" << endl;
        cout << three.x << " y: " << three.y << endl;

        ThisUsage tu(1,2);
        ThisUsage tu2 = tu;
        cout << "tu" << endl;
        tu.print();
        tu.setX(10);
        cout << "tu2 after chaning tu" << endl;
        tu2.print();
        ClassStaticVariable csv1 (ClassStaticVariable(1));
        ClassStaticVariable csv2[5] = {{1}, {2}, {3}, {4}, {5}};
        ClassStaticVariable *csv3 = new ClassStaticVariable(2);
        cout << "count after 2 elements and 1 array of length 5" << endl;
        cout << csv1.n << endl;
        delete csv3;
        cout << "count after del" << endl;
        cout << ClassStaticVariable::n << endl;
    }
    FriendExample friendExample(2,3);
    FriendExample friendExample2 (duplicate(friendExample));
    cout << "friend area" << endl;
    cout << friendExample2.area() << endl;
}


