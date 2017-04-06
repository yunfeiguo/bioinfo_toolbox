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
#include "FriendClassExample.h"
#include "InheritanceExample.h"
#include "AbstractBaseClassExample.h"
#include "TemplateExample.h"
#include "vectorExample.h"
#include "templateSpecialization.h"
#include "stdException.h"
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <fcntl.h>
#include <exception>

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
        std::cout << "area: " << x.area() << std::endl;

        SecondClass *y = new SecondClass(3,4);
        std::cout << "area: " << y->area() << std::endl;
        delete y;
        SecondClass z(3,4);
        std::cout << "area: " << z.area() << std::endl;
//desctructor will be automatically called once z is out of scope.
        OperatorOverloading one(1,2);
        OperatorOverloading two(2,3);
        OperatorOverloading three;
        three = one + two;
        std::cout << "add one and two using overloaded +" << std::endl;
        std::cout << three.x << " y: " << three.y << std::endl;

        ThisUsage tu(1,2);
        ThisUsage tu2 = tu;
        std::cout << "tu" << std::endl;
        tu.print();
        tu.setX(10);
        std::cout << "tu2 after chaning tu" << std::endl;
        tu2.print();
        ClassStaticVariable csv1 (ClassStaticVariable(1));
        ClassStaticVariable csv2[5] = {{1}, {2}, {3}, {4}, {5}};
        ClassStaticVariable *csv3 = new ClassStaticVariable(2);
        std::cout << "count after 2 elements and 1 array of length 5" << std::endl;
        std::cout << csv1.n << std::endl;
        delete csv3;
        std::cout << "count after del" << std::endl;
        std::cout << ClassStaticVariable::n << std::endl;
        FriendExample friendExample(2,3);
        FriendExample friendExample2 (duplicate(friendExample));
        std::cout << "friend area" << std::endl;
        std::cout << friendExample2.area() << std::endl;
        MySquare mySquare;
        mySquare.setL(5);
        MyRectangle myRectangle;
        myRectangle.convert(mySquare);
        std::cout << "myRectangle after conversion" << std::endl;
        std::cout << myRectangle.area() << std::endl;

        Process p(1);
        LinuxProcess lp(2);
        std::cout << "regular process" << std::endl;
        std::cout << p.processGetPid() << std::endl;
        std::cout << "linux process" << std::endl;
        std::cout << lp.linuxProcessGetPid() << std::endl; //cannot call processGetPid() directly because it is a protected member
        Bird bird(1);
        Insect insect(2);
        Animal *a1 = &bird;
        Animal *a2 = &insect;
        a1->move();
        a2->move();
        std::cout << "max of 1 and 2 " << getMax<int>(1,2) << std::endl;
        std::cout << "min of 1 and 2L " << getMin<int,long>(1,2L) << std::endl;
        testVector();
        testTemplateSpecialization();
        testMyException();
    }
    seqan::CharString id;
    seqan::Dna5String seq;
    seqan::CharString seqFileName = "/Users/guoy28/Downloads/bioinfo_toolbox/sandbox/cpp/the_cpp_language_tutorial/resources/example.fa";

    try {
        seqan::SeqFileIn seqFileIn(toCString(seqFileName));
        readRecord(id, seq, seqFileIn);
        std::cout << id << '\t' << seq << '\n';
    } catch (std::exception const & e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
    }
}
