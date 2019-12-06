//
// Created by guoy28 on 3/12/17.
//

#include "templateSpecialization.h"
#include<iostream>
void testTemplateSpecialization() {
    myContainer<int> x(7);
    myContainer<char> y('c');
    std::cout << "int increase" << std::endl;
    std::cout << x.increase() << std::endl;

    std::cout << "char upper and increase" << std::endl;
    std::cout << y.uppercase() << std::endl;
    //y.increase() is not available
    //myContainer with template and myContainer<char> are actually 2 different classes
    //they do NOT share methods with each other
}
