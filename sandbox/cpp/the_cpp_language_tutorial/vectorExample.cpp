//
// Created by guoy28 on 2/28/17.
//

#include "vectorExample.h"
#include <vector>
#include <iostream>

void testVector() {
    size_t size = 10;
    int_vec_t v(size);//intialize to 0 or new objects
    std::vector<char> c;
    c.reserve(3); //no initilization
    int staticArray[size];
    int *dynamicArray = new int[size];

    for (int i = 0; i < size; i++) {
        staticArray[i] = i;
        dynamicArray[i] = i;
        v[i] = i;
        c.push_back(i + 'a');
    }
    std::cout << "check arrays" << std::endl;
    for (int i = 0; i < size; i++) {
        std::cout << "static: " << staticArray[i] << std::endl;
        std::cout << "dynamic: " << dynamicArray[i] << std::endl;
        std::cout << "vector: " << v[i] << std::endl;
        std::cout << "vector: " << c[i] << std::endl;
    }
    std::cout << v.size() << " elements remaining" << std::endl;

    delete [] dynamicArray;
    try {
        v.at(size + 1);
    } catch(std::out_of_range o) {
        std::cout << "exception: " << o.what() << std::endl;
    }
}