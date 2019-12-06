//
// Created by guoy28 on 3/12/17.
//

#include "stdException.h"
#include <exception>
#include <iostream>

void testMyException() {
    myException exp;
    try {
        throw exp;
    } catch(std::exception& e) {
        std::cout << e.what() << std::endl;
    }
}