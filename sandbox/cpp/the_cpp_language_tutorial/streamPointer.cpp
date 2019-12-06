//
// Created by guoy28 on 4/12/17.
//

#include "streamPointer.h"
#include <fstream>
#include <iostream>

void streamPointerExample() {
    std::fstream myfile;
    myfile.open(__FILE__);
    if (myfile.is_open()) {
        long begin = myfile.tellg();
        myfile.seekg(0, std::ios::end);
        long end = myfile.tellg();
        myfile.close();
        std::cout << "file size: " << end - begin << std::endl;
    } else {
        std::cout << "failed to open" << std::endl;
    }
}