//
// Created by guoy28 on 4/12/17.
//

#include "fileio.h"
#include <iostream>
#include <fstream>

void fileioExample() {
    std::ofstream myfile;
    myfile.open("/tmp/test");
    myfile << "testtesttest" << std::endl;
    myfile.close();
    std::cout << "written to /tmp/test" << std::endl;
}