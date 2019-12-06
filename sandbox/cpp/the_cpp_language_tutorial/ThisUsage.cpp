//
// Created by guoy28 on 1/30/17.
//
#include <iostream>
using namespace std;
#include "ThisUsage.h"
//ThisUsage ThisUsage::operator= (const ThisUsage& param) {
ThisUsage& ThisUsage::operator= (const ThisUsage& param) {
    x = param.x;
    y = param.y;
    return *this; //note here we return *this, not this itself (which is a pointer), why? because return type is ThisUsage&
                    //needs an actual object
}

ThisUsage::ThisUsage(int x, int y) {
    this->x = x;
    this->y = y;
}


void ThisUsage::print() {
    cout << "ThisUsage class" << endl;
    cout << x << endl;
    cout << y << endl;
}

void ThisUsage::setX(int newX) {
    x = newX;
}