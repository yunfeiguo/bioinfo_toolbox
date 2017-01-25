//
// Created by guoy28 on 1/22/17.
//

#include "OperatorOverloading.h"

OperatorOverloading OperatorOverloading::operator+(OperatorOverloading two) {
    OperatorOverloading sum;
    sum.x = x + two.x;
    sum.y = y + two.y;
    return(sum);
}

OperatorOverloading::OperatorOverloading(int a, int b) {
    x = a;
    y = b;
}