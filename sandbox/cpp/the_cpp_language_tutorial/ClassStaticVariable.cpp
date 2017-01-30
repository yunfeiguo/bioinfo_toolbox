//
// Created by guoy28 on 1/30/17.
//

#include "ClassStaticVariable.h"

int ClassStaticVariable::n = 0;

ClassStaticVariable::ClassStaticVariable(int x) {
    n++;
}
ClassStaticVariable::~ClassStaticVariable() {
    n--;
}