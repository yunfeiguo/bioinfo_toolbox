//
// Created by guoy28 on 2/2/17.
//

#include "AbstractBaseClassExample.h"
#include <iostream>
using namespace std;

Animal::Animal(int i) {
    cout << "I am " << i << " years/months old" << endl;
}
Bird::Bird(int i) : Animal(i) {}
Insect::Insect(int i) : Animal(i) {}

void Bird::move() {
    cout << "fly" << endl;
}

void Insect::move() {
    cout << "fly or crawl" << endl;
}
