//
// Created by guoy28 on 1/30/17.
//

#include "FriendClassExample.h"

int MyRectangle::area() {
    return w*h;
}
void MyRectangle::setH(int h) {this->h = h;}
void MyRectangle::setW(int w) {this->w = w;}
void MySquare::setL(int l) {this->l = l;}

int MySquare::area() {
    return l*l;
}

void MyRectangle::convert(MySquare m) {
    this->w = this->h = m.l;
}
