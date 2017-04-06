//
// Created by guoy28 on 2/3/17.
//

//#include "TemplateExample.h"
template <class T>
T getMax(T a, T b) {
    return a > b ? a : b;
}

template <class X, class Y>
X getMin(X a, Y b) {
    return a > b? (X) b : a;
};