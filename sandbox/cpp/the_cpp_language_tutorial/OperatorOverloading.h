//
// Created by guoy28 on 1/22/17.
//

#ifndef THE_CPP_LANGUAGE_TUTORIAL_OPERATOROVERLOADING_H
#define THE_CPP_LANGUAGE_TUTORIAL_OPERATOROVERLOADING_H


class OperatorOverloading {
public:
    int x, y;
    OperatorOverloading() {};
    OperatorOverloading(int a, int b);
    OperatorOverloading operator +(OperatorOverloading);
};


#endif //THE_CPP_LANGUAGE_TUTORIAL_OPERATOROVERLOADING_H
