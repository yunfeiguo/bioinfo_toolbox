//
// Created by guoy28 on 1/18/17.
//

#ifndef THE_CPP_LANGUAGE_TUTORIAL_FIRSTCLASS_H
#define THE_CPP_LANGUAGE_TUTORIAL_FIRSTCLASS_H


class FirstClass {
    int x,y;
public:
    void setValues(int x0, int y0) {
        x = x0;
        y = y0;
    }
    int area() {
        return x*y;
    }

};

class SecondClass {
    int x,y;
public:
    SecondClass(int a, int b);
    ~SecondClass();
    int area();
};

#endif //THE_CPP_LANGUAGE_TUTORIAL_FIRSTCLASS_H
