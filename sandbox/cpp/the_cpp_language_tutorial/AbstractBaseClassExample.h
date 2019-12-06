//
// Created by guoy28 on 2/2/17.
//

#ifndef THE_CPP_LANGUAGE_TUTORIAL_ABSTRACTBASECLASSEXAMPLE_H
#define THE_CPP_LANGUAGE_TUTORIAL_ABSTRACTBASECLASSEXAMPLE_H


class Animal {
protected:
    int age;
public:
    Animal(int);
    virtual void move() =0;
};

class Bird:public Animal {
public:
    Bird(int);
    void move();
};

class Insect:public Animal {
public:
    Insect(int);
    void move();
};

#endif //THE_CPP_LANGUAGE_TUTORIAL_ABSTRACTBASECLASSEXAMPLE_H
