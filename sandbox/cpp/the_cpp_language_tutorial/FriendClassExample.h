//
// Created by guoy28 on 1/30/17.
//

#ifndef THE_CPP_LANGUAGE_TUTORIAL_FRIENDCLASSEXAMPLE_H
#define THE_CPP_LANGUAGE_TUTORIAL_FRIENDCLASSEXAMPLE_H

class MySquare;
class MyRectangle {
    int w,h;
public:
    int area();
    void setW(int);
    void setH(int);
    /**
     * change current MyRectangle to same
     * dimension as an MySquare
     */
    void convert(MySquare);
};

class MySquare {
    int l;
public:
    int area();
    void setL(int);
    friend class MyRectangle;
};


#endif //THE_CPP_LANGUAGE_TUTORIAL_FRIENDCLASSEXAMPLE_H
