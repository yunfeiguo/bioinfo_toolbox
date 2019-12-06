//
// Created by guoy28 on 1/30/17.
//

#ifndef THE_CPP_LANGUAGE_TUTORIAL_THISUSAGE_H
#define THE_CPP_LANGUAGE_TUTORIAL_THISUSAGE_H


class ThisUsage {
    int x;
    int y;
public:
    ThisUsage(int,int);
    void print();
    void setX(int);
    ThisUsage& operator= (const ThisUsage&);

};


#endif //THE_CPP_LANGUAGE_TUTORIAL_THISUSAGE_H
