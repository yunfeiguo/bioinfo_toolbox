//
// Created by guoy28 on 1/30/17.
//

#ifndef THE_CPP_LANGUAGE_TUTORIAL_FRIENDEXAMPLE_H
#define THE_CPP_LANGUAGE_TUTORIAL_FRIENDEXAMPLE_H


class FriendExample {
    int x;
    int y;
public:
    FriendExample(int, int);
    int area();
    friend FriendExample duplicate(const FriendExample&);
    friend FriendExample& duplicate2(const FriendExample&);
};


#endif //THE_CPP_LANGUAGE_TUTORIAL_FRIENDEXAMPLE_H
