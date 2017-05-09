//
// Created by guoy28 on 1/30/17.
//

#include "FriendExample.h"
FriendExample::FriendExample(int a, int b){
    this->x = a;
    this->y = b;
}
int FriendExample::area(){
    return x*y;
}
FriendExample duplicate(const FriendExample& old){
    FriendExample t (old.x,old.y);
    return t;
}
/*
 * don't return reference of a local variable, which
 * dissappears after the function returns.
FriendExample& duplicate2(const FriendExample& old){
    FriendExample t (old.x,old.y);
    return t;
}
 */
