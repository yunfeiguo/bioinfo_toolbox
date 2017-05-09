//
// Created by guoy28 on 3/29/17.
//

#ifndef THE_CPP_LANGUAGE_TUTORIAL_CASTING_H
#define THE_CPP_LANGUAGE_TUTORIAL_CASTING_H


class CastingBase {
    virtual void dummy() {}
};
class CastingBaseDerived : public CastingBase {
    int a;
};

void tryCasting();

#endif //THE_CPP_LANGUAGE_TUTORIAL_CASTING_H
