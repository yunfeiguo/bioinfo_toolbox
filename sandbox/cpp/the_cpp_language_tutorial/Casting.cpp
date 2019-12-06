//
// Created by guoy28 on 3/29/17.
//

#include "Casting.h"
#include<iostream>
#include<typeinfo>
#include<exception>

void tryCasting() {
    CastingBase* pba = new CastingBaseDerived;//implicit casting under the hood
    CastingBase* pbb = new CastingBase;
    CastingBaseDerived* pd;

    try {
        pd = dynamic_cast<CastingBaseDerived*>(pba);
        if (pd == 0) std::cout << "null pointer on castingbasederived" << std::endl;
        pd = dynamic_cast<CastingBaseDerived*>(pbb);
        if (pd == 0) std::cout << "null pointer on castingbase" << std::endl;
    } catch (std::exception& e) {
        std::cout << e.what() << std::endl;
    }

    typeid(pba).name();

}
