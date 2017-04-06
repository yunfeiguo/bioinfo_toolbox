//
// Created by guoy28 on 3/12/17.
//

#ifndef THE_CPP_LANGUAGE_TUTORIAL_STDEXCEPTION_H
#define THE_CPP_LANGUAGE_TUTORIAL_STDEXCEPTION_H

#include <exception>
class stdException {

};

class myException : public std::exception {
    virtual const char* what() const throw() {
        return "my excpetion here";
    }
};

void testMyException();
#endif //THE_CPP_LANGUAGE_TUTORIAL_STDEXCEPTION_H
