//
// Created by guoy28 on 3/12/17.
//

#ifndef THE_CPP_LANGUAGE_TUTORIAL_TEMPLATESPECIALIZATION_H
#define THE_CPP_LANGUAGE_TUTORIAL_TEMPLATESPECIALIZATION_H


class templateSpecialization {

};

template <class T>
class myContainer {
    T element;
public:
    myContainer (T arg) {
        element = arg;
    }
    T increase() {
        return ++element;
    }
};

template <>
class myContainer<char> {
    char element;
public:
    myContainer(char arg) {
        element = arg;
    }
    char uppercase() {
        if ((element >= 'a') && (element <= 'z')) {
            element += 'A' - 'a';
        }
        return element;
    }
};

void testTemplateSpecialization();

#endif //THE_CPP_LANGUAGE_TUTORIAL_TEMPLATESPECIALIZATION_H
