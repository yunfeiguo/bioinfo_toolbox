//
// Created by guoy28 on 2/2/17.
//

#ifndef THE_CPP_LANGUAGE_TUTORIAL_INHERITANCEEXAMPLE_H
#define THE_CPP_LANGUAGE_TUTORIAL_INHERITANCEEXAMPLE_H

class Runnable {
public:
    int pid;
    Runnable(int);
    int getPid();
};

//all public and protected members of Runnable will become private in Process
class Process : private Runnable {
public:
    Process(int);
    int processGetPid();
};

//LinuxProcess can't inherit any public and protected members of Process
//because they are private members of Process
class LinuxProcess : protected Process {
public:
    LinuxProcess(int);
    int linuxProcessGetPid();
};



#endif //THE_CPP_LANGUAGE_TUTORIAL_INHERITANCEEXAMPLE_H
