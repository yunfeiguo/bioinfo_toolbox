//
// Created by guoy28 on 2/2/17.
//

#include "InheritanceExample.h"

Runnable::Runnable(int p) {
    pid = p;
}
int Runnable::getPid() {
    return pid;
}
Process::Process(int p) : Runnable(p) {
    //pid = p; //this won't work as pid is private member of Runnable
}
int Process::processGetPid() {
    return pid;
}

LinuxProcess::LinuxProcess(int p) : Process(p) {
}

int LinuxProcess::linuxProcessGetPid() {
    return processGetPid();
}
