#ifndef WORK_PRELIMDM_TAGGERISR_H_
#define WORK_PRELIMDM_TAGGERISR_H_

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"

#include <vector>

class taggerISR() {

taggerISR();
taggerISR(std::vector<Jet*>& jets);


};

#endif  
