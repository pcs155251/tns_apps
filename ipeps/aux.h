#ifndef __AUX_H_INCLUDED__
#define __AUX_H_INCLUDED__

#include "uni10.hpp"
#include "ipeps.h"
#include <cstdio>

template<typename T>
void meaSimpTimQuants( ipeps<T>& ansatz, UniTensor<T>& ham, double field, string &dataPath );

template<typename T>
void meaFullTimQuants( ipeps<T>& ansatz, UniTensor<T>& ham, double field, string &dataPath );

#include "aux.cpp"

#endif
