#ifndef _OPERATOR_H_
#define _OPERATOR_H_

#include "uni10.hpp"

void spin_check(uni10_float32 spin);

uni10::Matrix<uni10_double64> matSx(uni10_float32 spin=0.5);

uni10::Matrix<std::complex<double>> matSy(uni10_float32 spin=0.5);

uni10::Matrix<uni10_double64> matSp(uni10_float32 spin=0.5);

uni10::Matrix<uni10_double64> matSm(uni10_float32 spin=0.5);

uni10::Matrix<uni10_double64> matSz(uni10_float32 spin=0.5);

uni10::Matrix<std::complex<double>> twoSiteTRStransform();

#endif
