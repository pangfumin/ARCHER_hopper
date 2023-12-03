//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: functionDynamicForwardTest.h
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 11:20:34
//

#ifndef FUNCTIONDYNAMICFORWARDTEST_H
#define FUNCTIONDYNAMICFORWARDTEST_H

// Include Files
// #include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
extern void DynamicForward(
    double Mw, double Mc, double Mp, 
    double Iwx, double Iwy, double Iwz,
    double Icx, double Icy, double Icz, 
    double Ipx, double Ipy, double Ipz,
    double Rw, double Lc, double Lcp,  double g, 
    double Tw, double Tp,
    double q1, double q2,  double q5, 
    double q1_d, double q2_d, double q3_d, double q4_d, double q5_d, 
    double M_nonlin[25],
    double RHS_nonlin[5]);

#endif
//
// File trailer for functionDynamicForwardTest.h
//
// [EOF]
//
