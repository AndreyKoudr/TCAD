/**
BSD 2-Clause License

Copyright (c) 2025, Andrey Kudryavtsev (andrewkoudr@hotmail.com)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*******************************************************************************

  Templated CAD

  tdata.h

  Data on aifoils etc.

*******************************************************************************/

#pragma once

#include "tpoint.h"

namespace tcad {

/** NACA0012 points. */
template <class T> const std::vector<std::pair<T,T>> NACA0012xy = {
  {0.000000, 0.000000},
  {0.000685, 0.004611},
  {0.002739, 0.009114},
  {0.006156, 0.013503},
  {0.010926, 0.017770},
  {0.017037, 0.021904},
  {0.024472, 0.025893},
  {0.033210, 0.029726},
  {0.043227, 0.033389},
  {0.054497, 0.036867},
  {0.066987, 0.040145},
  {0.080665, 0.043211},
  {0.095492, 0.046049},
  {0.111427, 0.048648},
  {0.128428, 0.050996},
  {0.146447, 0.053083},
  {0.165435, 0.054902},
  {0.185340, 0.056447},
  {0.206107, 0.057714},
  {0.227680, 0.058702},
  {0.250000, 0.059412},
  {0.273005, 0.059848},
  {0.296632, 0.060015},
  {0.320816, 0.059921},
  {0.345492, 0.059575},
  {0.370590, 0.058989},
  {0.396044, 0.058175},
  {0.421783, 0.057148},
  {0.447736, 0.055923},
  {0.473832, 0.054515},
  {0.500000, 0.052940},
  {0.526168, 0.051216},
  {0.552264, 0.049358},
  {0.578217, 0.047383},
  {0.603956, 0.045307},
  {0.629410, 0.043147},
  {0.654508, 0.040917},
  {0.679184, 0.038634},
  {0.703368, 0.036311},
  {0.726995, 0.033962},
  {0.750000, 0.031603},
  {0.772320, 0.029246},
  {0.793893, 0.026905},
  {0.814660, 0.024593},
  {0.834565, 0.022323},
  {0.853553, 0.020107},
  {0.871572, 0.017959},
  {0.888573, 0.015891},
  {0.904508, 0.013914},
  {0.919335, 0.012042},
  {0.933013, 0.010286},
  {0.945503, 0.008658},
  {0.956773, 0.007168},
  {0.966790, 0.005826},
  {0.975528, 0.004642},
  {0.982963, 0.003626},
  {0.989074, 0.002783},
  {0.993844, 0.002120},
  {0.997261, 0.001644},
  {0.999315, 0.001356},
  {1.000000, 0.0}
//!!!!!!  {1.000000, 0.001260}
};

/** E178 points, these points start from TE, go around LE along 
  the upper surface and then back to TE. */
template <class T> std::vector<TPoint<T>> E178 = {
  TPoint<T>(1.00000, 0.00000),
  TPoint<T>(0.99678, 0.00019),
  TPoint<T>(0.98722, 0.00086),
  TPoint<T>(0.97156, 0.00218),
  TPoint<T>(0.95011, 0.00428),
  TPoint<T>(0.92328, 0.00726),
  TPoint<T>(0.89158, 0.01111),
  TPoint<T>(0.85553, 0.01581),
  TPoint<T>(0.81571, 0.02124),
  TPoint<T>(0.77271, 0.02731),
  TPoint<T>(0.72716, 0.03385),
  TPoint<T>(0.67971, 0.04067),
  TPoint<T>(0.63103, 0.04753),
  TPoint<T>(0.58180, 0.05414),
  TPoint<T>(0.53261, 0.06005),
  TPoint<T>(0.48378, 0.06484),
  TPoint<T>(0.43563, 0.06834),
  TPoint<T>(0.38845, 0.07036),
  TPoint<T>(0.34250, 0.07086),
  TPoint<T>(0.29801, 0.06984),
  TPoint<T>(0.25523, 0.06750),
  TPoint<T>(0.21461, 0.06409),
  TPoint<T>(0.17657, 0.05969),
  TPoint<T>(0.14150, 0.05444),
  TPoint<T>(0.10977, 0.04845),
  TPoint<T>(0.08166, 0.04185),
  TPoint<T>(0.05745, 0.03476),
  TPoint<T>(0.03731, 0.02735),
  TPoint<T>(0.02143, 0.01980),
  TPoint<T>(0.00985, 0.01236),
  TPoint<T>(0.00263, 0.00545),
  TPoint<T>(0.00000, 0.00000),
  //TPoint<T>(0.00000, 0.00018),
  TPoint<T>(0.00302,-0.00470),
  TPoint<T>(0.01190,-0.00876),
  TPoint<T>(0.02598,-0.01222),
  TPoint<T>(0.04524,-0.01497),
  TPoint<T>(0.06953,-0.01698),
  TPoint<T>(0.09870,-0.01829),
  TPoint<T>(0.13248,-0.01896),
  TPoint<T>(0.17054,-0.01905),
  TPoint<T>(0.21249,-0.01867),
  TPoint<T>(0.25788,-0.01788),
  TPoint<T>(0.30619,-0.01678),
  TPoint<T>(0.35687,-0.01546),
  TPoint<T>(0.40934,-0.01398),
  TPoint<T>(0.46296,-0.01243),
  TPoint<T>(0.51712,-0.01086),
  TPoint<T>(0.57114,-0.00933),
  TPoint<T>(0.62440,-0.00787),
  TPoint<T>(0.67625,-0.00651),
  TPoint<T>(0.72607,-0.00529),
  TPoint<T>(0.77326,-0.00421),
  TPoint<T>(0.81725,-0.00328),
  TPoint<T>(0.85749,-0.00250),
  TPoint<T>(0.89350,-0.00187),
  TPoint<T>(0.92482,-0.00136),
  TPoint<T>(0.95106,-0.00092),
  TPoint<T>(0.97194,-0.00049),
  TPoint<T>(0.98729,-0.00014),
  TPoint<T>(0.99677,-0.00001),
  TPoint<T>(1.00000,-0.00000)
};

/** NACA1410 points, these points start from TE, go around LE along 
  the upper surface and then back to TE. */
template <class T> std::vector<TPoint<T>> NACA1410 = {
//  {1.00000, 0.00105},
  TPoint<T>(1.00000, 0.00000),  // TE
  TPoint<T>(0.95021, 0.00832),
  TPoint<T>(0.90034, 0.01513),
  TPoint<T>(0.80049, 0.02741),
  TPoint<T>(0.70051, 0.03804),
  TPoint<T>(0.60042, 0.04692),
  TPoint<T>(0.50025, 0.05385),
  TPoint<T>(0.40000, 0.05836),
  TPoint<T>(0.29937, 0.05940),
  TPoint<T>(0.24907, 0.05809),
  TPoint<T>(0.19880, 0.05531),
  TPoint<T>(0.14861, 0.05062),
  TPoint<T>(0.09854, 0.04338),
  TPoint<T>(0.07358, 0.03837),
  TPoint<T>(0.04870, 0.03194),
  TPoint<T>(0.02398, 0.02297),
  TPoint<T>(0.01174, 0.01639),
  TPoint<T>(0.00000, 0.00000),  // LE
  TPoint<T>(0.01326, -0.01515),
  TPoint<T>(0.02602, -0.02055),
  TPoint<T>(0.05130, -0.02726),
  TPoint<T>(0.07642, -0.03157),
  TPoint<T>(0.10146, -0.03462),
  TPoint<T>(0.15139, -0.03844),
  TPoint<T>(0.20120, -0.04031),
  TPoint<T>(0.25093, -0.04091),
  TPoint<T>(0.30063, -0.04064),
  TPoint<T>(0.40000, -0.03936),
  TPoint<T>(0.49975, -0.03439),
  TPoint<T>(0.59958, -0.02914),
  TPoint<T>(0.69949, -0.02304),
  TPoint<T>(0.79951, -0.01629),
  TPoint<T>(0.89966, -0.00901),
  TPoint<T>(0.94979, -0.00512),
  TPoint<T>(1.00000,  0.00000)  // TE
//  {1.00000, -0.00105}
};

/** Middle line of contour of Kilo propeller blade in metres. Every point contains :
  X - Z along span
  Y - X(Z) of middle line 
  Z - blade HALF-width from middle line,
  W - twist angle around middle curve in degrees. */
template <class T> std::vector<TPoint<T>> KiloBlade = {
  TPoint<T>(0.0,    0.0,    0.2,    -60.0),
  TPoint<T>(0.1,  -0.02,    0.2,    -56.25),
  TPoint<T>(0.2,  -0.05,    0.2,    -52.5),
  TPoint<T>(0.3,  -0.08,    0.21,   -48.75),
  TPoint<T>(0.4,  -0.13,    0.23,   -45.0),
  TPoint<T>(0.5,  -0.15,    0.25,   -41.25),
  TPoint<T>(0.6,  -0.16,    0.28,   -37.5),
  TPoint<T>(0.7,  -0.17,    0.30,   -33.75),
  TPoint<T>(0.8,  -0.16,    0.31,   -30.0),
  TPoint<T>(0.9,  -0.14,    0.30,   -26.25),
  TPoint<T>(1.0,  -0.11,    0.28,   -22.5),
  TPoint<T>(1.1,  -0.08,    0.27,   -18.75),
  TPoint<T>(1.2,   0.00,    0.26,   -15.0),
  TPoint<T>(1.3,   0.08,    0.25,   -11.25),
  TPoint<T>(1.4,   0.14,    0.22,   -7.5),
  TPoint<T>(1.5,   0.24,    0.15,   -3.75),
  TPoint<T>(1.57,  0.33,    0.0,    -0.375)
};

/** Kilo propeller hub in Z-X coordinates. */
template <class T> std::vector<TPoint<T>> KiloPropHub = {
  TPoint<T>(0.0,    0.0,    -1.1),
  TPoint<T>(0.13,   0.0,    -1.0),
  TPoint<T>(0.178,  0.0,    -0.9),
  TPoint<T>(0.205,  0.0,    -0.8),
  TPoint<T>(0.222,  0.0,    -0.7),
  TPoint<T>(0.236,  0.0,    -0.6),
  TPoint<T>(0.24,   0.0,    -0.5),
  TPoint<T>(0.24,   0.0,    -0.4),
  TPoint<T>(0.24,   0.0,    -0.3),
  TPoint<T>(0.24,   0.0,    -0.2),
  TPoint<T>(0.24,   0.0,    -0.1),
  TPoint<T>(0.24,   0.0,     0.0)
};

/** Kilo propeller hub end in Z-X coordinates. */
template <class T> std::vector<TPoint<T>> KiloPropHubEnd = {
  TPoint<T>(0.24,   0.0,     0.0),
  TPoint<T>(0.0,    0.0,     0.0)
};

/** Kilo hull in Z-X coordinates. */
#define KHCOEF 0.33

template <class T> std::vector<std::vector<TPoint<T>>> KiloHull = {
  {
    TPoint<T>(0.0 * KHCOEF,    0.0,    -29.7),
    TPoint<T>(0.24,            0.0,    -29.7)
  },
  {
    TPoint<T>(0.24,            0.0,    -29.7),
    TPoint<T>(7.5 * KHCOEF,    0.0,    -60.0 * KHCOEF),
    TPoint<T>(12.5 * KHCOEF,   0.0,    -30.0 * KHCOEF),
    TPoint<T>(15.0 * KHCOEF,   0.0,    -15.0 * KHCOEF),
    TPoint<T>(15.0 * KHCOEF,   0.0,      0.0 * KHCOEF),
  },
  {
    TPoint<T>(15.0 * KHCOEF,   0.0,      0.0 * KHCOEF),
    TPoint<T>(15.0 * KHCOEF,   0.0,    104.0 * KHCOEF),
  },
  {
    TPoint<T>(15.0 * KHCOEF,   0.0,    104.0 * KHCOEF),
    TPoint<T>(15.0 * KHCOEF,   0.0,    110.0 * KHCOEF),
    TPoint<T>(12.5 * KHCOEF,   0.0,    120.0 * KHCOEF),
    TPoint<T>( 8.0 * KHCOEF,   0.0,    130.0 * KHCOEF),
    TPoint<T>( 0.1 * KHCOEF,   0.0,    134.0 * KHCOEF),
    TPoint<T>( 0.0 * KHCOEF,   0.0,    134.0 * KHCOEF)
  }
};

}
