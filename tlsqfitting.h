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

  tlsqfitting.h

  Fitting with regular poly of power up to 4

*******************************************************************************/

#pragma once

#include "tpoints.h"

namespace tcad {

/** Least square fitting */
template <class T> class TLSQFitting {
public:
  enum {
    LSFITTING_MAXPOWER = 4
  };

  // polynomial coefficients
  T coefs[LSFITTING_MAXPOWER] = {0.0};
  // must be 0,1,2,3 ... LSFITTING_MAXPOWER-1 
  int power = -1;

  /** Constructor(s) */
  TLSQFitting()
  {
    power = NOT_DEFINED;
  }

  /** Make fitting. */
  bool fit(int ppower, int numpoints, T x[], T y[])
  {
                                // assume unsuccessful output
    power = NOT_DEFINED;
                                // limit the range
    LIMIT(ppower,1,LSFITTING_MAXPOWER - 1);
                                // allocate matrices for linear system
    int N = ppower + 1;
    std::vector<T> A(N * N,0.0);
                                // fill matrices
    T sum = 0.0;
    for (int i = 0; i < N; i++)
    {
	                            // system matrix
      for (int j = 0; j < N; j++)
      {
        sum = 0.0;
        for (int k = 0; k < numpoints; k++)
	      {
          sum += pow(x[k],i) * pow(x[k],j);
	      }
        A[i*N+j] = sum;
      }

                                // right-hand side
		  sum = 0.0;
      for (int k = 0; k < numpoints; k++)
      {
        sum += y[k] * pow(x[k],i);
      }
      coefs[i] = sum;
    }
                                // solve the system
    if (!solveSystem(N,&A[0],&coefs[0],TOLERANCE(T)))
    {	  
		  return false;
    }	
                                // set sign that it is OK	  
    power = ppower; 

    return true;
  } 

  /** Get value on original dimensional x */
  T GetValue(T x)
  {
    if (power != NOT_DEFINED)
    {
      T sum = 0.0;
      for (int i = 0; i < (power + 1); i++)
      {
	      sum += coefs[i] * pow(x,i);
	    }

	    return sum;
    } else
    {
      return 0.0;
    }
  }
                            
};

}
