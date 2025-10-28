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

  tbeziersurface.h

  Bezier surface (regular composite of Bezier patches)

  dimensions : 2 (U,V parameters)

*******************************************************************************/

#pragma once

#include "tbasesurface.h"
#include "tbezierpatch.h"

namespace tcad {

template <class T> class TBezierSurface : public TBaseSurface<T> {
public:

  //===== Construction =========================================================

  /** Constructor. */
  TBezierSurface() : TBaseSurface<T>()
  {
  }

  /** Constructor. */
  TBezierSurface() : TBaseSurface<T>()
  {
    this->K1 = 3;
    this->K2 = 3;

    this->cpoints.clear();

    this->cpoints.insert(this->cpoints.end(),SU0->controlPoints().begin(),SU0->controlPoints().end());
    this->cpoints.insert(this->cpoints.end(),S1V->controlPoints().begin(),S1V->controlPoints().end());
    this->cpoints.insert(this->cpoints.end(),SU1->controlPoints().begin(),SU1->controlPoints().end());
    this->cpoints.insert(this->cpoints.end(),S0V->controlPoints().begin(),S0V->controlPoints().end());

    update();
  }

  /** Destructor. */
  virtual ~TBezierSurface() {}

  //===== Abstract =============================================================

  /** Get k-th derivative on parameter U,V [0..1]. 0-derivative is position(U,V) (see below). 
    No UV derivatives. */
  virtual TPoint<T> derivative(T U, T V, Parameter onparameter, int k)
  {
    if (k == 0)
    {
    } else if (k == 1)
    {
      if (onparameter == PARAMETER_U)
      {
      } else if (onparameter == PARAMETER_V)
      {
        Vvector[3][0] = 3.0 * V * V;
      }
    } else if (k == 2)
    {
      if (onparameter == PARAMETER_U)
      {
        Uvector[0][0] = 0.0; 
        Uvector[0][1] = 0.0; 
        Uvector[0][2] = 2.0; 
        Uvector[0][3] = 6.0 * U; 
        Vvector[0][0] = 1.0;
        Vvector[1][0] = V;
        Vvector[2][0] = V * V;
        Vvector[3][0] = V * V * V;
      } else if (onparameter == PARAMETER_V)
      {
        Uvector[0][0] = 1.0; 
        Uvector[0][1] = U; 
        Uvector[0][2] = U * U; 
        Uvector[0][3] = U * U * U; 
        Vvector[0][0] = 0.0;
        Vvector[1][0] = 0.0;
        Vvector[2][0] = 2.0;
        Vvector[3][0] = 6.0 * V;
      }
    }

    Matrix<TPoint<T>,4,1> Bright = this->Q * Vvector;
    Matrix<TPoint<T>,1,1> result = Bright.transpose() * Uvector.transpose();

    return result[0][0];
  }

  /** Update after any change in control points. */
  virtual void update()
  {
  }

  /** Access to U parameters. */
  std::vector<T> &parametersU()
  {
    return parmsU;
  }

  /** Access to V parameters. */
  std::vector<T> &parametersV()
  {
    return parmsV;
  }

  /** Find a patch for these U,V. */
  bool findBezierPatch(T U, T V, TBezierPatch &patch, T &u, T &v)
  {
    LIMIT(U,0.0,1.0);
    LIMIT(V,0.0,1.0);

    // get two columns with patch between
    T u = 0.0;
    int indexU = findParametricInterval(parmsU,U,&u);
    assert(indexU >= 0);

    std::vector<TPoint<T>> col0,col1;
    this->getColumn(indexU,col0);
    this->getColumn(indexU + 1,col1);

    // get two rows with patch between
    T v = 0.0;
    int indexV = findParametricInterval(parmsV,V,&v);
    assert(indexV >= 0);

    std::vector<TPoint<T>> row0,row1;
    this->getRow(indexV,row0);
    this->getRow(indexV + 1,row1);

    TBezierSegment SU0,S1V,SU1,S0V;


  }

private:

  // every row/columns makes a Bezier curve of K1/K2 points;
  // every curve contains (K1/K2 + 1) / 4 Bezier segments

  //// number of columns minus 1
  //int K1 = 0;
  //// number of rows munus 1
  //int K2 = 0;

  //// (K1 + 1) * (K2 + 1) control points, call update() after every change
  //std::vector<TPoint<T>> cpoints;

  // parameters [0..1] parameterised by length
  std::vector<T> parmsU;
  std::vector<T> parmsV;
};

}
