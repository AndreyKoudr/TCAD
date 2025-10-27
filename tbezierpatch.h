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

  tbezierpatch.h

  Bezier patch

  dimensions : 2 (U,V parameters)

*******************************************************************************/

#pragma once

#include "tmatrix.h"
#include "tbeziersegment.h"
#include "tbasesurface.h"

namespace tcad {

/** This is Ferguson matrix. */
template <class T> constexpr std::array<std::array<T,4>,4> Mferg = {
{
  {{ T(1.0), T(0.0), T(0.0), T(0.0) }},
  {{ T(0.0), T(0.0), T(1.0), T(0.0) }},
  {{T(-3.0), T(3.0),T(-2.0),T(-1.0) }},
  {{ T(2.0),T(-2.0), T(1.0), T(1.0) }}
}
};

template <class T> class TBezierPatch : public TBaseSurface<T> {
public:

  //===== Construction =========================================================

  /** Constructor. */
  TBezierPatch() : TBaseSurface<T>()
  {
    this->K1 = 3;
    this->K2 = 3;
  }

  /** Constructor. */
  TBezierPatch(TBezierSegment<T> *SU0, TBezierSegment<T> *S1V, TBezierSegment<T> *SU1, 
    TBezierSegment<T> *S0V) : TBaseSurface<T>()
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
  virtual ~TBezierPatch() {}

  //===== Abstract =============================================================

  /** Get k-th derivative on parameter U,V [0..1]. 0-derivative is position(U,V) (see below). */
  virtual TPoint<T> derivative(T U, T V, Parameter onparameter, int k)
  {
    Matrix<T,4,4> M(Mferg<T>);
    Matrix<T,4,4> Mt = M.transpose();

    if (k == 0)
    {
      Matrix<T,1,4> Uvector;
      Uvector[0][0] = 1.0; 
      Uvector[0][1] = U; 
      Uvector[0][2] = U * U; 
      Uvector[0][3] = U * U * U; 
      Matrix<T,4,1> Vvector;
      Vvector[0][0] = 1.0;
      Vvector[1][0] = V;
      Vvector[2][0] = V * V;
      Vvector[3][0] = V * V * V;

      Matrix<T,1,4> left = Uvector * M;
      Matrix<T,4,1> right = Mt * Vvector;
      Matrix<TPoint<T>,4,1> Bright = this->B * right;
      Matrix<TPoint<T>,1,1> result = Bright.transpose() * left.transpose();

      return result[0][0];
    } else
    {
//!!!!!!!

return TPoint<T>();
    }
  }

  /** Update after any change in control points. */
  virtual void update()
  {
    // update B matrix
    updateB();
  }

//protected:
//
//  // number of columns minus 1
//  int K1 = 0;
//  // number of rows munus 1
//  int K2 = 0;
//
//  // (K1 + 1) * (K2 + 1) control points, call update() after every change
//  std::vector<TPoint<T>> cpoints;

/** cpoints make 4 Bezier segment along 4 boundaries :
                      cpoints 8..11
                 ------------>-------------
                |                          |
                |                          |
cpoints 12..15  ^                          ^ cpoints 4..7
                |                          |
                |                          |
                 ------------>-------------
                      cpoints 0..3
*/
  
private:

  /** See "Computational geometry for design and manufactute by Faux and Pratt. */

  /** this matrix is updated in update() */
  Matrix<TPoint<T>,4,4> B;

  /** Update B matrix, call after a change in update(). */
  void updateB()
  {
    TBezierSegment<T> s0(this->cpoints[0],this->cpoints[1],this->cpoints[2],this->cpoints[3]);
    TBezierSegment<T> s1(this->cpoints[4],this->cpoints[5],this->cpoints[6],this->cpoints[7]);
    TBezierSegment<T> s2(this->cpoints[8],this->cpoints[9],this->cpoints[10],this->cpoints[11]);
    TBezierSegment<T> s3(this->cpoints[12],this->cpoints[13],this->cpoints[14],this->cpoints[15]);

		this->B[0][0] = s0.derivative(0.0,0);
		this->B[1][0] = s0.derivative(1.0,0);
		this->B[2][0] = s0.derivative(0.0,1);
		this->B[3][0] = s0.derivative(1.0,1);

		this->B[1][2] = s1.derivative(0.0,1);
		this->B[1][3] = s1.derivative(1.0,1);

		this->B[0][1] = s2.derivative(0.0,0);
		this->B[1][1] = s2.derivative(1.0,0);
		this->B[2][1] = s2.derivative(0.0,1);
		this->B[3][1] = s2.derivative(1.0,1);

		this->B[0][2] = s3.derivative(0.0,1);
		this->B[0][3] = s3.derivative(1.0,1);

    // zero twist, uv ders zero
	  this->B[2][2] = TPoint<T>();
	  this->B[2][3] = TPoint<T>();
	  this->B[3][2] = TPoint<T>();
	  this->B[3][3] = TPoint<T>();
  }

};

}
