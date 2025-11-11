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
#include "tbeziercurve.h"
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
  TBezierSurface(const TBezierSurface &other) : TBaseSurface<T>()
  {
    this->K1 = other.K1;
    this->K2 = other.K2;
    this->cpoints = other.cpoints;
    this->parmsU = other.parmsU;
    this->parmsV = other.parmsV;

    update();
  }

  /** Assignment operator. */
  TBezierSurface &operator = (const TBezierSurface &other)  
  {
    this->K1 = other.K1;
    this->K2 = other.K2;
    this->cpoints = other.cpoints;
    this->parmsU = other.parmsU;
    this->parmsV = other.parmsV;

    update();

    return *this;
  }

  /** Constructor. Every points[i] is a row of points (U-changing). The number of
    cols/rows is arbitrary : every row/col is approximated as TBezierCurve<T>. */
  TBezierSurface(std::vector<std::vector<TPoint<T>>> &points, int numsegmentsU, int numsegmentsV,
    CurveEndType startU = END_CLAMPED, CurveEndType endU = END_CLAMPED, 
    CurveEndType startV = END_CLAMPED, CurveEndType endV = END_CLAMPED, 
    bool keepmeshrefinementU = true, bool keepmeshrefinementV = true) : TBaseSurface<T>()
  {
    this->K1 = numsegmentsU * 4 - 1;
    this->K2 = numsegmentsV * 4 - 1;

    this->cpoints.clear();
    this->cpoints.resize((this->K1 + 1) * (this->K2 + 1));

    std::vector<TPoint<T>> temp;
    int tempK1 = 0;
    int tempK2 = int(points.size()) - 1;

    // loop by rows
    for (int i = 0; i < int(points.size()); i++)
    {
      TBezierCurve<T> crow(points[i],numsegmentsU,startU,endU,keepmeshrefinementU);
      temp.insert(temp.end(),crow.controlPoints().begin(),crow.controlPoints().end());
      if (i == 0)
      {
        this->parmsU = crow.parameters();
        tempK1 = int(crow.controlPoints().size()) - 1;
      }
    }

    // now we need to redivide temp columns in V direction
    for (int i = 0; i <= tempK1; i++)
    {
      std::vector<TPoint<T>> col;
      tcad::getColumn(temp,tempK1,tempK2,i,col);
      TBezierCurve<T> ccol(col,numsegmentsV,startV,endV,keepmeshrefinementV);
      if (i == 0)
        this->parmsV = ccol.parameters();

      // set this column
      this->setColumn(i,ccol.controlPoints());
    }

    update();
  }

  /** Destructor. */
  virtual ~TBezierSurface() {}

  //===== Abstract =============================================================

  /** Get k-th derivative on parameter U,V [0..1]. 0-derivative is position(U,V) (see below). 
    No UV derivatives. */
  virtual TPoint<T> derivative(T U, T V, Parameter onparameter, int k)
  {
    TBezierPatch<T> patch;
    T u = 0.0;
    T v = 0.0;
    T DU = 1.0;
    T DV = 1.0;
    if (findBezierPatch(U,V,patch,u,v,DU,DV))
    {
      TPoint<T> r = patch.derivative(u,v,onparameter,k);

      if (k == 0)
      {
      } else if (k == 1)
      {
        if (onparameter == PARAMETER_U)
        {
          r /= DU;
        } else if (onparameter == PARAMETER_V)
        {
          r /= DV;
        }
      } else if (k == 2)
      {
        if (onparameter == PARAMETER_UU)
        {
          r /= (DU * DU);
        } else if (onparameter == PARAMETER_VV)
        {
          r /= (DV * DV);
        }
      }

      return r;
    } else
    {
      return TPoint<T>();
    }
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
  bool findBezierPatch(T U, T V, TBezierPatch<T> &patch, T &u, T &v, T &DU, T &DV)
  {
    LIMIT(U,0.0,1.0);
    LIMIT(V,0.0,1.0);

    // get two columns with patch between
    u = 0.0;
    int indexU = findParametricInterval(parmsU,U,&u);
    assert(indexU >= 0);
    if (indexU < 0)
      return false;
    int indexU4 = indexU * 4;

    std::vector<TPoint<T>> col0,col1;
    this->getColumn(indexU4,col0);
    this->getColumn(indexU4 + 3,col1);

    // get two rows with patch between
    v = 0.0;
    int indexV = findParametricInterval(parmsV,V,&v);
    assert(indexV >= 0);
    if (indexV < 0)
      return false;
    int indexV4 = indexV * 4;

    std::vector<TPoint<T>> row0,row1;
    this->getRow(indexV4,row0);
    this->getRow(indexV4 + 3,row1);

    TBezierSegment<T> SU0,S1V,SU1,S0V;
    SU0.init(row0[indexU4],row0[indexU4 + 1],row0[indexU4 + 2],row0[indexU4 + 3]);
    SU1.init(row1[indexU4],row1[indexU4 + 1],row1[indexU4 + 2],row1[indexU4 + 3]);
    S1V.init(col1[indexV4],col1[indexV4 + 1],col1[indexV4 + 2],col1[indexV4 + 3]);
    S0V.init(col0[indexV4],col0[indexV4 + 1],col0[indexV4 + 2],col0[indexV4 + 3]);

    patch.init(&SU0,&S1V,&SU1,&S0V);

    DU = parmsU[indexU + 1] - parmsU[indexU];
    DV = parmsV[indexV + 1] - parmsV[indexV];

    return true;
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
