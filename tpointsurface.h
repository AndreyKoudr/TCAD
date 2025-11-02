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

  tpointsurface.h

  Point surface (a regular net of straight-line segments)

  dimensions : 2 (U,V parameters)

*******************************************************************************/

#pragma once

#include "tbasesurface.h"

namespace tcad {

template <class T> class TPointSurface : public TBaseSurface<T> {
public:

  //===== Construction =========================================================

  /** Constructor. */
  TPointSurface() : TBaseSurface<T>()
  {
  }

  /** Constructor. Every points[i] is a row of points (U-changing). */
  TPointSurface(std::vector<std::vector<TPoint<T>>> &points,
    bool parametersbynumbersU = false, bool parametersbynumbersV = false) : TBaseSurface<T>()
  {
    parmsbynumbersU = parametersbynumbersU;
    parmsbynumbersV = parametersbynumbersV;

    this->cpoints.clear();

    this->K1 = int(points[0].size()) - 1;
    this->K2 = int(points.size()) - 1;

    for (int i = 0; i < int(points.size()); i++)
    {
      this->cpoints.insert(this->cpoints.end(),points[i].begin(),points[i].end());
    }

    update();
  }

  /** Destructor. */
  virtual ~TPointSurface() {}

  //===== Abstract =============================================================

  /** Get k-th derivative on parameter U,V [0..1]. 0-derivative is position(U,V) (see below). 
    No UV derivatives. */
  virtual TPoint<T> derivative(T U, T V, Parameter onparameter, int k)
  {
    std::array<TPoint<T>,4> corners;
    T u = 0.0;
    T v = 0.0;
    T DU = 1.0;
    T DV = 1.0;
    if (findPatch(U,V,corners,u,v,DU,DV))
    {
      TPoint<T> r;

      if (k == 0)
      {
        TPoint<T> func;
        rectShapeFunc01(u,v,func);
        r = corners[0] * func.X + corners[1] * func.Y + corners[2] * func.Z + corners[3] * func.W;
      } else if (k == 1)
      {
        if (onparameter == PARAMETER_U)
        {
          TPoint<T> func;
          rectShapeFuncDerU01(u,v,func);
          r = corners[0] * func.X + corners[1] * func.Y + corners[2] * func.Z + corners[3] * func.W;
          r /= DU;
        } else if (onparameter == PARAMETER_V)
        {
          TPoint<T> func;
          rectShapeFuncDerV01(u,v,func);
          r = corners[0] * func.X + corners[1] * func.Y + corners[2] * func.Z + corners[3] * func.W;
          r /= DV;
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
    std::vector<TPoint<T>> row;
    getRow(this->cpoints,this->K1,this->K2,0,row);
    prepareParameters(row,parmsU,true,parmsbynumbersU);

    std::vector<TPoint<T>> col;
    getColumn(this->cpoints,this->K1,this->K2,0,col);
    prepareParameters(col,parmsV,true,parmsbynumbersV);
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

  /** Find a patch for these U,V. corners are numbered counter clockwise from 
    the lower left corner. */
  bool findPatch(T U, T V, std::array<TPoint<T>,4> &corners, T &u, T &v, T &DU, T &DV)
  {
    LIMIT(U,0.0,1.0);
    LIMIT(V,0.0,1.0);

    // get two columns with patch between
    u = 0.0;
    int indexU = findParametricInterval(parmsU,U,&u);
    assert(indexU >= 0);
    if (indexU < 0)
      return false;

    // get two rows with patch between
    v = 0.0;
    int indexV = findParametricInterval(parmsV,V,&v);
    assert(indexV >= 0);
    if (indexV < 0)
      return false;

    std::vector<TPoint<T>> row0,row1;
    this->getRow(indexV,row0);
    this->getRow(indexV + 1,row1);

    corners[0] = row0[indexU];
    corners[1] = row0[indexU + 1];
    corners[2] = row1[indexU + 1];
    corners[3] = row1[indexU];

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

  // parameters [0..1] parameterised by length or numbers
  bool parmsbynumbersU = false;
  bool parmsbynumbersV = false;
  std::vector<T> parmsU;
  std::vector<T> parmsV;
};

}
