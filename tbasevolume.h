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

  tbasevolume.h

  Basic abstract class for volumes

  dimensions : 3 (U,V,W parameters)

*******************************************************************************/

#pragma once

#include "tbasics.h"
#include "tpoint.h"
#include "ttransform.h"

namespace tcad {

template <class T> class TBaseVolume {
public:

  //===== Construction =========================================================

  /** Constructor. */
  TBaseVolume() {}

  /** Destructor. */
  virtual ~TBaseVolume() {}

  //===== Abstract =============================================================

  /** Get k-th derivative on parameter U,V [0..1]. 0-derivative is position(U,V) (see below). */
  virtual TPoint<T> derivative(T U, T V, T W,Parameter onparameter, int k) = 0;

  /** Update after any change in control points. */
  virtual void update() = 0;

  //===== Operations ===========================================================

  /** Same as 0-th derivative*/
  virtual TPoint<T> position(T U, T V, T W)
  {
    return this->derivative(U,V,W,PARAMETER_ANY,0);
  }

  /** Get control points. */
  virtual std::vector<TPoint<T>> &controlPoints()
  {
    return cpoints;
  }

  /** Make transform. */
  virtual void makeTransform(TTransform<T> *transform)
  {
    for (auto &p : cpoints)
    {
      p = transform->applyTransform(p);
    }

    // call virtual update()
    this->update();
  }

  /** Are equal? */
  virtual bool equal(TBaseVolume &other, T tolerance)
  {
    T diff = difference(cpoints,other.cpoints);
    return (diff >= 0.0 && diff < tolerance);
  }

  /** Find values of parameters U,V for a point on (or close to) the curve. Returns -1 in failure. */
  virtual TPoint<T> findUVforPoint(TPoint<T> p, int numpointsU = MANY_POINTS, int numpointsV = MANY_POINTS)
  {
    // create fine mesh of points
    std::vector<TPoint<T>> points;
    int k1 = 0;
    int k2 = 0;
    createPoints(points,nullptr,&k1,&k2,numpointsU,numpointsV);

    // temp
    TPoint<T> proj;
    int seg = 0;
    T u = 0.0;

    // look for closest point on rows
    T U = -1.0;
    T mindist = std::numeric_limits<T>::max();
    for (int i = 0; i < k2 + 1; i++)
    {
      std::vector<TPoint<T>> row;
      tcad::getRow(points,k1,k2,i,row);
      if (projectPointOnPoints(row,p,proj,&seg,&u))
      {
        T dist = !(p - proj);
        if (dist < mindist)
        {
          TPoint<T> p0 = row[seg];
          TPoint<T> p1 = row[seg + 1];
          U = p0.W + (p1.W - p0.W) * u;
          mindist = dist;
        }
      }
    }

    // look for closest point on columns
    T V = -1.0;
    mindist = std::numeric_limits<T>::max();
    for (int i = 0; i < k1 + 1; i++)
    {
      std::vector<TPoint<T>> col;
      tcad::getColumn(points,k1,k2,i,col);
      if (projectPointOnPoints(col,p,proj,&seg,&u))
      {
        T dist = !(p - proj);
        if (dist < mindist)
        {
          TPoint<T> p0 = col[seg];
          TPoint<T> p1 = col[seg + 1];
          V = p0.W + (p1.W - p0.W) * u;
          mindist = dist;
        }
      }
    }

    return TPoint<T>(U,V);
  }

  /** Calculate min/max. */
  virtual bool calculateMinMax(TPoint<T> *min, TPoint<T> *max, TPoint<T> *imin = nullptr, 
    TPoint<T> *imax = nullptr)
  {
    return tcad::calculateMinMax(cpoints,min,max,imin,imax);
  }

public: //!!!!!!!

  // number of columns minus 1
  int K1 = 0;
  // number of rows munus 1
  int K2 = 0;

  int K3 = 0;

protected:

  // (K1 + 1) * (K2 + 1) * (K3 + 1) control points, call update() after every change
  std::vector<TPoint<T>> cpoints;

};

}
