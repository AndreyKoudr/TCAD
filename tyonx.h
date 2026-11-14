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

  tyonx.h

  Simple Y on X function. Only X,Y are involved. X,Y are monotonically
  increasing, no checks.

*******************************************************************************/

#pragma once

#include "tbasics.h"
#include "tpoint.h"
#include "tpoints.h"

namespace tcad {

/** Simple Y on X function. Only X,Y are involved. X,Y are monotonically
  increasing, no checks. */
template <class T> class TYonX {
public:
  std::vector<TPoint<T>> points;

  /** Constructor. */
  TYonX()
  {
  }

  /** Constructor. */
  TYonX(const std::vector<TPoint<T>> &ppoints)
  {
    points = ppoints;
  }

  /** Copy constructor. */
  TYonX(const TYonX &other)  
  {
    this->points = other.points;
  }

  /** Assignment operator. */
  TYonX &operator = (const TYonX &other)  
  {
    this->points = other.points;

    return *this;
  }

  /** Swap X and Y. */
  void swapXY()
  {
    for (auto &p : points)
    {
      auto temp = p.X; p.X = p.Y; p.Y = p.temp;
    }
  }

  /** Calculate length. */
  T length()
  {
    return calculateLength<T>(points);
  }

  /** Calculate tolerance. */
  T tolerance()
  {
    return length() * std::numeric_limits<T>::epsilon() * static_cast<T>(10.0);
  }

  /* Find Y on X. */
  T findYonX(T X)
  {
    LIMIT(X,points.front().X,points.back().X);

    T tol = tolerance();

    for (int i = 0; i < int(points.size()) - 1; i++)
    {
      TPoint<T> p0 = points[i];
      TPoint<T> p1 = points[i + 1];

      if (X >= p0.X - tol && X <= p1.X + tol)
      {
        T seglen = std::abs(p1.X - p0.X);
        T U = (seglen > tol) ? ((X - p0.X) / seglen) : 0.5;
        LIMIT(U,0.0,1.0);
        T Y = p0.Y + (p1.Y - p0.Y) * U;

        return Y;
      }
    }

    return -1.0;
  }
  
  /* Find X on Y. */
  T findXonY(T Y)
  {
    swapXY();

    T res = findYonX(Y);

    swapXY();

    return res;
  }
};

}