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

  tplane.h

  TPlane

*******************************************************************************/

#pragma once

#include "tpoint.h"
#include "tmisc.h"

namespace tcad {

/** Intersection of segment with plane */
enum {
  SPINTR_DEGONPLANE = -1,	///< segment is degenerate and both vectors are on plane
  SPINTR_NONE,				    ///< no intersection at all
  SPINTR_TOUCHING1,			  ///< segment has finite length within tolerance, one TOUCHING intersection, by vector V1
  SPINTR_TOUCHING2,			  ///< TOUCHING intersection, by vector V2
  SPINTR_NORMAL1,			    ///< segment has finite length, one intersection
  SPINTR_NORMAL2			    ///< segment has finite length, both vectors lie on plane
};

/** Intersection bw segment and plane. Simple, one intersection. */
template <class T> bool segPlaneIntersect(TPoint<T> V1, TPoint<T> V2, TPoint<T> N, T D,
  TPoint<T> *I, T *U, T tolerance = 0.0)
{
  T K1 = V1 * N;
  T K2 = V2 * N;
  T DK = K2 - K1;
  if (std::abs(DK) <= TOLERANCE(T))
  {
    return false;
  } else
  {
    T len = !(V2 - V1);
    double parmtolerance = parmTolerance(len,tolerance);

    *U = (-K1 - D) / DK;
    if ((*U >= -parmtolerance) && (*U <= 1.0 + parmtolerance))
    {
      *I = V1 + (V2 - V1) * (*U);
      return true;
    } else
    {
      return false;
    }
  }
}

/** Intersection bw segment and plane. All cases considered. */
template <class T> int segPlaneIntersect(TPoint<T> V1, TPoint<T> V2, TPoint<T> N, T D,
  TPoint<T> *I1, T *U1, TPoint<T> *I2, T *U2, T *quality = nullptr)
{
  T K1 = V1 * N;
  T K2 = V2 * N;
  T DK = K2 - K1;

  // V1-V2 length
  TPoint<T> DV = V2 - V1;
  T len = !DV;
  T parmtolerance = parmTolerance(len);

  // Ax + By + Cz + D = ...
  T equation1 = K1 + D;
  T equation2 = K2 + D;
  bool onplane1 = (std::abs(equation1) <= TOLERANCE(T));
  bool onplane2 = (std::abs(equation2) <= TOLERANCE(T));

  // two close points within tolerance
  if (len <= TOLERANCE(T))
  {
    if (onplane1)
    {
      // two close points within tolerance which lie on plane
      *I1 = V1;
      *U1 = 0.0;
      if (quality != nullptr) *quality = 0.0;
      return SPINTR_DEGONPLANE;
    } else
    {
      // two close points far from the plane
      if (quality != nullptr) *quality = 0.0;
      return SPINTR_NONE;
    }
  } else
  {
    // the points are apart and lie on the plane
    if (onplane1 && onplane2)
    {
      *I1 = V1;
      *U1 = 0.0;
      *I2 = V2;
      *U2 = 1.0;
      if (quality != nullptr) *quality = 1.0 - std::abs((+DV) * (+N));
      return SPINTR_NORMAL2;
    } else
    {
      // solution
      *U1 = -equation1 / DK;

      // always report intersection point
      *I1 = V1 + (V2 - V1) * (*U1);

      // compare parameter by tolerance
      if ((*U1 >= 0.0 - parmtolerance) && (*U1 <= 1.0 + parmtolerance))
      {
        if (*U1 >= -parmtolerance && *U1 <= parmtolerance)
        {
          if (quality != nullptr) *quality = 1.0 - std::abs((+DV) * (+N));
          return SPINTR_TOUCHING1;
        } else if (*U1 >= 1.0 - parmtolerance && *U1 <= 1.0 + parmtolerance)
        {
          if (quality != nullptr) *quality = 1.0 - std::abs((+DV) * (+N));
          return SPINTR_TOUCHING2;
        } else
        {
          if (quality != nullptr) *quality = 1.0 - std::abs((+DV) * (+N));
          return SPINTR_NORMAL1;
        }
      } else
      {
        if (quality != nullptr) *quality = 0.0;
        return SPINTR_NONE;
      }
    }
  }
}

template <typename T> struct TPlane {

  typedef TPoint<T> vector3;

  // plane normal (normalised) - normalised together with plane constant
  vector3 normal = vector3(1.0,0.0,0.0);
  T constant = 0.0;

  /** Constructors */
  TPlane() = default;

  /** Copy constructor */
  TPlane(const TPlane& copy)
  {
    normal = copy.normal;
    constant = copy.constant;
  }

  /** Assignment operator. */
  TPlane& operator = (const TPlane& copy)
  {
    normal = copy.normal;
    constant = copy.constant;
    return *this;
  }

  /** Cconstructor from 4 scalars. */
  TPlane(const T A, const T B, const T C, const T D)
  {
    normal = vector3(A,B,C);
    constant = D;
  }

  /** Constructor from normal and scalar. */
  TPlane(const vector3 &N, const T D)
  {
    normal = N;
    constant = D;
  }

  /** Constructor form normal and point on plane. */
  TPlane(const vector3 &N, const vector3 &P)
  {
    normal = N;
    constant = - (P * N);
  }

  /** Constructor from 3 points */
  TPlane(const vector3 &V1, const vector3 &V2, const vector3 &V3, bool& OK)
  {
    OK = makePlaneOf3Vectors<T>(V1,V2,V3,normal,constant);
  }

  /** Tolerant equality */
  bool operator == (TPlane<T> &other) const
  {
    return (normal == other.normal && std::abs(constant - other.constant) < normal.tolerance());
  }

  /** Get distance (signed) from point to plane */
  T distance(const vector3 &V) const
  {
    return ((V * normal) + constant) / (normal * normal);
  }

  //* Get distance (signed) to another plane. */
  T distance(const TPlane& other) const
  {
    // normals must be normalised
    assert((!normal - static_cast<T>(1.0)) < normal.tolerance());
    assert((!other.normal - static_cast<T>(1.0)) < normal.tolerance());

    return (other.constant - constant);
  }

  /** Parallel to another plane? */
  bool is_parallel(const TPlane& other) const
  {
    // normals must be normalised
    assert((!normal - static_cast<T>(1.0)) < normal.tolerance());
    assert((!other.normal - static_cast<T>(1.0)) < normal.tolerance());

    T angular_tolerance = std::numeric_limits<T>::epsilon() * static_cast<T>(100.0);

    return (!(normal ^ other.normal) < angular_tolerance);
  }

  /** Test intersection of a straight-line segment 1-2 and a plane; 
    two intersections are possible only if segment lies on the plane; 
    returned are SPINTR_... */
  int segPlaneIntersect(TPoint<T> V1, TPoint<T> V2, TPoint<T> *I1, T *U1, TPoint<T> *I2, 
    T *U2, T *quality = nullptr)
  {
    return tcad::segPlaneIntersect(V1,V2,normal,constant,I1,U1,I2,U2,quality);
  }

  /** Just one intersection is considered */
  bool segmentIntersect(TPoint<T> V1, TPoint<T> V2, TPoint<T> *I1, T *U1, T tolerance) const
  {
    int result = tcad::segPlaneIntersect(V1,V2,normal,constant,I1,U1,tolerance);

    if (result == SPINTR_TOUCHING1 || result == SPINTR_TOUCHING2 || result == SPINTR_NORMAL1)
    {
      return true;
    } else
    {
      return false;
    }
  }
};
 
}
