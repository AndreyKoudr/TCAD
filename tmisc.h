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

  tmisc.h

  Miscellaneous functions

*******************************************************************************/

#pragma once

#include "tbasics.h"
#include "tsystems.h"
#include "tpoint.h"

namespace tcad {

/** How to straighten three vectors. */
typedef enum {By13,By12,By23,None} StraightenTypes;

/** Parametric tolerance on length */
template <class T> T parmTolerance(T len, T tolerance = 0.0)
{
  double parmtol = (len < TOLERANCE(T) || tolerance <= 0.0) ? PARM_TOLERANCE : (tolerance / len);
  return parmtol;
}

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
  {1.000000, 0.001260}
};

/** Get sign. */
template <class T> T sign(T value)
{
  return (value >= 0.0) ? +1 : -1; 
}

/** All busy? */
bool allBusy(std::vector<bool> &ibusy);

/** Find first not busy, return index, -1 in case of failure. */
int findFirstNotBusy(std::vector<bool> &ibusy);

/** Divide point indices (into sharp corners). size is the number of points. */
void divideByIndices(std::vector<int> &indices, int size, std::vector<std::pair<int,int>> &division);

/** Get next increasing element of an array. */
template <class TI> TI Tnext(TI size, TI &index)
{
  if (++index >= TI(size))
  {
    index = 0;
  }
  return index;
}

/** Derivative of (x^n)' = n * x ^ (n - 1) */
template <class T> T derivativeXn(T x, int power)
{
  return (power > 0) ? T(power) * pow(x,power - 1) : 0.0;
}

/** Build a plane defined by normal N (normalised) and scalar D. Vectors 01 and 02  should 
  not be collinear. */
template <class T> bool makePlaneOf3Vectors(const TPoint<T> &V0, const TPoint<T> &V1, const TPoint<T> &V2, 
  TPoint<T> &N, T &D)
{
  N = (V1 - V0) ^ (V2 - V0);
  T len = !N;
  if (len > TOLERANCE(T))
  {
    // normal normalised
    N = N * (static_cast<T>(1.0) / len); D = -(N * V0); 
    return true;
  }
  else
  {
    return false;
  }
}

/** Intersection of segments in XY. */
template <typename T> bool intersectSegmentsXY(T X1, T Y1, T X2, T Y2, T X3, T Y3, T X4, T Y4,
  T *t1, T *t2, T *Xi, T *Yi)
{
  T a11, a12, a21, a22, b1, b2;
  a11 = X2 - X1; a12 = X3 - X4;
  a21 = Y2 - Y1; a22 = Y3 - Y4;
  b1 = X3 - X1; b2 = Y3 - Y1;
  *t1 = *t2 = 0;
  if (solveSystem2x2(a11,a12,a21,a22,b1,b2,t1,t2))
  {
    *Xi = X1 + (X2 - X1) * (*t1); 
    *Yi = Y1 + (Y2 - Y1) * (*t1);
    return true;
  } else
  {
    return false;
  }
}

/** Find closest points between two straight-line segments in 3D. */
template <class T> bool intersectSegments(TPoint<T> p0, TPoint<T> p1, TPoint<T> v0, TPoint<T> v1,
  T &l, T &m, T &dist, TPoint<T> *ip = nullptr, TPoint<T> *iv = nullptr)
{
  l = m = dist = 0.0;

  TPoint<T> a = p0;
  TPoint<T> b = p1 - p0;
  TPoint<T> c = v0;
  TPoint<T> d = v1 - v0;

  std::vector<T> A(4);
  std::vector<T> B(2);

  TPoint<T> ac = a - c;
  A[0] = b * b;
  A[1] = - (d * b);
  A[2] = b * d;
  A[3] = - (d * d);

  B[0] = - (ac * b);
  B[1] = - (ac * d);

  bool ok = solveSystem2x2(A[0],A[1],A[2],A[3],B[0],B[1],&l,&m);

  if (ok)
  {
    TPoint<T> p = p0 + b * l;
    TPoint<T> v = v0 + d * m;

    dist = !(p - v);

    if (ip)
      *ip = p;
    if (iv)
      *iv = v;
  }

  return ok;
}

/** Calculate porjection of a point on straight-line segment. */
template <typename T> bool projectPointOnSegment(TPoint<T> p, TPoint<T> p0, TPoint<T> p1, 
  TPoint<T> *intr, T *t, T parmtolerance = PARM_TOLERANCE)
{
  TPoint<T> v10 = p1 - p0;
  T len = !v10;
  TPoint<T> dir = v10 / len;

  TPoint<T> v00 = p - p0;
  *t = v00 * dir;
  *intr = p0 + dir * (*t);
  *t /= len;

  if (*t >= 0.0 - parmtolerance && *t <= 1.0 + parmtolerance)
  {
    return true;
  } else
  {
    return false;
  }
}

/** Find interval in a MONOTONICALLY INCREASING table [0..1], returns -1 as failure. 
  Returned index is never end point, always one less. */
template <class T> int findParametricInterval(const std::vector<T> &table, const T value, T *u = nullptr)
{
  if (table.size() < 2)
    return -1;

  for (int i = 0; i < int(table.size()) - 1; i++)
  {
    T tolerance = (i == int(table.size()) - 2) ? TOLERANCE(T) : 0.0;
    if (table[i + 1] + tolerance >= value)
    {
      if (u)
      {
        T dt = table[i + 1] - table[i];
        *u = (dt > TOLERANCE(T)) ? (value - table[i]) / dt : 0.0;
        LIMIT(*u,0.0,1.0);
      }
      return i;
    }
  }

  return -1;
}

/** Find interval in a table. */
template <class T> int findInterval(const std::vector<T> &table, const T value)
{
  // size, signed integer to avoid problems with size() - 2 (last interval)
  int size = static_cast<int>(table.size());

  // increasing?
  bool increasing = ((table[size - 1] - table[0]) >= 0.0);

  // table must be monotonic
#ifdef _DEBUG
  if (increasing)
  {
    for (int i = 0; i < size - 1; i++)
    {
      if (table[i + 1] < table[i])
        assert(false && "Table not monotonic");
    }
  }
  else
  {
    for (int i = 0; i < size - 1; i++)
    {
      if (table[i + 1] > table[i])
        assert(false && "Table not monotonic");
    }
  }
#endif

  // lower and upper, originally out of scope
  int lower = -1;
  int upper = size - 1;

  // min/max
  T min = increasing ? table[0] : table[size - 1];
  T max = increasing ? table[size - 1] : table[0];

  // tolerance
  T tolerance = TOLERANCE(T);

  // obvious outcome
  if (value < min - tolerance)
    return -1;
  if (value > max + tolerance)
    return -1;

  while ((upper - lower) > 1)
  {
    // middle point
    int middle = (lower + upper) >> 1;
    if ((value >= table[middle]) == increasing)
    {
      lower = middle;
    }
    else
    {
      upper = middle;
    }
  }

  // we are going to return lower, correct it just in case
  if (lower < 0) lower = 0;
  if (lower > size - 2) lower = size - 2;

  // done!
  return lower;
}

/** Spline basis. */
template <class T> void splineBasis(int order, T U, std::vector<T> &knots, std::vector<T> &basis)
{
  T d, e;
  int i;

  // just 0.0 and 1.0
  T zero = static_cast<T>(0.0);
  T unit = static_cast<T>(1.0);

  // number of points
  int numpoints = static_cast<int>(knots.size());

  // clear
  std::fill(basis.begin(),basis.end(),0.0);

  // degree
  int degree = order - 1;

  // tolerance
  T tolerance = TOLERANCE(T);

  // if we calculate basis functions at all points
  int point0 = 0;
  int point1 = numpoints;

  // B-spline basis has local support, so most of basis function values are zero
  // find knots interval for U by bisection
  int interval = findInterval<T>(knots,U);
  assert(interval != -1 && "Unable to find parametric interval for parameter");

  // take only points around i
  point0 = interval - degree - 1;
  point1 = interval + degree + 2;
  if (point0 < 0) point0 = 0;
  if (point1 > numpoints) point1 = numpoints;

  // calculate first order basis functions
  for (i = point0; i < (point1 - 1); i++)
  {
    T d = knots[i + 1] - knots[i];
    if ((U >= knots[i]) && (U < knots[i + 1]) && d > tolerance)
      basis[i] = unit;
  }

  if (std::abs(U - knots[numpoints - 1]) <= tolerance)
  {
    basis[numpoints - order - 1] = unit;
  }

  // calculate higher order basis functions
  for (int k = 2; k <= order; k++)
  {
    for (i = point0; i < point1 - k; i++)
    {
      // if lower order basis function is zero bail out
      T kd = knots[i + k - 1] - knots[i];

      if (std::abs(kd) > tolerance) d = ((U - knots[i]) * basis[i]) / kd;
      else d = zero;

      kd = knots[i + k] - knots[i + 1];
      if (std::abs(kd) > tolerance) e = ((knots[i + k] - U) * basis[i + 1]) / kd;
      else e = zero;

      basis[i] = d + e;
    }
  }

#ifdef _DEBUG
  T sum = 0.0;
  T sumtolerance = 0.00000001;

  for (int i = 0; i < numpoints; i++)
  {
    sum += basis[i];
  }
  if (fabs(sum - unit) > sumtolerance)
  {
    assert((fabs(sum - unit) <= sumtolerance) && "Error in spline basis");
  }
#endif

}
 
template <class T> void straightenThreePoints(TPoint<T> &p1, TPoint<T> &p2, TPoint<T> &p3,
  StraightenTypes straightenType = By13)
{
  if (straightenType == None)
    return;

  TPoint<T> p13 = p3 - p1; T len13 = !p13;
  TPoint<T> p12 = p2 - p1; T len12 = !p12;
  TPoint<T> p23 = p3 - p2; T len23 = !p23;

  T tolerance = TOLERANCE(T);

  if ((len13 > tolerance) && (len12 > tolerance) && (len23 > tolerance))
  {
    switch (straightenType) {
      case By13: { p13 = +p13; break; }
      case By12: { p13 = +p12; break; }
      case By23: { p13 = +p23; break; }
      default : return;
    }

    p12 = p13; p12 = p12 * len12; p1 = p2 - p12;
    p23 = p13; p23 = p23 * len23; p3 = p2 + p23;
  }
}

/** Linear rectangle, U,V [-1..+1] */
template <class T> void rectShapeFunc(T U, T V, TPoint<T> &func)
{
  func.X = (1.0 - U) * (1.0 - V) * 0.25;
  func.Y = (1.0 + U) * (1.0 - V) * 0.25;
  func.Z = (1.0 + U) * (1.0 + V) * 0.25;
  func.W = (1.0 - U) * (1.0 + V) * 0.25;
}

/** Linear rectangle, U,V [-1..+1] */
template <class T> void rectShapeFuncDerU(T U, T V, TPoint<T> &func)
{
  func.X = -(1.0 - V) * 0.25;
  func.Y = +(1.0 - V) * 0.25;
  func.Z = +(1.0 + V) * 0.25;
  func.W = -(1.0 + V) * 0.25;
}

/** Linear rectangle, U,V [-1..+1] */
template <class T> void rectShapeFuncDerV(T U, T V, TPoint<T> &func)
{
  func.X = -(1.0 - U) * 0.25;
  func.Y = -(1.0 + U) * 0.25;
  func.Z = +(1.0 + U) * 0.25;
  func.W = +(1.0 - U) * 0.25;
}

/** Linear rectangle, U,V [0..1] */
template <class T> void rectShapeFunc01(T U, T V, TPoint<T> &func)
{
  T u = U * 2.0 - 1.0;
  T v = V * 2.0 - 1.0;
  rectShapeFunc(u,v,func);
}

/** Linear rectangle, U,V [0..1] */
template <class T> void rectShapeFuncDerU01(T U, T V, TPoint<T> &func)
{
  T u = U * 2.0 - 1.0;
  T v = V * 2.0 - 1.0;
  rectShapeFuncDerU(u,v,func);

  // from -1..+1 to 0..1
  func *= 2.0;
}

/** Linear rectangle, U,V [0..1] */
template <class T> void rectShapeFuncDerV01(T U, T V, TPoint<T> &func)
{
  T u = U * 2.0 - 1.0;
  T v = V * 2.0 - 1.0;
  rectShapeFuncDerV(u,v,func);

  // from -1..+1 to 0..1
  func *= 2.0;
}

/** Linear rectangle, U,V [-1..+1] */
template <class T> TPoint<T> rectCoord(TPoint<T> corners[4], T U, T V)
{
  T func[4];
  rectShapeFunc<T>(U,V,func);

  TPoint<T> sum = corners[0] * func[0] + corners[1] * func[1] + corners[2] * func[2] + corners[3] * func[3];
  
  return sum;
}

/** Get unifom list of parameters for a number of points. */
template <class T> void getUniParms(int numpoints, std::vector<T> &parms)
{
  parms.clear();

  for (int i = 0; i < numpoints; i++)
  {
    T U = T(i) / T(numpoints - 1);
    parms.push_back(U);
  }
}

/** Get row of control points. Each point contains U parameter in W. */
template <class T> void getRow(std::vector<TPoint<T>> &cpoints, int K1, int K2, int index, std::vector<TPoint<T>> &row)
{
  assert(index >= 0 && index <= K2);

  row.clear();

  int n = index  * (K1 + 1);
  for (int i = 0; i < K1 + 1; i++)
  {
    row.push_back(cpoints[n + i]);

    T U = T(i) / T(K1);
    row.back().W = U;
  }
}

/** Set row of control points. */
template <class T> void setRow(std::vector<TPoint<T>> &cpoints, int K1, int K2, int index, const std::vector<TPoint<T>> &row)
{
  assert(row.size() == K1 + 1);

  int n = index  * (K1 + 1);
  for (int i = 0; i < K1 + 1; i++)
  {
    cpoints[n + i] = row[i];
  }
}

/** Get column of control points. Each point contains V parameter in W. */
template <class T> void getColumn(std::vector<TPoint<T>> &cpoints, int K1, int K2, int index, std::vector<TPoint<T>> &column)
{
  assert(index >= 0 && index <= K1);

  column.clear();

  for (int i = 0; i < K2 + 1; i++)
  {
    column.push_back(cpoints[index + (i * (K1 + 1))]);

    T V = T(i) / T(K2);
    column.back().W = V;
  }
}

/** Set column of control points. */
template <class T> void setColumn(std::vector<TPoint<T>> &cpoints, int K1, int K2, int index, 
  const std::vector<TPoint<T>> &column)
{
  assert(column.size() == K2 + 1);

  for (int i = 0; i < K2 + 1; i++)
  {
    cpoints[index + (i * (K1 + 1))] = column[i];
  }
}

/** Get row of points of 2D array. */
template <class T> void getRow(std::vector<std::vector<TPoint<T>>> &points, int index, 
  std::vector<TPoint<T>> &row)
{
  assert(index >= 0 && index < int(points.size()));

  row = points[index];
}

/** Set row of points of 2D array. */
template <class T> void setRow(std::vector<std::vector<TPoint<T>>> &points, int index, 
  std::vector<TPoint<T>> &row)
{
  assert(index >= 0 && index < int(points.size()));
  assert(row.size() == points[index].size());

  points[index] = row;
}

/** Get column of points of 2D array. */
template <class T> void getColumn(std::vector<std::vector<TPoint<T>>> &points, int index, 
  std::vector<TPoint<T>> &col)
{
  assert(index >= 0 && index < int(points[0].size()));

  col.clear();
  for (int i = 0; i < int(points.size()); i++)
  {
    col.push_back(points[i][index]);
  }
}

/** Set column of points of 2D array. */
template <class T> void setColumn(std::vector<std::vector<TPoint<T>>> &points, int index, 
  std::vector<TPoint<T>> &col)
{
  assert(index >= 0 && index < int(points[index].size()));
  assert(col.size() == points.size());

  for (int i = 0; i < int(points.size()); i++)
  {
    points[i][index] = col[i];
  }
}

/** Get derivatives (normalised) of points of 2D array at U = 0 boundary. */
template <class T> void getDerivativesU0(std::vector<std::vector<TPoint<T>>> &points,  
  std::vector<TPoint<T>> &der)
{
  der.clear();

  std::vector<TPoint<T>> p0,p1;
  getColumn(points,0,p0);
  getColumn(points,1,p1);
  
  std::vector<TPoint<T>> d = p1 - p0;
  der = +d;
}

/** Get derivatives (normalised) of points of 2D array at U = 1 boundary. */
template <class T> void getDerivativesU1(std::vector<std::vector<TPoint<T>>> &points,  
  std::vector<TPoint<T>> &der)
{
  der.clear();

  std::vector<TPoint<T>> p0,p1;
  getColumn(points,int(points[0].size() - 2),p0);
  getColumn(points,int(points[0].size() - 1),p1);
  
  std::vector<TPoint<T>> d = p1 - p0;
  der = +d;
}

/** Get derivatives (normalised) of points of 2D array at V = 0 boundary. */
template <class T> void getDerivativesV0(std::vector<std::vector<TPoint<T>>> &points,  
  std::vector<TPoint<T>> &der)
{
  der.clear();

  std::vector<TPoint<T>> p0,p1;
  getRow(points,0,p0);
  getRow(points,1,p1);
  
  std::vector<TPoint<T>> d = p1 - p0;
  der = +d;
}

/** Get derivatives (normalised) of points of 2D array at V = 1 boundary. */
template <class T> void getDerivativesV1(std::vector<std::vector<TPoint<T>>> &points,  
  std::vector<TPoint<T>> &der)
{
  der.clear();

  std::vector<TPoint<T>> p0,p1;
  getRow(points,int(points.size()) - 2,p0);
  getRow(points,int(points.size()) - 1,p1);
  
  std::vector<TPoint<T>> d = p1 - p0;
  der = +d;
}

/* Test piercing triangle by vector point0->point1; if pierced, returns U - 
  parameter to mark point between point0 and point1. U check for [0..1] is here. */
template <typename T> bool segTriIntersect(const TPoint<T> &point0, const TPoint<T> &point1,
  const std::array<TPoint<T>,3> &coords, T &U, TPoint<T> &intersection, 
  const T tolerance, const T parmtolerance = PARM_TOLERANCE)
{
                                      // get vector point0->point1 
  TPoint<T> v01 = point1 - point0;
                                      // get facet plane 
  TPoint<T> N; T D;
  if (!makePlaneOf3Vectors(coords[0],coords[1],coords[2],N,D))
    return false;
                                      // N * (p1 - p0) 
  T r = N * v01;
  U = std::abs(r);
                                      // view direction parallel to facet plane 
  if (U < tolerance)  
    return false;
                                      // get intersection between point0  &point1 
  U = (-D - point0 * N) / r;
  if (U < -parmtolerance || U > 1.0 + parmtolerance)
    return false;
                                      // test if intersection inside facet contour 
  TPoint<T> c = v01 * U; 
  intersection = point0 + c;
                                      // scan each facet side 
  TPoint<T> olddir;
  bool inside = true;
  for (int i = 0; i < 3; i++)
  {
    int i1 = (i < 2) ? (i + 1) : 0;

    TPoint<T> v = intersection - coords[i];
    TPoint<T> dv = coords[i1] - coords[i];
    TPoint<T> dir = dv ^ v;

    if (i > 0) 
    {
      if ((dir * olddir) < T(0.0)) 
      {
        inside = false;
        break;
      }
    }

    olddir = dir;
  }
                                      // additional check by projecting intersection on every edge
  bool onedge = false;
  if (!inside)
  {
    for (int i = 0; i < 3; i++)
    {
      int i1 = (i < 2) ? (i + 1) : 0;
      TPoint<T> p0 = coords[i];
      TPoint<T> p1 = coords[i1];

      // project intersection on every edge with tolerance
      TPoint<T> intr;
      T t = 0.0;
      if (projectPointOnSegment(intersection,p0,p1,&intr,&t,parmtolerance))
      {
        T dist = !(intersection - intr);
        if (dist < tolerance)
        {
          onedge = true;
          break;
        }
      }
    }
  }

  return inside || onedge;
}

/** Refine parameter U [0..1] value near ends, with powerstart and powerend.
  power 1.0 means no refinement, use 0.5 -> U^0.5 for rounded edge. */
template <typename T> T refineParameter(T U, T startpower = 1.0, T endpower = 1.0)
{
  T u = U * 2.0 - 1.0;
  if (u < 0.0)
  {
    u = sign(u) * pow(std::abs(u),startpower);
  } else
  {
    u = sign(u) * pow(std::abs(u),endpower);
  }
  U = (u + 1.0) * 0.5;

  LIMIT(U,0.0,1.0);

  return U;
}

/** Make B-spline knots. */
template <class T> void makeKnots(int K, int M, std::vector<T> &knots)
{
                              // fill knots, n is number of knots
  int n = K + M + 2;
                              // allocate
  knots.resize(n,0.0);

  knots[0] = 0.0;
  T d = 0.0;
  for (int i = 1; i < n; i++)
  {
    if ((i < (M + 1)) || (i > n - (M + 1)))
    {
      d = 0.0;
    } else
    {
      d = 1.0;
    }

    knots[i] = knots[i - 1] + d;
  }

  T Uend = knots[K + M + 1];

  for (int i = 0; i < n; i++)
  {
    knots[i] /= Uend;
  }
}

}
