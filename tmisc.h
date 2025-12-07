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
  TPoint<T> dp = p1 - p0;
  TPoint<T> dv = v1 - v0;

  if (!dp < TOLERANCE(T) || !dv < TOLERANCE(T))
    return false;

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
  assert(index >= 0 && index < int(points[0].size()));
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
  if (U < TOLERANCE(T))
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
//!!!!!!!        if (dist < tolerance)
        if (dist < TOLERANCE(T))
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

/** Triangle area. Always non-negative. */
template <class T> T triangleArea(TPoint<T> p1, TPoint<T> p2, TPoint<T> p3)
{
  T A,B,C,P;
  TPoint<T> p12,p23,p31;

  p12 = p2 - p1; p23 = p3 - p2; p31 = p3 - p1;
  A = !p12; B = !p23; C = !p31;

  // semiperimeter
  P = (A + B + C) * 0.5;

  // Geron's formula
  T area = P * (P - A) * (P - B) * (P - C); 
  if (area < 0.0) 
    area = 0.0;

  return sqrt(area);
}

/** Get three barythencric coordinate of a point inside triangle. The point MUST BE inside or close
  to the edge. */
template <class T> TPoint<T> barycentricCoord(std::array<TPoint<T>,3> &coords, TPoint<T> &intersection)
{
  // whole triangle area
  T area = triangleArea(coords[0],coords[1],coords[2]);

  if (area > TOLERANCE(T))
  {
    // first node 0 coordinate, zero at side 1-2
    T area0 = triangleArea(intersection,coords[1],coords[2]);

    // node 1 coordinate, zero at side 2-0
    T area1 = triangleArea(intersection,coords[2],coords[0]);

    T e0 = area0 / area;
    T e1 = area1 / area;

    LIMIT(e0,0.0,1.0);
    LIMIT(e1,0.0,1.0);
    T e2 = 1.0 - e0 - e1;
    LIMIT(e2,0.0,1.0);

    TPoint<T> ecoord(e0,e1,e2);

    return ecoord;
  } else
  {
    // what is not correct it designates failure
    return TPoint<T>();
  }
}

/** Convert 1D points size (K1 + 1) x (K2 + 1) into 2D points[K2 + 1][K1 + 1]. */
template <class T> void points1Dto2D(std::vector<TPoint<T>> &points1D, int K1, int K2,
  std::vector<std::vector<TPoint<T>>> &points2D)
{
  assert(int(points1D.size()) == (K1 + 1) * (K2 + 1));

  points2D.clear();
  int count = 0;
  for (int i = 0; i <= K2; i++)
  {
    points2D.push_back(std::vector<TPoint<T>>());
    for (int j = 0; j <= K1; j++)
    {
      points2D.back().push_back(points1D[count++]);
    }
  }
}

/** Convert 2D points[K2 + 1][K1 + 1] into 1D points size (K1 + 1) x (K2 + 1). */
template <class T> void points2Dto1D(std::vector<std::vector<TPoint<T>>> &points2D,
  std::vector<TPoint<T>> &points1D, int *K1 = nullptr, int *K2 = nullptr)
{
  points1D.clear();
  for (int i = 0; i < int(points2D.size()); i++)
  {
    for (int j = 0; j < int(points2D[i].size()); j++)
    {
      points1D.push_back(points2D[i][j]);
    }
  }

  if (K1)
    *K1 = int(points2D[0].size()) - 1;
  if (K2)
    *K2 = int(points2D.size()) - 1;
}

/** Boundary parameters for surface are measured from node 0 (U = V = 0.0)
  counter-clockwise, max value being 4.0; index being boundary piece number 0..3 
  correspondingly. 

        parm = 3.0      index = 2       parm = 2.0
            3------------------------------2
            |                              |
            |                              |
            |                              |
index = 3   |                              | index = 1
            ^ V                            |
            |    U                         |
parm = 4.0  0---->-------------------------1
        parm = 0.0      index = 0      parm = 1.0

*/
template <class T> bool UVToBoundaryParm(T U, T V, T &parm, T parmtolerance = PARM_TOLERANCE)
{
  if (std::abs(V - 0.0) < parmtolerance)
  {
    parm = U;
    return true;
  } else if (std::abs(U - 1.0) < parmtolerance)
  {
    parm = 1.0 + V;
    return true;
  } else if (std::abs(V - 1.0) < parmtolerance)
  {
    parm = 2.0 + (1.0 - U);
    return true;
  } else if (std::abs(U - 0.0) < parmtolerance)
  {
    parm = 3.0 + (1.0 - V);
    return true;
  } else
  {
    // unreal
    parm = -1.0;
    return false;
  }
}

/** Convert boundary parameter [0..4] to U,V [0..1]. */
template <class T> TPoint<T> boundaryParmToUV(T parm)
{
  // integer part is side number
  int side = int(parm);

  // fractional part is parameter inside the side
  T r = parm - T(side);
  LIMIT(r,0.0,1.0);

  // side can only be 0..3
  side = side % 4;
  if (side == 0)
  {
    return TPoint<T>(r,0.0);
  } else if (side == 1)
  {
    return TPoint<T>(1.0,r);
  } else if (side == 2)
  {
    return TPoint<T>(1.0 - r,1.0);
  } else if (side == 3)
  {
    return TPoint<T>(0.0,1.0 - r);
  } else
  {
    return TPoint<T>();
  }
}

/** Get parameter on the opposite side. parm being [0..4] */
template <class T> T oppositeParm(T parm)
{
  int i = (int) parm;
  T r = parm - T(i);
  T oparm = T(i) + 2.0;
  oparm += (1.0 - r);
  if (oparm >= 4.0) oparm -= 4.0;
  return oparm;
}

/** Get next parameter value when closing the intersection curve. */
template <class T> bool nextBoundaryParm(T parm, T endparm, T &nextparm,
  T parmtolerance = PARM_TOLERANCE)
{
  // ended?
  T dist = std::abs(endparm - parm);
  if (dist < parmtolerance)
    return false;

  int side = int(parm);
  T r = parm - T(side);

  if (std::abs(r - 0.0) < parmtolerance || std::abs(r - 1.0) < parmtolerance)
  {
    nextparm = parm + 1.0;
  } else
  {
    nextparm = parm + (1.0 - r);
  }

  if (endparm >= parm && endparm <= nextparm)
  {
    nextparm = endparm;
    return true;
  }

  if (nextparm >= 4.0)
    nextparm -= 4.0;

  if (endparm >= parm && endparm <= nextparm)
    nextparm = endparm;

  return true;
}

/** Extend min/max box by coef. */
template <class T> void extendMinMax(TPoint<T> &min, TPoint<T> &max, T coef = 1.0)
{
  TPoint<T> c = (min + max) * 0.5;
  TPoint<T> d = (max - min) * 0.5 * coef;
  min = c - d;
  max = c + d;
}

/** Extend min/max box by coef. */
template <class T> void extendBox(std::pair<TPoint<T>,TPoint<T>> &box, T coef = 1.0)
{
  extendMinMax(box.first,box.second,coef);
}

/** Boxes intersect? */
template <class T> bool boxesOverlap(std::pair<TPoint<T>,TPoint<T>> &box0,
  std::pair<TPoint<T>,TPoint<T>> &box1)
{
  bool overlap = 
    box0.first.X < box1.second.X &&
    box1.first.X < box0.second.X &&
    box0.first.Y < box1.second.Y &&
    box1.first.Y < box0.second.Y &&
    box0.first.Z < box1.second.Z &&
    box1.first.Z < box0.second.Z;

  return overlap;
}

/** Generates ranges array of numranges + 1 elements. */
template <class T> void getRanges(T max, int numranges, std::vector<T> &ranges)
{
  T step = max / numranges;
  LIMIT_MIN(step, 1);

  ranges.resize(numranges + 1);
  T value = T(0);

  for (size_t i = 0; i <= numranges; i++)
  {
    if (value >= max)
      value = max;

    ranges[i] = value;

    value += step;
  }

  ranges[numranges] = max;
}

/** Bezier basis. */
template <class T> void BezierBasis(T U, std::vector<T> &knots, std::vector<T> &basis, 
  int &segment) 
{
  assert(knots.size() > 0);
  static const int degree = 3;

  // number of knots (num segments plus 1)
  int numKnots = static_cast<int>(knots.size());

  // allocate and clear space for basis functions
  int numPoints = (numKnots - 1) * degree + 1;
  basis.resize(numPoints,0.0);

  // tolerance
  T tolerance = TOLERANCE(T);

  // find knots interval (segment) for U by bisection
  segment = findInterval<T>(knots,U);
  assert(segment != -1 && "Unable to find parametric interval for parameter");

  // get local parameter within segment [0..1]
  T d = knots[segment + 1] - knots[segment];
  assert(d > tolerance);

  T u = 0;
  if (d > tolerance)
  {
    u = (U - knots[segment]) / d;
  }

  // starting non-zero basis function point
  int i0 = segment * degree;

  T u1 = 1.0 - u;
  basis[i0] = u1 * u1 * u1;
  basis[i0 + 1] = 3.0 * u * u1 * u1;
  basis[i0 + 2] = 3.0 * u * u * u1;
  basis[i0 + 3] = u * u * u;
}

/** Bezier basis, first derivative. Derivatives are scaled by the whole length. */
template <class T> void BezierBasisDer1(T U, std::vector<T> &knots, std::vector<T> &basis, 
  int &segment) 
{
  assert(knots.size() > 0);
  static const int degree = 3;

  // number of knots (num segments plus 1)
  int numKnots = static_cast<int>(knots.size());

  // allocate and clear space for basis functions
  int numPoints = (numKnots - 1) * degree + 1;
  basis.resize(numPoints,0.0);

  // tolerance
  T tolerance = TOLERANCE(T);

  // find knots interval (segment) for U by bisection
  segment = findInterval<T>(knots,U);
  assert(segment != -1 && "Unable to find parametric interval for parameter");

  // get local parameter within segment [0..1]
  T d = knots[segment + 1] - knots[segment];
  assert(d > tolerance);

  T u = 0;
  if (d > tolerance)
  {
    u = (U - knots[segment]) / d;
  }

  // starting non-zero basis function point
  int i0 = segment * degree;

  T u1 = 1.0 - u;
  T u2 = u * u;

  // multiply by d to scale into whole length
  basis[i0] = -3.0 * u1 * u1;
  basis[i0 + 1] = 3.0 - 12.0 * u + 9.0 * u2;
  basis[i0 + 2] = 6.0 * u - 9.0 * u2;
  basis[i0 + 3] = 3.0 * u2;

  if (d > tolerance)
  {
    basis[i0] /= d;
    basis[i0 + 1] /= d;
    basis[i0 + 2] /= d;
    basis[i0 + 3] /= d;
  }
}

/** Bezier basis, second derivative. Derivatives are scaled by the whole length. */
template <class T> void BezierBasisDer2(T U, std::vector<T> &knots, std::vector<T> &basis, 
  int &segment) 
{
  assert(knots.size() > 0);
  static const int degree = 3;

  // number of knots (num segments plus 1)
  int numKnots = static_cast<int>(knots.size());

  // allocate and clear space for basis functions
  int numPoints = (numKnots - 1) * degree + 1;
  basis.resize(numPoints,0.0);

  // tolerance
  T tolerance = TOLERANCE(T);

  // find knots interval (segment) for U by bisection
  segment = findInterval<T>(knots,U);
  assert(segment != -1 && "Unable to find parametric interval for parameter");

  // get local parameter within segment [0..1]
  T d = knots[segment + 1] - knots[segment];
  assert(d > tolerance);

  T u = 0;
  if (d > tolerance)
  {
    u = (U - knots[segment]) / d;
  }

  // starting non-zero basis function point
  int i0 = segment * degree;

  T u1 = 1.0 - u;
  T u2 = u * u;

  // multiply by d to scale into whole length
  basis[i0] = 6.0 * u1;
  basis[i0 + 1] = -12.0 + 18.0 * u;
  basis[i0 + 2] = 6.0 - 18.0 * u;
  basis[i0 + 3] = 6.0 * u;

  T d2 = d * d;

  // not sure that the vector length is right
  if (d2 > tolerance)
  {
    basis[i0] /= d2;
    basis[i0 + 1] /= d2;
    basis[i0 + 2] /= d2;
    basis[i0 + 3] /= d2;
  }
}

/** Location index for 2D regular mesh of points. */
template <class T> int getIndex(int K1, int K2, int i, int j)
{
  assert (i < (K1 + 1));
  assert (j < (K2 + 1));

  return j * (K1 + 1) + i;
}

/** Location index for 3D regular mesh of points. */
template <class T> int getIndex(int K1, int K2, int K3, int i, int j, int k)
{
  assert (i < (K1 + 1));
  assert (j < (K2 + 1));
  assert (k < (K3 + 1));

  return k * (K2 + 1) * (K1 + 1) + j * (K1 + 1) + i;
}

/** Get row of control points (numbers) for 3D regular mesh of points. */
template <class T> void getULinePoints(int K1, int K2, int K3, int j, int k, std::vector<int> &points)
{
  int index = tcad::getIndex<T>(K1,K2,K3,0,j,k);
  for (int l = 0; l < (K1 + 1); l++)
  {
    points.push_back(index);
    index += 1;
  }
}

/** Get column of control points (numbers) for 3D regular mesh of points. */
template <class T> void getVLinePoints(int K1, int K2, int K3, int i, int k, std::vector<int> &points)
{
  int index = tcad::getIndex<T>(K1,K2,K3,i,0,k);
  for (int l = 0; l < (K2 + 1); l++)
  {
    points.push_back(index);
    index += (K1 + 1);
  }
}

/** Get layer of control points (numbers) for 3D regular mesh of points. */
template <class T> void getWLinePoints(int K1, int K2, int K3, int i, int j, std::vector<int> &points)
{
  int numUV = (K1 + 1) * (K2 + 1);
  int index = tcad::getIndex<T>(K1,K2,K3,i,j,0);
  for (int l = 0; l < (K3 + 1); l++)
  {
    points.push_back(index);
    index += numUV;
  }
}

/** Get row of control points for 3D regular mesh of points. 
  cpoints have size (K1 + 1) * (K2 + 1) * (K3 + 1). */
template <class T> void getRow(std::vector<TPoint<T>> &cpoints, int K1, int K2, int K3, 
  int j, int k, std::vector<TPoint<T>> &points)
{
  std::vector<int> ipoints;
  tcad::getULinePoints<T>(K1,K2,K3,j,k,ipoints);

  points.clear();
  for (int i = 0; i < int(ipoints.size()); i++)
  {
    points.push_back(cpoints[ipoints[i]]);
  }

  assert(int(points.size()) == K1 + 1);
}

/** Set row of control points for 3D regular mesh of points.
  cpoints have size (K1 + 1) * (K2 + 1) * (K3 + 1). */
template <class T> void setRow(std::vector<TPoint<T>> &cpoints, int K1, int K2, int K3, 
  int j, int k, std::vector<TPoint<T>> &points)
{
  assert(int(points.size()) == K1 + 1);

  std::vector<int> ipoints;
  tcad::getULinePoints<T>(K1,K2,K3,j,k,ipoints);

  for (int i = 0; i < int(ipoints.size()); i++)
  {
    cpoints[ipoints[i]] = points[i];
  }
}

/** Get column of control points for 3D regular mesh of points.
  cpoints have size (K1 + 1) * (K2 + 1) * (K3 + 1). */
template <class T> void getColumn(std::vector<TPoint<T>> &cpoints, int K1, int K2, int K3, 
  int i, int k, std::vector<TPoint<T>> &points)
{
  std::vector<int> ipoints;
  tcad::getVLinePoints<T>(K1,K2,K3,i,k,ipoints);

  points.clear();
  for (int i = 0; i < int(ipoints.size()); i++)
  {
    points.push_back(cpoints[ipoints[i]]);
  }

  assert(int(points.size()) == K2 + 1);
}

/** Set column of control points for 3D regular mesh of points. 
  cpoints have size (K1 + 1) * (K2 + 1) * (K3 + 1). */
template <class T> void setColumn(std::vector<TPoint<T>> &cpoints, int K1, int K2, int K3, 
  int i, int k, std::vector<TPoint<T>> &points)
{
  assert(int(points.size()) == K2 + 1);

  std::vector<int> ipoints;
  tcad::getVLinePoints<T>(K1,K2,K3,i,k,ipoints);

  for (int i = 0; i < int(ipoints.size()); i++)
  {
    cpoints[ipoints[i]] = points[i];
  }
}

/** Get layer of control points for 3D regular mesh of points.
  cpoints have size (K1 + 1) * (K2 + 1) * (K3 + 1). */
template <class T> void getLayer(std::vector<TPoint<T>> &cpoints, int K1, int K2, int K3, 
  int i, int j, std::vector<TPoint<T>> &points)
{
  std::vector<int> ipoints;
  tcad::getWLinePoints<T>(K1,K2,K3,i,j,ipoints);

  points.clear();
  for (int i = 0; i < int(ipoints.size()); i++)
  {
    points.push_back(cpoints[ipoints[i]]);
  }

  assert(int(points.size()) == K3 + 1);
}

/** Set layer of control points for 3D regular mesh of points.
  cpoints have size (K1 + 1) * (K2 + 1) * (K3 + 1). */
template <class T> void setLayer(std::vector<TPoint<T>> &cpoints, int K1, int K2, int K3, 
  int i, int j, std::vector<TPoint<T>> &points)
{
  assert(int(points.size()) == K3 + 1);

  std::vector<int> ipoints;
  tcad::getWLinePoints<T>(K1,K2,K3,i,j,ipoints);

  for (int i = 0; i < int(ipoints.size()); i++)
  {
    cpoints[ipoints[i]] = points[i];
  }
}

/** Remove first and last element. */
template <class T> void removeFirstLast(std::vector<T> &knots, std::vector<T> &reduced)
{
  reduced = knots;

  reduced.erase(reduced.begin());
  reduced.erase(reduced.end() - 1);
}

/**
  Faces (all normals point outside) :
              7-------------------6 
             /|        /         /|    
            / |      (5)->      / |     
           /  |        |       /  |          
          4---------<-(3)-----5   |    
          |   |               |   |     
          |(0)|     ^         |(1)|    
          | | |     |         ||/ |
          |/  3----(2)->------|---2 
          |  /                |  / 
  WZk     ^ /VYj       (4)->  | /
  layers  |/ columns   /      |/
          0-->----------------1
          UXi rows

  cpoints are ALL control points, (K1 + 1) * (K2 + 1) * (K3 + 1)
*/

template <class T> void getFace(std::vector<TPoint<T>> &cpoints, int faceno, 
  int K1, int K2, int K3, std::vector<std::vector<TPoint<T>>> &points)
{
  points.clear();

  if (faceno == 0)
  {
    for (int k = 0; k <= K3; k++)
    {
      std::vector<TPoint<T>> temp;
      getColumn<T>(cpoints,K1,K2,K3,0,k,temp);
      std::reverse(temp.begin(),temp.end());
      points.push_back(temp);
    }
  } else if (faceno == 1)
  {
    for (int k = 0; k <= K3; k++)
    {
      std::vector<TPoint<T>> temp;
      getColumn<T>(cpoints,K1,K2,K3,K1,k,temp);
      points.push_back(temp);
    }
  } else if (faceno == 2)
  {
    for (int k = 0; k <= K3; k++)
    {
      std::vector<TPoint<T>> temp;
      getRow<T>(cpoints,K1,K2,K3,0,k,temp);
      points.push_back(temp);
    }
  } else if (faceno == 3)
  {
    for (int k = 0; k <= K3; k++)
    {
      std::vector<TPoint<T>> temp;
      getRow<T>(cpoints,K1,K2,K3,K2,k,temp);
      std::reverse(temp.begin(),temp.end());
      points.push_back(temp);
    }
  } else if (faceno == 4)
  {
    for (int j = K2; j >= 0; j--)
    {
      std::vector<TPoint<T>> temp;
      getRow<T>(cpoints,K1,K2,K3,j,0,temp);
      points.push_back(temp);
    }
  } else if (faceno == 5)
  {
    for (int j = 0; j <= K2; j++)
    {
      std::vector<TPoint<T>> temp;
      getRow<T>(cpoints,K1,K2,K3,j,K3,temp);
      points.push_back(temp);
    }
  }
}

/** Approximate size along U. */
template <class T> T Usize(std::vector<std::vector<TPoint<T>>> &points)
{
  int K1 = int(points[0].size() - 1);
  int K2 = int(points.size() - 1);

  T len = 0.0;
  for (int i = 0; i <= K2; i++)
  {
    std::vector<TPoint<T>> row;
    getRow(points,i,row);
    len += calculateLength(row);
  }

  len /= T(K2 + 1);

  return len;
}

/** Approximate size along V. */
template <class T> T Vsize(std::vector<std::vector<TPoint<T>>> &points)
{
  int K1 = int(points[0].size() - 1);
  int K2 = int(points.size() - 1);

  T len = 0.0;
  for (int i = 0; i <= K1; i++)
  {
    std::vector<TPoint<T>> col;
    getColumn(points,i,col);
    len += calculateLength(col);
  }

  len /= T(K1 + 1);

  return len;
}

/** Reverse points rows. */
template <class T> void reverseRows(std::vector<std::vector<TPoint<T>>> &points)
{
  int K1 = int(points[0].size() - 1);
  int K2 = int(points.size() - 1);

  for (int i = 0; i <= K2; i++)
  {
    std::vector<TPoint<T>> temp;
    getRow(points,i,temp);
    std::reverse(temp.begin(),temp.end());
    setRow(points,i,temp);
  }
}

/** Reverse points columns. */
template <class T> void reverseColumns(std::vector<std::vector<TPoint<T>>> &points)
{
  int K1 = int(points[0].size() - 1);
  int K2 = int(points.size() - 1);

  for (int i = 0; i <= K1; i++)
  {
    std::vector<TPoint<T>> temp;
    getColumn(points,i,temp);
    std::reverse(temp.begin(),temp.end());
    setColumn(points,i,temp);
  }
}

}
