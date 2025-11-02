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

  toperations.h

  Simplified operations on curves (for easy usage)

*******************************************************************************/

#pragma once

#include "tpoints.h"

#include "tbasecurve.h"

#include "tlsqsegment.h"
#include "torthosegment.h"
#include "tbeziersegment.h"

#include "tpointcurve.h"
#include "tbeziercurve.h"
#include "tsplinecurve.h"

namespace tcad {

//===== Auxiliary ==============================================================

/** Make coefficients with maximum 1.0 at maxposition [0..1], and decreasing to zero at ends.
  zero coefs are returned if power == 0.0. */
template <class T> void makeHat(int n, std::vector<T> &coefs, T power = 0.5, T maxposition = 0.5)
{
  coefs.clear();
  coefs.resize(n,0.0);

  int middle = int(n * maxposition);
  int middle2 = n - middle;

  // do not touch first and last points
  if (power != 0.0)
  {
    for (int i = 1; i < n - 1; i++)
    {
      T coef = (i <= middle) ? T(std::abs(i - middle)) / T(middle) :
        T(std::abs(i - middle)) / T(middle2);

      coef = 1.0 - coef;
      LIMIT(coef,0.0,1.0);

      coef = pow(coef,power);

      coefs[i] = coef;
    }
  }
}

/** f(0) = 1 and f(1) = 0 with zero derivatives at 0 and 1. */
template <class T> T makeHalfHat(T U)
{
  LIMIT(U,0.0,1.0);
  static TBezierSegment<T> s(TPoint<T>(0,1),TPoint<T>(0.666666,1),TPoint<T>(0.666666,0),TPoint<T>(1,0));
  return s.derivative(U,0).Y;
}

//===== Operations =============================================================

/** Remove duplicate points in a curve, only neighbour duplicates are excluded. */
template <class T> bool removeDuplicateNeighbours(std::vector<TPoint<T>> &points, T tolerance)
{
  // do not sort coordinates (false), exclude only neighbours
  bool ok = removeDuplicates(points,false,tolerance);
  return ok;
}

/** Remove ALL duplicate points like mesh nodes, points are sorted by XYZ. */
template <class T> bool removeDuplicateNodes(std::vector<TPoint<T>> &points, T tolerance)
{
  // sort coordinates, then exclude
  bool ok = removeDuplicates(points,true,tolerance);
  return ok;
}

/** Smooth curve points by orthogonal polynomials (very hard, very smooth). Do not use poly power
  higher than 10 as gamma function grows very quickly like factorial. You can create any curve 
  from the smoothed points. */
template <class T> void smoothPointsByOrtho(std::vector<TPoint<T>> &points, 
  CurveEndType start, CurveEndType end, int power = 8, int integration = GAUSSINT_8)
{
  // (1) prepare parameters by length to place new smoothed points to original 
  // parameteric positions
  std::vector<T> parms;
  prepareParameters(points,parms,true,false);

  // (2) create orthogonal polynomial from points, it smoothes points by poly, makes many new points
  TOrthoSegment<T> ssegment(points,start,end,power,integration);
  std::vector<TPoint<T>> spoints;
  ssegment.createPoints(spoints);

  // (3) make smoothed point curve, it has many points
  TPointCurve<T> smoothed(spoints);

  // (4) reparameterise (redivide) into a desired number of points, here we
  // take the original number of points
  std::vector<TPoint<T>> newpoints;
  for (int i = 0; i < int(parms.size()); i++)
  {
    TPoint<T> p = smoothed.derivative(parms[i],0);
    newpoints.push_back(p);
  }

  // (5) return new smoothed points
  points = newpoints;
}

/** Smooth curve points by a qubic Bezier segment (very hard, very smooth). You can create any curve 
  from the smoothed points. */
template <class T> void smoothPointsByBezier(std::vector<TPoint<T>> &points, 
  CurveEndType start, CurveEndType end)
{
  // (1) prepare parameters by length to place new smoothed points to original 
  // parameteric positions
  std::vector<T> parms;
  prepareParameters(points,parms,true,false);

  // (2) create a cubic Bezier segment from points
  TBezierSegment<T> ssegment(points,start,end);
  std::vector<TPoint<T>> spoints;
  ssegment.createPoints(spoints);

  // (3) make smoothed point curve, it has many points
  TPointCurve<T> smoothed(spoints);

  // (4) reparameterise (redivide) into a desired number of points, here we
  // take the original number of points
  std::vector<TPoint<T>> newpoints;
  for (int i = 0; i < int(parms.size()); i++)
  {
    TPoint<T> p = smoothed.derivative(parms[i],0);
    newpoints.push_back(p);
  }

  // (5) return new smoothed points
  points = newpoints;
}

/** Smooth curve points by a Bezier curve of numsegments. You can create any curve 
  from the smoothed points. */
template <class T> void smoothPointsByBezierCurve(std::vector<TPoint<T>> &points, int numsegments,
  CurveEndType start, CurveEndType end)
{
  assert(numsegments >= 1);

  // (1) prepare parameters by length to place new smoothed points to original 
  // parameteric positions
  std::vector<T> parms;
  prepareParameters(points,parms,true,false);

  // (2) create a cubic Bezier segment from points
  TBezierCurve<T> scurve(points,numsegments,start,end);
  std::vector<TPoint<T>> spoints;
  scurve.createPoints(spoints);

  // (3) make smoothed point curve, it has many points
  TPointCurve<T> smoothed(spoints);

  // (4) reparameterise (redivide) into a desired number of points, here we
  // take the original number of points
  std::vector<TPoint<T>> newpoints;
  for (int i = 0; i < int(parms.size()); i++)
  {
    TPoint<T> p = smoothed.derivative(parms[i],0);
    newpoints.push_back(p);
  }

  // (5) return new smoothed points
  points = newpoints;
}

/** Smooth curve points by spline of k intervals. You can create any curve 
  from the smoothed points. If clamped... spline start and end directions are taken from
  points, otherwise a natural spline is built (second derivative of XYZ on parameter
  is zero) at the end. */
template <class T> void smoothPointsBySplineCurve(std::vector<TPoint<T>> &points, int k, 
  CurveEndType start, CurveEndType end, int degree = SPLINE_DEGREE)
{
  assert(k >= 1);

  // (1) prepare parameters by length to place new smoothed points to original 
  // parameteric positions
  std::vector<T> parms;
  prepareParameters(points,parms,true,false);

  // (2) create a spline from points
  TSplineCurve<T> scurve(points,k,start,end,degree);
  std::vector<TPoint<T>> spoints;
  scurve.createPoints(spoints);

  // (3) make smoothed point curve, it has many points
  TPointCurve<T> smoothed(spoints);

  // (4) reparameterise (redivide) into a desired number of points, here we
  // take the original number of points
  std::vector<TPoint<T>> newpoints;
  for (int i = 0; i < int(parms.size()); i++)
  {
    TPoint<T> p = smoothed.derivative(parms[i],0);
    newpoints.push_back(p);
  }

  // (5) return new smoothed points
  points = newpoints;
}

/** Make NACA airfoil wing surface. */
template <class T> void makeNACASurface(std::vector<std::vector<TPoint<T>>> &points,
  int numspans, T resizeX = 0.99, T resizeY = 0.99, T moveZ = 0.02, T dadegZ = -0.5, T dadegY = -0.5)
{
  points.clear();

  // rescale to -0.5, +0.5]
  std::vector<TPoint<T>> airfoil;
  for (auto p : NACA0012xy<T>)
  {
    airfoil.push_back(TPoint<T>(p.first - 0.5,p.second));
  }

  points.push_back(airfoil);

  for (int i = 0; i < numspans - 1; i++)
  {
    TTransform<T> t;
    t.Resize(TPoint<T>(resizeX,resizeY,1.0));
    t.Rotate(TPoint<T>(0.0,0.0,1.0),dadegZ * CPI);
    t.Translate(TPoint<T>(0.0,0.0,moveZ));
    t.Rotate(TPoint<T>(0.0,1.0,0.0),dadegY * CPI);

    for (int j = 0; j < int(airfoil.size()); j++)
    {
      airfoil[j] = t.applyTransform(airfoil[j]);
    }

    points.push_back(airfoil);
  }
}

/** Make cylinder points with elliptical cross-sections in XY plane with axis along Z. */
template <class T> void makeCylinder(int numsections, T Zmin, T Zmax, int numpoints, T a, T b, 
  std::vector<std::vector<TPoint<T>>> &points, T adegfrom = 0.0, T adegto = 180.0)
{
  T dZ = (Zmax - Zmin) / T(numsections - 1);
  for (int i = 0; i < numsections; i++)
  {
    T Z = Zmin + dZ * T(i);

    std::vector<TPoint<T>> section;
    makeEllipseXY(numpoints,TPoint<T>(),a,b,section,adegfrom,adegto);

    TTransform<T> t;
    t.Translate(TPoint<T>(0.0,0.0,Z));
    makeTransform(section,&t);

    points.push_back(section);
  }
}

}