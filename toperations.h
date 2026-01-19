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

#include "tpointsurface.h"
#include "tsplinesurface.h"

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
  CurveEndType start, CurveEndType end, int power = 8, int integration = GAUSSINT_8, 
  bool keependpoints = true, int numnewpoints = -1, T refinestartU = 1.0, T refineendU = 1.0)
{
  // (1) prepare parameters by length to place new smoothed points to original 
  // parameteric positions
  std::vector<T> parms;
  if (numnewpoints > 1)
  {
    prepareUniformParameters(numnewpoints,parms,refinestartU,refineendU);
  } else
  {
    prepareParameters(points,parms,true,false);
  }

  // (2) create orthogonal polynomial from points, it smoothes points by poly, makes many new points
  TOrthoSegment<T> ssegment(points,start,end,power,integration);
  std::vector<TPoint<T>> spoints;
  ssegment.createPoints(spoints,(numnewpoints > 1) ? numnewpoints : MANY_POINTS);

  std::vector<TPoint<T>> newpoints;

  // (3) make smoothed point curve, it has many points
  TPointCurve<T> smoothed(spoints);

  // (4) reparameterise (redivide) into a desired number of points, here we
  // take the original number of points
  for (int i = 0; i < int(parms.size()); i++)
  {
    TPoint<T> p = smoothed.derivative(parms[i],0);
    newpoints.push_back(p);
  }

  if (keependpoints)
  {
    newpoints.front() = points.front();
    newpoints.back() = points.back();
  }

  // (5) return new smoothed points
  points = newpoints;
}

/** Smooth curve points by a qubic Bezier segment (very hard, very smooth). You can create any curve 
  from the smoothed points. */
template <class T> void smoothPointsByBezier(std::vector<TPoint<T>> &points, 
  CurveEndType start, CurveEndType end, int numnewpoints = -1)
{
  // (1) prepare parameters by length to place new smoothed points to original 
  // parameteric positions
  std::vector<T> parms;
  if (numnewpoints > 1)
  {
    prepareUniformParameters(numnewpoints,parms);
  } else
  {
    prepareParameters(points,parms,true,false);
  }

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

//!!!!!!
  //if (start == END_CLAMPED)
  //{
  //  TPoint<T> dir = +(startDirection(points));
  //  T len = !(newpoints[1] - newpoints[0]);
  //  newpoints[1] = newpoints[0] + dir * len;
  //}

  //if (end == END_CLAMPED)
  //{
  //  TPoint<T> dir = +(endDirection(points));
  //  T len = !(newpoints[newpoints.size() - 1] - newpoints[newpoints.size() - 2]);
  //  newpoints[newpoints.size() - 2] = newpoints[newpoints.size() - 1] + dir * len;
  //}

  // (5) return new smoothed points
  points = newpoints;
}

/** Smooth curve points by a Bezier curve of numsegments. You can create any curve 
  from the smoothed points. */
template <class T> void smoothPointsByBezierCurve(std::vector<TPoint<T>> &points, int numsegments,
  CurveEndType start, CurveEndType end, int numnewpoints = -1)
{
  assert(numsegments >= 1);

  // (1) prepare parameters by length to place new smoothed points to original 
  // parameteric positions
  std::vector<T> parms;
  if (numnewpoints > 1)
  {
    prepareUniformParameters(numnewpoints,parms);
  } else
  {
    prepareParameters(points,parms,true,false);
  }

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

//!!!!!!
  //if (start == END_CLAMPED)
  //{
  //  TPoint<T> dir = +(startDirection(points));
  //  T len = !(newpoints[1] - newpoints[0]);
  //  newpoints[1] = newpoints[0] + dir * len;
  //}

  //if (end == END_CLAMPED)
  //{
  //  TPoint<T> dir = +(endDirection(points));
  //  T len = !(newpoints[newpoints.size() - 1] - newpoints[newpoints.size() - 2]);
  //  newpoints[newpoints.size() - 2] = newpoints[newpoints.size() - 1] + dir * len;
  //}

  // (5) return new smoothed points
  points = newpoints;
}

/** Smooth curve points by spline of k intervals. You can create any curve 
  from the smoothed points. If clamped... spline start and end directions are taken from
  points, otherwise a natural spline is built (second derivative of XYZ on parameter
  is zero) at the end. */
template <class T> void smoothPointsBySplineCurve(std::vector<TPoint<T>> &points, int k, 
  CurveEndType start, CurveEndType end, int degree = SPLINE_DEGREE, int numnewpoints = -1)
{
  assert(k >= 1);

  // (1) prepare parameters by length to place new smoothed points to original 
  // parameteric positions
  std::vector<T> parms;
  if (numnewpoints > 1)
  {
    prepareUniformParameters(numnewpoints,parms);
  } else
  {
    prepareParameters(points,parms,true,false);
  }

  // (2) create a spline from points
  TSplineCurve<T> scurve(points,k,degree,start,end);
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

/** Transform list of surfaces. */
template <class T> void makeTransform(std::vector<TSplineSurface<T> *> &surfaces, TTransform<T> *t)
{
  for (int i = 0; i < int(surfaces.size()); i++)
  {
    surfaces[i]->makeTransform(t);
  }
}

/** Delete surfaces. */
template <class T> void deleteSurfaces(std::vector<TSplineSurface<T> *> &surfaces)
{
  for (int i = 0; i < int(surfaces.size()); i++)
  {
    DELETE_CLASS(surfaces[i]);
  }

  surfaces.clear();
}

/** Flip (reverse) surfaces. */
template <class T> void reverseSurfaces(std::vector<TSplineSurface<T> *> &surfaces)
{
  for (int i = 0; i < int(surfaces.size()); i++)
  {
    surfaces[i]->reverseU();
  }
}

/** Copy and transform list of surfaces. */
template <class T> void copySurfaces(
  // input
  std::vector<TSplineSurface<T> *> &surfaces,
  // output 
  std::vector<TSplineSurface<T> *> &newsurfaces,
  // parms
  TTransform<T> *t = nullptr, bool reverseU = false, bool reverseV = false)
{
  for (int i = 0; i < int(surfaces.size()); i++)
  {
    // copy
    TSplineSurface<T> *newsurface = new TSplineSurface<T>(*surfaces[i]);

    // transform
    if (t)
      newsurface->makeTransform(t);

    // reverseU/V
    if (reverseU)
      newsurface->reverseU();
    if (reverseV)
      newsurface->reverseV();

    newsurfaces.push_back(newsurface);
  }
}

/** Make surface duplicates to make no more than one loop to keep Rhino happy. */
template <class T> void makeSingleLoop(std::vector<TSplineSurface<T> *> &surfaces, 
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV)
{
  int n = int(surfaces.size());
  for (int i = 0; i < n; i++)
  {
    if (boundariesUV[i].size() > 1)
    {
      for (int j = 1; j < int(boundariesUV[i].size()); j++)
      {
        TSplineSurface<T> *copy = new TSplineSurface<T>(*surfaces[i]);
        surfaces.push_back(copy);
        boundariesUV.push_back(std::vector<std::vector<std::vector<tcad::TPoint<T>>>>());
        boundariesUV.back().push_back(boundariesUV[i][j]);
      }
      
      boundariesUV[i].resize(1);
    }
  }
}

/** Get piece of boundary curve in XYZ from parametric boundariesUV. */
template <class T> void getBoundaryPartXYZ(std::vector<tcad::TSplineSurface<T> *> &surfaces, 
//surface     loop        part        points
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV,
  int surface, int loop, int part, 
  std::vector<TPoint<T>> &points)
{
  points.clear();

  // this curve should be XYZ, not UV
  for (auto &UV : boundariesUV[surface][loop][part])
  {
    points.push_back(surfaces[surface]->position(UV.X,UV.Y));
  }
}

/** Get loop index in XYZ. */
template <class T> void getLoopXYZ(std::vector<tcad::TSplineSurface<T> *> &surfaces, 
//surface     loop        part        points
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV,
  int surface, int loopindex, std::vector<std::vector<TPoint<T>>> &loop)
{
  loop.clear();

  for (int i = 0; i < int(boundariesUV[surface][loopindex].size()); i++)
  {
    std::vector<TPoint<T>> points;
    getBoundaryPartXYZ(surfaces,boundariesUV,surface,loopindex,i,points);

    loop.push_back(points);
  }
}

/** Prepare outer loops (index 0) in XYZ coorinates. If a loop is empty, it generates
  a default outer loop. */
template <class T> void prepareLoopsXYZ(
  std::vector<TSplineSurface<T> *> &surfaces,
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV,
  std::vector<std::vector<std::vector<TPoint<T>>>> &loops)
{
  loops.clear();

  assert(surfaces.size() == boundariesUV.size());

  for (int i = 0; i < int(boundariesUV.size()); i++)
  {
    if (boundariesUV[i].empty())
    {
      // create default outer boundary
      closeOuterBoundary(boundariesUV[i]);
      std::vector<std::vector<TPoint<T>>> loop;
      getLoopXYZ(surfaces,boundariesUV,i,0,loop); // 0 - outer loop

      loops.push_back(loop);

      // restore status quo
      boundariesUV[i].clear();
    } else
    {
      for (int k = 0; k < int(boundariesUV[i].size()); k++) // all loops
      {
        std::vector<std::vector<TPoint<T>>> loop;
        getLoopXYZ(surfaces,boundariesUV,i,k,loop); 

        loops.push_back(loop);
      }
    }
  }
}

/** Prepare outer loops (index 0) in XYZ coorinates. If a loop is empty, it generates
  a default outer loop. The same as previous but output as a list of boundary curves in 
  TPointCurve<T>. */
template <class T> void prepareLoopsXYZ(
  std::vector<TSplineSurface<T> *> &surfaces,
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV,
  std::vector<std::vector<TPointCurve<T>>> &curves)
{
  curves.clear();

  assert(surfaces.size() == boundariesUV.size());

  for (int i = 0; i < int(boundariesUV.size()); i++)
  {
    curves.push_back(std::vector<TPointCurve<T>>());

    if (boundariesUV[i].empty())
    {
      // create default outer boundary
      closeOuterBoundary(boundariesUV[i]);
      std::vector<std::vector<TPoint<T>>> loop;
      getLoopXYZ(surfaces,boundariesUV,i,0,loop); // 0 - outer loop

      for (int j = 0; j < int(loop.size()); j++)
      {
        TPointCurve<T> curve(loop[j]);
        curves.back().push_back(curve);
      }

      // restore status quo
      boundariesUV[i].clear();
    } else
    {
      for (int k = 0; k < int(boundariesUV[i].size()); k++) // all loops
      {
        std::vector<std::vector<TPoint<T>>> loop;
        getLoopXYZ(surfaces,boundariesUV,i,k,loop); 

        for (int j = 0; j < int(loop.size()); j++)
        {
          TPointCurve<T> curve(loop[j]);
          curves.back().push_back(curve);
        }
      }
    }
  }
}

/** Approximate surfaces min/max. */
template <class T> std::pair<TPoint<T>,TPoint<T>> getMinMax(std::vector<TSplineSurface<T> *> &surfaces)
{
  std::pair<TPoint<T>,TPoint<T>> result;

  for (int i = 0; i < int(surfaces.size()); i++)
  {
    std::pair<TPoint<T>,TPoint<T>> minmax = surfaces[i]->getMinMax();
    if (i == 0)
    {
      result = minmax;
    } else
    {
      result.first = pointMin<T>(result.first,minmax.first);
      result.second = pointMax<T>(result.second,minmax.second);
    }
  }

  return result;
}

/** Approximate surfaces size. */
template <class T> T surfacesSize(std::vector<TSplineSurface<T> *> &surfaces)
{
  std::pair<TPoint<T>,TPoint<T>> minmax = getMinMax(surfaces);
  return !(minmax.second - minmax.first);
}

/** Close face loops for intersected faces by boundaries. */
template <class T> void closeFaceLoops(std::vector<TSplineSurface<T> *> &surfaces,
  // surface  // loop 0   // 4 pieces // piece contents
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV,
  std::vector<std::vector<std::vector<std::vector<TPoint<T>>>>> &boundaries,
  std::vector<int> &exlist,
  T maxseglen, T parmtolerance = PARM_TOLERANCE, int manypoints = MANY_POINTS2D)
{
  exlist.clear();

  for (int i = 0; i < int(surfaces.size()); i++)
  {
    if (!boundaries[i].empty())
    {
//!!! do not do      boundariesUV[i].clear();

      std::vector<std::vector<TPoint<T>>> boundary;
      for (int j = 0; j < int(boundaries[i].size()); j++)
      {
        for (int k = 0; k < int(boundaries[i][j].size()); k++)
        {
          //std::vector<int> list;
          //if (!findOverlapping(boundary,boundaries[i][j][k],true,list,tolerance,parmtolerance))
          if (boundaries[i][j][k].size() > 1)
            boundary.push_back(boundaries[i][j][k]);
        }
      }

      if (!boundary.empty())
      {
        bool reversed = false;
        int index = longestSizePiece(boundary);
        int boundaryside = boundaryLine<T>(boundary[index],parmtolerance,&reversed);
        if (boundaryside >= 0)
        {
          if (reversed)
          {
            exlist.push_back(i);
          } else
          {
            closeOuterBoundary(boundariesUV[i]);
          }
        } else
        {
          bool ok = surfaces[i]->closeBoundaryLoop(boundary,boundariesUV[i],
            maxseglen * 2.0,parmtolerance,manypoints - 1);
        }
      }
    }
  }
}

/** Close connected faces with same boundary direction. */
template <class T> void closeConnectedFaces(std::vector<TSplineSurface<T> *> &surfaces,
  // surface  // loop 0   // 4 pieces // piece contents
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV,
  std::vector<int> &exlist,
  T tolerance, T bigtolerance, int numattempts = 500)
{
  for (int a = 0; a < numattempts; a++)
  {
    std::vector<int> list;

    std::vector<std::vector<TPointCurve<T>>> curves;

    prepareLoopsXYZ<T>(surfaces,boundariesUV,curves);

//outputDebugString(""); 

    for (int i = 0; i < int(boundariesUV.size()); i++)
    {
      // we need to make a correct outer loop taking its orientation from cut surfaces
      if (boundariesUV[i].empty())
      {
        // find an edge connected to a outerloop edge from cut surfaces which already
        // have a proper orientation after cuts

        std::vector<bool> busy(boundariesUV.size(),false);
        for (int j = 0; j < int(boundariesUV.size()); j++)
        {
          // this orientation is correct after cuts
          busy[j] = ((i != j) && !boundariesUV[j].empty());
        }

        std::pair<int,int> res;
        bool reversed = false;
        T mindist = 0.0;
        if (findClosest<T>(i,curves,busy,res,tolerance,bigtolerance,&reversed,&mindist))
        {

//outputDebugString("i+ " + to_string(i) + " j " + to_string(res.first) + " k " + to_string(res.second) + 
//  " mindist " + to_string(mindist) + 
//  " tolerance " + to_string(tolerance,12) + " bigtolerance " + to_string(bigtolerance,12));

          // it must be empty
          bool keepempty = (std::find(exlist.begin(),exlist.end(),i) != exlist.end());

          if (reversed ^ keepempty)
            list.push_back(i);
        } else
        {
//outputDebugString("i- " + to_string(i) + " mindist " + to_string(mindist) + 
//  " tolerance " + to_string(tolerance,12) + " bigtolerance " + to_string(bigtolerance,12)); 
        }
      }
    }

    if (list.empty())
      break;

    for (int i : list)
    {
      closeOuterBoundary(boundariesUV[i]);
    }
  }
}


}