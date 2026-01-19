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

  tblocks.h

  Making (solid) blocks like a surface of revolution, airfoil etc.

*******************************************************************************/

#pragma once

#include "tpoints.h"
#include "tmisc.h"
#include "tjacobipoly.h"
#include "tbasesurface.h"
#include "tsplinesurface.h"
#include "toperations.h"
#include "tedge.h"

#include <map>

//#define DEBUG_BLOCKS
//#ifdef NDEBUG
//  #undef DEBUG_BLOCKS
//#endif

#ifdef DEBUG_BLOCKS
  #include "export.h"
#endif

namespace tcad {

/** Make NACA airfoil wing surface. */
template <class T> void makeNACASurface(std::vector<std::vector<TPoint<T>>> &points,
  int numspans, TPoint<T> resizecoef = TPoint<T>(1,1,1), T resizeX = 0.99, T resizeY = 0.99, T moveZ = -0.02, 
  T dadegZ = -0.5, T dadegY = -0.5)
{
  points.clear();

  // rescale to -0.5, +0.5]
  std::vector<TPoint<T>> airfoil;
  for (auto p : NACA0012xy<T>)
  {
    TPoint<T> pos(p.first - 0.5,p.second);
    pos.X *= resizecoef.X;
    pos.Y *= resizecoef.Y;
    pos.Z *= resizecoef.Z;

    airfoil.push_back(pos);
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

/** Make 3D airfoil (blade) from camber surface and thickness. Camber surface and thickness must have 
  the same number of points in both directions, numbered first along U (chord), then along V.
  The tickness must be zero at U = 0 and U = 1. */
template <class T> void makeAirfoil(std::vector<std::vector<TPoint<T>>> &camberpoints,
  std::vector<std::vector<T>> &thickness, std::vector<TSplineSurface<T> *> &surfaces, 
  int M1 = SPLINE_DEGREE, int M2 = SPLINE_DEGREE, 
  CurveEndType startU = END_CLAMPED, CurveEndType endU = END_FREE,
  CurveEndType startV = END_FREE, CurveEndType endV = END_FREE)
{
  assert(camberpoints.size() == thickness.size());
  assert(camberpoints[0].size() == thickness[0].size());

//!!!!!!  surfaces.clear();

  int K1 = int(camberpoints[0].size()) - 1;
  int K2 = int(camberpoints.size()) - 1;

  // make camber surface to calculate normals, the surfaces parametrisation is by numbers
  TPointSurface<T> cambersurface(camberpoints,true,true);

  std::vector<std::vector<TPoint<T>>> upperpoints,lowerpoints;

  for (int i = 0; i <= K2; i++)
  {
    T V = T(i) / T(K2);

    upperpoints.push_back(std::vector<TPoint<T>>());
    lowerpoints.push_back(std::vector<TPoint<T>>());
    for (int j = 0; j <= K1; j++)
    {
      T U = T(j) / T(K1);

      TPoint<T> pos = cambersurface.position(U,V);
      TPoint<T> derU = cambersurface.derivative(U,V,PARAMETER_U,1);
      TPoint<T> derV = cambersurface.derivative(U,V,PARAMETER_V,1);
      TPoint<T> normal = (+(derU ^ derV)) * (thickness[i][j] * 0.5); //!!! half thickness

      upperpoints.back().push_back(pos + normal);
      lowerpoints.back().push_back(pos - normal);
    }
  }

  reverseRows(lowerpoints);

  // approximate
#if 0
  TSplineSurface<T> *upper = new TSplineSurface<T>(upperpoints,M1,M2,startU,endU,startV,endV);
  TSplineSurface<T> *lower = new TSplineSurface<T>(lowerpoints,M1,M2,startU,endU,startV,endV);
#else
  TSplineSurface<T> *upper = new TSplineSurface<T>(upperpoints,K1,M1,K2,M2,startU,endU,startV,endV);
  TSplineSurface<T> *lower = new TSplineSurface<T>(lowerpoints,K1,M1,K2,M2,endU,startU,startV,endV); // it was reversed
#endif

  surfaces.push_back(upper);
  surfaces.push_back(lower);
}

/** Make airfoil from upper and lower surface coordinates. Points must be numbered from 
  TE along upper surface, around LE, then along the lower surface back to TE, like in
  http://airfoiltools.com/airfoil/details?airfoil=naca1410-il.
  TE is supposed to have max X coordinate, LE - min X coordinate.
  Output is single upperlower numpoints array with X from LE to TE contained as 
  X of every point, Y as upper Y coordinate and Z(!) as lower(!) Y coordinate. 
  Points near edges can be refined by setting endrefinementcoef < 1.0.
  Only X,Y point components are valid. 
  Points are approximated with orthogonal polynomials in x,y, not periodic, with
  alpha,beta to treat infinite derivatives at edges.
  Return is approximation accuracy for upper and lower surfaces. */
template <class T> std::pair<T,T> makeAirfoilPointsXY(std::vector<TPoint<T>> &points,
  bool roundLE, bool roundTE, int numpoints, std::vector<TPoint<T>> &upperlower,
  bool reverse = false, int degree = 5, int integration = OTHER_INTEGRATION, T endrefinementcoef = 0.5)
{
  // divide points into two parts
  TPoint<T> min,max,imin,imax;

  if (!calculateMinMax<T>(points,&min,&max,&imin,&imax))
  {
    return std::pair<T,T>(-1.0,-1.0);
  }

  int LEindex = ROUND(imin.X);
  T xmin = min.X;
  T xmax = max.X;
  T dx = xmax - xmin;

  std::vector<TPoint<T>> tempupper,templower;

  for (int i = LEindex; i >= 0; i--)
  {
    tempupper.push_back(points[i]);
  }

  for (int i = LEindex; i < int(points.size()); i++)
  {
    templower.push_back(points[i]);

    // lower surface has negative Y, Jacoby poly does not like negative values
    templower.back().Y *= -1.0;
  }

  // approximate both surfaces with ortho poly
  TJacobiPoly<T> fu(roundTE ? 0.5 : 0.0,roundLE ? 0.5 : 0.0);

  std::vector<T> x,y,z;
  splitXYZ<T>(tempupper,x,y,z);

  bool ok = fu.fit(degree,x,y,integration);
  if (!ok)
  {
    return std::pair<T,T>(-1.0,-1.0);
  }

  T uaccuracy = fu.accuracy(x,y);

  // now the lower surface
  TJacobiPoly<T> fl(roundTE ? 0.5 : 0.0,roundLE ? 0.5 : 0.0);

  splitXYZ<T>(templower,x,y,z);

  ok = fl.fit(degree,x,y,integration);
  if (!ok)
  {
    return std::pair<T,T>(-1.0,-1.0);
  }

  T laccuracy = fl.accuracy(x,y);

  // now make actual points
  for (int i = 0; i < numpoints; i++)
  {
    T U = T(i) / T(numpoints - 1);
    T u = U * 2.0 - 1.0;
    T un = pow(std::abs(u),endrefinementcoef) * sign(u);
    U = (un + 1.0) * 0.5;

    // this is X
    T x = xmin + dx * U;

    // y for upper and lower
    T yu = fu.getValue(un);
    T yl = fl.getValue(un);

    // lower surface has negative Y, Jacoby poly does not like negative values
    upperlower.push_back(TPoint<T>(x,yu,-yl));
  }

  // correct start/end just in case
  upperlower.front() = TPoint<T>(points[LEindex].X,points[LEindex].Y,points[LEindex].Y);
  upperlower.back() = TPoint<T>(points[0].X,points[0].Y,points[0].Y);

  if (reverse)
  {
    TPoint<T> min,max;
    calculateMinMax(upperlower,&min,&max);
    T dx = max.X - min.X;

    for (auto &p : upperlower)
    {
      p.X = dx - p.X;
    }
  }

  return std::pair<T,T>(uaccuracy,laccuracy);
}

/** Generate regular net in X,Y,Z (Z - span) of camber/thickness points for use in makeAirfoil(). 
  The upperlower points can be generated by makeAirfoilPointsXY(). Every upperlower point contains 
  X in X coordinate, upper Y in Y coordinate and lower Y in Z coordinate. The points go from LE 
  (left, min X) to TE (max X). 
  Middle line of contour of propeller blade in metres. Every point contains :
  X - Z along span
  Y - X(Z) of middle line 
  Z - blade HALF-width from middle line
  W - twist angle around middle curve in degrees. */
template <class T> void makeBladeCamberThickness(std::vector<TPoint<T>> &upperlower, 
  std::vector<TPoint<T>> &blade, int numchordpoints, int numspanpoints,
  std::vector<std::vector<TPoint<T>>> &camberpoints, std::vector<std::vector<T>> &thickness,
  T tolerance, bool smoothcontour = true)
{
  camberpoints.clear();
  thickness.clear();

  // blade tip is rounded
  bool roundtip = (std::abs(blade.back().Z) < tolerance);

  // middle points for middle line
  std::vector<TPoint<T>> middlepoints;
  // width on Z
  std::vector<TPoint<T>> width;
  // twist on Z
  std::vector<TPoint<T>> twist;
  for (auto &p : blade)
  {
    middlepoints.push_back(TPoint<T>(p.Y,0.0,p.X));
    width.push_back(TPoint<T>(p.Z,0.0,p.X)); // half-width on Z
    twist.push_back(TPoint<T>(p.W,0.0,p.X)); // twist on Z
  }

  // smooth 
  if (smoothcontour)
  {
#if 1 //!!!!!!!
    smoothPointsBySplineCurve(middlepoints,5,END_FIXED,END_FREE);
    smoothPointsBySplineCurve(width,5,END_FIXED,END_FREE);
    smoothPointsBySplineCurve(twist,5,END_FIXED,END_FREE);
#else
    smoothPointsByBezier(middlepoints,END_FIXED,END_FREE);
    smoothPointsByBezier(width,END_FIXED,END_FREE);
    smoothPointsByBezier(twist,END_FIXED,END_FREE);
#endif
  }

#ifdef DEBUG_BLOCKS
  {
    std::vector<std::vector<TPoint<T>>> cpoints;
    cpoints.push_back(middlepoints);
    cpoints.push_back(width);
    cpoints.push_back(twist);
  
    saveLinesIges(cpoints,"twist.iges");
  }
#endif

  // now make camber and thickness
  TSplineCurve<T> middlespline(middlepoints,SPLINE_DEGREE);
  TSplineCurve<T> widthspline(width,SPLINE_DEGREE);
  TSplineCurve<T> twistspline(twist,SPLINE_DEGREE);

  // airfoil coordinates
  std::vector<TPoint<T>> upper,lower;
  for (auto &p : upperlower)
  {
    upper.push_back(TPoint<T>(p.X,p.Y));
    lower.push_back(TPoint<T>(p.X,p.Z));
  }

  // we need dx to rescale
  TPoint<T> umin,umax;
  calculateMinMax<T>(upper,&umin,&umax);
  TPoint<T> lmin,lmax;
  calculateMinMax<T>(lower,&lmin,&lmax);
  TPoint<T> min = pointMin<T>(umin,lmin);
  TPoint<T> max = pointMax<T>(umax,lmax);

  TPoint<T> centre = (min + max) * 0.5;
  T dx = max.X - min.X;
  T DW = 1.0 / T(numspanpoints - 1);

  for (int k = 0; k < numspanpoints; k++)
  {
    if (roundtip && k == numspanpoints - 1)
    {
      // get actual contours from camber surface
      std::vector<TPoint<T>> LE,TE;
      getColumn<T>(camberpoints,0,LE);
      getColumn<T>(camberpoints,int(camberpoints[0].size()) - 1,TE);

      TPoint<T> dir0 = -(+(endDirection(LE)));
      TPoint<T> dir1 = +(endDirection(TE));
      TPoint<T> p0 = LE.back();
      TPoint<T> p1 = TE.back();
      T len = !(p1 - p0) * 1.5; //!!!!!!!

      TBezierSegment<T> segment(p0,p1,dir0,dir1,len); 

      std::vector<TPoint<T>> lpoints;

      T margin = 0.01;
      T Umin = margin;
      T Umax = 1.0 - margin;
      int numpoints = 10;
      T DU = (Umax - Umin) / T(numpoints - 1);
      for (int j = 0; j < numpoints; j++)
      {
        T U = Umin + T(j) * DU;
        lpoints.push_back(segment.position(U));
      }

#ifdef DEBUG_BLOCKS
      {
        std::vector<std::vector<TPoint<T>>> cpoints;
        cpoints.push_back(lpoints);
  
        saveLinesIges(cpoints,std::string("roundtip.iges"));
      }
#endif


      TPointCurve<T> pupper(lpoints); 
      TPointCurve<T> plower(lpoints);

      camberpoints.push_back(std::vector<TPoint<T>>());
      thickness.push_back(std::vector<T>());

      for (int i = 0; i < numchordpoints; i++)
      {
        T U = T(i) / T(numchordpoints - 1);
        TPoint<T> pu = pupper.position(U);
        TPoint<T> pl = plower.position(U);

        TPoint<T> p = (pu + pl) * 0.5;
        camberpoints.back().push_back(p);
        thickness.back().push_back(0.0);
      }
    } else
    {
      T W = T(k) / T(numspanpoints - 1);

      // left and right LE/TE points
      TPoint<T> p0,p1,dp;

      // middle line direction
      TPoint<T> mdir = +(middlespline.derivative(W,1));
      // centre
      TPoint<T> mc = middlespline.position(W);
      // width in X
      T halfwidth = widthspline.position(W).X;
      // twist in X
      T twistdeg = twistspline.position(W).X;

//!!! the first line must be horizontal
      if (k == 0)
      {
        p0 = mc - TPoint<T>(halfwidth);
        p1 = mc + TPoint<T>(halfwidth);
        dp = p1 - p0;
        // direction of middle line
        mdir = +(dp ^ TPoint<T>(0.0,1.0,0.0));
      } else
      {
        TPoint<T> move = +(mdir ^ TPoint<T>(0.0,1.0,0.0)) * halfwidth;
        p0 = mc - move;
        p1 = mc + move;
        dp = p1 - p0;
      }

      // centre between lines
      TPoint<T> c = (p0 + p1) * 0.5;
      // rescale coefficient
      T len = !dp;
      T rescalecoef = len / dx;

      // we now need to stretch upperlower points between p0 and p1
      std::vector<TPoint<T>> tempupper = upper;
      std::vector<TPoint<T>> templower = lower;

      // apply transform to airfoil XY points
      TTransform<T> t;
      t.Translate(-centre);
      t.Resize(TPoint<T>(rescalecoef,rescalecoef,rescalecoef));

      // rotate ~around Z

      // current normal in XY plane
      TPoint<T> normal(0.0,0.0,1.0);

      TPoint<T> rotaxisZ = +(normal ^ mdir);
      T rotangleZ = normal < mdir;
      t.Rotate(rotaxisZ,rotangleZ); 

      makeTransform<T>(tempupper,&t);
      makeTransform<T>(templower,&t);

      // rotate ~around middle line by twist angle
      t.LoadIdentity();
      t.Rotate(mdir,twistdeg * CPI);

      t.Translate(c);

      makeTransform<T>(tempupper,&t);
      makeTransform<T>(templower,&t);

      TPointCurve<T> pupper(tempupper); 
      TPointCurve<T> plower(templower);

#ifdef DEBUG_BLOCKS
      {
        std::vector<std::vector<TPoint<T>>> cpoints;
        cpoints.push_back(tempupper);
        cpoints.push_back(templower);
  
        saveLinesIges(cpoints,std::string("tempupperlower") + to_string(k) + ".iges");
      }
#endif

      camberpoints.push_back(std::vector<TPoint<T>>());
      thickness.push_back(std::vector<T>());

      if (numchordpoints == int(tempupper.size()))
      {
        for (int i = 0; i < numchordpoints; i++)
        {
          TPoint<T> pu = tempupper[i];
          TPoint<T> pl = templower[i];

          TPoint<T> p = (pu + pl) * 0.5;
          T d = !(pu - pl);
          camberpoints.back().push_back(p);
          thickness.back().push_back(d);
        }
      } else
      {
        for (int i = 0; i < numchordpoints; i++)
        {
          T U = T(i) / T(numchordpoints - 1);
          TPoint<T> pu = pupper.position(U);
          TPoint<T> pl = plower.position(U);

          TPoint<T> p = (pu + pl) * 0.5;
          T d = !(pu - pl);
          camberpoints.back().push_back(p);
          thickness.back().push_back(d);
        }
      }
    }
  }
}

/** Generate regular net in X,Y,Z (Z - span) of camber/thickness points for use in makeAirfoil(). 
  The upperlower points can be generated by makeAirfoilPointsXY(). Every upperlower point contains 
  X in X coordinate, upper Y in Y coordinate and lower Y in Z coordinate. The points go from LE 
  (left, min X) to TE (max X). 
  leftright represent left/right blade contour lines like those in KiloWing1Contour. */
template <class T> void makeBladeCamberThickness(std::vector<TPoint<T>> &upperlower, 
  std::vector<std::array<TPoint<T>,2>> &leftright, int numchordpoints, int numspanpoints,
  std::vector<std::vector<TPoint<T>>> &camberpoints, std::vector<std::vector<T>> &thickness,
  T tolerance, bool smoothcontour = true)
{
  camberpoints.clear();
  thickness.clear();

  std::vector<TPoint<T>> left,right;
  for (auto &p : leftright)
  {
    left.push_back(p[0]);
    right.push_back(p[1]);
  }

  // smooth 
  if (smoothcontour)
  {
#if 1 //!!!!!!!
    smoothPointsBySplineCurve(left,5,END_FIXED,END_FREE);
    smoothPointsBySplineCurve(right,5,END_FIXED,END_FREE);
#else
    smoothPointsByBezier(middlepoints,END_FIXED,END_FREE);
    smoothPointsByBezier(width,END_FIXED,END_FREE);
    smoothPointsByBezier(twist,END_FIXED,END_FREE);
#endif
  }

#ifdef DEBUG_BLOCKS
  {
    std::vector<std::vector<TPoint<T>>> cpoints;
    cpoints.push_back(left);
    cpoints.push_back(right);
  
    saveLinesIges(cpoints,"leftright.iges");
  }
#endif

  // now make camber and thickness
  TSplineCurve<T> leftspline(left,SPLINE_DEGREE);
  TSplineCurve<T> rightspline(right,SPLINE_DEGREE);

  // airfoil coordinates
  std::vector<TPoint<T>> upper,lower;
  for (auto &p : upperlower)
  {
    upper.push_back(TPoint<T>(p.X,p.Y));
    lower.push_back(TPoint<T>(p.X,p.Z));
  }

  // we need dx to rescale
  TPoint<T> umin,umax;
  calculateMinMax<T>(upper,&umin,&umax);
  TPoint<T> lmin,lmax;
  calculateMinMax<T>(lower,&lmin,&lmax);
  TPoint<T> min = pointMin<T>(umin,lmin);
  TPoint<T> max = pointMax<T>(umax,lmax);

  TPoint<T> centre = (min + max) * 0.5;
  T dx = max.X - min.X;
  T DW = 1.0 / T(numspanpoints - 1);

  for (int k = 0; k < numspanpoints; k++)
  {
    T W = T(k) / T(numspanpoints - 1);

    // left and right LE/TE points
    TPoint<T> p0 = leftspline.position(W);
    TPoint<T> p1 = rightspline.position(W);
    TPoint<T> dp = p1 - p0;
    // centre between lines
    TPoint<T> c = (p0 + p1) * 0.5;

    // rescale coefficient
    T len = !dp;
    T rescalecoef = len / dx;

    // we now need to stretch upperlower points between p0 and p1
    std::vector<TPoint<T>> tempupper = upper;
    std::vector<TPoint<T>> templower = lower;

    // apply transform to airfoil XY points
    TTransform<T> t;
    t.Translate(-centre);
    t.Resize(TPoint<T>(rescalecoef,rescalecoef,rescalecoef));

    // rotate
    TPoint<T> rotaxis = +(dp ^ TPoint<T>(1.0,0.0,0.0));
    T rotangle = dp < TPoint<T>(1.0,0.0,0.0);
    t.Rotate(rotaxis,-rotangle); 

    makeTransform<T>(tempupper,&t);
    makeTransform<T>(templower,&t);

    // rotate ~around middle line by twist angle
    t.LoadIdentity();
    t.Translate(c);

    makeTransform<T>(tempupper,&t);
    makeTransform<T>(templower,&t);

    TPointCurve<T> pupper(tempupper); 
    TPointCurve<T> plower(templower);

#ifdef DEBUG_BLOCKS
    {
      std::vector<std::vector<TPoint<T>>> cpoints;
      cpoints.push_back(tempupper);
      cpoints.push_back(templower);
  
      saveLinesIges(cpoints,std::string("tempupperlower") + to_string(k) + ".iges");
    }
#endif

    camberpoints.push_back(std::vector<TPoint<T>>());
    thickness.push_back(std::vector<T>());

    if (numchordpoints == int(tempupper.size()))
    {
      for (int i = 0; i < numchordpoints; i++)
      {
        TPoint<T> pu = tempupper[i];
        TPoint<T> pl = templower[i];

        TPoint<T> p = (pu + pl) * 0.5;
        T d = !(pu - pl);
        camberpoints.back().push_back(p);
        thickness.back().push_back(d);
      }
    } else
    {
      for (int i = 0; i < numchordpoints; i++)
      {
        T U = T(i) / T(numchordpoints - 1);
        TPoint<T> pu = pupper.position(U);
        TPoint<T> pl = plower.position(U);

        TPoint<T> p = (pu + pl) * 0.5;
        T d = !(pu - pl);
        camberpoints.back().push_back(p);
        thickness.back().push_back(d);
      }
    }
  }
}

/** Make points for surfaces of revolution around Z. 
  The contour has axial coordinate in countourAxisCoord and radial coordinate in countourRadialCoord.
  numpoints is that along circumference. */
template <class T> void makeSurfaceOfRevolution(std::vector<TPoint<T>> &contour, 
  Axes countourAxialCoord, Axes countourRadialCoord,
  int numpoints, T adegfrom, T adegto, std::vector<std::vector<TPoint<T>>> &points)
{
  int numsections = int(contour.size());

  for (int i = 0; i < numsections; i++)
  {
    T r = contour[i][countourRadialCoord];

    std::vector<TPoint<T>> section;
    makeEllipseXY(numpoints,TPoint<T>(),r,r,section,adegfrom,adegto);

    TTransform<T> t;
    t.Translate(TPoint<T>(0.0,0.0,contour[i][countourAxialCoord]));
    makeTransform(section,&t);

    points.push_back(section);
  }
}

/** Make single face of a surface of revolution around Z. 
  The contour has axial coordinate in countourAxisCoord and radial coordinate in countourRadialCoord.
  numfaces and pointsperface are those along circumference. 
  Resulting faces have U as circumferential and V - parametric axial coordinate. */
template <class T> TSplineSurface<T> *makeSurfaceOfRevolution(std::vector<TPoint<T>> &contour,
  Axes countourAxialCoord, Axes countourRadialCoord,
  int numpoints, T adegfrom, T adegto, int K1, int K2, 
  int M1 = SPLINE_DEGREE, int M2 = SPLINE_DEGREE,
  CurveEndType startU = END_CLAMPED, CurveEndType endU = END_CLAMPED, // must be round
  CurveEndType startV = END_FREE, CurveEndType endV = END_FREE) 
{
  // make points
  std::vector<std::vector<TPoint<T>>> points;
  makeSurfaceOfRevolution(contour,countourAxialCoord,countourRadialCoord,numpoints,adegfrom,adegto,points);

  // create cylindrical surface
  TSplineSurface<T> *surface = new TSplineSurface<T>(points,K1,M1,K2,M2,
    startU,endU,startV,endV);

  return surface;
}

/** Make surfaces of revolution around Z. 
  The contour has axial coordinate in countourAxisCoord and radial 
  coordinate in countourRadialCoord.
  numfaces and pointsperface are those along circumference. 
  Resulting faces have U as circumferential and V - parametric axial coordinate. */
template <class T> void makeSurfacesOfRevolution(std::vector<TPoint<T>> &contour,
  Axes countourAxialCoord, Axes countourRadialCoord,
  int numfaces, int pointsperface, int K1, int K2, 
  std::vector<TSplineSurface<T> *> &surfaces, 
  int M1 = SPLINE_DEGREE, int M2 = SPLINE_DEGREE,
  CurveEndType startU = END_CLAMPED, CurveEndType endU = END_CLAMPED, // must be round
  CurveEndType startV = END_FREE, CurveEndType endV = END_FREE) 
{
  T da = 360.0 / T(numfaces);
  for (int i = 0; i < numfaces; i++)
  {
    T a0 = T(i) * da;
    T a1 = T(i + 1) * da;

    TSplineSurface<T> *surface = makeSurfaceOfRevolution(contour,countourAxialCoord,countourRadialCoord,
      pointsperface,a0,a1,K1,K2,M1,M2,startU,endU,startV,endV);

    surfaces.push_back(surface);
  }
}

/** Make surfaces of revolution around Z (multiple faces around Z from 0 t0 360 deg). 
  The contour has axial coordinate in countourAxisCoord and radial coordinate in countourRadialCoord.
  numfaces and pointsperface are those along circumference. 
  Resulting faces have U as circumferential and V - parametric axial coordinate. */
template <class T> void makeSurfacesOfRevolution(std::vector<std::vector<TPoint<T>>> &contour,
  Axes countourAxialCoord, Axes countourRadialCoord,
  int numfaces, int pointsperface, int K1, int K2, 
  std::vector<TSplineSurface<T> *> &surfaces, 
  int M1 = SPLINE_DEGREE, int M2 = SPLINE_DEGREE,
  CurveEndType startU = END_CLAMPED, CurveEndType endU = END_CLAMPED, // must be round
  CurveEndType startV = END_FREE, CurveEndType endV = END_FREE) 
{
  for (int i = 0; i < int(contour.size()); i++)
  {
    makeSurfacesOfRevolution(contour[i],countourAxialCoord,countourRadialCoord,
      numfaces,pointsperface,K1,K2,surfaces,M1,M2,startU,endU,startV,endV);
  }
}

/** Intersect surfaces, make trimming curves to make a solid. It is a raw function to be used inside
  TBrep<T>. It creates surfaces and boundary copies. Surfaces should be deleted by deleteSurfaces(). 
  Always returns true now. */
template <class T> bool makeTrimming(
  // input
  std::vector<TSplineSurface<T> *> &surfaces0,
  std::vector<TSplineSurface<T> *> &surfaces1,
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV0,
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV1,
  // boolean operation
  Boolean operation,
  // result
  std::vector<TSplineSurface<T> *> &surfaces,
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV,
  // the same edge on another surface (defined by intersections), 100% reliable, no tolerances
  T tolerance, T parmtolerance = PARM_TOLERANCE, 
  bool clearboundaries = false, bool bodyleft = true, int manypoints = MANY_POINTS2D,
  bool improveintersection = false, int maxiter = 20, T relaxcoef = 0.5,
#if 0 //!!!!!!!
  T refinestartU = 0.5, T refineendU = 0.5, 
  T refinestartV = 0.5, T refineendV = 0.5,
  int maxcoef = 1) 
#else
  T refinestartU = 1.0, T refineendU = 1.0, //!!!!!!!
  T refinestartV = 1.0, T refineendV = 1.0,
  int maxcoef = 1)
#endif
{
  // temporary copies of surfaces and boundaries
  std::vector<TSplineSurface<T> *> tsurfaces0,tsurfaces1;
  copySurfaces<T>(surfaces0,tsurfaces0);
  copySurfaces<T>(surfaces1,tsurfaces1);
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> tboundariesUV0 = boundariesUV0;
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> tboundariesUV1 = boundariesUV1;

  // normally normal must be outside, it is forced inside
  if (operation == UNITE)
  {
  } else if (operation == SUBTRACT)
  {
    reverseSurfaces(tsurfaces1);
  } else if (operation == INTERSECT)
  {
  }

  if (clearboundaries)
  {
    tboundariesUV0.clear();
    tboundariesUV1.clear();

    prepareOuterBoundary(tsurfaces0,tboundariesUV0);
    prepareOuterBoundary(tsurfaces1,tboundariesUV1);
  }

  // prepare triangles
  // set good surface subdivisions for triangles (they better be uniform for threading)
  std::vector<TSplineSurface<T> *> allsurfaces = tsurfaces0;
  allsurfaces.insert(allsurfaces.end(),tsurfaces1.begin(),tsurfaces1.end());

  T surfacessize0 = surfacesSize(tsurfaces0);
  T surfacessize1 = surfacesSize(tsurfaces1);
  T minsurfacessize = std::min<T>(surfacessize0,surfacessize1);

  // find min surfaces size
  T minsize = std::numeric_limits<T>::max();
  for (int i = 0; i < int(allsurfaces.size()); i++)
  {
    T usize = allsurfaces[i]->Usize();
    T vsize = allsurfaces[i]->Vsize();
    T size = (usize + vsize) * 0.5;

    minsize = std::min(minsize,size);
  }

  // minsize will correspond to manypoints division,
  // make triangles for all surfaces together with bounding boxes
  std::vector<TTriangles<T> *> tris0;
  std::vector<TTriangles<T> *> tris1;
  std::vector<std::vector<std::array<TPoint<T>,8>>> tboxes0;
  std::vector<std::vector<std::array<TPoint<T>,8>>> tboxes1;

  for (int i = 0; i < int(tsurfaces0.size()); i++)
  {
    TTriangles<T> *tris = new TTriangles<T>();

    T usize = tsurfaces0[i]->Usize();
    T vsize = tsurfaces0[i]->Vsize();

    int coef0 = int(usize / minsize);
    int coef1 = int(vsize / minsize);

    LIMIT(coef0,1,maxcoef);
    LIMIT(coef1,1,maxcoef);

    int numpointsU = manypoints * coef0;
    int numpointsV = manypoints * coef1;

    tsurfaces0[i]->createTriangles(*tris,numpointsU,numpointsV,
      refinestartU,refineendU,refinestartV,refineendV);

    std::vector<std::array<TPoint<T>,8>> boxes;
    tris->makeFaceBoxes(boxes);

    tboxes0.push_back(boxes);
    tris0.push_back(tris);
  }

  for (int i = 0; i < int(tsurfaces1.size()); i++)
  {
    TTriangles<T> *tris = new TTriangles<T>();

    T usize = tsurfaces1[i]->Usize();
    T vsize = tsurfaces1[i]->Vsize();

    int coef0 = int(usize / minsize);
    int coef1 = int(vsize / minsize);

    LIMIT(coef0,1,maxcoef);
    LIMIT(coef1,1,maxcoef);

    int numpointsU = manypoints * coef0;
    int numpointsV = manypoints * coef1;
    tsurfaces1[i]->createTriangles(*tris,numpointsU,numpointsV,
      refinestartU,refineendU,refinestartV,refineendV);

    std::vector<std::array<TPoint<T>,8>> boxes;
    tris->makeFaceBoxes(boxes);

    tboxes1.push_back(boxes);
    tris1.push_back(tris);
  }

  // "delayed" boundary pieces
  std::vector<std::vector<std::vector<std::vector<TPoint<T>>>>> boundaries0(tsurfaces0.size());
  std::vector<std::vector<std::vector<std::vector<TPoint<T>>>>> boundaries1(tsurfaces1.size());

  // min/max of all surfaces
  TPoint<T> min,max;
  T maxedge = 0.0;

  std::vector<std::pair<TPoint<T>,TPoint<T>>> boxes0;

  for (int i = 0; i < int(tsurfaces0.size()); i++)
  {
    std::pair<TPoint<T>,TPoint<T>> mm = tsurfaces0[i]->getMinMax();
    T d = !(mm.second - mm.first);
    maxedge = std::max<T>(maxedge,d);

    if (i == 0)
    {
      min = mm.first;
      max = mm.second;
    } else
    {
      min = pointMin<T>(min,mm.first);
      max = pointMax<T>(max,mm.second);
    }

    extendBox(mm,1.1); //!!!
    boxes0.push_back(mm);
  }

  std::vector<std::pair<TPoint<T>,TPoint<T>>> boxes1;

  for (int i = 0; i < int(tsurfaces1.size()); i++)
  {
    std::pair<TPoint<T>,TPoint<T>> mm = tsurfaces1[i]->getMinMax();
    T d = !(mm.second - mm.first);
    maxedge = std::max<T>(maxedge,d);

    min = pointMin<T>(min,mm.first);
    max = pointMax<T>(max,mm.second);

    extendBox(mm,1.1); //!!!
    boxes1.push_back(mm);
  }

  maxedge *= 1.1; //!!!!!!
  // makeTrimming() always returns true (used in operators)
#if 1
  extendMinMax<T>(min,max,1.0,tolerance);
#endif

  TPoint<T> dmm = max - min;

  // makeTrimming() always returns true (used in operators)
#if 0
  // no intersection possible
  if (dmm.X < tolerance || dmm.Y < tolerance || dmm.Z < tolerance)
  {
    deleteTriangles<T>(tris0);
    deleteTriangles<T>(tris1);
    return false;
  }
#endif

  // cells big enough for spacial partitioning to identify if two triangles
  // may intersect, they can if they have a node in the same cell only
  OBackground<T> cells(min,max,(LINT) (dmm.X / maxedge),(LINT) (dmm.Y / maxedge),(LINT) (dmm.Z / maxedge),8);

  // total number of cells
  LINT numcells = cells.numBackgroundCells();

  // now distribute all surfaces over cells
  std::vector<std::set<LINT>> cellfaces0(numcells);
  for (LINT i = 0; i < int(tsurfaces0.size()); i++)
  {
    std::pair<TPoint<T>,TPoint<T>> mm = tsurfaces0[i]->getMinMax();

    std::array<TPoint<T>,8> corners;
    makeBox<T>(mm.first,mm.second,corners);
   
    for (int j = 0; j < 8; j++)
    {
      LINT index = cells.findCell(corners[j]);
      assert(index >= 0);

      cellfaces0[index].insert(i);
    }
  }

  // other faces in cells
  std::vector<std::set<LINT>> cellfaces1(numcells);
  for (LINT i = 0; i < int(tsurfaces1.size()); i++)
  {
    std::pair<TPoint<T>,TPoint<T>> mm = tsurfaces1[i]->getMinMax();

    std::array<TPoint<T>,8> corners;
    makeBox<T>(mm.first,mm.second,corners);
   
    for (int j = 0; j < 8; j++)
    {
      LINT index = cells.findCell(corners[j]);
      assert(index >= 0);

      cellfaces1[index].insert(i);
    }
  }

  // now distinguish cells which are (1) not empty (2) contain faces from
  // this and other surfaces (tsurfaces0 and tsurfaces1)
  std::vector<LINT> activecells;
  for (LINT i = 0; i < numcells; i++)
  {
    if (!cellfaces0[i].empty() && !cellfaces1[i].empty())
    {
      activecells.push_back(i);
    }
  }
  int numactive = int(activecells.size());

  // for tolerance
  T maxseglen = 0.0;
  T bigtolerance = 0.0;

  // now make intersections
  for (int cl = 0; cl < numactive; cl++)
  {
    LINT acell = activecells[cl];

    for (auto i : cellfaces0[acell])
    {
#ifdef DEBUG_BLOCKS
      outputDebugString(std::string("i ") + to_string(i));
#endif

      TSplineSurface<T> *surface0 = tsurfaces0[i];

      for (auto j : cellfaces1[acell])
      {
#ifdef DEBUG_BLOCKS
        outputDebugString(std::string("j ") + to_string(j));
#endif

        // simple measure
        if (!boxesOverlap(boxes0[i],boxes1[j]))
          continue;

        TSplineSurface<T> *surface1 = tsurfaces1[j];

        // this may switch off in this loop at the slightest failure
        bool improve = improveintersection;

        std::vector<std::vector<TPoint<T>>> intersections; 
        std::vector<std::vector<TPoint<T>>> boundary0,boundary1,iboundary0,iboundary1;

  #ifdef DEBUG_BLOCKS
        outputDebugString("");
  #endif

        // most important part in the whole play
        bool ok = tris0[i]->intersect(*tris1[j],int(bodyleft),tboxes0[i],tboxes1[j],
          intersections,parmtolerance,&boundary0,&boundary1,NUM_THREADS); 

        if (ok)
        {
#if 1
          iboundary0 = boundary0;
          iboundary1 = boundary1;
#else
          for (int k = 0; k < int(boundary0.size()); k++) 
          {
            if (boundaryLine<T>(boundary0[k],parmtolerance))
            {
              TPoint<T> dir = midDirection<T>(boundary0[k]);
              TPoint<T> middle = midPoint<T>(boundary0[k]);
              int side = boundarySide<T>(middle,parmtolerance);
              TPoint<T> sidedir = boundarySideDir<T>(side);

              if (!(dir > sidedir)) 
              {
                iboundary0.push_back(boundary0[k]);
              }
            } else
            {
              iboundary0.push_back(boundary0[k]);
            }

            if (boundaryLine<T>(boundary1[k],parmtolerance))
            {
              TPoint<T> dir = midDirection<T>(boundary1[k]);
              TPoint<T> middle = midPoint<T>(boundary1[k]);
              int side = boundarySide<T>(middle,parmtolerance);
              TPoint<T> sidedir = boundarySideDir<T>(side);

              if (!(dir > sidedir)) 
              {
                iboundary1.push_back(boundary1[k]);
              }
            } else
            {
              iboundary1.push_back(boundary1[k]);
            }
          }
#endif

          if (improve)
          {
            for (int k = 0; k < iboundary0.size(); k++)
            {
              for (int l = 0; l < int(iboundary0[k].size()); l++)
              {
                TPoint<T> parms(iboundary0[k][l].X,iboundary0[k][l].Y,iboundary1[k][l].X,iboundary1[k][l].Y);

                //!!! this may spoil everything, planar faces? are problematic
                if (improveIntersection<T>(tsurfaces0[i],tsurfaces1[j],parms,maxiter,relaxcoef,parmtolerance))
                {
                  iboundary0[k][l].X = parms.X;
                  iboundary0[k][l].Y = parms.Y;
                  iboundary1[k][l].X = parms.Z;
                  iboundary1[k][l].Y = parms.W;
                } else
                {
                  // one failure stops all
                  improve = false;
                  break;
                }
              }

              // one failure stops all
              if (!improve)
                break;

              // this improvement near boundaries tends to go to it and
              // make duplicate points
              correctParmsParallel(iboundary0[k],iboundary1[k],parmtolerance);
            }
          }

          // if still improve
          boundary0 = iboundary0;
          boundary1 = iboundary1;

          assert(boundary0.size() == boundary1.size());

          // calculate difference in boundary lines
          T maxdiff = 0.0;

          for (int k = 0; k < boundary0.size(); k++)
          {
            std::vector<TPoint<T>> intr0,intr1;
            for (int l = 0; l < int(boundary0[k].size()); l++)
            {
              TPoint<T> p0 = surfaces0[i]->position(boundary0[k][l].X,boundary0[k][l].Y);
              TPoint<T> p1 = surfaces1[j]->position(boundary1[k][l].X,boundary1[k][l].Y);
              intr0.push_back(p0);
              intr1.push_back(p1);

              T d = !(p1 - p0);
              maxdiff = std::max<T>(maxdiff,d);
            }
          }

          bigtolerance = std::max<T>(bigtolerance,maxdiff);

          for (int k = 0; k < int(boundary0.size()); k++)
          {
            T segmin = 0.0;
            T segmax = 0.0;
            segmentLenMinMax<T>(boundary0[k],segmin,segmax);

            maxseglen = std::max(maxseglen,segmax);
          }

          for (int k = 0; k < int(boundary1.size()); k++)
          {
            T segmin = 0.0;
            T segmax = 0.0;
            segmentLenMinMax<T>(boundary1[k],segmin,segmax);

            maxseglen = std::max(maxseglen,segmax);
          }

  #ifdef DEBUG_BLOCKS
          outputDebugString(std::string(" maxseglen = ") + to_string(maxseglen));
  #endif
        }

        if (ok)
        {
#if 0
          // if operation is union, check intersection direction : the tsurfaces0
          // cut must be directed opposite to tsurfaces0 normal
          int testpoint = int(boundary0[0].size()) / 2;
          TPoint<T> testpoint0 = boundary0[0][testpoint];
          TPoint<T> testpoint1 = boundary1[0][testpoint];

          // two normals at test point
          TPoint<T> normal0 = +surface0->normal(testpoint0.X,testpoint0.Y);
          TPoint<T> normal1 = +surface1->normal(testpoint1.X,testpoint1.Y);

          // intersection direction in model space
          TPoint<T> idir = direction<T>(intersections[0],testpoint);

          TPoint<T> n0i = normal0 ^ idir;
          bool dirok = n0i > normal1;

          bool reverse0 = !dirok;
#endif

          bool reverse0 = (operation == UNITE);

          if (reverse0)
          {
            for (auto &b : boundary0)
            {
              std::reverse(b.begin(),b.end());
            }
            std::reverse(boundary0.begin(),boundary0.end());
          }

          bool reverse1 = !reverse0; //???????

          if (reverse1)
          {
            for (auto &b : boundary1)
            {
              std::reverse(b.begin(),b.end());
            }
            std::reverse(boundary1.begin(),boundary1.end());
          }

          boundaries0[i].push_back(boundary0);
          boundaries1[j].push_back(boundary1);
        }
      }
    }
  }

  // close loops of faces intersected by boundaries
  std::vector<int> exlist0;
  closeFaceLoops(tsurfaces0,tboundariesUV0,boundaries0,exlist0,maxseglen,parmtolerance,manypoints);

  // close loops of faces intersected by boundaries
  std::vector<int> exlist1;
  closeFaceLoops(tsurfaces1,tboundariesUV1,boundaries1,exlist1,maxseglen,parmtolerance,manypoints); 

  // extend faces from cut faces keeping orientation
  closeConnectedFaces(tsurfaces0,tboundariesUV0,exlist0,tolerance,tolerance * 100.0); //!!!!!!!

  // extend faces from cut faces keeping orientation
  closeConnectedFaces(tsurfaces1,tboundariesUV1,exlist1,tolerance,tolerance * 100.0); //!!!!!!!

//for (int i = 0; i < int(surfaces1.size()); i++)
//{
//  saveTrimmedSurfaceIges<T>(surfaces1[i],tboundariesUV1[i],to_string(i) + ".iges"); 
//}

  deleteTriangles<T>(tris0);
  deleteTriangles<T>(tris1);

  // combine all into one list
  deleteSurfaces(surfaces);
  surfaces = tsurfaces0;
  surfaces.insert(surfaces.end(),tsurfaces1.begin(),tsurfaces1.end());

  boundariesUV = tboundariesUV0;
  boundariesUV.insert(boundariesUV.end(),tboundariesUV1.begin(),tboundariesUV1.end());

  ASSERT(surfaces.size() == boundariesUV.size());

  // make single loop everywhere to keep Rhino happy
  //!!!!!!! makeSingleLoop(surfaces,boundariesUV); //!!!!!!! not correct for a single hole

  return true;
}

/** Make rounding points on top of wall. */
template <class T> void makeRoundingPoints(TSplineSurface<T> *wall, 
  T radius, T insidesign, int K2, std::vector<std::vector<TPoint<T>>> &points)
{
  assert(radius > TOLERANCE(T));

  std::vector<TPoint<T>> toppoints;
  wall->getRow(wall->K2,toppoints);

  points.resize(K2 + 1,std::vector<TPoint<T>>(toppoints.size()));

  for (int i = 0; i < int(toppoints.size()); i++)
  {
    std::vector<TPoint<T>> col;
    wall->getColumn(i,col);

    // column direction pointing "up"
    TPoint<T> coldir = +(col[col.size() - 1] - col[col.size() - 2]);

    // direction along points to the right
    TPoint<T> pointsdir = direction(toppoints,i);

    // "inside" direction
    TPoint<T> insidedir = (pointsdir ^ coldir) * insidesign;

    // Bezier segment
    TPoint<T> p0 = toppoints[i];
    TPoint<T> p3 = p0 + insidedir * radius + coldir * radius;
    TPoint<T> p1 = p0 + coldir * (radius * BEZIER_CIRCLECOEF);
    TPoint<T> p2 = p3 - insidedir * (radius * BEZIER_CIRCLECOEF);

    TBezierSegment<T> seg(p0,p1,p2,p3);

    for (int j = 0; j <= K2; j++)
    {
      T V = T(j) / T(K2);
      TPoint<T> p = seg.position(V); 

      points[j][i] = p;
    }
  }
}

/** Make flat top with rounded corners between two curves. The flat part U = 0 and U = 1 has degenerated
  edges. points0 are on V = 0, points1 - on V = 1. points0 and points1 batter have same length. */
template <class T> void makeFlatRoundedTop(TSplineSurface<T> *wall0, TSplineSurface<T> *wall1,
  T radius, std::vector<TSplineSurface<T> *> &surfaces, int K2top = 8, int K2round = 8, int degree = SPLINE_DEGREE,
  // for curve intersections
  T tolerancecoef = 0.01, T parmtolerance = PARM_TOLERANCE, int manypoints = MANY_POINTS)
{
  int K1 = wall0->K1;
  CurveEndType startU = wall0->clampedStartU() ? END_CLAMPED : END_FREE;
  CurveEndType endU = wall0->clampedEndU() ? END_CLAMPED : END_FREE;
  CurveEndType startV = wall0->clampedStartV() ? END_CLAMPED : END_FREE;
  CurveEndType endV = wall0->clampedEndV() ? END_CLAMPED : END_FREE;

  std::vector<TPoint<T>> points0,points1;
  wall0->getRow(wall0->K2,points0);
  wall1->getRow(wall1->K2,points1);

  // this is not desirable
  if (points1.size() > points0.size())
  {
    TPointCurve<T> curve(points0);
    curve.redivide(int(points1.size()),points0);
  } else if (points0.size() > points1.size())
  {
    TPointCurve<T> curve(points1);
    curve.redivide(int(points0.size()),points1);
  } 

  if (radius < TOLERANCE(T))
  {
    // make flat top
    std::vector<std::vector<TPoint<T>>> flatpoints;
    flatpoints.push_back(points0);
    flatpoints.push_back(points1);

    TSplineSurface<T> *flat = new TSplineSurface<T>(flatpoints,K1,degree,K2top,degree,
      startU,endU,startV,endV); 

    surfaces.push_back(flat);
  } else
  {
    std::vector<std::vector<TPoint<T>>> wall0points,wall1points;

    makeRoundingPoints(wall0,radius,-1.0,K2round,wall0points);
    makeRoundingPoints(wall1,radius,-1.0,K2round,wall1points);

    for (int i = 0; i < int(wall0points.size()); i++)
    {
      std::reverse(wall0points[i].begin(),wall0points[i].end());
    }

    for (int i = 1; i < int(wall0points.size()); i++)
    {
      std::vector<TPoint<T>> UV;

      TPointCurve<T> curve0(wall0points[i]);
      TPointCurve<T> curve1(wall1points[i]);

      T tolerance = curve0.length() * tolerancecoef;

      int numintrs = curve0.findIntersections(curve1,UV,tolerance,parmtolerance,manypoints); 

      if (numintrs == 2)
      {
        curve0.cutPiece(int(wall0points[i].size()),UV[0].X,UV[1].X,wall0points[i]);
        curve1.cutPiece(int(wall1points[i].size()),UV[0].Y,UV[1].Y,wall1points[i]);

        TPointCurve<T> ccurve0(wall0points[i]);
        TPointCurve<T> ccurve1(wall1points[i]);

        ccurve0.redivide(int(wall0points[i].size()),wall0points[i],false,45.0,0.5);
        ccurve1.redivide(int(wall1points[i].size()),wall1points[i],false,45.0,0.5);

        // connect starts/ends, just in case
        wall0points[i].front() = wall1points[i].front() = (wall0points[i].front() + wall1points[i].front()) * 0.5;
        wall0points[i].back() = wall1points[i].back() = (wall0points[i].back() + wall1points[i].back()) * 0.5;
      }
    }

    // flat part
    makeTwoWallsAndFlatBottomTopRound(wall0points,wall1points,surfaces,false,true, 
      K1,K2round,K2top,degree,startU,endU,startV,endV);
  }
}

/** Make flat top. */
template <class T> TSplineSurface<T> *makeFlatTop(std::vector<TPoint<T>> &points0, std::vector<TPoint<T>> &points1, 
  int K1 = 40, int K2 = 8, int degree = SPLINE_DEGREE,
  CurveEndType startU = END_FREE, CurveEndType endU = END_FREE,
  CurveEndType startV = END_FREE, CurveEndType endV = END_FREE)
{
  std::vector<std::vector<TPoint<T>>> flatpoints;
  flatpoints.push_back(points0);
  flatpoints.push_back(points1);

  TSplineSurface<T> *flat = new TSplineSurface<T>(flatpoints,K1,degree,K2,degree,
    startU,endU,startV,endV); 

  return flat;
}

/** Make two walls defined by contours of points from bottom to top along Z (points.back() 
  has the highest Z). All the contours in points must have the same size. wall0 is "upper" 
  airfoil surcae, wall1 is lower. */
template <class T> void makeTwoWallsAndFlatBottomTopRound(
  std::vector<std::vector<TPoint<T>>> &wall0points,
  std::vector<std::vector<TPoint<T>>> &wall1points,
  std::vector<TSplineSurface<T> *> &surfaces, bool coverbottom, bool covertop, 
  int K1 = 40, int K2 = 40, int K2top = 8, int degree = SPLINE_DEGREE,
  CurveEndType startU = END_FREE, CurveEndType endU = END_FREE,
  CurveEndType startV = END_FREE, CurveEndType endV = END_FREE)
{
  //===== make first wall 0 (upper surface) =====
  TSplineSurface<T> *wall0 = new TSplineSurface<T>(wall0points,K1,degree,K2,degree,
    startU,endU,startV,endV);

  // one boundary for flat top
  std::vector<TPoint<T>> points0 = wall0points.back();
  std::vector<TPoint<T>> bpoints0 = wall0points.front();

  // make correct normal
  wall0->reverseU();

  //===== make second wall 1 (lower surface) =====
  TSplineSurface<T> *wall1 = new TSplineSurface<T>(wall1points,K1,degree,K2,degree,
    startU,endU,startV,endV);

  // one boundary for flat top
  std::vector<TPoint<T>> points1 = wall1points.back();
  std::vector<TPoint<T>> bpoints1 = wall1points.front();

  //===== flat bottom/top =====
  TSplineSurface<T> *flattop = nullptr;
  TSplineSurface<T> *flatbottom = nullptr;

  if (coverbottom)
  {
    // make flat bottom
    flatbottom = makeFlatTop(bpoints0,bpoints1,K1,K2top,degree,
      startU,endU,startV,endV);
  }

  if (covertop)
  {
    // make flat top
    std::reverse(points0.begin(),points0.end());
    std::reverse(points1.begin(),points1.end());

    flattop = makeFlatTop(points0,points1,K1,K2top,degree,
      startU,endU,startV,endV);
  }

  if (flatbottom)
    surfaces.push_back(flatbottom);
  surfaces.push_back(wall0);
  surfaces.push_back(wall1);
  if (flattop)
    surfaces.push_back(flattop);
}

/** Make two walls (symmetric around XZ plane) defined by contours of points from 
  bottom to top along Z (points.back() has the highest Z). All the contours in points 
  must have the same size and be symmetric around XZ plane. */
template <class T> void makeTwoWallsAndFlatBottomTop(std::vector<std::vector<TPoint<T>>> &points,
  std::vector<TSplineSurface<T> *> &surfaces, bool coverbottom, bool covertop, 
  int K1 = 40, int K2 = 40, int K2top = 8, int degree = SPLINE_DEGREE,
  CurveEndType startU = END_FREE, CurveEndType endU = END_FREE,
  CurveEndType startV = END_FREE, CurveEndType endV = END_FREE)
{
  // another symmetric wall
  std::vector<std::vector<TPoint<T>>> points1 = points;

  // transform for symmetric part arounf XZ plane
  TTransform<T> t;
  t.LoadIdentity();
  t.Resize(TPoint<T>(1.0,-1.0,1.0));

  makeTransform<T>(points1,&t);

  makeTwoWallsAndFlatBottomTopRound(points,points1,surfaces,coverbottom,covertop, 
    K1,K2,K2top,degree,startU,endU,startV,endV);
}

/** Convert blade contour of left-right XY points into middle line. Every output point contains :
  X - Z along span
  Y - X(Z) of middle line 
  Z - blade HALF-width from middle line,
  W - twist angle around middle curve in degrees. */
template <class T> void contourToBlade(std::vector<std::array<TPoint<T>,2>> &contour, std::vector<TPoint<T>> &blade, T coef)
{
  blade.clear();

  TPoint<T> offset = (contour[0][0] + contour[0][1]) * 0.5;

  for (auto c : contour)
  {
    TPoint<T> c0 = c[0] - offset;
    TPoint<T> c1 = c[1] - offset;

    TPoint<T> middle = (c0 + c1) * 0.5;
    TPoint<T> d = c1 - c0;
    TPoint<T> p;
    p.X = middle.Y;
    p.Y = middle.X;
    p.Z = (!d) * 0.5;
    p.W = 0.0;

    blade.push_back(p * coef);
  }
}
  

}