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
#include "tsplinesurface.h"
#include "toperations.h"


//#define DEBUG_BLOCKS
#ifdef NDEBUG
  #undef DEBUG_BLOCKS
#endif

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
  int degree = 5, int integration = OTHER_INTEGRATION, T endrefinementcoef = 0.5)
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

/** Make surface of revolution points with elliptical cross-sections in XY plane 
  with axis along Z. The contour must be defined in Z-X coordinates, normally
  with decreasing Z. numpoints is that along circumference. 
  U is along circumference from adegfrom to adegto, V along axis of symmetry Z. */
template <class T> void makeSurfaceOfRevolution(std::vector<TPoint<T>> &contour, 
  int numpoints, T adegfrom, T adegto, std::vector<std::vector<TPoint<T>>> &points)
{
  int numsections = int(contour.size());

  for (int i = 0; i < numsections; i++)
  {
    T r = contour[i].X;

    std::vector<TPoint<T>> section;
    makeEllipseXY(numpoints,TPoint<T>(),r,r,section,adegfrom,adegto);

    TTransform<T> t;
    t.Translate(TPoint<T>(0.0,0.0,contour[i].Z));
    makeTransform(section,&t);

    points.push_back(section);
  }
}

/** Make surface of revolution around Z (single face). 
  The contour must be defined in Z-X coordinates,
  normally with decreasing Z. numpoints is that along circumference. 
  U is circumferential, V - along Z axis. */
template <class T> TSplineSurface<T> *makeSurfaceOfRevolution(std::vector<TPoint<T>> &contour,
  int numpoints, T adegfrom, T adegto, int K1, int K2, 
  int M1 = SPLINE_DEGREE, int M2 = SPLINE_DEGREE,
  CurveEndType startU = END_CLAMPED, CurveEndType endU = END_CLAMPED, // must be round
  CurveEndType startV = END_FREE, CurveEndType endV = END_FREE) 
{
  std::vector<std::vector<TPoint<T>>> points;
  makeSurfaceOfRevolution(contour,numpoints,adegfrom,adegto,points);

  // create cylindrical surface
  TSplineSurface<T> *surface = new TSplineSurface<T>(points,K1,M1,K2,M2,
    startU,endU,startV,endV);

  return surface;
}

/** Make surfaces of revolution around Z (multiple faces around Z from 0 t0 360 deg). 
  The contour must be defined in Z-X coordinates,
  normally with decreasing Z. numfaces is that along circumference. 
  U is circumferential, V - along Z axis. */
template <class T> void makeSurfacesOfRevolution(std::vector<TPoint<T>> &contour,
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

    TSplineSurface<T> *surface = makeSurfaceOfRevolution(contour,
      pointsperface,a0,a1,K1,K2,M1,M2,startU,endU,startV,endV);

    surfaces.push_back(surface);
  }
}

/** Intersect surfaces to make a solid. */
template <class T> void makeSolid(
  std::vector<TSplineSurface<T> *> &surfaces0,
  std::vector<TSplineSurface<T> *> &surfaces1,
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV0,
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV1,
  T tolerance, T &bigtolerance, T parmtolerance = PARM_TOLERANCE, int manypoints = MANY_POINTS2D,
  bool improveintersection = true, int maxiter = 20, T relaxcoef = 0.5)
{
  boundariesUV0.clear();
  boundariesUV1.clear();

  // step 1 : make outer boundaries
  for (int i = 0; i < int(surfaces0.size()); i++)
  {
    std::vector<std::vector<TPoint<T>>> loop;
    surfaces0[i]->closeOuterBoundaryLoop(loop,manypoints);

    boundariesUV0.push_back(std::vector<std::vector<std::vector<tcad::TPoint<T>>>>());
    boundariesUV0.back().push_back(loop);
  }
  for (int i = 0; i < int(surfaces1.size()); i++)
  {
    std::vector<std::vector<TPoint<T>>> loop;
    surfaces1[i]->closeOuterBoundaryLoop(loop,manypoints);

    boundariesUV1.push_back(std::vector<std::vector<std::vector<tcad::TPoint<T>>>>());
    boundariesUV1.back().push_back(loop);
  }

  // "delayed" boundary pieces
  std::vector<std::vector<std::vector<std::vector<TPoint<T>>>>> boundaries0(surfaces0.size());
  std::vector<std::vector<std::vector<std::vector<TPoint<T>>>>> boundaries1(surfaces1.size());

  // estmate big tolerance as max difference between boundary lines to be used in making
  // soldi downstream
  bigtolerance = 0.0;
  T maxseglen = 0.0;

  // step 2 : intersect with other surfaces
  for (int i = 0; i < int(surfaces0.size()); i++)
  {
    for (int j = 0; j < int(surfaces1.size()); j++)
    {
      std::vector<std::vector<TPoint<T>>> intersections; 
      std::vector<std::vector<TPoint<T>>> boundary0,boundary1;

// outputDebugString("================================================");

      bool ok = surfaces0[i]->intersect(*surfaces1[j],intersections,boundary0,boundary1,
        tolerance,parmtolerance,
        manypoints,manypoints,1.0,1.0,1.0,1.0, 
#ifdef _DEBUG
        manypoints,manypoints,1.0,1.0,1.0,1.0,false);
//!!!!!!        manypoints,manypoints,1.0,1.0,1.0,1.0,true);
#else
        manypoints,manypoints,1.0,1.0,1.0,1.0,false);
#endif

      if (ok)
      {
        if (improveintersection)
        {
          for (int k = 0; k < boundary0.size(); k++)
          {
            for (int l = 0; l < int(boundary0[k].size()); l++)
            {
              TPoint<T> parms(boundary0[k][l].X,boundary0[k][l].Y,boundary1[k][l].X,boundary1[k][l].Y);

              if (improveIntersection<T>(surfaces0[i],surfaces1[j],parms,maxiter,relaxcoef,parmtolerance))
              {
                boundary0[k][l].X = parms.X;
                boundary0[k][l].Y = parms.Y;
                boundary1[k][l].X = parms.Z;
                boundary1[k][l].Y = parms.W;
              }
            }
          }
        }

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

          T segmin = 0.0;
          T segmax = 0.0;
          segmentLenMinMax<T>(boundary0[k],segmin,segmax);

          maxseglen = std::max(maxseglen,segmax);
        }

        bigtolerance = std::max<T>(bigtolerance,maxdiff);

#ifdef DEBUG_BLOCKS
        outputDebugString(std::string("maxdiff = ") + to_string(maxdiff) + 
          std::string(" bigtolerance = ") + to_string(bigtolerance) +
          std::string(" maxseglen = ") + to_string(maxseglen));
#endif
      }

      for (auto &b : boundary1)
      {
        std::reverse(b.begin(),b.end());
      }
      std::reverse(boundary1.begin(),boundary1.end());

      if (ok)
      {
        // close boundary on first surface
        bool ok0 = surfaces0[i]->closeBoundaryLoop(boundary0,boundariesUV0[i],tolerance,
          maxseglen * 2.0,PARM_TOLERANCE,manypoints);
        if (!ok0)
        {
          boundaries0[i].push_back(boundary0);
        }
        bool ok1 = surfaces1[j]->closeBoundaryLoop(boundary1,boundariesUV1[j],tolerance,
          maxseglen * 2.0,PARM_TOLERANCE,manypoints);
        if (!ok1)
        {
          boundaries1[j].push_back(boundary1);
        }
      }
    }
  }

  // combine delayed boundaries
  for (int i = 0; i < int(surfaces0.size()); i++)
  {
    if (!boundaries0[i].empty())
    {
      std::vector<std::vector<TPoint<T>>> boundary;
      for (int j = 0; j < int(boundaries0[i].size()); j++)
      {
        for (int k = 0; k < int(boundaries0[i][j].size()); k++)
        {
          boundary.push_back(boundaries0[i][j][k]);
        }
      }

      bool ok = surfaces0[i]->closeBoundaryLoop(boundary,boundariesUV0[i],tolerance,
        maxseglen * 2.0,PARM_TOLERANCE,manypoints);
    }
  }

  for (int i = 0; i < int(surfaces1.size()); i++)
  {
    if (!boundaries1[i].empty())
    {
      std::vector<std::vector<TPoint<T>>> boundary;
      for (int j = 0; j < int(boundaries1[i].size()); j++)
      {
        for (int k = 0; k < int(boundaries1[i][j].size()); k++)
        {
          boundary.push_back(boundaries1[i][j][k]);
        }
      }

      bool ok = surfaces1[i]->closeBoundaryLoop(boundary,boundariesUV1[i],tolerance,
        maxseglen * 2.0,PARM_TOLERANCE,manypoints);
    }
  }
}

}