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

  tbrep.h

  Common class for breps (surfaces + trimming curves).

*******************************************************************************/

#pragma once

#include "tpoints.h"
#include "tmisc.h"
#include "tsplinesurface.h"
#include "tsolid.h"
#include "tblocks.h"
#include "export.h"

namespace tcad {

/** Common class for breps (surfaces + trimming curves). For boolean operations 
  and ease of use. */
template <class T> struct TBrep {
public:

  //===== Data ======

  // untrimmed surfaces
  std::vector<TSplineSurface<T> *> surfaces;

  // trimming curves for surfaces
  //surface   loop        piece      UV points in X,Y  
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> boundariesUV;

  // tolerance for the whole model
  T tolerance = 0.0;

  //===== Contruction/destruction =====

  /** Constructor. */
  TBrep()
  {
  }

  /** Destructor. */
  ~TBrep()
  {
    deleteSurfaces(surfaces);
  }

  /** Copy constructor. */
  TBrep(const TBrep &other)  
  {
    tolerance = other.tolerance;
  
    deleteSurfaces(surfaces);

    for (int i = 0; i < int(other.surfaces.size()); i++)
    {
      TSplineSurface<T> *surface = new TSplineSurface<T>(*(other.surfaces[i]));
      surfaces.push_back(surface);
    }

    boundariesUV = other.boundariesUV;
  }

  /** Assignment operator. */
  TBrep &operator = (const TBrep &other)  
  {
    tolerance = other.tolerance;

    deleteSurfaces(surfaces);

    for (int i = 0; i < int(other.surfaces.size()); i++)
    {
      TSplineSurface<T> *surface = new TSplineSurface<T>(*(other.surfaces[i]));
      surfaces.push_back(surface);
    }

    boundariesUV = other.boundariesUV;

    return *this;
  }

  /** Constructor with single tolerance. */
  TBrep(T ptolerance)
  {
    tolerance = ptolerance;
  }

  /** Constructor, takes ownership of psurfaces pointers, do not delete them. */
  TBrep(std::vector<TSplineSurface<T> *> &psurfaces, 
    std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &pboundariesUV,
    T ptolerance = -1.0, T parmtolerance = PARM_TOLERANCE)
  {
    surfaces = psurfaces;
    boundariesUV = pboundariesUV;

    if (ptolerance > 0.0)
    {
      tolerance = ptolerance;
    } else
    {
      std::pair<TPoint<T>,TPoint<T>> minmax = getMinMax(surfaces);
      tolerance = !(minmax.second - minmax.first) * parmtolerance;
    }
  }

  /** Constructor, takes ownership of psurfaces pointers, do not delete them. */
  TBrep(std::vector<TSplineSurface<T> *> &psurfaces, 
    T ptolerance = -1.0, T parmtolerance = PARM_TOLERANCE)
  {
    surfaces = psurfaces;

    // make boundaries
    prepareOuterBoundary<T>(surfaces,boundariesUV);

    if (ptolerance > 0.0)
    {
      tolerance = ptolerance;
    } else
    {
      std::pair<TPoint<T>,TPoint<T>> minmax = getMinMax(surfaces);
      tolerance = !(minmax.second - minmax.first) * parmtolerance;
    }
  }

  /** Clear all. */
  void clear()
  {
    deleteSurfaces(surfaces);
    boundariesUV.clear();
    //!!! incorrect tolerance = 0.0;
  }

  /** Set tolerance. */
  void setTolerance(T ptolerance)
  {
    tolerance = ptolerance;
  }

  void closeOuterBoundary()
  {
    boundariesUV.clear();
    tcad::closeOuterBoundary<T>(surfaces,boundariesUV);

    assert(surfaces.size() == boundariesUV.size());
  }

  void clearOuterBoundary()
  {
    // make boundaries
    tcad::clearOuterBoundary<T>(surfaces,boundariesUV);
  }

  //===== Booleans =====

  /** Booleans : union. */
  TBrep operator + (TBrep &other) 
  {
    assert(tolerance > TOLERANCE(T));
    assert(other.tolerance > TOLERANCE(T));

    // result
    std::vector<TSplineSurface<T> *> rsurfaces;
    std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> rboundariesUV;

    makeTrimming<T>(surfaces,other.surfaces,boundariesUV,other.boundariesUV,
      UNITE,rsurfaces,rboundariesUV,tolerance);

    TBrep result(rsurfaces,rboundariesUV,std::max<T>(tolerance,other.tolerance));

    return result;
  }

  /** Booleans : subtraction. */
  TBrep operator - (TBrep &other)
  {
    assert(tolerance > TOLERANCE(T));
    assert(other.tolerance > TOLERANCE(T));

    // result
    std::vector<TSplineSurface<T> *> rsurfaces;
    std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> rboundariesUV;

    makeTrimming<T>(surfaces,other.surfaces,boundariesUV,other.boundariesUV,
      SUBTRACT,rsurfaces,rboundariesUV,tolerance);

    TBrep result(rsurfaces,rboundariesUV,std::max<T>(tolerance,other.tolerance));

    return result;
  }

  /** Booleans : intersection. */
  TBrep operator ^ (TBrep &other)
  {
    assert(tolerance > TOLERANCE(T));
    assert(other.tolerance > TOLERANCE(T));

    // result
    std::vector<TSplineSurface<T> *> rsurfaces;
    std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> rboundariesUV;

    makeTrimming<T>(surfaces,other.surfaces,boundariesUV,other.boundariesUV,
      INTERSECT,rsurfaces,rboundariesUV,tolerance);

    TBrep result(rsurfaces,rboundariesUV,std::max<T>(tolerance,other.tolerance));

    return result;
  }

  //==== Transforms =====

  /** Transform. */
  void makeTransform(TTransform<T> *t)
  {
    tcad::makeTransform<T>(surfaces,t);
  }

  //===== Solid? =====

  /** Solid? */
  bool isSolid()
  {
    std::vector<TPoint<T>> vertices;  
    std::vector<std::array<LINT,11>> edges;

    return createSolidEdges(surfaces,boundariesUV,vertices,edges,tolerance);
  }

  //===== Faces =====

  /** Delete face from surfaces. */
  void deleteFace(int index)
  {
    DELETE_CLASS(surfaces[index]);
    surfaces.erase(surfaces.begin() + index);
    boundariesUV.erase(boundariesUV.begin() + index);
  }

  /** Add faces from another brep. */
  void addFaces(TBrep<T> &other)
  {
    for (int i = 0; i < int(other.surfaces.size()); i++)
    {
      surfaces.push_back(new TSplineSurface<T>(*other.surfaces[i]));
    }
    boundariesUV.insert(boundariesUV.end(),other.boundariesUV.begin(),other.boundariesUV.end());
  }

  //===== Blocks =====

  /** Max box between min and max. */
  void makeBox(TPoint<T> min, TPoint<T> max, 
    std::string name = "",
    int m1 = SPLINE_DEGREE, int m2 = SPLINE_DEGREE, int m3 = SPLINE_DEGREE,
    int numpointsU = 11, int numpointsV = 11, int numpointsW = 11)
  {
    // clear all
    //!!! do not clear();

    // make 3D box
    TSplineVolume<T> box(min,max,m1,m2,m3,numpointsU,numpointsV,numpointsW);

    // extract faces
    for (int i = 0; i < 6; i++)
    {
      TSplineSurface<T> *face = box.getFace(i);
      surfaces.push_back(face);
    }

    // name for debugging
    nameSurfaces(name);

    // make boundaries
    prepareOuterBoundary<T>(surfaces,boundariesUV);
  }

  /** Make surface of revolution around Z (multiple faces around Z from 0 t0 360 deg). 
    The contour has axial coordinate in countourAxisCoord and radial coordinate in countourRadialCoord.
    numfaces is that along circumference. 
    U is circumferential, V - axial parameteric directions. */
  void makeAxisymmetricBody(std::vector<std::vector<TPoint<T>>> &contour,
    Axes countourAxialCoord, Axes countourRadialCoord,
    std::string name = "",
    T angdegfrom = 0.0, T angdegto = 360.0, 
    int numfaces = 16, int pointsperface = 17, int K1 = 16, int K2 = 16, 
    int M1 = SPLINE_DEGREE, int M2 = SPLINE_DEGREE,
    CurveEndType startU = END_CLAMPED, CurveEndType endU = END_CLAMPED, // must be round
    CurveEndType startV = END_FREE, CurveEndType endV = END_FREE)
  {
    // clear all
    //!!! do not clear();

    // make surfaces
    makeSurfacesOfRevolution<T>(contour,countourAxialCoord,countourRadialCoord,numfaces,pointsperface,K1,K2, 
      surfaces,angdegfrom,angdegto,M1,M2,startU,endU,startV,endV);

    // name for debugging
    nameSurfaces(name);

    // make boundaries
    prepareOuterBoundary<T>(surfaces,boundariesUV);
  }

  /** Make sphere of radius R. */
  void makeSphere(T R,
    std::string name = "",
    // for faces around
    int numfaces = 2, int pointsperface = 64, 
    // for a countour along
    int numcontourpoints = MANY_POINTS2D, 
    T angdegfrom = 0.0, T angdegto = 360.0,
    int K1 = 32, int K2 = 32, 
    int M1 = SPLINE_DEGREE, int M2 = SPLINE_DEGREE,
    CurveEndType startU = END_CLAMPED, CurveEndType endU = END_CLAMPED, // must be round
    CurveEndType startV = END_CLAMPED, CurveEndType endV = END_CLAMPED)
  {
    std::vector<std::vector<TPoint<T>>> contour;

    // in XY, around (0,0)
    makeEllipseXY<T>(numcontourpoints,numfaces,TPoint<T>(),R,R,contour,0.0,180.0);

    // contour goes with decreasing Z (ellipse X from 0 to 180), so
    reverse(contour); //!!!

    // contour in XY
    makeAxisymmetricBody(contour,AxisX,AxisY,
      name,angdegfrom,angdegto,numfaces * 2,pointsperface,K1,K2,M1,M2,startU,endU,startV,endV);
  }

  /** Make ellipsoid with a,b,c as semiaxes. */
  void makeEllipsoid(T a, T b, T c,
    std::string name = "",
    int numfaces = 7, int pointsperface = 8, 
    int numcontourpoints = MANY_POINTS2D, 
    int K1 = 32, int K2 = 32, 
    int M1 = SPLINE_DEGREE, int M2 = SPLINE_DEGREE,
    CurveEndType startU = END_CLAMPED, CurveEndType endU = END_CLAMPED, // must be round
    CurveEndType startV = END_CLAMPED, CurveEndType endV = END_CLAMPED)
  {
    makeSphere(1.0,name,numfaces,pointsperface,numcontourpoints,
      K1,K2,M1,M2,startU,endU,startV,endV);

    TTransform<T> t;
    t.Resize(TPoint<T>(a,b,c));

    makeTransform(surfaces,&t);
  }

  /** Make cylinder of radius R of length L around Z with centre at 0,0,0. */
  void makeCylinder(T L, T Rbottom, T Rtop,
    std::string name = "",
    // num faces around
    int numfaces = 4, 
    // num faces along L
    int numfacesL = 4,
    // we go from bottom (min Z)
    int numfacesbottom = 0,
    // to top (max Z)
    int numfacestop = 0,
    // points to construct spline edge
    int pointsperface = 64, 
    // for a countour along
    int numcontourpoints = MANY_POINTS2D, 
    T angdegfrom = 0.0, T angdegto = 360.0,    
    int K1 = 32, int K2 = 32, 
    int M1 = SPLINE_DEGREE, int M2 = SPLINE_DEGREE,
    CurveEndType startU = END_CLAMPED, CurveEndType endU = END_CLAMPED, // must be round
    CurveEndType startV = END_FIXED, CurveEndType endV = END_FIXED)
  {
    assert(numfaces >= 2);
    assert(numfacesL >= 1);

    std::vector<std::vector<TPoint<T>>> contour;
    std::vector<TPoint<T>> points;

    if (numfacesbottom)
    {
      makeStraightLine<T>(TPoint<T>(-L * 0.5,0.0),TPoint<T>(-L * 0.5,Rbottom),numfacesbottom + 1,points);
      contour.push_back(points);
    }

    makeStraightLine<T>(TPoint<T>(-L * 0.5,Rbottom),TPoint<T>(L * 0.5,Rtop),numfacesL + 1,points);
    contour.push_back(points);

    if (numfacestop)
    {
      makeStraightLine<T>(TPoint<T>(L * 0.5,Rtop),TPoint<T>(L * 0.5,0.0),numfacestop + 1,points);
      contour.push_back(points);
    }

    // contour in XY
    makeAxisymmetricBody(contour,AxisX,AxisY,
      name,angdegfrom,angdegto,numfaces,pointsperface,K1,K2,M1,M2,startU,endU,startV,endV);
  }

  //===== Export =====

  /** Clear names to display surface number in IGES. */
  void clearNames()
  {
    tcad::clearNames<T>(surfaces);
  }

  /** Name surfaces for debugging and IGES. */
  void nameSurfaces(std::string prefix)
  {
    tcad::nameSurfaces(surfaces,prefix);
  }

  /** Save trimmed surfaces. */
  bool saveSurfacesIges(const std::string &filename, int splinedegree = SPLINE_DEGREE)
  {
    return ::saveTrimmedSurfacesIges<T>(surfaces,boundariesUV,filename);
  }

  /** Save trimmed surfaces as solid in IGES. All surfaces must have closed boundaries. */
  bool saveSolidIges(const std::string &filename,
    T tolerance, T parmtolerance = PARM_TOLERANCE, int splinedegree = SPLINE_DEGREE, int numdigits = 18, 
    std::vector<std::vector<TPoint<T>>> *badedges = nullptr)
  {
    return ::saveSolidIges<T>(surfaces,boundariesUV,filename,
      tolerance,parmtolerance,SPLINE_DEGREE,18,badedges);
  }

};

}