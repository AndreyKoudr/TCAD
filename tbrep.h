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
  void makeTransform(TTransform<T> *t, int index0 = -1, int index1 = -1)
  {
    tcad::makeTransform<T>(surfaces,t,index0,index1);
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
    if (surfaces.size() == boundariesUV.size())
      boundariesUV.erase(boundariesUV.begin() + index);
    DELETE_CLASS(surfaces[index]);
    surfaces.erase(surfaces.begin() + index);
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

  /** Exclude identical faces. n^2. */
  void excludeInnerFaces(T tolerance, 
    int numpointsU = MANY_POINTS2D, int numpointsV = MANY_POINTS2D)
  {
    int attempts = int(surfaces.size());

    for (int a = 0; a < attempts; a++)
    {
      bool found = false;
      for (int i = 0; i < int(surfaces.size()); i++)
      {
        for (int j = i + 1; j < int(surfaces.size()); j++)
        {
          bool reversedU,reversedV;

          if (surfaces[j]->equal(*surfaces[i],tolerance,reversedU,reversedV,
            numpointsU,numpointsV))
          {
            DELETE_CLASS(surfaces[j]);
            surfaces.erase(surfaces.begin() + j);
            // i is less, delete after
            DELETE_CLASS(surfaces[i]);
            surfaces.erase(surfaces.begin() + i);

            found = true;
          }
        }
      }

      if (!found)
        break;
    }

    // make empty boundaries
    clearOuterBoundary();
    prepareOuterBoundary<T>(surfaces,boundariesUV);
  }

  //===== Blocks =====

  /** Max box between min and max. bottom is face with Zmin, top - with Zmax. */
  void makeBox(TPoint<T> min, TPoint<T> max, 
    std::string name = "",
    int facemask = 0xFF,
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

    // delete unmasked faces
    int n = int(surfaces.size());
    for (int i = 5; i >= 0; i--)
    {
      if (!bitSet<T>(facemask,i))
      {
        deleteFace(n - 6 + i); 
      }
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

  /** Make cylinder of radius R of length L around Z with centre at 0,0,0. if numfacesbottom = 0 or
    numfacestop = 0, there is no bottom/top correspondingly. */
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

  /** Make a wing.

      airfoil. 
    points in XY, these points start from TE, go around LE along 
    the upper surface and then back to TE, like this

    template <class T> std::vector<TPoint<T>> NACA0009 = {
      TPoint<T>(1.00000, 0.0),
      TPoint<T>(0.99572, 0.00057),
      TPoint<T>(0.98296, 0.00218),
      ...
      TPoint<T>(0.00428, 0.00767),
      TPoint<T>(0.00107, 0.00349),
      TPoint<T>(0.0,     0.0),
      TPoint<T>(0.00107, -0.00349),
      TPoint<T>(0.00428, -0.00767),
      ...
      TPoint<T>(0.98296, -0.00218),
      TPoint<T>(0.99572, -0.00057),
      TPoint<T>(1.00000, 0.0)
    };
    points must be closed! The chord must be close to 1.0.

      contour. 
    Wing contour of left-right points in XZ plane. Z coordinate is along span, like this:

    template <class T> std::vector<std::array<TPoint<T>,2>> KiloWing1Contour = {
      { TPoint<T>(0.0 * KHCOEF, 0.0,  0.0 * KHCOEF), TPoint<T>(22.0 * KHCOEF, 0.0, 5.5 * KHCOEF) },
      { TPoint<T>(0.0 * KHCOEF, 0.0,  2.5 * KHCOEF), TPoint<T>(18.0 * KHCOEF, 0.0, 7.5 * KHCOEF) },
      { TPoint<T>(0.0 * KHCOEF, 0.0,  4.0 * KHCOEF), TPoint<T>(17.0 * KHCOEF, 0.0, 9.0 * KHCOEF) },
      { TPoint<T>(0.0 * KHCOEF, 0.0,  8.5 * KHCOEF), TPoint<T>(16.0 * KHCOEF, 0.0, 11.5 * KHCOEF) },
      { TPoint<T>(0.0 * KHCOEF, 0.0, 12.5 * KHCOEF), TPoint<T>(15.5 * KHCOEF, 0.0, 13.0 * KHCOEF) },
      { TPoint<T>(0.0 * KHCOEF, 0.0, 14.0 * KHCOEF), TPoint<T>(15.0 * KHCOEF, 0.0, 14.0 * KHCOEF) }};
    contour is in real model coordinates, not [0..1].

    closetoptip - true if to close top tip with topradius, no fillet if topradius == 0.0
    closebottomtip - true if to close bottom tip with bottomradius, no fillet if bottomradius == 0.0
    K1 - number of spline intervals along U (along chord)
    K2 - number of spline intervals along V (along span)
  */
  bool makeWing(std::vector<TPoint<T>> &airfoil, bool roundLE, bool roundTE,
    std::vector<std::array<TPoint<T>,2>> &contour, bool closetoptip, bool closebottomtip,
    T topradius, T bottomradius,
    T tolerance, T parmtolerance = PARM_TOLERANCE, int K1 = 40, int K2 = 20)
  {
    // step 1 : make airfoil in XY plane

    if (!closed<T>(airfoil,parmtolerance))
      return false;

    this->tolerance = tolerance;

    // make upper and lower curves, still in XY, make LE lowest X (reverse)
    std::vector<TPoint<T>> upperlower; 
    std::pair<T,T> res = makeAirfoilPointsXY<T>(airfoil,roundLE,roundTE,K1,upperlower,true);

    // step 2 : make 3D camber surface and thickness from upperlower and
    // contour as a regular set of points

    std::vector<std::vector<TPoint<T>>> camberpoints; 
    std::vector<std::vector<T>> thickness;
    makeBladeCamberThickness<T>(upperlower,contour,K1,K2,camberpoints,thickness,tolerance); 

    // create camber surface
    TSplineSurface<T> cambersurface(camberpoints,K1,SPLINE_DEGREE,K2,SPLINE_DEGREE,
      END_CLAMPED,END_CLAMPED,END_FREE,END_FREE);

    //#define DEBUG_DIR std::string("C:/AndrewK/MyProjects2/temp/")
    //::saveSurfaceIges(&cambersurface,DEBUG_DIR + "Camber surface.iges");

    // step 3 : make 3D wing of two upper and lower spline surfaces

    clear(); // this is important as the first two surfaces are upper/lower airfoil

    makeAirfoil<T>(camberpoints,thickness,surfaces,SPLINE_DEGREE,SPLINE_DEGREE,
      END_CLAMPED,END_CLAMPED,END_FREE,END_FREE);

    //::saveSurfacesIges(surfaces,DEBUG_DIR + "Surfaces.iges"); 

    // step 4 : close and round tips

    if (closetoptip)
    {
      makeFlatRoundedTop(surfaces[0],surfaces[1],topradius,surfaces);
    }

    if (closebottomtip)
    {
      surfaces[0]->reverseU();
      surfaces[1]->reverseU();
      surfaces[0]->reverseV();
      surfaces[1]->reverseV();

      makeFlatRoundedTop(surfaces[0],surfaces[1],bottomradius,surfaces);

      surfaces[0]->reverseU();
      surfaces[1]->reverseU();
      surfaces[0]->reverseV();
      surfaces[1]->reverseV();
    }

    // make boundaries
    prepareOuterBoundary<T>(surfaces,boundariesUV);

    return true;
  }

  /** Make box with with four rounded corners (int XY plane) of three boxes 
    and four quarter-cylinders. bottom is faces with Zmin, top - with Zmax. */
  bool makeBoxRounded4(T dx, T dy, T dz, T radius,
    bool makebottom, bool maketop,
    int numfacesaround = 4, T parmtolerance = PARM_TOLERANCE,     
    int m1 = SPLINE_DEGREE, int m2 = SPLINE_DEGREE, int m3 = SPLINE_DEGREE,
    int numpointsU = 11, int numpointsV = 11, int numpointsW = 11,
    int numcontourpoints = MANY_POINTS2D)
  {
    // tolerance and radius must be non-zero
    assert(tolerance >= TOLERANCE(T));
    assert(radius >= tolerance);

    if (tolerance < TOLERANCE(T) || radius < tolerance)
      return false;

    // boundaries of 9 XY faces
    T x[4] = {-dx * 0.5, -dx * 0.5 + radius, +dx * 0.5 - radius, +dx * 0.5};
    T y[4] = {-dy * 0.5, -dy * 0.5 + radius, +dy * 0.5 - radius, +dy * 0.5};

    // make 5 boxes
    int count = 0;
    for (int i = 0; i < 3; i++)
    {
      T y0 = y[i];
      T y1 = y[i + 1];

      for (int j = 0; j < 3; j++)
      {
        if (count == 0 || count == 2 || count == 6 || count == 8)
        {
          count++;
          continue;
        }

        T x0 = x[j];
        T x1 = x[j + 1];

        int mask = 0xFF;

        if (!makebottom)
          mask = clearBit<T>(mask,4);
        if (!maketop)
          mask = clearBit<T>(mask,5);

        if (count == 1 || count == 7)
        {
          mask = clearBit<T>(mask,0);
          mask = clearBit<T>(mask,1);
        } else if (count == 3 || count == 5)
        {
          mask = clearBit<T>(mask,2);
          mask = clearBit<T>(mask,3);
        }

        makeBox(
          TPoint<T>(x0,y0,-dz * 0.5),
          TPoint<T>(x1,y1,+dz * 0.5), 
          "",mask,1,m2,m3,numpointsU,numpointsV,numpointsW); //!!! m1 = 1

        count++;
      }
    }

    TTransform<T> t;

    // add 4 quarter-cylinders at corners

    int last = int(surfaces.size());
    makeCylinder(dz,radius,radius,"",numfacesaround,1,int(makebottom),int(maketop),64,
      numcontourpoints,180.0,270.0,32,32,m1,m2);
    t.LoadIdentity();
    t.Translate(TPoint<T>(x[1],y[1]));
    makeTransform(&t,last,int(surfaces.size() - 1));

    last = int(surfaces.size());
    makeCylinder(dz,radius,radius,"",numfacesaround,1,int(makebottom),int(maketop),64,
      numcontourpoints,270.0,360.0,32,32,m1,m2);
    t.LoadIdentity();
    t.Translate(TPoint<T>(x[2],y[1]));
    makeTransform(&t,last,int(surfaces.size() - 1));

    last = int(surfaces.size());
    makeCylinder(dz,radius,radius,"",numfacesaround,1,int(makebottom),int(maketop),64,
      numcontourpoints,0.0,90.0,32,32,m1,m2);
    t.LoadIdentity();
    t.Translate(TPoint<T>(x[2],y[2]));
    makeTransform(&t,last,int(surfaces.size() - 1));

    last = int(surfaces.size());
    makeCylinder(dz,radius,radius,"",numfacesaround,1,int(makebottom),int(maketop),64,
      numcontourpoints,90.0,180.0,32,32,m1,m2);
    t.LoadIdentity();
    t.Translate(TPoint<T>(x[1],y[2]));
    makeTransform(&t,last,int(surfaces.size() - 1));

    excludeInnerFaces(tolerance,3,3); //!!!

    return true;
  }

  /** Make box with with two rounded ends (int XY plane) of two boxes 
    and two half-cylinders. bottom is faces with Zmin, top - with Zmax. */
  bool makeBoxRounded2(T dx, T radius, T dz,
    bool makebottom, bool maketop,
    bool roundleft = true, bool roundright = true,
    int numfacesaround = 4, T parmtolerance = PARM_TOLERANCE,     
    int m1 = SPLINE_DEGREE, int m2 = SPLINE_DEGREE, int m3 = SPLINE_DEGREE,
    int numpointsU = 11, int numpointsV = 11, int numpointsW = 11,
    int numcontourpoints = MANY_POINTS2D)
  {
    // tolerance and radius must be non-zero
    assert(tolerance >= TOLERANCE(T));
    assert(radius >= tolerance);

    if (tolerance < TOLERANCE(T) || radius < tolerance)
      return false;

    // boundaries of 6 XY faces
    T x[4] = {-dx * 0.5, -dx * 0.5 + radius, +dx * 0.5 - radius, +dx * 0.5};
    T y[3] = {-radius, 0.0, +radius};

    // make 2 boxes
    int count = 0;
    for (int i = 0; i < 2; i++)
    {
      T y0 = y[i];
      T y1 = y[i + 1];

      for (int j = 0; j < 3; j++)
      {
        if (count == 0 || count == 2 || count == 3 || count == 5)
        {
          count++;
          continue;
        }

        T x0 = x[j];
        T x1 = x[j + 1];

        int mask = 0xFF;

        if (!makebottom)
          mask = clearBit<T>(mask,4);
        if (!maketop)
          mask = clearBit<T>(mask,5);
        if (roundleft)
          mask = clearBit<T>(mask,0);
        if (roundright)
          mask = clearBit<T>(mask,1);

        if (count == 1)
        {
          mask = clearBit<T>(mask,3);
        } else if (count == 4)
        {
          mask = clearBit<T>(mask,2);
        }

        makeBox(
          TPoint<T>(x0,y0,-dz * 0.5),
          TPoint<T>(x1,y1,+dz * 0.5), 
          "",mask,1,m2,m3,numpointsU,numpointsV,numpointsW); //!!! m1 = 1

        count++;
      }
    }

    TTransform<T> t;

    // add 4 quarter-cylinders at corners

    int last = int(surfaces.size());

    if (roundleft)
    {
      makeCylinder(dz,radius,radius,"",numfacesaround,1,int(makebottom),int(maketop),64,
        numcontourpoints,180.0,270.0,32,32,m1,m2);
      t.LoadIdentity();
      t.Translate(TPoint<T>(x[1],y[1]));
      makeTransform(&t,last,int(surfaces.size() - 1));

      last = int(surfaces.size());
      makeCylinder(dz,radius,radius,"",numfacesaround,1,int(makebottom),int(maketop),64,
        numcontourpoints,90.0,180.0,32,32,m1,m2);
      t.LoadIdentity();
      t.Translate(TPoint<T>(x[1],y[1]));
      makeTransform(&t,last,int(surfaces.size() - 1));
    }

    if (roundright)
    {
      last = int(surfaces.size());
      makeCylinder(dz,radius,radius,"",numfacesaround,1,int(makebottom),int(maketop),64,
        numcontourpoints,270.0,360.0,32,32,m1,m2);
      t.LoadIdentity();
      t.Translate(TPoint<T>(x[2],y[1]));
      makeTransform(&t,last,int(surfaces.size() - 1));

      last = int(surfaces.size());
      makeCylinder(dz,radius,radius,"",numfacesaround,1,int(makebottom),int(maketop),64,
        numcontourpoints,0.0,90.0,32,32,m1,m2);
      t.LoadIdentity();
      t.Translate(TPoint<T>(x[2],y[1]));
      makeTransform(&t,last,int(surfaces.size() - 1));
    }

  //  excludeInnerFaces(tolerance,3,3); //!!! not needed

    return true;
  }

  /** Make a letter from 8x14 font, every pixel size dx x dy x dz is a box,
    bottom, top - in Z direction. LOWER!!! left corner is (0,0,0). */
  bool makeLetter8x14(unsigned char letter, T dx, T dy, T dz,
    bool makebottom, bool maketop, T tolerance, 
    T parmtolerance = PARM_TOLERANCE,     
    int m1 = SPLINE_DEGREE, int m2 = SPLINE_DEGREE, int m3 = SPLINE_DEGREE,
    int numpointsU = 11, int numpointsV = 11, int numpointsW = 11,
    int numcontourpoints = MANY_POINTS2D)
  {
    std::array<unsigned char,14> &bits = Font8x14<T>[letter];

    // must be 14
    int bytesperchar = sizeof(Font8x14<T>[letter]);
    assert(bytesperchar == 14);

    // rows, upside down, LOWER!!! left corner is (0,0,0)
    for (int i = 0; i < bytesperchar; i++)
    {
      // pixels
      for (int j = 0; j < 8; j++)
      {
        if (bitSet<T>(bits[i],7 - j)) // bits are reverse
        {
          int mask = 0xFF;

          if (!makebottom)
            mask = clearBit<T>(mask,4);
          if (!maketop)
            mask = clearBit<T>(mask,5);

          T x0 = T(j) * dx;
          T x1 = x0 + dx;
          T y0 = T(i) * dy;
          T y1 = y0 + dy;
          makeBox(
            TPoint<T>(x0,y0,-dz * 0.5),
            TPoint<T>(x1,y1,+dz * 0.5), 
            "",mask,1,1,1,numpointsU,numpointsV,numpointsW); //!!! m1 = 1
        }
      }
    }

    excludeInnerFaces(tolerance,3,3);

    return true;
  }

  /** Make text with 8x14 font, every pixel size dx x dy x dz is a box,
    bottom, top - in Z direction. LOWER!!! left corner is (0,0,0). */
  bool makeText8x14(const std::string &text, T dx, T dy, T dz,
    bool makebottom, bool maketop, T tolerance, 
    T parmtolerance = PARM_TOLERANCE,     
    int m1 = SPLINE_DEGREE, int m2 = SPLINE_DEGREE, int m3 = SPLINE_DEGREE,
    int numpointsU = 11, int numpointsV = 11, int numpointsW = 11,
    int numcontourpoints = MANY_POINTS2D)
  {
    for (int i = 0; i < int(text.size()); i++)
    {
      TBrep<T> letter(tolerance);
      letter.makeLetter8x14(text[i],dx,dy,dz,makebottom,maketop,tolerance, 
        parmtolerance,m1,m2,m3,numpointsU,numpointsV,numpointsW,numcontourpoints);

      TTransform<T> t;
      t.Translate(TPoint<T>(T(i) * dx * 8));

      letter.makeTransform(&t);

      addFaces(letter);
    }

    return true;
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