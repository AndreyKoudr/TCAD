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

  tbasesurface.h

  Basic abstract class for surfaces and surface patches

  dimensions : 2 (U,V parameters)

*******************************************************************************/

#pragma once

#include "tbasics.h"
#include "tpoint.h"
#include "tmisc.h"
#include "ttransform.h"
#include "tbasecurve.h"
#include "ttriangles.h"
#include "tsystems.h"
#include "tboundary.h"

namespace tcad {

// U or V or W for e.g. set derivatives on parameter
enum Parameter {
  PARAMETER_U,      // for surfaces and volumes
  PARAMETER_UU = PARAMETER_U,    
  PARAMETER_V,      // for surfaces and volumes
  PARAMETER_VV = PARAMETER_V,    
  PARAMETER_W,      // for volumes
  PARAMETER_WW = PARAMETER_W,    

  PARAMETER_ANY,
                    // second cross derivatives
  PARAMETER_UV,     
  PARAMETER_UW,  
   
  PARAMETER_VU,  
  PARAMETER_VW,   
  
  PARAMETER_WU, 
  PARAMETER_WV     
};

// booleans
enum Boolean {
  UNITE,
  SUBTRACT,
  INTERSECT
};

/** Outer loop, not hole. */
template <class T> bool outerLoop(std::vector<std::vector<TPoint<T>>> &loop, T parmtolerance = PARM_TOLERANCE)
{
  TPoint<T> loopnormal = getLoopNormal(loop,parmtolerance);

  // this is a hole, we need an outer loop as loop[0]
  // sure it is not possible to have both holes and patches at the same time
  bool outer = (loopnormal > TPoint<T>(0.0,0.0,1.0));

  return outer;
}

template <class T> class TBaseSurface {
public:

  //===== Construction =========================================================

  /** Constructor. */
  TBaseSurface() {}

  /** Destructor. */
  virtual ~TBaseSurface() {}

  //===== Abstract =============================================================

  /** Get k-th derivative on parameter U,V [0..1]. 0-derivative is position(U,V) (see below). */
  virtual TPoint<T> derivative(T U, T V, Parameter onparameter, int k) = 0;

  /** Update after any change in control points. */
  virtual void update() = 0;

  //===== Operations ===========================================================

  /** Same as 0-th derivative. */
  virtual TPoint<T> position(T U, T V)
  {
    return this->derivative(U,V,PARAMETER_ANY,0);
  }

  /** Normal (not normalised). */
  virtual TPoint<T> normal(T U, T V)
  {
    return this->derivative(U,V,PARAMETER_U,1) ^ this->derivative(U,V,PARAMETER_V,1);
  }

  /** Get control points. */
  virtual std::vector<TPoint<T>> &controlPoints()
  {
    return cpoints;
  }

  /** Get row of control points. */
  void getRow(int index, std::vector<TPoint<T>> &row)
  {
    tcad::getRow(this->cpoints,K1,K2,index,row);
  }

  /** Set row of control points. */
  void setRow(int index, const std::vector<TPoint<T>> &row)
  {
    tcad::setRow(this->cpoints,K1,K2,index,row);
  }

  /** Get column of control points. */
  void getColumn(int index, std::vector<TPoint<T>> &column)
  {
    tcad::getColumn(this->cpoints,K1,K2,index,column);
  }

  /** Set column of control points. */
  void setColumn(int index, const std::vector<TPoint<T>> &column)
  {
    tcad::setColumn(this->cpoints,K1,K2,index,column);
  }

  /** Generate a regular set of actual points on the surface to
    (1) create another type of surface from these points (new K1 = numpointsU - 1),
        (new K2 = numpointsV - 1) or 
    (2) compare/involve another type of geometry to compare/interact 

    A set of according parameter values is generated as well.
    Set refine... to 0.5 at a corresponding end to refine,
    1.0 has no effect.

    UVpoints, if not null, contain U,V parameter values in X,Y
  */
  virtual void createPoints(std::vector<TPoint<T>> &points, 
    std::vector<TPoint<T>> *UVpoints,
    int *k1 = nullptr, int *k2 = nullptr,
    int numpointsU = MANY_POINTS2D, int numpointsV = MANY_POINTS2D,
    T refinestartU = 1.0, T refineendU = 1.0, 
    T refinestartV = 1.0, T refineendV = 1.0)
  {
    points.clear();

    for (int i = 0; i < numpointsV; i++)
    {
      T V = T(i) / T(numpointsV - 1);

      V = refineParameter(V,refinestartV,refineendV);

      for (int j = 0; j < numpointsU; j++)
      {
        T U = T(j) / T(numpointsU - 1);

        U = refineParameter(U,refinestartU,refineendU);

        TPoint<T> p = this->position(U,V);
        points.push_back(p);

        if (UVpoints)
          UVpoints->push_back(TPoint<T>(U,V));
      }
    }

    if (k1)
      *k1 = numpointsU - 1;
    if (k2)
      *k2 = numpointsV - 1;
  }

  /** Generate a uniform set of actual points on the surface to
    (1) create another type of surface from these points (new K1 = numpointsU - 1),
        (new K2 = numpointsV - 1) or 
    (2) compare/involve another type of geometry to compare/interact 

    The mesh is made as unform as possible, no guarantee that the number of points
    is exactly as specified.

    A set of according parameter values is generated as well. XYZstepU, XYZstepV are 
    real model coordinates.

    UVpoints, if not null, contain U,V parameter values in X,Y
  */
  virtual void createPointsUniform(std::vector<TPoint<T>> &points, 
    std::vector<TPoint<T>> *UVpoints,
    T XYZstepU, T XYZstepV, 
    int *k1 = nullptr, int *k2 = nullptr, int manypoints = MANY_POINTS2D, 
    T Umin = 0.0, T Umax = 1.0, T Vmin = 0.0, T Vmax = 1.0,
    int minpoints = 11, int maxpoints = 301) //!!!!!!!
  {
    assert(XYZstepU > TOLERANCE(T));
    assert(XYZstepV > TOLERANCE(T));
    assert((Umax - Umin) > TOLERANCE(T));
    assert((Vmax - Vmin) > TOLERANCE(T));

    points.clear();

    // maybe surface is small, we need to decrease step
    // approximate # divisions
    T usize = Usize() * (Umax - Umin);
    T vsize = Vsize() * (Vmax - Vmin);

    int Udivs = int(usize / XYZstepU);
    int Vdivs = int(vsize / XYZstepV);

    // refine/derefine step
    if (Udivs < (minpoints - 1))
      XYZstepU = usize / T(minpoints - 1);
    if (Vdivs < (minpoints - 1))
      XYZstepV = vsize / T(minpoints - 1);
    if (Udivs > (maxpoints - 1))
      XYZstepU = usize / T(maxpoints - 1);
    if (Vdivs > (maxpoints - 1))
      XYZstepV = vsize / T(maxpoints - 1);

    // just in case, default step
    T defaultstep = 1.0 / T(manypoints - 1);

    // U,V to make an almost regular mesh
    std::vector<T> Us,Vs;

    // make U divisions along middle line V = (Vmin + Vmax) * 0.5
    // last used U increment
    T Ustep = defaultstep;
    // current U
    T U = Umin;
    Us.push_back(U);

    while (U < Umax)
    {
      T dleft = Umax - U;

      if (dleft < Ustep * 0.1) //!!!
      {
        U = Umax;
        Us.push_back(U);
      } else
      {
        TPoint<T> der = this->derivative(U,(Vmin + Vmax) * 0.5,PARAMETER_U,1);
        T len = !der;
        Ustep = (len > TOLERANCE(T)) ? (XYZstepU / len) : defaultstep;

        U += Ustep;
        LIMIT_MAX(U,Umax);
        Us.push_back(U);
      }
    }

    // make V divisions along middle line U = (Umin + Umax) * 0.5
    // last used V increment
    T Vstep = defaultstep;
    // current V
    T V = Vmin;
    Vs.push_back(V);

    while (V < Vmax)
    {
      T dleft = Vmax - V;

      if (dleft < Vstep * 0.1) //!!!
      {
        V = Vmax;
        Vs.push_back(V);
      } else
      {
        TPoint<T> der = this->derivative((Umin + Umax) * 0.5,V,PARAMETER_V,1);
        T len = !der;
        Vstep = (len > TOLERANCE(T)) ? (XYZstepV / len) : defaultstep;

        V += Vstep;
        LIMIT_MAX(V,Vmax);
        Vs.push_back(V);
      }
    }

    for (int i = 0; i < int(Vs.size()); i++)
    {
      T V = Vs[i];

      for (int j = 0; j < int(Us.size()); j++)
      {
        T U = Us[j];

        TPoint<T> p = this->position(U,V);
        points.push_back(p);

        if (UVpoints)
          UVpoints->push_back(TPoint<T>(U,V));
      }
    }

    if (k1)
      *k1 = int(Us.size()) - 1;
    if (k2)
      *k2 = int(Vs.size()) - 1;
  }

  /** Create triangles on a regular set points.
    Set refine... to 0.5 at a corresponding end to refine,
    1.0 has no effect. */
  virtual bool createTriangles(TTriangles<T> &tris,
    int numpointsU = MANY_POINTS2D, int numpointsV = MANY_POINTS2D,
    T refinestartU = 1.0, T refineendU = 1.0, 
    T refinestartV = 1.0, T refineendV = 1.0,
    T stepU = 0.0, T stepV = 0.0, bool uniform = false,
    T Umin = 0.0, T Umax = 1.0, T Vmin = 0.0, T Vmax = 1.0,
    int minpoints = 11, int maxpoints = 301) 
  {
    tris.clear();

    std::vector<TPoint<T>> points;
    std::vector<TPoint<T>> UVpoints;
    int k1 = 0;
    int k2 = 0;

    if (uniform)
    {
      assert(stepU > TOLERANCE(T));
      assert(stepV > TOLERANCE(T));

      int numpoints = std::max<int>(numpointsU,numpointsV);
      createPointsUniform(points,&UVpoints,stepU,stepV,&k1,&k2,numpoints,Umin,Umax,Vmin,Vmax,minpoints,maxpoints);
    } else
    {
      createPoints(points,&UVpoints,&k1,&k2,numpointsU,numpointsV,refinestartU,refineendU,refinestartV,refineendV);
    }

    for (int i = 0; i < k2; i++)
    {
      std::vector<TPoint<T>> row0,row1;
      tcad::getRow(points,k1,k2,i,row0);
      tcad::getRow(points,k1,k2,i + 1,row1);

      std::vector<TPoint<T>> UVrow0,UVrow1;
      tcad::getRow(UVpoints,k1,k2,i,UVrow0);
      tcad::getRow(UVpoints,k1,k2,i + 1,UVrow1);

      for (int j = 0; j < int(row0.size() - 1); j++)
      {
        TPoint<T> p0 = row0[j];
        TPoint<T> p1 = row0[j + 1];
        TPoint<T> p2 = row1[j + 1];
        TPoint<T> p3 = row1[j];

        TPoint<T> UVp0 = UVrow0[j];
        TPoint<T> UVp1 = UVrow0[j + 1];
        TPoint<T> UVp2 = UVrow1[j + 1];
        TPoint<T> UVp3 = UVrow1[j];

        T d0 = !(p2 - p0);
        T d1 = !(p3 - p1);

        if (d0 < d1)
        {
          tris.addTri(p2,p0,p1,0.0);
          tris.addTri(p0,p2,p3,0.0);

          tris.UVcorners.push_back(std::array<TPoint<T>,3> {UVp2,UVp0,UVp1});
          tris.UVcorners.push_back(std::array<TPoint<T>,3> {UVp0,UVp2,UVp3});
        } else
        {
          tris.addTri(p1,p3,p0,0.0);
          tris.addTri(p3,p1,p2,0.0);

          tris.UVcorners.push_back(std::array<TPoint<T>,3> {UVp1,UVp3,UVp0});
          tris.UVcorners.push_back(std::array<TPoint<T>,3> {UVp3,UVp1,UVp2});
        }
      }
    }

    std::pair<TPoint<T>,TPoint<T>> mm = tris.minmax();
    TPoint<T> d = mm.second - mm.first;
    T tolerance = !d * PARM_TOLERANCE;

    // remove duplicate nodes and renumber corners

#if 1
    bool ok = true;
#else
    // this can be very very slow for 200x200 quads and 240000 nodes
    bool ok = tris.buildConnectivityArray(tolerance);
#endif
    return ok;
  }

  /** Reverse U. Normal is changed to opposite. */
  virtual void reverseU()
  {
    // reverse rows of control poins
    for (int i = 0; i <= K2; i++)
    {
      std::vector<TPoint<T>> row;
      this->getRow(i,row);
      std::reverse(row.begin(),row.end());
      this->setRow(i,row);
    }
  }

  /** Reverse V. Normal is changed to opposite. */
  virtual void reverseV()
  {
    // reverse columns of control poins
    for (int i = 0; i <= K1; i++)
    {
      std::vector<TPoint<T>> col;
      getColumn(i,col);
      std::reverse(col.begin(),col.end());
      setColumn(i,col);
    }
  }

  /** Approximate size along U. */
  T Usize()
  {
    T len = 0.0;
    for (int i = 0; i <= K2; i++)
    {
      std::vector<TPoint<T>> row;
      this->getRow(i,row);
      len += calculateLength(row);
    }

    len /= T(K2 + 1);

    return len;
  }

  /** Approximate size along V. */
  T Vsize()
  {
    T len = 0.0;
    for (int i = 0; i <= K1; i++)
    {
      std::vector<TPoint<T>> col;
      this->getColumn(i,col);
      len += calculateLength(col);
    }

    len /= T(K1 + 1);

    return len;
  }

  /** Approximate size. */
  T maxSize()
  {
    return std::max<T>(Usize(),Vsize());
  }

  /** Approximate size. */
  T minSize()
  {
    return std::min<T>(Usize(),Vsize());
  }

  /** Get XYZ step along U to make approximately numpoints. */
  T Ustep(int numpoints = MANY_POINTS2D, T Umin = 0.0, T Umax = 1.0)
  {
    assert(Umax > Umin);
    T step = Usize() * (Umax - Umin) / T(numpoints - 1);
    return step;
  }

  /** Get XYZ step along V to make approximately numpoints. */
  T Vstep(int numpoints = MANY_POINTS2D, T Vmin = 0.0, T Vmax = 1.0)
  {
    assert(Vmax > Vmin);
    T step = Vsize() * (Vmax - Vmin) / T(numpoints - 1);
    return step;
  }

  /** Get curvature. */
  virtual std::pair<T,T> curvature(T U, T V)
  {
    TPoint<T> Fu = this->derivative(U,V,PARAMETER_U,1);
    TPoint<T> Fv = this->derivative(U,V,PARAMETER_V,1);
    TPoint<T> Fuu = this->derivative(U,V,PARAMETER_U,2);
    TPoint<T> Fvv = this->derivative(U,V,PARAMETER_V,2);

	  T Ulen = !Fu;
	  T curvatureU = (Ulen < TOLERANCE(T)) ? 0.0 : (!(Fu ^ Fuu) / (Ulen * Ulen * Ulen));

	  T Vlen = !Fv;
	  T curvatureV = (Vlen < TOLERANCE(T)) ? 0.0 : (!(Fv ^ Fvv) / (Vlen * Vlen * Vlen));

    return std::pair<T,T>(curvatureU,curvatureV);
  }

  /** Get U curvature. */
  virtual T curvatureU(T U, T V)
  {
    TPoint<T> Fu = this->derivative(U,V,PARAMETER_U,1);
    TPoint<T> Fuu = this->derivative(U,V,PARAMETER_U,2);

	  T Ulen = !Fu;
	  T curvatureU = (Ulen < TOLERANCE(T)) ? 0.0 : (!(Fu ^ Fuu) / (Ulen * Ulen * Ulen));

    return curvatureU;
  }

  /** Get V curvature. */
  virtual T curvatureV(T U, T V)
  {
    TPoint<T> Fv = this->derivative(U,V,PARAMETER_V,1);
    TPoint<T> Fvv = this->derivative(U,V,PARAMETER_V,2);

	  T Vlen = !Fv;
	  T curvatureV = (Vlen < TOLERANCE(T)) ? 0.0 : (!(Fv ^ Fvv) / (Vlen * Vlen * Vlen));

    return curvatureV;
  }

  /** Get XYZ step along U to achieve distance t to straight-line edge. Full U size in case 
    of failure (curvature radius infinite). */
  virtual T curvatureStepU(T t, T U, T V)
  {
    T curvature = curvatureU(U,V);
    // straight line
    if (curvature < TOLERANCE(T))
    {
      return this->Usize();
    } else
    {
      T r = 1.0 / curvature;
      T step = 2.0 * sqrt(2.0 * t * r);
      return step;
    }
  }

  /** Get XYZ step along V to achieve distance t to straight-line edge. Full V size in case 
    of failure (curvature radius infinite). */
  virtual T curvatureStepV(T t, T U, T V)
  {
    T curvature = curvatureV(U,V);
    // straight line
    if (curvature < TOLERANCE(T))
    {
      return this->Vsize();
    } else
    {
      T r = 1.0 / curvature;
      T step = 2.0 * sqrt(2.0 * t * r);
      return step;
    }
  }

  /** Make transform. */
  virtual void makeTransform(TTransform<T> *transform)
  {
    for (auto &p : cpoints)
    {
      p = transform->applyTransform(p);
    }

    // call virtual update()
    this->update();
  }

  /** Are equal? */
  virtual bool equal(TBaseSurface &other, T tolerance, int numpointsU = MANY_POINTS2D, 
    int numpointsV = MANY_POINTS2D)
  {
    std::vector<TPoint<T>> points,otherpoints;
    createPoints(points,nullptr,nullptr,nullptr,numpointsU,numpointsV);
    other.createPoints(otherpoints,nullptr,nullptr,nullptr,numpointsU,numpointsV);

    T diff = difference(points,otherpoints);
    return (diff >= 0.0 && diff < tolerance);
  }

  /** Are equal? Report if reversed in U or V. */
  bool equal(TBaseSurface &other, T tolerance, bool &reversedU, bool &reversedV,
    int numpointsU = MANY_POINTS2D, int numpointsV = MANY_POINTS2D)
  {
    reversedU = reversedV = false;

    if (this->equal(other,tolerance,numpointsU,numpointsV))
    {
      return true;
    } 

    other.reverseU();
    bool ok = this->equal(other,tolerance,numpointsU,numpointsV);
    other.reverseU();

    if (ok)
    {
      reversedU = true;
      reversedV = false;
      return true;
    }

    other.reverseV();
    ok = this->equal(other,tolerance,numpointsU,numpointsV);
    other.reverseV();

    if (ok)
    {
      reversedU = false;
      reversedV = true;
      return true;
    }

    other.reverseU();
    other.reverseV();
    ok = this->equal(other,tolerance,numpointsU,numpointsV);
    other.reverseV();
    other.reverseU();

    if (ok)
    {
      reversedU = true;
      reversedV = true;
      return true;
    }

    return false;
  }

  /** Find values of parameters U,V for a point on (or close to) the surface. 
    Create finer mesh by createPoints() for more accurate results. */
  static TPoint<T> findUVforPoint(std::vector<TPoint<T>> &points,
    std::vector<TPoint<T>> &UVpoints, int k1, int k2, TPoint<T> p) 
  {
    // temp
    TPoint<T> proj;
    int seg = 0;
    T u = 0.0;

    // look for closest point on rows
    T U = -1.0;
    T mindist = std::numeric_limits<T>::max();
    for (int j = 0; j < k2 + 1; j++)
    {
      std::vector<TPoint<T>> row;
      tcad::getRow(points,k1,k2,j,row);
      if (projectPointOnPoints(row,p,proj,&seg,&u))
      {
        T dist = !(p - proj);
        if (dist < mindist)
        {
          TPoint<T> p0 = row[seg];
          TPoint<T> p1 = row[seg + 1];
          int index0 = tcad::getIndex<T>(k1,k2,seg,j);
          int index1 = tcad::getIndex<T>(k1,k2,seg + 1,j);
          U = UVpoints[index0].X + (UVpoints[index1].X - UVpoints[index0].X) * u;
          mindist = dist;
        }
      }
    }

    // look for closest point on columns
    T V = -1.0;
    mindist = std::numeric_limits<T>::max();
    for (int i = 0; i < k1 + 1; i++)
    {
      std::vector<TPoint<T>> col;
      tcad::getColumn(points,k1,k2,i,col);
      if (projectPointOnPoints(col,p,proj,&seg,&u))
      {
        T dist = !(p - proj);
        if (dist < mindist)
        {
          TPoint<T> p0 = col[seg];
          TPoint<T> p1 = col[seg + 1];
          int index0 = tcad::getIndex<T>(k1,k2,i,seg);
          int index1 = tcad::getIndex<T>(k1,k2,i,seg + 1);
          V = UVpoints[index0].Y + (UVpoints[index1].Y - UVpoints[index0].Y) * u;
          mindist = dist;
        }
      }
    }

    return TPoint<T>(U,V);
  }

  /** Find approximate (linearised) UV for points. */
  void findUVForPoints(std::vector<TPoint<T>> &points, std::vector<TPoint<T>> &UVs,
    int numpointsU = MANY_POINTS2D, int numpointsV = MANY_POINTS2D)
  {
    std::vector<TPoint<T>> spoints,sUVpoints;
    int k1 = 0;
    int k2 = 0;
    createPoints(spoints,&sUVpoints,&k1,&k2,numpointsU,numpointsV);

    for (int i = 0; i < int(points.size()); i++)
    {
      TPoint<T> UV = findUVforPoint(spoints,sUVpoints,k1,k2,points[i]);
      UVs.push_back(UV);
    }
  }

  /** Calculate min/max. */
  virtual bool calculateMinMax(TPoint<T> *min, TPoint<T> *max, TPoint<T> *imin = nullptr, 
    TPoint<T> *imax = nullptr, int numpointsU = MANY_POINTS2D, int numpointsV = MANY_POINTS2D)
  {
    std::vector<TPoint<T>> points;
    createPoints(points,nullptr,nullptr,nullptr,numpointsU,numpointsV);

    return tcad::calculateMinMax(points,min,max,imin,imax);
  }

  /** Cut by plane. */
  bool intersectByPlane(TPlane<T> &plane, std::vector<std::vector<TPoint<T>>> &lines, 
    std::vector<std::vector<TPoint<T>>> &boundary,
    T tolerance, T parmtolerance = TOLERANCE(T), 
    int numpointsU = MANY_POINTS2D, int numpointsV = MANY_POINTS2D,
    T refinestartU = 1.0, T refineendU = 1.0, 
    T refinestartV = 1.0, T refineendV = 1.0)
  {
    TTriangles<T> tris;
    if (createTriangles(tris,numpointsU,numpointsV,refinestartU,refineendU,refinestartV,refineendV))
    {
      return tris.intersectByPlane(plane,lines,tolerance,parmtolerance,&boundary);
    } else
    {
      return false;
    }
  }

  /** Cut out a part of surface from U0 to U1 and from V0 to V1 into list of points. */
  template <class T> void cutPiece(T Ufrom, T Uto, T Vfrom, T Vto, std::vector<TPoint<T>> &points,
    int numpointsU = MANY_POINTS2D, int numpointsV = MANY_POINTS2D)
  {
    LIMIT(Ufrom,0.0,1.0);
    LIMIT(Uto,0.0,1.0);
    LIMIT(Vfrom,0.0,1.0);
    LIMIT(Vto,0.0,1.0);

    points.clear();
    T DU = (Uto - Ufrom) / T(numpointsU - 1);
    T DV = (Vto - Vfrom) / T(numpointsV - 1);

    for (int i = 0; i < numpointsV; i++)
    {
      T V = Vfrom + T(i) * DV;
      for (int j = 0; j < numpointsU; j++)
      {
        T U = Ufrom + T(j) * DU;

        TPoint<T> p = position(U,V);
        points.push_back(p);
      }
    }
  }

  /** Get XYZ difference along two boundary (cut) lines. */
  template <class T> T boundaryDifference(TBaseSurface<T> &other,
    std::vector<std::vector<TPoint<T>>> &boundary0,
    std::vector<std::vector<TPoint<T>>> &boundary1)
  {
    assert(boundary0.size() == boundary1.size());

    // calculate difference in boundary lines
    T maxdiff = 0.0;

    for (int k = 0; k < boundary0.size(); k++)
    {
      //std::vector<TPoint<T>> intr0,intr1;
      for (int l = 0; l < int(boundary0[k].size()); l++)
      {
        TPoint<T> p0 = position(boundary0[k][l].X,boundary0[k][l].Y);
        TPoint<T> p1 = other.position(boundary1[k][l].X,boundary1[k][l].Y);
        //intr0.push_back(p0);
        //intr1.push_back(p1);

        T d = !(p1 - p0);
        maxdiff = std::max<T>(maxdiff,d);
      }
    }

    return maxdiff;
  }

  /** Prepare triangles from UV min/max. This is a cyclic part of intersect(). */
  template <class T> bool makeIntersectionTriangles(TBaseSurface<T> &other, 
    TTriangles<T> &tris, TTriangles<T> &othertris, 
    TPoint<T> &UVmin, TPoint<T> &UVmax, TPoint<T> &oUVmin, TPoint<T> &oUVmax, 
    int minpoints, int maxpoints, bool makesamestep, T tolerance,
    T parmtolerance = PARM_TOLERANCE)
  { 
    // free memory
    tris.clear();
    othertris.clear();

    TPoint<T> UVcentre = (UVmin + UVmax) * 0.5;
    TPoint<T> oUVcentre = (oUVmin + oUVmax) * 0.5;
    T stepU = curvatureStepU(tolerance,UVcentre.X,UVcentre.Y); 
    T stepV = curvatureStepV(tolerance,UVcentre.X,UVcentre.Y);
    T ostepU = other.curvatureStepU(tolerance,oUVcentre.X,oUVcentre.Y);
    T ostepV = other.curvatureStepV(tolerance,oUVcentre.X,oUVcentre.Y);

    if (makesamestep)
    {
      T stepmax = std::min<T>(stepU,stepV);
      T ostepmax = std::min<T>(ostepU,ostepV);

      //!!!!!!!stepU = stepV = ostepU = ostepV = (stepmax + ostepmax) * 0.5;
      stepU = stepV = ostepU = ostepV = std::min<T>(stepmax,ostepmax);
    }

//outputDebugString(std::string("stepU ") + to_string(stepU) + std::string(" stepV ") + to_string(stepV) +
//  std::string(" ostepU ") + to_string(ostepU) + std::string(" ostepV ") + to_string(ostepV)); 

    if (createTriangles(tris,MANY_POINTS2D,MANY_POINTS2D,1.0,1.0,1.0,1.0,stepU,stepV,true,
      UVmin.X,UVmax.X,UVmin.Y,UVmax.Y,minpoints,maxpoints) && 
      other.createTriangles(othertris,MANY_POINTS2D,MANY_POINTS2D,1.0,1.0,1.0,1.0,ostepU,ostepV,true,
      oUVmin.X,oUVmax.X,oUVmin.Y,oUVmax.Y,minpoints,maxpoints)) 
    {
      return true;
    } else
    {
      return false;
    }
  }

  /** Prepare triangles by refining triangles involved in active cells. */
  template <class T> bool makeIntersectionTriangles(TBaseSurface<T> &other, 
    TTriangles<T> &tris, TTriangles<T> &othertris, 
    OBackground<T> &cells, std::vector<std::vector<LINT>> &celltris, std::vector<std::vector<LINT>> &ocelltris,
    std::vector<LINT> &activecells)
  { 
    TTriangles<T> newtris;
    TTriangles<T> newothertris;

    for (int i = 0; i < int(activecells.size()); i++)
    {
      int index = int(activecells[i]);

      for (int j = 0; j < int(celltris[index].size()); j++)
      {
        int face = int(celltris[index][j]);
        std::array<TPoint<T>,3> corners = tris.threeCorners(face);

        newtris.addTri(corners[0],corners[1],corners[2],0.0);
        newtris.UVcorners.push_back(tris.UVcorners[face]);
      }
      for (int j = 0; j < int(ocelltris[index].size()); j++)
      {
        int face = int(ocelltris[index][j]);
        std::array<TPoint<T>,3> corners = othertris.threeCorners(face);

        newothertris.addTri(corners[0],corners[1],corners[2],0.0);
        newothertris.UVcorners.push_back(othertris.UVcorners[face]);
      }
    }

    tris = newtris;
    othertris = newothertris;

    return true;
  }

  /** Refine triangles. */
  template <class T> bool makeIntersectionTriangles(TBaseSurface<T> &other, 
    TTriangles<T> &tris, TTriangles<T> &othertris, T coef)
  { 
    int divisions = ROUND(coef);
    LIMIT(divisions,2,MAX_TRISUBDIVS); 

    if (divisions > 1)
    {
      // barycentric coordinates of small triangles
      std::vector<std::array<TPoint<T>,3>> barcoords;
      subdivideTriangle(divisions,barcoords);

      TTriangles<T> newtris;
      TTriangles<T> newothertris;

      for (int face = 0; face < tris.numFaces(); face++)
      {
        std::array<TPoint<T>,3> UVcorners = tris.UVcorners[face];

        for (int i = 0; i < int(barcoords.size()); i++)
        {
          std::array<TPoint<T>,3> uvcorners;
          std::array<TPoint<T>,3> corners;
          for (int k = 0; k < 3; k++)
          {
            TPoint<T> uv = barycentricValue<T>(UVcorners[0],UVcorners[1],UVcorners[2],barcoords[i][k]);
            uvcorners[k] = uv;
            corners[k] = this->position(uv.X,uv.Y);
          }

          newtris.addTri(corners[0],corners[1],corners[2],0.0);
          newtris.UVcorners.push_back(uvcorners);
        }
      }

      for (int face = 0; face < othertris.numFaces(); face++)
      {
        std::array<TPoint<T>,3> UVcorners = othertris.UVcorners[face];

        for (int i = 0; i < int(barcoords.size()); i++)
        {
          std::array<TPoint<T>,3> uvcorners;
          std::array<TPoint<T>,3> corners;
          for (int k = 0; k < 3; k++)
          {
            TPoint<T> uv = barycentricValue<T>(UVcorners[0],UVcorners[1],UVcorners[2],barcoords[i][k]);
            uvcorners[k] = uv;
            corners[k] = other.position(uv.X,uv.Y);
          }

          newothertris.addTri(corners[0],corners[1],corners[2],0.0);
          newothertris.UVcorners.push_back(uvcorners);
        }
      }

      tris = newtris;
      othertris = newothertris;
    }

    return true;
  }

  /** Prepare triangles and cells to cut out intersection part. This is a cyclic 
    part of intersect(). */
  template <class T> bool makeIntersectionCells(TBaseSurface<T> &other, 
    TTriangles<T> &tris, TTriangles<T> &othertris, 
    OBackground<T> &cells, std::vector<std::vector<LINT>> &celltris, std::vector<std::vector<LINT>> &ocelltris,
    std::vector<LINT> &activecells,
    TPoint<T> &UVmin, TPoint<T> &UVmax, TPoint<T> &oUVmin, TPoint<T> &oUVmax, 
    T parmtolerance = PARM_TOLERANCE, T maxedgecoef = MAXEDGE_COEF, T increaseUVcoef = MINMAX_COEF)
  { 
    // cutout intersecting piece and make refined intersection on them
    T minedge,maxedge;
    tris.getEdgeMinMax(minedge,maxedge);
    T ominedge,omaxedge;
    othertris.getEdgeMinMax(ominedge,omaxedge);

    // max edge length
    maxedge = std::max(maxedge,omaxedge) * maxedgecoef; 

    // failure if no active cells
    return tris.prepareCells(othertris,maxedge,cells,celltris,ocelltris,activecells,
      UVmin,UVmax,oUVmin,oUVmax,parmtolerance,increaseUVcoef,print);
  }

  /** Find intersection curve(s) with another surface. Returns number of intersection curves. 
    boundary0/1 contain U,V in X,Y for every intersection curve for both surfaces. */
  template <class T> int intersect(TBaseSurface<T> &other, bool bodyleft, std::vector<std::vector<TPoint<T>>> &intersections, 
    std::vector<std::vector<TPoint<T>>> &boundary0, std::vector<std::vector<TPoint<T>>> &boundary1,
    T tolerance,
    T parmtolerance = PARM_TOLERANCE, 
    int numpointsU0 = MANY_POINTS2D, int numpointsV0 = MANY_POINTS2D,
    T refinestartU0 = 1.0, T refineendU0 = 1.0, 
    T refinestartV0 = 1.0, T refineendV0 = 1.0,
    int numpointsU1 = MANY_POINTS2D, int numpointsV1 = MANY_POINTS2D,
    T refinestartU1 = 1.0, T refineendU1 = 1.0, 
    T refinestartV1 = 1.0, T refineendV1 = 1.0,
    bool debug = false, int preintrcycles = PREINTR_CYCLES, int intrcycles = INTR_CYCLES, 
    int maxboundarypoints = MAX_BOUNDPOINTS)
  {
    // generated triangles
    TTriangles<T> tris,othertris;

    // cells for space partitioning
    OBackground<T> cells; 
    std::vector<std::vector<LINT>> celltris,ocelltris;
    std::vector<LINT> activecells;

    // parameter limits in active cells (where tris intersect)
    TPoint<T> UVmin(0.0,0.0);
    TPoint<T> UVmax(1.0,1.0);
    TPoint<T> oUVmin(0.0,0.0);
    TPoint<T> oUVmax(1.0,1.0);

    T tol = tolerance * pow(10.0,preintrcycles - 1); 

    // shrink UVmin/max around intersection line
    for (int i = 0; i < preintrcycles; i++)
    {
//outputDebugString("coarsening " + to_string(i)); 

      if (!makeIntersectionTriangles(other,tris,othertris,UVmin,UVmax,oUVmin,oUVmax,11,301,(i >= preintrcycles - 1), //!!!!!!!
        tol,parmtolerance))
        return 0;

      if (!makeIntersectionCells(other,tris,othertris,
        cells,celltris,ocelltris,activecells,UVmin,UVmax,oUVmin,oUVmax,parmtolerance))
        return 0;

      // decrease tolerance, increase refinement
      tol *= 0.1; 
    }

    // now try intersection and check accuracy
    for (int i = 0; i < intrcycles; i++)
    {
//outputDebugString("intersecting " + to_string(i));

      // make tris which only used in active cells
      if (!makeIntersectionTriangles(other,tris,othertris,cells,celltris,ocelltris,activecells))
        return 0;

      // make intersection
      if (!tris.intersect(othertris,bodyleft,intersections,parmtolerance,&boundary0,&boundary1))
        return 0;

      for (int k = 0; k < int(boundary0.size()); k++)
      {
        decimatePoints(boundary0[k],maxboundarypoints);
        decimatePoints(boundary1[k],maxboundarypoints);
      }

      // calculate difference in boundary lines
      T maxdiff = boundaryDifference(other,boundary0,boundary1);

//outputDebugString( 
//  std::string("i = ") + to_string(i) +
//  std::string(" tris ") + to_string(tris.numFaces()) + " " + to_string(othertris.numFaces()) + " " + 
//  std::string(" cells ") + to_string(activecells.size()) +
//  std::string(" maxdiff ") + to_string(maxdiff,12) + 
//  std::string(" bigtolerance ") + to_string(tolerance * BIGTOLERANCE_COEF,12));

      // success
      bool ok = maxdiff < tolerance * BIGTOLERANCE_COEF;
      if (ok || i == intrcycles - 1 ||  //!!!!!!!
        tris.numFaces() > MAX_TRIS || othertris.numFaces() > MAX_TRIS)
      {
//if (ok) 
//{
//  outputDebugString(std::string("ok i = ") + to_string(i)); 
//} else
//{
//  outputDebugString(std::string("exit i = ") + to_string(i));
//}
        return int(intersections.size());
      } else
      {
        // increase refinement
        T coef = maxdiff * 2.0 / (tolerance * BIGTOLERANCE_COEF);

//outputDebugString(std::string("coef ") + to_string(coef,12) + " i = " + to_string(i)); 

        // refine triangles
        if (!makeIntersectionTriangles(other,tris,othertris,coef))
          return 0;

        // make new cells
        if (!makeIntersectionCells(other,tris,othertris,
          cells,celltris,ocelltris,activecells,UVmin,UVmax,oUVmin,oUVmax,parmtolerance))
          return 0;
      }
    }

    // failure
    return 0;
  }

  /** Find intersection curve(s) by a plane. Returns number of intersection curves. 
    boundary0/1 contain U,V in X,Y for every intersection curve for both surfaces. */
  template <class T> int intersectByPlane(TPlane<T> &plane, std::vector<std::vector<TPoint<T>>> &intersections, 
    std::vector<std::vector<TPoint<T>>> &boundary,
    T tolerance, T parmtolerance = TOLERANCE(T), 
    int numpointsU = MANY_POINTS2D, int numpointsV = MANY_POINTS2D,
    T refinestartU = 1.0, T refineendU = 1.0, 
    T refinestartV = 1.0, T refineendV = 1.0)
  {
    TTriangles<T> tris,othertris;

    if (createTriangles(tris,numpointsU,numpointsV,refinestartU,refineendU,refinestartV,refineendV))
    {
      if (tris.intersectByPlane(plane,intersections,tolerance,parmtolerance,&boundary))
      {
        return int(intersections.size());
      } else
      {
        return 0;
      }
    } else
    {
      return 0;
    }
  }

  /** Convert U/V boundary into XYZ points. */
  template <class T> void boundaryIntoPoints(std::vector<std::vector<TPoint<T>>> &UVboundary,
    std::vector<std::vector<TPoint<T>>> &points)
  {
    points.clear();

    for (int i = 0; i < int(UVboundary.size()); i++)
    {
      points.push_back(std::vector<TPoint<T>>());
      for (int j = 0; j < int(UVboundary[i].size()); j++)
      {
        TPoint<T> p = this->position(UVboundary[i][j].X,UVboundary[i][j].Y);
        points.back().push_back(p);
      }
    }
  }

  /** Convert U/V into XYZ points. */
  template <class T> void UVIntoPoints(std::vector<TPoint<T>> &UV,
    std::vector<TPoint<T>> &points)
  {
    points.clear();

    for (int i = 0; i < int(UV.size()); i++)
    {
      TPoint<T> p = this->position(UV[i].X,UV[i].Y);
      points.push_back(p);
    }
  }

  /** Extend segment to boundary by moving point p1 in the direction p0->p1. 
    intr is intersection point, distance to boundary is returned, -1.0 in 
    case of failure. */
  template <class T> T extendToBoundary(TPoint<T> p0, TPoint<T> p1, TPoint<T> &intr,
    T parmtolerance = PARM_TOLERANCE)
  {
    T mindist = std::numeric_limits<T>::max();
    bool found = false;
    for (int i = 0; i < 4; i++)
    {
      int i1 = (i < 3) ? (i + 1) : 0;
      TPoint<T> c0 = cornerUV<T>[i];
      TPoint<T> c1 = cornerUV<T>[i1];

      T t1,t2,Xi,Yi;
      if (intersectSegmentsXY<T>(p0.X,p0.Y,p1.X,p1.Y,c0.X,c0.Y,c1.X,c1.Y,
        &t1,&t2,&Xi,&Yi) && t1 > 1.0 - parmtolerance) // forward
      {
        T dx = Xi - p1.X;
        T dy = Yi - p1.Y;
        T d = sqrt(dx * dx + dy * dy);
        if (d < mindist)
        {
          mindist = d;
          intr.X = Xi;
          intr.Y = Yi;
        }
        found = true;
      }
    }

    return found ? mindist : -1.0;
  }

  /** Find a cut piece (except busy) closest by its start to boundary. */
  template <class T> int findClosestStart(std::vector<std::vector<TPoint<T>>> &cut,
    T &mindist, TPoint<T> &intr, std::vector<bool> *busy = nullptr)
  {
    // find closest piece to the current loop end among not busy
    int closest = -1;

    mindist = std::numeric_limits<T>::max();
    for (int i = 0; i < int(cut.size()); i++)
    {
      if (busy && (*busy)[i])
        continue;
      
      T dist = extendToBoundary<T>(cut[i][1],cut[i][0],intr);
      if (dist >= 0.0 && dist < mindist)
      {
        mindist = dist;
        closest = i;
      }
    }

    if (busy && closest >= 0)
      (*busy)[closest] = true;

    return closest;
  }

  /** Find a cut piece (except busy) closest by its end to boundary. */
  template <class T> int findClosestEnd(std::vector<std::vector<TPoint<T>>> &cut,
    T &mindist, TPoint<T> &intr, std::vector<bool> *busy = nullptr)
  {
    // find closest piece to the current loop end among not busy
    int closest = -1;

    mindist = std::numeric_limits<T>::max();
    for (int i = 0; i < int(cut.size()); i++)
    {
      if (busy && (*busy)[i])
        continue;
      
      T dist = extendToBoundary<T>(cut[i][cut[i].size() - 2],cut[i][cut[i].size() - 1],intr);
      if (dist >= 0.0 && dist < mindist)
      {
        mindist = dist;
        closest = i;
      }
    }

    if (busy && closest >= 0)
      (*busy)[closest] = true;

    return closest;
  }

  /** Order cuts from starting piece. Pieces all have a correct orientation. */
  template <class T> bool orderCuts(int starting, std::vector<std::vector<TPoint<T>>> &cut, 
    std::vector<std::vector<TPoint<T>>> &ordered, T tolerance, T parmtolerance = PARM_TOLERANCE)
  {
    std::vector<bool> busy(cut.size(),false);

    ordered.push_back(cut[starting]);
    busy[starting] = true;

    // all pieces must have a correct direction set during intersection
    while (!allBusy(busy))
    {
      // success : delete this ordered piece from cuts
      if (boundaryPoint(ordered.back().back(),parmtolerance))
      {
        // cleanup cut
        for (int i = int(busy.size()) - 1; i >= 0; i--)
        {
          if (busy[i])
            cut.erase(cut.begin() + i);
        }

        return true;
      }

      // find closest piece to the current ordered end among not busy
      int closest = findClosestFront(cut,busy,ordered.back().back(),tolerance);

      if (closest < 0)
      {
        return false;
      }

      // attach next piece
      ordered.push_back(cut[closest]);

      // remove gap between the two if any
      TPoint<T> p = (ordered[ordered.size() - 2].back() + ordered[ordered.size() - 1].front()) * 0.5;
      ordered[ordered.size() - 2].back() = ordered[ordered.size() - 1].front() = p;
      busy[closest] = true;

      // success : delete this ordered piece from cuts
      if (boundaryPoint(ordered.back().back(),parmtolerance))
      {
        // cleanup cut
        for (int i = int(busy.size()) - 1; i >= 0; i--)
        {
          if (busy[i])
            cut.erase(cut.begin() + i);
        }

        return true;
      }
    }

    if (!ordered.empty() && boundaryPoint(ordered.back().back(),parmtolerance))
    {
      // cleanup cut
      for (int i = int(busy.size()) - 1; i >= 0; i--)
      {
        if (busy[i])
          cut.erase(cut.begin() + i);
      }

      return true;
    } else
    {
      return false;
    }
  }

  /** Cut contains U,V in X,Y. We need to start from a boundaryPoint() and end with a 
    boundaryPoint(), keeping the same direction with no reversals. Normally there are
    no more than 2 pieces in the cut. tolerance is used to compare piece ends when connecting,
    it may be much bigger than parmtolerance used for boundary identification. */
  template <class T> int orderCutPieces(std::vector<std::vector<TPoint<T>>> &cut, 
    std::vector<std::vector<std::vector<TPoint<T>>>> &allordered, T tolerance, T parmtolerance = PARM_TOLERANCE)
  {
    allordered.clear();

    // trivial
    if (cut.size() == 0)
      return 0;

    int count = 0;
    bool found = false;
    do {

      found = false;
      // find starting piece
      int starting = -1;
      for (int i = 0; i < int(cut.size()); i++)
      {
        if (boundaryPoint(cut[i].front(),parmtolerance))
        {
          starting = i;
          break;
        }
      }

      if (starting < 0)
        break;

      // order pieces from starting
      std::vector<std::vector<TPoint<T>>> ordered;
      // cut is modified
      if (orderCuts(starting,cut,ordered,tolerance,parmtolerance))
      {
        if (boundaryPoint<T>(ordered.front().front(),parmtolerance))
          correctBoundaryPoint<T>(ordered.front().front());
        if (boundaryPoint<T>(ordered.back().back(),parmtolerance))
          correctBoundaryPoint<T>(ordered.back().back());

        allordered.push_back(ordered);
        count++;
        found = true;
      }
    } while (found);


    return count;
  }

  /** Cut contains U,V in X,Y. Extract a closed loop if exists. cut is modified, 
    may not be empty at exit. All pieces are supposed to have a correct direction set
    during intersections. */
  template <class T> bool extractLoop(int starting, std::vector<std::vector<TPoint<T>>> &cut, 
    std::vector<std::vector<TPoint<T>>> &loop, T tolerance, T parmtolerance = PARM_TOLERANCE)
  {
    std::vector<bool> busy(cut.size(),false);

    loop.push_back(cut[starting]);
    busy[starting] = true;

    if (boundaryPoint(loop.front().front(),parmtolerance) ||
      boundaryPoint(loop.back().back(),parmtolerance))
    {
      return false;
    }

    // contour closed?
    T dist = !(loop.front().front() - loop.back().back());
    if (dist < tolerance)
    {
      // cleanup cut
      for (int i = int(busy.size()) - 1; i >= 0; i--)
      {
        if (busy[i])
          cut.erase(cut.begin() + i);
      }

      return true;
    }

    // all pieces must have a correct direction set during intersection
    while (!allBusy(busy))
    {
      // find closest piece to the current ordered end among not busy
      int closest = findClosestFront(cut,busy,loop.back().back(),tolerance);

      if (closest < 0)
      {
        return false;
      }

      // attach next piece
      loop.push_back(cut[closest]);
      busy[closest] = true;

      if (boundaryPoint(loop.front().front(),parmtolerance) ||
        boundaryPoint(loop.back().back(),parmtolerance))
      {
        return false;
      }

      // contour closed?
      T dist = !(loop.front().front() - loop.back().back());
      if (dist < tolerance)
      {
        // cleanup cut
        for (int i = int(busy.size()) - 1; i >= 0; i--)
        {
          if (busy[i])
            cut.erase(cut.begin() + i);
        }

        return true;
      }
    }

    return false;
  }

  /** Cut contains U,V in X,Y. Extract all closed loops. cut is modified, 
    may not be empty at exit. All pieces are supposed to have a correct direction set
    during intersections. */
  template <class T> int extractLoops(std::vector<std::vector<TPoint<T>>> &cut, 
    std::vector<std::vector<std::vector<TPoint<T>>>> &loops, T tolerance, 
    T parmtolerance = PARM_TOLERANCE, bool reverse = false)
  {
    int count = 0;

    bool found = false;
    do {
      found = false;
      for (int i = 0; i < int(cut.size()); i++)
      {
        std::vector<std::vector<TPoint<T>>> loop;
        if (extractLoop(i,cut,loop,tolerance,parmtolerance))
        {
          if (reverse)
          {
            std::reverse(loop.begin(),loop.end());
            for (auto &l : loop)
            {
              std::reverse(l.begin(),l.end());
            }
          }
          loops.push_back(loop);
          found = true;
          count++;
          break;
        }
      }
    } while (found);

    return count;
  }

  /** Close UV boundary. cutUV is a cut across surface in UV coordinates. cutFromOuter is true
    for unions. */
  template <class T> bool closeBoundaryLoop(std::vector<std::vector<TPoint<T>>> &cutUV,
    std::vector<std::vector<std::vector<TPoint<T>>>> &loops, 
    T bigtolerance = 0.01, T parmtolerance = PARM_TOLERANCE, int numdivisions = 100,
    T maxparmgap = 0.001)
 //!!!!!!!   T maxparmgap = 1.1 / T(MANY_POINTS2D - 1))
  {
    if (cutUV.empty())
      return false;

    // step 1 : prepare outer loop
    std::vector<std::vector<TPoint<T>>> outerloop;

    if (!loops.empty()) //!!!
    {
      outerloop = loops[0];
    } else
    {
      closeOuterBoundaryLoop(outerloop,numdivisions);
    }

    // step 2 : combine cut pieces into a single line
    std::vector<std::vector<TPoint<T>>> cut = cutUV;

    // inner loops
    std::vector<std::vector<std::vector<TPoint<T>>>> innerloops;
#if 0
    int n = extractLoops(cut,innerloops,parmtolerance,parmtolerance); //!!!!!!!
#else
    int n = extractLoops(cut,innerloops,parmtolerance * 1000.0,parmtolerance); //!!!!!!!
#endif
  
    // all done
    if (n)
    {
      // get loop direction, is it a hole or space around the hole?
      int numholes = 0;
      int numpatches = 0;
      for (int i = 0; i < n; i++)
      {
        // this is a hole, we need an outer loop as loop[0]
        //!!! incorrect - sure it is not possible to have both holes and patches at the same time
        // it is possible if a hole is inside another hole, it has the "outer" direction
        bool hole = !outerLoop(innerloops[i],parmtolerance);
        if (hole)
        {
          numholes++;
        } else
        {
          numpatches++;
        }
      }

      // there are holes, we need an outer loop as loop[0]
      if (numholes)
      {
        if (loops.empty()) //!!!
          loops.insert(loops.end(),outerloop); 
      }

      loops.insert(loops.end(),innerloops.begin(),innerloops.end());

      if (cut.empty()) 
        return true;
    }

    // step 3 : order pieces to connect start to ends
    std::vector<std::vector<std::vector<TPoint<T>>>> allordered;
    if (!orderCutPieces(cut,allordered,bigtolerance,parmtolerance))
    {
      // ... second attempt
      std::vector<bool> busy(cut.size(),false);

      // try to extend cut to reach the boundary from both ends
      T mindist0 = 0.0;
      T mindist1 = 0.0;
      TPoint<T> intr0, intr1;
      int closest0 = findClosestStart(cut,mindist0,intr0,&busy);
      int closest1 = findClosestEnd(cut,mindist1,intr1,&busy);

      if (closest0 >= 0 && mindist0 < maxparmgap)
      {
        cut[closest0].front() = intr0;
      }

      if (closest1 >= 0 && mindist1 < maxparmgap)
      {
        cut[closest1].back() = intr1;
      }

      // again
      if (!orderCutPieces(cut,allordered,bigtolerance,parmtolerance))
      {
        return false;
      }
    }

//if (print)
//{
//  TPoint<T> p0 = this->position(allordered[0][0].front().X,allordered[0][0].front().Y);
//  TPoint<T> p1 = this->position(allordered[1][0].back().X,allordered[1][0].back().Y);
//
//  T UVdist = !(allordered[0][0].front() - allordered[1][0].back());
//  T dist = !(p1 - p0);
//
//  outputDebugString("dist " + to_string(dist,12) + " UVdist " + to_string(UVdist,12));
//}

    // two cuts from boundary to boundary with one same point
    // combine two boundary pieces into one to make a single loop
    bool twocutsfromboundarytoboundary = false;
    if (allordered.size() == 2)
    {
      if (boundaryPoint<T>(allordered[0].front().front(),parmtolerance) &&
        boundaryPoint<T>(allordered[0].front().back(),parmtolerance) &&
        boundaryPoint<T>(allordered[1].front().front(),parmtolerance) &&
        boundaryPoint<T>(allordered[1].front().back(),parmtolerance))
      {
        T dist0 = !(allordered[0].front().back() - allordered[1].front().front());
        T dist1 = !(allordered[1].front().back() - allordered[0].front().front());

        if (dist0 < parmtolerance * 1000.0) //!!!!!!!
        {
          allordered[0].front().back() = allordered[1].front().front() = 
            (allordered[0].front().back() + allordered[1].front().front()) * 0.5;

          std::vector<std::vector<std::vector<TPoint<T>>>> newallordered;
          newallordered.push_back(allordered[0]);
          newallordered[0].insert(newallordered[0].end(),allordered[1].begin(),allordered[1].end());
          allordered = newallordered;

          twocutsfromboundarytoboundary = true;
        } else if (dist1 < parmtolerance * 1000.0) //!!!!!!!
        {
          allordered[1].front().back() = allordered[0].front().front() = 
            (allordered[1].front().back() + allordered[0].front().front()) * 0.5;

          std::vector<std::vector<std::vector<TPoint<T>>>> newallordered;
          newallordered.push_back(allordered[1]);
          newallordered[0].insert(newallordered[0].end(),allordered[0].begin(),allordered[0].end());
          allordered = newallordered;

          twocutsfromboundarytoboundary = true;
        }
      }
    }

    // step 4 : embed cut into outer loop or close cuts
    // outerloop contains parts; cut is a single curve
    for (int i = 0; i < int(allordered.size()); i++)
    {
      // special case
      bool newloop = false;

redo:

      std::vector<std::vector<TPoint<T>>> cut = allordered[i];

      // make a single continuous curve from outer loop with marking sharp corners :
      // they mark every START of new curve piece
      std::vector<TPoint<T>> allpoints;
      std::vector<TPoint<T>> cutpoints;

      // find intersections, construct new points
      std::vector<TPoint<T>> UV;

      int numintrs = intersectLoopByCut(outerloop,cut,allpoints,cutpoints,UV,parmtolerance * 1000.0);

      // remove middle point
      if (twocutsfromboundarytoboundary && numintrs == 3)
      {
        UV.erase(UV.begin() + 1);
        numintrs--;
      }

      // new loop all inside single face
      if (numintrs < 2)
      {
        if (!newloop)
        {
          // it maybe another loop, recreate full outer loop
          closeOuterBoundaryLoop(outerloop,numdivisions);

          newloop = true;
          goto redo;
        }
      } else if (numintrs == 2)
      {
        // these will be cut points from intersection 1 to intersection 2
        std::vector<TPoint<T>> newpoints;

        // combine cut parts from cut and loop
        if (!makeNewLoopPoints<T>(cutpoints,allpoints,UV,newpoints,true,parmtolerance))
          return false;

        // divide outer loop back into parts
        std::vector<std::vector<TPoint<T>>> newouterloop; 
        divideByDuplicates<T>(newpoints,newouterloop,parmtolerance);

        // correct points on the boundary
        correctBoundaryPoints(newouterloop);

        // update outer loops[0]
        if (loops.empty() || newloop)
        {
          loops.push_back(newouterloop);
        } else
        {
          loops[0] = newouterloop;
        }
        outerloop = newouterloop;
      } else
      {
        continue;
      }
    }

    return true;
  }

  /** Get min/max from control points. */
  virtual std::pair<TPoint<T>,TPoint<T>> getMinMax()
  {
    std::pair<TPoint<T>,TPoint<T>> minmax;
    tcad::calculateMinMax<T>(this->cpoints,&minmax.first,&minmax.second);
    return minmax;
  }

public: //!!!

  // number of columns minus 1
  int K1 = 0;
  // number of rows munus 1
  int K2 = 0;

  // name for debugging purposes
  std::string name;

protected:

  // (K1 + 1) * (K2 + 1) control points, call update() after every change
  std::vector<TPoint<T>> cpoints;

};

/** Close outer UV boundary by 4 pieces. */
template <class T> void closeOuterBoundaryLoop(std::vector<std::vector<TPoint<T>>> &closedboundary, 
  int numdivisions = 100)
{
  closedboundary.clear();

  for (int i = 0; i < 4; i++)
  {
    int i1 = i + 1;
    if (i1 > 3)
      i1 = 0;

    TPoint<T> UV = cornerUV<T>[i];
    TPoint<T> nextUV = cornerUV<T>[i1];

    TPointCurve<T> line(UV,nextUV,numdivisions);

    // we shall keep a checksum at fronts of boundary pieces
    line.controlPoints().front().W = 0.0;

    closedboundary.push_back(line.controlPoints());
  }
}

/** Calculate min/max from control points. */
template <class T> bool calculateMinMax(std::vector<TBaseSurface<T> *> &surfaces,
  TPoint<T> &min, TPoint<T> &max)
{
  bool ok = false;

  for (int i = 0; i < int(surfaces.size()); i++)
  {
    TPoint<T> mi,ma;
    
    if (tcad::calculateMinMax(surfaces[i]->controlPoints(),&mi,&ma))
    {
      if (!ok)
      {
        min = mi;
        max = ma;
      } else
      {
        min = pointMin(min,mi);
        max = pointMax(max,ma);
      }

      ok = true;
    }
  }

  return ok;
}

/** Is parameter on edge? This stuff is for improvement of surface intersections. */
template <class T> bool parmOnEdge(T u, T v, bool onedge[4], T t[4], T tolerance = PARM_TOLERANCE)
{
                              // on any edge?
  bool result = false;

                              // zero all
  for (int i = 0; i < 4; i++)
  {
    onedge[i] = false;
    t[i] = -1;
  }

                              // close to edges?
  bool u0 = (std::abs(u) < tolerance);
  bool u1 = (std::abs(u - 1) < tolerance);
  bool v0 = (std::abs(v) < tolerance);
  bool v1 = (std::abs(v - 1) < tolerance);

  if (u0)
  {
    onedge[3] = true;
    t[3] = v;
    result = true;
  } else if (u1)
  {
    onedge[1] = true;
    t[1] = v;
    result = true;
  } 

  if (v0)
  {
    onedge[0] = true;
    t[0] = u;
    result = true;
  } else if (v1)
  {
    onedge[2] = true;
    t[2] = u;
    result = true;
  }

#if 0
DebugOutput(CString("u ") + CString(u,10) + " v " +CString(v,10) + " res = " + 
  CString((int) result));
#endif

  return result;
}

//        U
//    3--->------2----------2
//    |                     |
//    3                     |
//    |                     1
//  V ^                     ^ V
//    |                     |
//    0--->------0----------1
//        U

/** Is parameter on edge? This stuff is for improvement of surface intersections. */
template <class T> bool parmOnEdge(T u, T v, bool onedge[4], T t[4], bool oncorner[4], T tolerance = PARM_TOLERANCE)
{
                              // on any edge?
  bool result = parmOnEdge(u,v,onedge,t,tolerance);

  oncorner[0] = onedge[3] && onedge[0];
  oncorner[1] = onedge[0] && onedge[1];
  oncorner[2] = onedge[1] && onedge[2];
  oncorner[3] = onedge[2] && onedge[3];

  return result;
}

/** This stuff is for improvement of surface intersections. */
template <class T> bool solveSystemOneParmFixed(TBaseSurface<T> *F, TBaseSurface<T> *G, TPoint<T> &parms, 
  T A[16], T B[4], int type, T relaxcoef)
{
  TPoint<T> F0 = F->position(parms.X,parms.Y);
  TPoint<T> G0 = G->position(parms.Z,parms.W);
  TPoint<T> Fu = F->derivative(parms.X,parms.Y,PARAMETER_U,1);
  TPoint<T> Fv = F->derivative(parms.X,parms.Y,PARAMETER_V,1);
  TPoint<T> Gs = G->derivative(parms.Z,parms.W,PARAMETER_U,1);
  TPoint<T> Gt = G->derivative(parms.Z,parms.W,PARAMETER_V,1);

  makeSystemOneParmFixed(F0,G0,Fu,Fv,Gs,Gt,A,B);
  makeSystemOneParmFixedLastEquation(type,A,B);

  bool ok = solveSystem4x4<T>(A,B,TOLERANCE(T));

  if (ok)
  {
    parms.X += B[0] * relaxcoef;
    parms.Y += B[1] * relaxcoef;
    parms.Z += B[2] * relaxcoef;
    parms.W += B[3] * relaxcoef;

    LIMIT(parms.X,0,1);
    LIMIT(parms.Y,0,1);
    LIMIT(parms.Z,0,1);
    LIMIT(parms.W,0,1);
  }

  return ok;
}

//        U
//    3--->------2----------2
//    |                     |
//    3                     |
//    |                     1
//  V ^                     ^ V
//    |                     |
//    0--->------0----------1
//        U

/** This stuff is for improvement of surface intersections. */
template <class T> bool solveSystemBoundary(TBaseSurface<T> *F, TBaseSurface<T> *G, TPoint<T> &parms, 
  T A[9], T B[3], int edge, T relaxcoef)
{
                              // force G boundary coordinates
  switch (edge) {
    case 0 : parms.W = 0; break;
    case 1 : parms.Z = 1; break;
    case 2 : parms.W = 1; break;
    case 3 : parms.Z = 0; break;
    default : assert(false); break;
  }

  TPoint<T> F0 = F->position(parms.X,parms.Y);
  TPoint<T> G0 = G->position(parms.Z,parms.W);
  TPoint<T> Fu = F->derivative(parms.X,parms.Y,PARAMETER_U,1);
  TPoint<T> Fv = F->derivative(parms.X,parms.Y,PARAMETER_V,1);
  TPoint<T> Gs = G->derivative(parms.Z,parms.W,PARAMETER_U,1);
  TPoint<T> Gt = G->derivative(parms.Z,parms.W,PARAMETER_V,1);

  makeSystemBoundary(F0,G0,Fu,Fv,Gs,Gt,A,B,(edge == 0 || edge == 2));

  bool ok = solveSystemWithPivoting<T,int>(3,A,B,TOLERANCE(T));

  if (ok)
  {
    parms.X += B[0] * relaxcoef;
    parms.Y += B[1] * relaxcoef;

    switch (edge) {
      case 0 : case 2 : parms.Z += B[2] * relaxcoef; break;
      case 1 : case 3 : parms.W += B[2] * relaxcoef; break;
      default : assert(false); break;
    }

    LIMIT(parms.X,0,1);
    LIMIT(parms.Y,0,1);
    LIMIT(parms.Z,0,1);
    LIMIT(parms.W,0,1);
  }

  return ok;
}

/** This stuff is for improvement of surface intersections. */
template <class T> void makeSystemAllCases(TBaseSurface<T> *F, TBaseSurface<T> *G, TPoint<T> &parms,  
  T A[16], T B[4], TPoint<T> &inc)
{
                              // get derivatives
  TPoint<T> F0 = F->position(parms.X,parms.Y);
  TPoint<T> G0 = G->position(parms.Z,parms.W);
  TPoint<T> Fu = F->derivative(parms.X,parms.Y,PARAMETER_U,1);
  TPoint<T> Fv = F->derivative(parms.X,parms.Y,PARAMETER_V,1);
  TPoint<T> Gs = G->derivative(parms.Z,parms.W,PARAMETER_U,1);
  TPoint<T> Gt = G->derivative(parms.Z,parms.W,PARAMETER_V,1);

                              // first three equations
                              // they may be linearly dependent, e.g. for
                              // two coplanar surfaces, but leave it for now
  makeSystemOneParmFixed(F0,G0,Fu,Fv,Gs,Gt,A,B);

                              // fouth equation is relationship between 
                              // parameter increments
  A[3 * 4 + 0] = inc.X;
  A[3 * 4 + 1] = inc.Y;
  A[3 * 4 + 2] = inc.Z;
  A[3 * 4 + 3] = inc.W;
  B[3] = 0;
}

/** This stuff is for improvement of surface intersections. */
template <class T> bool improveIntersection(TBaseSurface<T> *F, TBaseSurface<T> *G, TPoint<T> &parms, 
  int maxiter = 100, T relaxcoef = 0.5, T tolerance = PARM_TOLERANCE, T maxparmchange = 0.1)
{
  TPoint<T> bestparms;
  TPoint<T> initparms = parms;

#ifdef _DEBUG
  T AA[16] = {0};
#endif

                              // original values
  parms = initparms;
  bestparms = parms;
  T bestdist = !(F->position(parms.X,parms.Y) - G->position(parms.Z,parms.W));

  for (int i = 0; i < maxiter; i++)
  {
                              // form a system to get new parameters
    TPoint<T> F0 = F->position(parms.X,parms.Y);
    TPoint<T> G0 = G->position(parms.Z,parms.W);
    TPoint<T> Fu = F->derivative(parms.X,parms.Y,PARAMETER_U,1);
    TPoint<T> Fv = F->derivative(parms.X,parms.Y,PARAMETER_V,1);
    TPoint<T> Gs = G->derivative(parms.Z,parms.W,PARAMETER_U,1);
    TPoint<T> Gt = G->derivative(parms.Z,parms.W,PARAMETER_V,1);

                              // right-hand size and solution
    T B[4] = {0};
    bool systemsolved = true;

    bool onedge1[4]; T t1[4]; 
    bool onedge2[4]; T t2[4];
    bool edge1 = parmOnEdge(parms.X,parms.Y,onedge1,t1,tolerance);
    bool edge2 = parmOnEdge(parms.Z,parms.W,onedge2,t2,tolerance);

                              // which system to solve?
    if (edge1 || edge2)
    {
                               // matrix 4 x 4
      T A[16] = {0};
      makeSystemOneParmFixed(F0,G0,Fu,Fv,Gs,Gt,A,B);

      if (edge1)
      {
        if (onedge1[0])
        {
                              // dv = 0
          A[3 * 4 + 1] = 1;
        } 
        if (onedge1[1])
        {
                              // du = 0
          A[3 * 4 + 0] = 1;
        } 
        if (onedge1[2])
        {
                              // dv = 0
          A[3 * 4 + 1] = 1;
        } 
        if (onedge1[3])
        {
                              // du = 0
          A[3 * 4 + 0] = 1;
        }
      }

      if (edge2)
      {
        if (onedge2[0])
        {
                              // dt = 0
          A[3 * 4 + 3] = 1;
        } 
        if (onedge2[1])
        {
                              // ds = 0
          A[3 * 4 + 2] = 1;
        } 
        if (onedge2[2])
        {
                              // dt = 0
          A[3 * 4 + 3] = 1;
        } 
        if (onedge2[3])
        {
                              // ds = 0
          A[3 * 4 + 2] = 1;
        }
      }
      B[3] = 0;

  #ifdef _DEBUG
      memmove(AA,A,sizeof(A));
  #endif

      if (!solveSystem4x4<T>(A,B,TOLERANCE(T)))
      {
        parms = bestparms;
        systemsolved = false;
      }
    } else
    {
                              // matrix 3 x 4
      T A[16] = {0};

                              // try to solve the system with one parameter fixed
                              // in all four directions
      bool ok4 = solveSystemOneParmFixed(F,G,parms,A,B,0,relaxcoef) &&
        solveSystemOneParmFixed(F,G,parms,A,B,1,relaxcoef) &&
        solveSystemOneParmFixed(F,G,parms,A,B,2,relaxcoef) &&
        solveSystemOneParmFixed(F,G,parms,A,B,3,relaxcoef);

      if (!ok4)
      {
        parms = bestparms;

        if (!solveSystemUnderdetermined3x4(F0,G0,Fu,Fv,Gs,Gt,A,B))
        {
          parms = bestparms;
          systemsolved = false;
        }
      }
    }

    if (systemsolved)
    {
      T dp = std::max<T>(
        std::max<T>(std::abs(B[0]),std::abs(B[1])),
        std::max<T>(std::abs(B[2]),std::abs(B[3])));

      if (dp > maxparmchange)
      {
        parms = bestparms;
        return false;
      }

      parms.X += B[0] * relaxcoef;
      parms.Y += B[1] * relaxcoef;
      parms.Z += B[2] * relaxcoef;
      parms.W += B[3] * relaxcoef;

      LIMIT(parms.X,0,1);
      LIMIT(parms.Y,0,1);
      LIMIT(parms.Z,0,1);
      LIMIT(parms.W,0,1);
    }

    T dist = !(F->position(parms.X,parms.Y) - G->position(parms.Z,parms.W));
    if (dist < tolerance)
    {
      return true;
    }

    if (dist < bestdist)
    {
      bestparms = parms;
      bestdist = dist;
    }
  }
  
  parms = bestparms;
  return false;
} 

/** This stuff is for improvement of surface intersections. */
template <class T> bool improveIntersectionSimple(TBaseSurface<T> *F, TBaseSurface<T> *G, TPoint<T> &parms, 
  TPoint<T> &inc, int maxiter = 100, T relaxcoef = 0.5, T tolerance = PARM_TOLERANCE)
{
  TPoint<T> bestparms;
  TPoint<T> initparms = parms;

#ifdef _DEBUG
  T AA[16] = {0};
#endif

                              // original values
  parms = initparms;
  bestparms = parms;
  T bestdist = !(F->position(parms.X,parms.Y) - G->position(parms.Z,parms.W));

  for (int i = 0; i < maxiter; i++)
  {
    bool systemsolved = true;

#if 1
    std::array<T,12> A; 
    std::array<T,4> B;

    TPoint<T> F0 = F->position(parms.X,parms.Y);
    TPoint<T> G0 = G->position(parms.Z,parms.W);
    TPoint<T> Fu = F->derivative(parms.X,parms.Y,PARAMETER_U,1);
    TPoint<T> Fv = F->derivative(parms.X,parms.Y,PARAMETER_V,1);
    TPoint<T> Gs = G->derivative(parms.Z,parms.W,PARAMETER_U,1);
    TPoint<T> Gt = G->derivative(parms.Z,parms.W,PARAMETER_V,1);

    if (!solveSystemUnderdetermined3x4(F0,G0,Fu,Fv,Gs,Gt,A,B))
    {
      parms = bestparms;
      systemsolved = false;
      break;
    }

    if (systemsolved)
    {
      parms.X += B[0] * relaxcoef;
      parms.Y += B[1] * relaxcoef;
      parms.Z += B[2] * relaxcoef;
      parms.W += B[3] * relaxcoef;

      LIMIT(parms.X,0,1);
      LIMIT(parms.Y,0,1);
      LIMIT(parms.Z,0,1);
      LIMIT(parms.W,0,1);
    }

#else
                              // right-hand size and solution
    T B[4] = {0};

                              // matrix 4 x 4
    T A[16] = {0};

                              // make system
    makeSystemAllCases(F,G,parms,A,B,inc);

 #ifdef _DEBUG
      memmove(AA,A,sizeof(A));
 #endif
                              // solve system
    if (!solveSystem4x4<T>(A,B,TOLERANCE(T)))
    {
      parms = bestparms;
      systemsolved = false;
    }

    if (systemsolved)
    {
      parms.X += B[0] * relaxcoef;
      parms.Y += B[1] * relaxcoef;
      parms.Z += B[2] * relaxcoef;
      parms.W += B[3] * relaxcoef;

      LIMIT(parms.X,0,1);
      LIMIT(parms.Y,0,1);
      LIMIT(parms.Z,0,1);
      LIMIT(parms.W,0,1);
    }

#endif

    T dist = !(F->position(parms.X,parms.Y) - G->position(parms.Z,parms.W));
    if (dist < tolerance)
    {
      return true;
    }

    if (dist < bestdist)
    {
      bestparms = parms;
      bestdist = dist;
    }
  }
  
  parms = bestparms;
  return false;
} 

/** Cut inside loop? */
template <class T> bool cutInsideLoop(std::vector<std::vector<TPoint<T>>> &cut,
  std::vector<std::vector<TPoint<T>>> &loop, T parmtolerance = PARM_TOLERANCE)
{
  //// test 1 : intersections
  //std::vector<TPoint<T>> newpoints;
  //std::vector<TPoint<T>> UV;

  // cut curve specifies a correct direction of the loop, it must go first here
  //for (int i = 0; i < int(cut.size()); i++)
  //{
  //  for (int j = 0; j < int(loop.size()); j++)
  //  {
  //    int numintrs = findIntersections(cut[i],loop[j],UV,parmtolerance); 

  //    if (numintrs > 0)
  //      return true;
  //  }
  //}

  // test 2 : inside?
  for (int i = 0; i < int(cut.size()); i++)
  {    
    for (int j = 0; j < int(cut[i].size()); j++)
    {
      bool inside = insideBoundary<T>(loop,cut[i][j],parmtolerance);
      if (!inside)
        return false;
    }
  }

  return true;
}

/** Cut inside loops? */
template <class T> bool cutInsideLoops(std::vector<std::vector<TPoint<T>>> &cut,
  std::vector<std::vector<std::vector<TPoint<T>>>> &loops, T parmtolerance = PARM_TOLERANCE)
{
  // loops not initialised
  if (loops.empty())
    return true;

  // test loop 0 (normally outer)
  bool loop0outer = outerLoop(loops[0],parmtolerance);

  if (loop0outer)
  {
    // full standard outer loop
    if (boundaryPoints(loops[0],parmtolerance))
      return true;

    //// one outer loop and holes //!!!
    //if (cutInsideLoop(cut,loops[0],parmtolerance)) 
    //  return true;

    // one outer loop and holes
    if (cutInsideLoop(cut,loops[0],parmtolerance)) 
    {
      // cut maybe inside a hole
      for (int i = 1; i < int(loops.size()); i++)
      {
        if (cutInsideLoop(cut,loops[i],parmtolerance))
          return false;
      }

      return true;
    } else
    {
      return false;
    }
  } else
  {
    // all loops are holes, test if cut is inside any loop
    for (int i = 0; i < int(loops.size()); i++)
    {
      if (cutInsideLoop(cut,loops[i],parmtolerance))
        return true;
    }
  }

  return false;
}

}
