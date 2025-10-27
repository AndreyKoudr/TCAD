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
#include "ttransform.h"
#include "tbasecurve.h"
#include "ttriangles.h"

namespace tcad {

// U or V or W
enum Parameter {
  PARAMETER_U,      // for surfaces and volumes
  PARAMETER_V,      // for surfaces and volumes
  PARAMETER_W,      // for volumes

  PARAMETER_ANY
};

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

  /** Same as 0-th derivative*/
  virtual TPoint<T> position(T U, T V)
  {
    return this->derivative(U,V,PARAMETER_ANY,0);
  }

  /** Get control points. */
  virtual std::vector<TPoint<T>> &controlPoints()
  {
    return cpoints;
  }

  /** Get row of control points. */
  void getRow(int index, std::vector<TPoint<T>> &row)
  {
    getRow(this->cpoints,K1,K2,row);
  }

  /** Set row of control points. */
  void setRow(int index, const std::vector<TPoint<T>> &row)
  {
    setRow(this->cpoints,K1,K2,row);
  }

  /** Get column of control points. */
  void getColumn(int index, std::vector<TPoint<T>> &column)
  {
    getColumn(this->cpoints,K1,K2,column);
  }

  /** Set column of control points. */
  void setColumn(int index, const std::vector<TPoint<T>> &column)
  {
    setColumn(this->cpoints,K1,K2,column);
  }

  /** Generate a uniform set of actual points on the surface to
    (1) create another type of surface from these points (new K1 = numpointsU - 1),
        (new K2 = numpointsV - 1) or 
    (2) compare/involve another type of geometry to compare/interact 

    A set of according parameter values is generated as well.
  */
  virtual void createPoints(std::vector<TPoint<T>> &points, 
    int *k1 = nullptr, int *k2 = nullptr,
    int numpointsU = MANY_POINTS, int numpointsV = MANY_POINTS)
  {
    points.clear();

    for (int i = 0; i < numpointsV; i++)
    {
      T V = T(i) / T(numpointsV - 1);
      for (int j = 0; j < numpointsU; j++)
      {
        T U = T(j) / T(numpointsU - 1);

        TPoint<T> p = this->position(U,V);
        points.push_back(p);
      }
    }

    if (k1)
      *k1 = numpointsU - 1;
    if (k2)
      *k2 = numpointsV - 1;
  }

  /** Create triangles on a regular set points. */
  virtual bool createTriangles(TTriangles<T> &tris,
    int numpointsU = MANY_POINTS, int numpointsV = MANY_POINTS)
  {
    tris.clear();

    std::vector<TPoint<T>> points;
    int k1 = 0;
    int k2 = 0;

    createPoints(points,&k1,&k2,numpointsU,numpointsV);

    for (int i = 0; i < k2; i++)
    {
      std::vector<TPoint<T>> row0,row1;
      tcad::getRow(points,k1,k2,i,row0);
      tcad::getRow(points,k1,k2,i + 1,row1);

      for (int j = 0; j < int(row0.size() - 1); j++)
      {
        TPoint<T> p0 = row0[j];
        TPoint<T> p1 = row0[j + 1];
        TPoint<T> p2 = row1[j + 1];
        TPoint<T> p3 = row1[j];

        T d0 = !(p2 - p0);
        T d1 = !(p3 - p1);

        if (d0 < d1)
        {
          tris.addTri(p2,p0,p1,0.0);
          tris.addTri(p0,p2,p3,0.0);
        } else
        {
          tris.addTri(p1,p3,p0,0.0);
          tris.addTri(p3,p1,p2,0.0);
        }
      }
    }

    std::pair<TPoint<T>,TPoint<T>> mm = tris.minmax();
    TPoint<T> d = mm.second - mm.first;
    T tolerance = !d * PARM_TOLERANCE;

    // remove duplicate nodes and renumber corners
    bool ok = tris.buildConnectivityArray(tolerance);

    return ok;
  }

  /** Reverse U. Normal is changed to opposite. */
  void reverseU()
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
  void reverseV()
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
  void Usize()
  {
    T len = 0.0;
    for (int i = 0; i <= K2; i++)
    {
      std::vector<TPoint<T>> row;
      this->getRow(i,row);
      len += calculateLength(row);
    }

    len /= T(K2 + 1);
  }

  /** Approximate size along V. */
  void Vsize()
  {
    T len = 0.0;
    for (int i = 0; i <= K1; i++)
    {
      std::vector<TPoint<T>> col;
      this->getColumn(i,col);
      len += calculateLength(col);
    }

    len /= T(K1 + 1);
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
  virtual bool equal(TBaseSurface &other, T tolerance, int numpointsU = MANY_POINTS, 
    int numpointsV = MANY_POINTS)
  {
    std::vector<TPoint<T>> points,otherpoints;
    createPoints(points,nullptr,nullptr,numpointsU,numpointsV);
    other.createPoints(otherpoints,nullptr,nullptr,numpointsU,numpointsV);

    T diff = difference(points,otherpoints);
    return (diff >= 0.0 && diff < tolerance);
  }

  /** Find values of parameters U,V for a point on (or close to) the curve. Returns -1 in failure. */
  virtual TPoint<T> findUVforPoint(TPoint<T> p, int numpointsU = MANY_POINTS, int numpointsV = MANY_POINTS)
  {
    // create fine mesh of points
    std::vector<TPoint<T>> points;
    int k1 = 0;
    int k2 = 0;
    createPoints(points,&k1,&k2,numpointsU,numpointsV);

    // temp
    TPoint<T> proj;
    int seg = 0;
    T u = 0.0;

    // look for closest point on rows
    T U = -1.0;
    T mindist = std::numeric_limits<T>::max();
    for (int i = 0; i < k2 + 1; i++)
    {
      std::vector<TPoint<T>> row;
      tcad::getRow(points,k1,k2,i,row);
      if (projectPointOnPoints(row,p,proj,&seg,&u))
      {
        T dist = !(p - proj);
        if (dist < mindist)
        {
          TPoint<T> p0 = row[seg];
          TPoint<T> p1 = row[seg + 1];
          U = p0.W + (p1.W - p0.W) * u;
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
          V = p0.W + (p1.W - p0.W) * u;
          mindist = dist;
        }
      }
    }

    return TPoint<T>(U,V);
  }

  /** Calculate min/max. */
  virtual bool calculateMinMax(TPoint<T> *min, TPoint<T> *max, TPoint<T> *imin = nullptr, 
    TPoint<T> *imax = nullptr, int numpointsU = MANY_POINTS, int numpointsV = MANY_POINTS)
  {
    std::vector<TPoint<T>> points;
    createPoints(points,nullptr,nullptr,numpointsU,numpointsV);

    return tcad::calculateMinMax(points,min,max,imin,imax);
  }

  /** Cut by plane. */
  bool cutByPlane(TPlane<T> &plane, std::vector<std::vector<TPoint<T>>> &lines, 
    T tolerance, int numpointsU = MANY_POINTS, int numpointsV = MANY_POINTS)
  {
    TTriangles tris;
    if (createTriangles(tris,numpointsU,numpointsV))
    {
      return tris.cutByPlane(plane,lines,tolerance);
    } else
    {
      return false;
    }
  }

  /** Cut out a part of surface from U0 to U1 and from V0 to V1 into list of points. */
  template <class T> void cutPiece(T Ufrom, T Uto, T Vfrom, T Vto, std::vector<TPoint<T>> &points,
    int numpointsU = MANY_POINTS, int numpointsV = MANY_POINTS)
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

  /** Find intersection curve(s) with another surface. Returns number od intersection curves. */
  template <class T> int intersect(TBaseSurface<T> &other, std::vector<std::vector<TPoint<T>>> &intersections, 
    T tolerance, T parmtolerance = TOLERANCE(T), 
    int numpointsU = MANY_POINTS, int numpointsV = MANY_POINTS)
  {
    TTriangles<T> tris,othertris;

    if (createTriangles(tris,numpointsU,numpointsV) &&
      other.createTriangles(othertris,numpointsU,numpointsV))
    {
      if (tris.intersect(othertris,intersections,tolerance,parmtolerance))
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

protected:

  // number of columns minus 1
  int K1 = 0;
  // number of rows munus 1
  int K2 = 0;

  // (K1 + 1) * (K2 + 1) control points, call update() after every change
  std::vector<TPoint<T>> cpoints;

};

}
