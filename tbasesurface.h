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

  /** Generate a uniform set of actual points on the surface to
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

  /** Create triangles on a regular set points.
    Set refine... to 0.5 at a corresponding end to refine,
    1.0 has no effect. */
  virtual bool createTriangles(TTriangles<T> &tris,
    int numpointsU = MANY_POINTS2D, int numpointsV = MANY_POINTS2D,
    T refinestartU = 1.0, T refineendU = 1.0, 
    T refinestartV = 1.0, T refineendV = 1.0)
  {
    tris.clear();

    std::vector<TPoint<T>> points;
    std::vector<TPoint<T>> UVpoints;
    int k1 = 0;
    int k2 = 0;

    createPoints(points,&UVpoints,&k1,&k2,numpointsU,numpointsV,refinestartU,refineendU,refinestartV,refineendV);

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
  virtual bool equal(TBaseSurface &other, T tolerance, int numpointsU = MANY_POINTS2D, 
    int numpointsV = MANY_POINTS2D)
  {
    std::vector<TPoint<T>> points,otherpoints;
    createPoints(points,nullptr,nullptr,nullptr,numpointsU,numpointsV);
    other.createPoints(otherpoints,nullptr,nullptr,nullptr,numpointsU,numpointsV);

    T diff = difference(points,otherpoints);
    return (diff >= 0.0 && diff < tolerance);
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

  /** Find intersection curve(s) with another surface. Returns number of intersection curves. 
    boundary0/1 contain U,V in X,Y for every intersection curve for both surfaces. */
  template <class T> int intersect(TBaseSurface<T> &other, bool bodyleft, std::vector<std::vector<TPoint<T>>> &intersections, 
    std::vector<std::vector<TPoint<T>>> &boundary0, std::vector<std::vector<TPoint<T>>> &boundary1,
    T parmtolerance = PARM_TOLERANCE, 
    int numpointsU0 = MANY_POINTS2D, int numpointsV0 = MANY_POINTS2D,
    T refinestartU0 = 1.0, T refineendU0 = 1.0, 
    T refinestartV0 = 1.0, T refineendV0 = 1.0,
    int numpointsU1 = MANY_POINTS2D, int numpointsV1 = MANY_POINTS2D,
    T refinestartU1 = 1.0, T refineendU1 = 1.0, 
    T refinestartV1 = 1.0, T refineendV1 = 1.0,
    bool debug = false)
  {
    TTriangles<T> tris,othertris;

    if (createTriangles(tris,numpointsU0,numpointsV0,
      refinestartU0,refineendU0,refinestartV0,refineendV0) &&
      other.createTriangles(othertris,numpointsU1,numpointsV1,
      refinestartU1,refineendU1,refinestartV1,refineendV1))
    {
#ifdef _DEBUG
      if (debug)
      {
        tris.tcad::TTriangles<T>::saveSTL("intersection0.stl","TCAD",true);
        othertris.tcad::TTriangles<T>::saveSTL("intersection1.stl","TCAD",true);
      }
#endif

      if (tris.intersect(othertris,bodyleft,intersections,parmtolerance,&boundary0,&boundary1))
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

  /** Find a cut piece (except busy) close to p. p is "tail", cut front is "head" to be
    connected to tail. */
  template <class T> int findClosest(std::vector<std::vector<TPoint<T>>> &cut,
    std::vector<bool> &busy, TPoint<T> p, T tolerance)
  {
    // find closest piece to the current loop end among not busy
    int closest = -1;

    T mindist = std::numeric_limits<T>::max();
    for (int i = 0; i < int(cut.size()); i++)
    {
      if (busy[i])
        continue;

      T dist = !(p - cut[i].front());
      if (dist < tolerance && dist < mindist)
      {
        mindist = dist;
        closest = i;
      }
    }

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
      int closest = findClosest(cut,busy,ordered.back().back(),tolerance);

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
      int closest = findClosest(cut,busy,loop.back().back(),tolerance);

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

  /** Prepare points before intersection. Divide peieces by duplicate points. */
  template <class T> void cutIntoPoints(std::vector<std::vector<TPoint<T>>> &cut, std::vector<TPoint<T>> &points)
  {
    points.clear();

    int count = 0;
    for (int i = 0; i < int(cut.size()); i++)
    {
      for (int j = 0; j < int(cut[i].size()); j++) // points may be duplicate, they mark new part
      {
        TPoint<T> p = cut[i][j];
        p.W = T(count++);
        points.push_back(p);
      }
    }
  }

  /** Close UV boundary. cutUV is a cut across surface in UV coordinates. cutFromOuter is true
    for unions. */
  template <class T> bool closeBoundaryLoop(std::vector<std::vector<TPoint<T>>> &cutUV,
    std::vector<std::vector<std::vector<TPoint<T>>>> &loops, 
    T bigtolerance = 0.01, T parmtolerance = PARM_TOLERANCE, int numdivisions = 100,
    T maxparmgap = 0.001)
  {
    if (cutUV.empty())
      return false;

    // step 1 : prepare outer loop
    std::vector<std::vector<TPoint<T>>> outerloop;
    closeOuterBoundaryLoop(outerloop,numdivisions);

    // step 2 : combine cut pieces into a single line
    std::vector<std::vector<TPoint<T>>> cut = cutUV;

    // inner loops
    std::vector<std::vector<std::vector<TPoint<T>>>> innerloops;
    int n = extractLoops(cut,innerloops,parmtolerance,parmtolerance);
 //!!!   int n = extractLoops(cut,innerloops,bigtolerance,parmtolerance);
  
    // all done
    if (n)
    {
      // get loop direction, is it a hole or space around the hole?
      int numholes = 0;
      int numpatches = 0;
      for (int i = 0; i < n; i++)
      {
        TPoint<T> loopnormal = getLoopNormal(innerloops[i],parmtolerance);

        // this is a hole, we need an outer loop as loop[0]
        // sure it is not possible to have both holes and patches at the same time
        bool hole = !(loopnormal > TPoint<T>(0.0,0.0,1.0));
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
      cutIntoPoints(outerloop,allpoints);

      // set W for cut
      std::vector<TPoint<T>> cutpoints;
      cutIntoPoints(cut,cutpoints);

      // find intersections, construct new points
      std::vector<TPoint<T>> newpoints;
      std::vector<TPoint<T>> UV;

      // cut curve specifies a correct direction of the loop, it must go first here
      int numintrs = findIntersections(cutpoints,allpoints,UV,parmtolerance); 

      // cleanup duplicate points
      if (UV.size() > 1)
      {
        std::vector<TPoint<T>> UVclean;

        for (int k = 0; k < int(UV.size()); k++)
        {
          bool found = false;
          for (int l = k + 1; l < int(UV.size()); l++)
          {
            T d0 = UV[l].X - UV[k].X;
            T d1 = UV[l].Y - UV[k].Y;
            T dist0 = std::abs(d0);
            T dist1 = std::abs(d1);

            found = (
              dist0 < parmtolerance ||
              dist1 < parmtolerance
            );

            if (found)
              break;
          }
          if (!found)
            UVclean.push_back(UV[k]);
        }

        UV = UVclean;

        numintrs = int(UV.size());
      }

      // new loop all inside single face
      if (numintrs == 0)
      {
        T dist = !(cut.front().front() - cut.back().back());
        if (dist < parmtolerance)
        {
          cut.front().front() = cut.back().back() = (cut.front().front() + cut.back().back()) * 0.5;

          // if this is a hole, we need an outer loop as loop[0]
          TPoint<T> loopnormal = getLoopNormal(cut,parmtolerance);
          bool hole = !(loopnormal > TPoint<T>(0.0,0.0,1.0));

          if (loops.empty())
            loops.push_back(outerloop);

          loops.push_back(cut);
        } else
        {
          if (!newloop)
          {
            // it maybe another loop, recreate full outer loop
            closeOuterBoundaryLoop(outerloop,numdivisions);

            newloop = true;
            goto redo;
          }
        }
      } else if (numintrs == 1)
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
        // first intersection point is from cutpoints start, cutpoints direction from
        // first to second intersection point defines cutting direction : void is to
        // the right, surface is to the left

        // all points
        int Useg0 = int(UV[0].X);
        T U0 = UV[0].X - T(Useg0);
        int Useg1 = int(UV[1].X);
        T U1 = UV[1].X - T(Useg1);

        // now go all cut points
        int Vseg0 = int(UV[0].Y);
        T V0 = UV[0].Y - T(Vseg0);
        int Vseg1 = int(UV[1].Y);
        T V1 = UV[1].Y - T(Vseg1);

        // just take all cut points void is to the right, surface is to the left
        cutOut(cutpoints,Useg0,U0,Useg1,U1,newpoints);

        if (newpoints.size() < 2) 
        {
          if (!newloop)
          {
            // it maybe another loop, recreate full outer loop
            closeOuterBoundaryLoop(outerloop,numdivisions);

            newloop = true;
            goto redo;
          }
        }

     //   std::vector<TPoint<T> corners;
     //   bool ok = goRoundBoundary<T>(newpoints.back(),newpoints.front(),corners,parmtolerance);

        // we need to proceed from newpoints last point to newpoints first
        // point; we always leave surface to the left, so there are TWO CASES here : 
        // point 1 has "greater" UV position when going along the boundary or "less"
        int n3 = 0;
        bool ok = cutOutGoingRound<T>(allpoints,Vseg1,V1,newpoints,n3,parmtolerance);
        if (n3 > 0)
        {
          ok = false;
        }

        // it must be always ok
        assert(ok);

        if (!ok)
          return false;

        // divide outer loop back into parts
        std::vector<std::vector<TPoint<T>>> newouterloop; 
        divideByDuplicates<T>(newpoints,newouterloop,parmtolerance);

        // correct points on the boundary
        for (int i = 0; i < int(newouterloop.size()); i++)
        {
          if (boundaryPoint(newouterloop[i].front(),parmtolerance))
          {
            correctBoundaryPoint(newouterloop[i].front(),parmtolerance);
          }
          if (boundaryPoint(newouterloop[i].back(),parmtolerance))
          {
            correctBoundaryPoint(newouterloop[i].back(),parmtolerance);
          }
        }

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

public: //!!!!!!!

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
                              // right-hand size and solution
    T B[4] = {0};
    bool systemsolved = true;

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

}
