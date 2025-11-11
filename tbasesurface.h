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

/** Surface corner UV values. 

                        side 2      
            3------------------------------2
            |                              |
            |                              |
            |                              |
   side 3   |                              | side 1
            ^ V                            |
            |    U                         |
            0---->-------------------------1
                       side 0      

*/

template <class T> const std::array<TPoint<T>,4> cornerUV =
{
  TPoint<T>(0.0,0.0,0.0),
  TPoint<T>(1.0,0.0,0.0),
  TPoint<T>(1.0,1.0,0.0),
  TPoint<T>(0.0,1.0,0.0)
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
  template <class T> int intersect(TBaseSurface<T> &other, std::vector<std::vector<TPoint<T>>> &intersections, 
    std::vector<std::vector<TPoint<T>>> &boundary0, std::vector<std::vector<TPoint<T>>> &boundary1,
    T tolerance, T parmtolerance = TOLERANCE(T), 
    int numpointsU0 = MANY_POINTS2D, int numpointsV0 = MANY_POINTS2D,
    T refinestartU0 = 1.0, T refineendU0 = 1.0, 
    T refinestartV0 = 1.0, T refineendV0 = 1.0,
    int numpointsU1 = MANY_POINTS2D, int numpointsV1 = MANY_POINTS2D,
    T refinestartU1 = 1.0, T refineendU1 = 1.0, 
    T refinestartV1 = 1.0, T refineendV1 = 1.0)
  {
    TTriangles<T> tris,othertris;

    if (createTriangles(tris,numpointsU0,numpointsV0,
      refinestartU0,refineendU0,refinestartV0,refineendV0) &&
      other.createTriangles(othertris,numpointsU1,numpointsV1,
      refinestartU1,refineendU1,refinestartV1,refineendV1))
    {
      if (tris.intersect(othertris,intersections,tolerance,parmtolerance,&boundary0,&boundary1))
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

  /** Close UV boundary. cutUV is a cut across surface in UV coordinates. */
  template <class T> bool closeBoundaryLoop(std::vector<std::vector<TPoint<T>>> &cutUV,
    std::vector<std::vector<TPoint<T>>> &closedboundary, 
    T tolerance, T bigtolerance = 0.1, T parmtolerance = PARM_TOLERANCE, int numdivisions = 100)
  {
    if (cutUV.empty())
      return false;

    // step 1 : combine cut pieces into a single line
    std::vector<TPoint<T>> cut;
    if (cutUV.size() == 1)
    {
      // all done
      cut = cutUV[0];
    } else
    {
      // this may be a number of intersection parts, combine them into one
      // with high tolerance (some segments may be missing due to coincident
      // triangle edges)
      //T bigtolerance = tolerance;
      //for (int i = 0; i < int(cutUV.size()); i++)
      //{
      //  T min,max;
      //  if (tcad::segmentLenMinMax(cutUV[i],min,max))
      //  {
      //    bigtolerance = std::max<T>(bigtolerance,max);
      //  }
      //}
      //bigtolerance *= 1.1;
      
      std::vector<std::vector<TPoint<T>>> lines;
      if (!curvesFromPieces(cutUV,lines,bigtolerance))
        return false;

      if (lines.size() == 1)
      {
        cut = lines[0];
      } else
      {
        // not continuous cut
        return false;
      }
    }

    closedboundary.push_back(cut);

    // e.g. a hole inside face
    T closetolerance = calculateLength(closedboundary[0]) * 0.01;
    if (closed(closedboundary[0],closetolerance))
    {
      if (!closed(closedboundary[0],tolerance))
      {
        closedboundary[0].push_back(closedboundary[0].front());
      }
      return true;
    }

    // step 2 : check if both cut ends lay on one of four surface boundaries
    T parm0 = 0.0;
    T parm1 = 0.0;
    if (!UVToBoundaryParm(cut.front().X,cut.front().Y,parm0,parmtolerance))
      return false;
    if (!UVToBoundaryParm(cut.back().X,cut.back().Y,parm1,parmtolerance))
      return false;

    // cut curve already must be correctly oriented to leave the 
    // remaining surface to the left

    // step 3 : starting from the curve end go along boundary to close it
    T startparm = parm1;
    T endparm = parm0;
    T parm = startparm;
    T nextparm = parm;

    while (nextBoundaryParm(parm,endparm,nextparm,parmtolerance))
    {
      assert(std::abs(nextparm - parm) > parmtolerance);

      TPoint<T> UV = boundaryParmToUV(parm);
      TPoint<T> nextUV = boundaryParmToUV(nextparm);

      // make parametric straight line between UV and nextUV
      T d = !(nextUV - UV);

      if (d > parmtolerance)
      {
        int numdivs = int(T(numdivisions) * d);
        LIMIT_MIN(numdivs,8);
        TPointCurve<T> line(UV,nextUV,numdivs);

        closedboundary.push_back(line.controlPoints());
      }

      parm = nextparm;
    }

    T dist = !(closedboundary.front().front() - closedboundary.back().back());
    bool ok = (dist < parmtolerance);
    assert(ok);

    return ok;
  }

  /** Close outer UV boundary by 4 pieces. */
  template <class T> void closeOuterBoundaryLoop(std::vector<std::vector<TPoint<T>>> &closedboundary, 
    int numdivisions = 100)
  {
    for (int i = 0; i < 4; i++)
    {
      int i1 = i + 1;
      if (i1 > 3)
        i1 = 0;

      TPoint<T> UV = cornerUV<T>[i];
      TPoint<T> nextUV = cornerUV<T>[i1];

      TPointCurve<T> line(UV,nextUV,numdivisions);

      closedboundary.push_back(line.controlPoints());
    }
  }

public: //!!!!!!!

  // number of columns minus 1
  int K1 = 0;
  // number of rows munus 1
  int K2 = 0;

protected:

  // (K1 + 1) * (K2 + 1) control points, call update() after every change
  std::vector<TPoint<T>> cpoints;

};

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

}
