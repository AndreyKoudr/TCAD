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

  tbasecurve.h

  Basic abstract class for curves

  dimensions : 1 (U parameter)

*******************************************************************************/

#pragma once

#include "tbasics.h"
#include "tpoint.h"
#include "ttransform.h"

namespace tcad {

// type of curve ends when we appriximate/interpolate a set of points representing a curve
enum CurveEndType {
  END_FREE,           // curve end left as it is after approximation (for all base segment/curve types)

  END_FIXED,          // end coordinate is coincident with the original point, working with
                      //  TBezierSegment
                      //  TBezierCurve  
                      //  TPointCurve (always holds automatically)
                      //  TSplineCurve (holds automatically)

  END_CLAMPED,        // both end coordinate and first XYZ derivatives on parameter are same as in original point,
                      // working with
                      //  TBezierCurve  
                      //  TPointCurve (always holds automatically)
                      //  TSplineCurve 

  END_ROUNDED,        // vertical tangent, working with
                      //  TOrthoSegment (only)

  END_TOTAL
};

template <class T> class TBaseCurve {
public:

  //===== Construction =========================================================

  /** Constructor. */
  TBaseCurve() {}

  /** Destructor. */
  virtual ~TBaseCurve() {}

  //===== Abstract =============================================================

  /** Get k-th derivative on on U [0..1]. 0-derivative is position. */
  virtual TPoint<T> derivative(T U, int k) = 0;

  /** Update after any change in control points. */
  virtual void update() = 0;

  //===== Operations ===========================================================

  /** Same as 0-th derivative*/
  virtual TPoint<T> position(T U)
  {
    return this->derivative(U,0);
  }

  /** Get control points. */
  virtual std::vector<TPoint<T>> &controlPoints()
  {
    return cpoints;
  }

  /** Get points along the curve to 
    (1) apply simple tpoints.h methods or
    (2) create another type of curve geometry from these points or 
    (3) compare/involve another type of curve geometry to compare/interact 

    Every point contains U curve parameter in W.
  */
  virtual void createPoints(std::vector<TPoint<T>> &points, int numpoints = MANY_POINTS)
  {
    points.clear();

    for (int i = 0; i < numpoints; i++)
    {
      T U = T(i) / T(numpoints - 1);
      TPoint<T> p = this->derivative(U,0);
      p.W = U; //!!! important for some procedures, e.g. intersectByPlane()
      points.push_back(p);
    }
  }

  /** Get curvature. */
  virtual T curvature(T U)
  {
	  TPoint<T> D1 = this->derivative(U,1);
	  TPoint<T> D2 = this->derivative(U,2);

	  T D1len = !D1;

    if (D1len < TOLERANCE(T))
    {
      return 0.0;
    } else
    {
      return !(D1 ^ D2) / (D1len * D1len * D1len);
    }
  }

  /** Reverse. */
  virtual void reverse()
  {
    std::reverse(cpoints.begin(),cpoints.end());

    // call virtual update()
    this->update();
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
  virtual bool equal(TBaseCurve &other, T tolerance, int numpoints = MANY_POINTS)
  {
    std::vector<TPoint<T>> points,otherpoints;
    createPoints(points,numpoints);
    other.createPoints(otherpoints,numpoints);

    T diff = difference(points,otherpoints);
    return (diff >= 0.0 && diff < tolerance);
  }

  /** Calculate length. */
  virtual T calculateLength(int numpoints = MANY_POINTS)
  {
    std::vector<TPoint<T>> points;
    createPoints(points,numpoints);

    return tcad::calculateLength(points);
  }

  /** Find value of parameter U for a point on (or close to) the curve. Returns -1 in failure. */
  virtual T findUforPoint(TPoint<T> p, int numpoints = MANY_POINTS)
  {
    std::vector<TPoint<T>> points;
    createPoints(points,numpoints);

    TPoint<T> proj;
    int seg = 0;
    T u = 0.0;
    if (projectPointOnPoints(points,p,proj,&seg,&u))
    {
      TPoint<T> p0 = points[seg];
      TPoint<T> p1 = points[seg + 1];
      T U = p0.W + (p1.W - p0.W) * u;
      return U;
    } else
    {
      return -1.0;
    }
  }

  /** Calculate min/max. */
  virtual bool calculateMinMax(TPoint<T> *min, TPoint<T> *max, TPoint<T> *imin = nullptr, 
    TPoint<T> *imax = nullptr, int numpoints = MANY_POINTS)
  {
    std::vector<TPoint<T>> points;
    createPoints(points,numpoints);

    return tcad::calculateMinMax(points,min,max,imin,imax);
  }

  /** Find all itersections by a plane. */
  virtual int intersectByPlane(TPlane<T> &plane, std::vector<T> &Upoints, 
    T tolerance, T parmtolerance = PARM_TOLERANCE, int numpoints = MANY_POINTS, std::vector<TPoint<T>> *intrs = nullptr)
  {
    std::vector<TPoint<T>> points;
    createPoints(points,numpoints);

    int n = tcad::intersectByPlane(points,plane,Upoints,tolerance,parmtolerance);

    if (n && intrs)
    {
      for (int i = 0; i < int(Upoints.size()); i++)
      {
        TPoint<T> p = this->derivative(Upoints[i],0);
        intrs->push_back(p);
      }
    }

    return n;
  } 

  /** Find all intersections with another curve. */
  template <class T> int findIntersections(TBaseCurve<T> &other,  std::vector<TPoint<T>> &UV, 
    T tolerance, T parmtolerance = PARM_TOLERANCE, int numpoints = MANY_POINTS)
  {
    std::vector<TPoint<T>> points,otherpoints;
    createPoints(points,numpoints);
    other.createPoints(otherpoints,numpoints);

    return tcad::findIntersections(points,otherpoints,UV,tolerance,parmtolerance);
  }

  /** Cut out a piece of curve from U0 to U1 into list of points. */
  template <class T> void cutPiece(int numpoints, T Ufrom, T Uto, std::vector<TPoint<T>> &points)
  {
    points.clear();
    T DU = (Uto - Ufrom) / T(numpoints - 1);

    for (int i = 0; i < numpoints; i++)
    {
      T U = Ufrom + T(i) * DU;
      TPoint<T> p = derivative(U,0);
      points.push_back(p);
    }
  }

protected:

  // control points, call update() after every change
  std::vector<TPoint<T>> cpoints;

};

/** Calculate min/max from control points. */
template <class T> bool calculateMinMax(std::vector<TBaseCurve<T> *> &curves,
  TPoint<T> &min, TPoint<T> &max)
{
  bool ok = false;

  for (int i = 0; i < int(curves->size()); i++)
  {
    TPoint<T> mi,ma;
    
    if (tcad::calculateMinMax(curves[i]->controlPoints(),&mi,&ma))
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
