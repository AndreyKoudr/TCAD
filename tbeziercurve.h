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

  tbeziercurve.h

  Curve of tbeziersegment-s.

*******************************************************************************/

#pragma once

#include "tbasecurve.h"
#include "tpointcurve.h"
#include "tbeziersegment.h"

namespace tcad {

template <class T> class TBezierCurve : public TBaseCurve<T> {
public:
  /** Constructor. */
  TBezierCurve() : TBaseCurve<T>()
  {
  }

  /** Constructor. */
  TBezierCurve(std::vector<TPoint<T>> &points, int numsegments, 
    bool exactstart = true, bool exactend = true, 
    bool clampedstart = false, bool clampedend = false,
    bool keepmeshrefinement = true, bool orthogonalLSQ = true) : TBaseCurve<T>()
  {
    TPointCurve<T> pointcurve(points,keepmeshrefinement);

    // number of points to build a bezier segment
    int numsegintervals = int(points.size()) / numsegments;
    LIMIT_MIN(numsegintervals,8);

    T DU = 1.0 / T(numsegments);
    T du = 1.0 / T(numsegintervals);
    for (int i = 0; i < numsegments; i++)
    {
      T U0 = DU * T(i);
      T U1 = DU * T(i + 1);

      std::vector<TPoint<T>> segpoints;
      for (int j = 0; j <= numsegintervals; j++)
      {
        T u = du * T(j);
        u = U0 + (U1 - U0) * u;
        TPoint<T> p = pointcurve.derivative(u,0);
        segpoints.push_back(p);
      }

      TBezierSegment<T> segment(segpoints,(i == 0 && exactstart),
        (i == (numsegments - 1) && exactend),orthogonalLSQ); 
      for (int j = 0; j < int(segment.controlPoints().size()); j++)
      {
        this->cpoints.push_back(segment.controlPoints()[j]);
      }
    }

    if (clampedstart)
    {
      int n0,n1,n2,n3;
      if (getSegment(0,n0,n1,n2,n3))
      {
        TPoint<T> dir = +(startDirection(points));
        T len = !(this->cpoints[n1] - this->cpoints[n0]);
        this->cpoints[n1] = this->cpoints[n0] + dir * len;
      }
    }

    if (clampedend)
    {
      int n0,n1,n2,n3;
      if (getSegment(numSegments() - 1,n0,n1,n2,n3))
      {
        TPoint<T> dir = +(endDirection(points));
        T len = !(this->cpoints[n3] - this->cpoints[n2]);
        this->cpoints[n2] = this->cpoints[n3] + dir * len;
      }
    }

    this->smooth();

    this->update();
  }

  /** Constructor. */
  TBezierCurve(std::vector<TPoint<T>> &points, int numsegments, 
    CurveEndType start, CurveEndType end,
    bool keepmeshrefinement = true, bool orthogonalLSQ = true) : 
  TBezierCurve(points,numsegments, 
    (start == END_FIXED || start == END_CLAMPED),(end == END_FIXED || end == END_CLAMPED), 
    (start == END_CLAMPED),(end == END_CLAMPED), 
    keepmeshrefinement,orthogonalLSQ)
  {
  }
  /** Copy constructor. */
  TBezierCurve(const TBezierCurve &other)  
  {
    this->cpoints = other.cpoints;

    this->update();
  }

  /** Assignment operator. */
  TBezierCurve &operator = (const TBezierCurve &other)  
  {
    this->cpoints = other.cpoints;

    this->update();

    return *this;
  }

  /** Initialise. */
  void init(std::vector<TPoint<T>> &controlpoints, std::vector<T> &parameters)
  {
    this->cpoints = controlpoints;
    this->parms = parameters;
  }

  /** Destructor. */
  virtual ~TBezierCurve() {}

  /** Get k-th derivative on on U [0..1]. 0-derivative is position. Parameterisation by
    length us used. */
  virtual TPoint<T> derivative(T U, int k)
  {
    LIMIT(U,0.0,1.0);

    TPoint<T> p;
    TBezierSegment<T> segment;
    T u = 0.0;
    if (findSegment(U,segment,u))
    {
      p = segment.derivative(u,k);
    }

    return p;
  }

  /** Update after any change in control points. */
  virtual void update()
  {
    this->parms.clear();
    len = 0.0;
    this->parms.push_back(len);

    for (int i = 0; i < numSegments(); i++)
    {
      TBezierSegment<T> segment;
      if (getSegment(i,segment))
      {
        T slen = segment.calculateLength();
        len += slen;
        this->parms.push_back(len);
      }
    }

    // normalise parameters to [0..1]
    for (T &l : parms)
    {
      l /= parms.back();
    }
  }

  /** Get number of Bezier segments in the curve. */
  int numSegments()
  {
    return int(this->cpoints.size()) / 4;
  }

  /** Find Bezier segment. */
  bool findSegment(T U, TBezierSegment<T> &segment, T &u)
  {
    int index = findParametricInterval(parms,U,&u);
    assert(index >= 0);
    if (index < 0)
      return false;

    int i4 = index * 4;
    segment.init(this->cpoints[i4],this->cpoints[i4 + 1],this->cpoints[i4 + 2],this->cpoints[i4 + 3]);

    return true;
  }

  /** i is segment number, get its 4 nodes. */
  bool getSegment(int i, int &n0, int &n1, int &n2, int &n3)
  {
    int numsegments = numSegments();
    if (numsegments < 1)
      return false;

    int i4 = i * 4;
    n0 = i4;
    n1 = i4 + 1;
    n2 = i4 + 2;
    n3 = i4 + 3;

    return true;
  }

  /** i is segment number. */
  bool getSegment(int i, TBezierSegment<T> &segment)
  {
    int n0,n1,n2,n3;
    if (getSegment(i,n0,n1,n2,n3))
    {
      TBezierSegment<T> s(this->cpoints[n0],this->cpoints[n1],this->cpoints[n2],this->cpoints[n3]);

      segment = s;

      return true;
    } else
    {
      return false;
    }
  }

  /** Smooth Bezier curve by C0 or C1 between segments using Bezier segment nodes. */
  void smooth(bool smooth0 = true, bool smooth1 = true, int index0 = -1, int index1 = -1)
  {
    int numsegments = numSegments();
    if (numsegments < 2)
      return;

    int i0 = (index0 == -1) ? 0 : index0;
    int i1 = (index1 == -1) ? numsegments - 1 : index1;
    LIMIT(i0,0,numsegments - 1);
    LIMIT(i1,0,numsegments - 1);

    if (smooth0)
    {
      for (int i = i0; i < i1; i++)
      {
        int n0,n1,n2,n3;
        int m0,m1,m2,m3;
        if (getSegment(i,n0,n1,n2,n3) && getSegment(i + 1,m0,m1,m2,m3))
        {
		      TPoint<T> v = (this->cpoints[n3] + this->cpoints[m0]) * 0.5;
          this->cpoints[n3] = this->cpoints[m0] = v;
        }
      }
    }

    if (smooth1)
    {

      for (int i = i0; i < i1; i++)
      {
        int n0,n1,n2,n3;
        int m0,m1,m2,m3;
        if (getSegment(i,n0,n1,n2,n3) && getSegment(i + 1,m0,m1,m2,m3))
        {
          straightenThreePoints(this->cpoints[n2],this->cpoints[n3],this->cpoints[m1]);
        }
      }
    }

    this->update();
  }

  /** Access to parameters. */
  std::vector<T> &parameters()
  {
    return parms;
  }

private:
  // length
  T len = 0.0;

  // segment ends parameters [0..1] parameterised by length; parameters are BETWEEN SEGMENTS,
  // every segment has its own length and parameterisation [0..1] inside
  std::vector<T> parms;
};

}