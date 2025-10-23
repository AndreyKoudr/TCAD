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

  tpointcurve.h

  Curve of many points with linear approximation between.

*******************************************************************************/

#pragma once

#include "tbasecurve.h"
#include "tpoints.h"

namespace tcad {

/** Curve of many points with linear approximation between. */
template <class T> class TPointCurve : public TBaseCurve<T> {
public:
  /** Constructor. */
  TPointCurve() : TBaseCurve<T>()
  {
  }

  /** Constructor. */
  TPointCurve(std::vector<TPoint<T>> &points, bool parametersbynumbers = false) : TBaseCurve<T>()
  {
    this->cpoints = points;
    parmsbynumbers = parametersbynumbers;

    this->update();
  }

  /** Constructor : make straight line. */
  TPointCurve(const TPoint<T> &v0, const TPoint<T> &v1, int numdivisions, bool parametersbynumbers = false) : TBaseCurve<T>()
  {
    this->cpoints.clear();
    parmsbynumbers = parametersbynumbers;

    TPoint<T> dv = (v1 - v0) / T(numdivisions);
    TPoint<T> v = v0;

    for (int i = 0; i <= numdivisions; i++)
    {
      this->cpoints.push_back(v);
      v += dv;
    }

    this->update();
  }

  /** Constructor : part of other line. */
  TPointCurve(TPointCurve &original, int i0, int i1, bool parametersbynumbers = false) : TBaseCurve<T>()
  {
    this->cpoints.clear();
    parmsbynumbers = parametersbynumbers;

    for (int i = i0; i <= i1; i++)
    {
      this->cpoints.push_back(original.controlPoints()[i]);
    }

    this->update();
  }

  /** Destructor. */
  virtual ~TPointCurve() {}

  /** Get k-th derivative on on U [0..1]. 0-derivative is poisition. Parameterisation by
    length us used. */
  virtual TPoint<T> derivative(T U, int k)
  {
    assert((this->cpoints.size() >= 2) && "Too few control points in point curve");

    LIMIT(U,0.0,1.0);

    T u = 0.0;
    int index = findParametricInterval(parms,U,&u);
    assert(index >= 0);

    if (k == 0)
    {
      TPoint<T> p0 = this->cpoints[index];
      TPoint<T> p1 = this->cpoints[index + 1];
      TPoint<T> p = p0 + (p1 - p0) * u;
      return p;
    } else if (k == 1)
    {
#if 1
      TPoint<T> der = derivative1(index);
#else
      T U0 = parms[index];
      T U1 = parms[index + 1];
      T DU = U1 - U0;

      TPoint<T> der = (std::abs(DU) > TOLERANCE(T)) ? ((p1 - p0) / DU) : 0.0;
#endif
			return der;
    } else if (k == 2)
    {
      TPoint<T> der = derivative2(index);
			return der;
    } else
    {
      return TPoint<T>();
    }
  }

  /** Update after any change in control points. */
  virtual void update()
  {
    prepareParameters(this->cpoints,parms,true,parmsbynumbers); 
    len = calculateLength(this->cpoints);
  }

  /** Redivide points. Set endrefinementcoef to < 1.0 (e.g. 0.5) to refine mesh near ends. */
  void redivide(int numPoints, std::vector<TPoint<T>> &newPoints,  
    bool treatSharpCorners = false, T sharpangledeg = 45.0, T endrefinementcoef = 1.0)
  {
    if (!this->cpoints.empty() && numPoints > 1)
    {
      // these are original sharp corners
      std::vector<TPoint<T>> corners;
      if (treatSharpCorners)
      {
        // do not exclude sharp corners
        std::vector<int> indices;
        findSharpCorners(this->cpoints,indices,TOLERANCE(T),sharpangledeg);
        // just in case
        corners.push_back(this->cpoints.front());
        corners.push_back(this->cpoints.back());
        // ...then all sharp
        for (int i = 0; i < int(indices.size()); i++)
          corners.push_back(this->cpoints[indices[i]]);
      }

      newPoints.clear();

      T DU = 1.0 / T(numPoints - 1);
      for (int i = 0; i < numPoints; i++)
      {
        T U = DU * T(i);

        T u = U * 2.0 - 1.0;
        T un = pow(std::abs(u),endrefinementcoef) * sign(u);
        U = (un + 1.0) * 0.5;

        TPoint<T> p = derivative(U,0);
        newPoints.push_back(p);
      }

      // just in case
      newPoints.front() = this->cpoints.front();
      newPoints.back() = this->cpoints.back();

      if (!corners.empty())
      {
        std::vector<bool> busy(newPoints.size(),false);
        for (int i = 0; i < int(corners.size()); i++)
        {
          T mindist = 0.0;
          int index = findClosest(newPoints,corners[i],&mindist);
          if (index >= 0 && !busy[index])
          {
            newPoints[index] = corners[i];
            busy[index] = true;
          }
        }
      }
    }
  }

  /** Redivide points. Set endrefinementcoef to < 1.0 (e.g. 0.5) to refine mesh near ends. */
  void redivide(int numPoints, bool treatSharpCorners = false, T sharpangledeg = 45.0, 
    T endrefinementcoef = 1.0)
  {
    std::vector<TPoint<T>> newPoints;
    redivide(numPoints,newPoints,treatSharpCorners,sharpangledeg,endrefinementcoef);
    this->cpoints = newPoints;

    this->update();
  }

  /** Shift this->cpoints to end (positive shift) or to start (negative shift). */
  void shift(int shift)
  { 
    for (int i = 0; i < std::abs(shift); i++)
    {
      if (shift > 0)
      {
        this->cpoints.push_back(this->cpoints.front());
        this->cpoints.erase(this->cpoints.begin());
      } else
      {
        this->cpoints.insert(this->cpoints.begin(),this->cpoints.back());
        this->cpoints.erase(this->cpoints.end() - 1);
      }
    }

    this->update();
  }

  /** Get max distance from other by projecting this this->cpoints on other (the other must be not shorter). 
    If exact, only exact projection (no distance to line ends) is considered. Returns -1 if
    such projection not found. */
  T diff(TPointCurve<T> &other, T *maxdeflection, T parmtolerance)
  {
    T maxdiff = 0.0;

    std::vector<T> diffs;
    for (int i = 0; i < int(this->cpoints.size()); i++)
    {
      TPoint<T> proj;
      bool ok = projectPointOnPoints(other.controlPoints(),this->cpoints[i],proj,nullptr,nullptr,parmtolerance);
      if (!ok)
        return -1.0;

      T diff = !(proj - this->cpoints[i]);
      diffs.push_back(diff);
      maxdiff = std::max<T>(diff,maxdiff);
    }

    if (maxdeflection)
    {
      *maxdeflection = 0.0;
      for (int i = 0; i < int(this->cpoints.size()); i++)
      {
        T ddiff = std::abs(diffs[i] - maxdiff);
        *maxdeflection = std::max<T>(ddiff,*maxdeflection);
      }
    }

    return maxdiff;
  }

  /** Get max distance from any other ends. -1.0 if failed. Return is ok if only a piece 
    of curves is coincident, not a single point if singlepointok == false.  */
  T endsDiff(TPointCurve<T>  &other, T tolerance, T parmtolerance = PARM_TOLERANCE, 
    int *countfirst = nullptr, int *countsecond = nullptr, bool singlepointok = true)
  {
    if (this->cpoints.empty() || other.controlPoints().empty())
      return -1.0;

    bool ok = false;
    T maxdiff = 0.0;
    T diff = 0.0;
    TPoint<T> proj;

    if (countfirst)
      *countfirst = 0;
    if (countsecond)
      *countsecond = 0;

    std::vector<TPoint<T>> projections;

    bool ok0 = projectPointOnPointsExact(other.controlPoints(),this->cpoints.front(),proj,nullptr,nullptr,parmtolerance);
    if (ok0)
    {
      diff = !(proj - this->cpoints.front());
      maxdiff = std::max<T>(diff,maxdiff);
      ok = true;
      if (countfirst)
        *countfirst++;

      addToProjections(projections,proj,tolerance);
    }

    bool ok1 = projectPointOnPoints(other.controlPoints(),this->cpoints.back(),proj,nullptr,nullptr,parmtolerance);
    if (ok1)
    {
      diff = !(proj - this->cpoints.back());
      maxdiff = std::max<T>(diff,maxdiff);
      ok = true;
      if (countfirst)
        *countfirst++;

      addToProjections(projections,proj,tolerance);
    }

    bool ok2 = projectPointOnPoints(this->cpoints,other.controlPoints().front(),proj,nullptr,nullptr,parmtolerance);
    if (ok2)
    {
      diff = !(proj - other.controlPoints().front());
      maxdiff = std::max<T>(diff,maxdiff);
      ok = true;
      if (countsecond)
        *countsecond++;

      addToProjections(projections,proj,tolerance);
    }

    bool ok3 = projectPointOnPoints(this->cpoints,other.controlPoints().back(),proj,nullptr,nullptr,parmtolerance);
    if (ok3)
    {
      diff = !(proj - other.controlPoints().back());
      maxdiff = std::max<T>(diff,maxdiff);
      ok = true;
      if (countsecond)
        *countsecond++;

      addToProjections(projections,proj,tolerance);
    }

    // not a single corner touch
    return (ok && (singlepointok || projections.size() > 1)) ? maxdiff : -1.0;
  }

  /** Unclose line : a fix after order(). */
  void unclose()
  {
    // find longest segment imax
    T min,dist1;
    int imin,imax;
    if (tcad::segmentLenMinMax(this->cpoints,min,dist1,&imin,&imax))
    {
      // distance between start and end
      T dist0 = !(this->cpoints.front() - this->cpoints.back());

      if (dist1 > dist0)
      {
        shift(imax + 1);
      }
    }
  }

  /** Order unordered points by finding closest edges to current line ends. 
    These edges are degenerated. */
  void order(T tolerance)
  {
    // order points
    std::vector<TPoint<T>> pieces = this->cpoints;
    this->cpoints.clear();

    tcad::curveFromPieces(pieces,this->cpoints,tolerance,true);

    // fix 
    unclose();

    this->update();
  }

  /** First derivative on parameter at node index. */
  TPoint<T> derivative1(int index)
  {
    int i0 = index - 1;
    int i1 = index + 1;
    if (i0 < 0)
    {
      i0++;
      i1++;
    }
    if (i1 > int(this->cpoints.size()) - 1)
    {
      i0--;
      i1--;
    }

    if (i0 >= 0 && i0 < int(this->cpoints.size()) &&
      i1 >= 0 && i1 < int(this->cpoints.size()))
    {
      TPoint<T> p0 = this->cpoints[i0];
      TPoint<T> p1 = this->cpoints[i1];
      TPoint<T> d = (p1 - p0) / (parms[i1] - parms[i0]);
      return d;
    } else
    {
      assert(false);
      return TPoint<T>();
    }
  }

  /** Second derivative on parameter at node index. */
  TPoint<T> derivative2(int index)
  {
    int i0 = index - 1;
    int i1 = index + 1;
    if (i0 < 0)
    {
      i0++;
      i1++;
    }
    if (i1 > int(this->cpoints.size()) - 1)
    {
      i0--;
      i1--;
    }

    if (i0 >= 0 && i0 < int(this->cpoints.size()) &&
      i1 >= 0 && i1 < int(this->cpoints.size()))
    {
      TPoint<T> d0 = derivative1(i0);
      TPoint<T> d1 = derivative1(i1);
      TPoint<T> d = (d1 - d0) / (parms[i1] - parms[i0]);
      return d;
    } else
    {
      assert(false);
      return TPoint<T>();
    }
  }


public:

  // parameterisation by numbers or by length
  bool parmsbynumbers = false;

private:
  // length
  T len = 0.0;

//!!!!!!!
public:
  // point parameters [0..1] parameterised by length
  std::vector<T> parms;
};

}