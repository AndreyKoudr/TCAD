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
    ok_ = prepareParameters(this->cpoints,parms,true,parmsbynumbers); 
    len = calculateLength(this->cpoints);
  }

  /** At least two non-identical points. */
  bool ok()
  {
    return ok_;
  }

  /** Access to parameters. */
  std::vector<T> &parameters()
  {
    return parms;
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

  /** Is closed? */
  bool closed(T tolerance = TOLERANCE(T))
  {
    T dist = !(this->cpoints.back() - this->cpoints.front());
    return (dist < tolerance);
  }

  /** Close curve if not closed. */
  bool close(T tolerance = TOLERANCE(T))
  {
    if (!closed(tolerance))
    {
      this->cpoints.push_back(this->cpoints.front());
      return true;
    } else
    {
      return false;
    }
  }

  /** Shift this->cpoints to end (positive shift) or to start (negative shift). 
    It amkes sense for closed curves. */
  bool shiftClosed(int shift, T tolerance = TOLERANCE(T))
  { 
    if (closed(tolerance))
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

      // remove duplicates
      removeDuplicates(this->cpoints,false,tolerance);

      // close
      close(tolerance);

      this->update();

      return true;
    } else
    {
      return false;
    }
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

  /** Order unordered points by finding closest edges to current line ends. 
    These edges are degenerated. It is supposed there is only ONE curve
    represented by pieces, not many. 
    maxedge ratio is a guessed max edge size to the model size. */
  void order(T tolerance, T maxedgeratio = 0.1)
  {
    // order points
    std::vector<std::vector<TPoint<T>>> pieces;
    pieces.push_back(this->cpoints);
    this->cpoints.clear();

    tcad::curveFromPieces(pieces,this->cpoints,tolerance,true,maxedgeratio);

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

  /** Find point with maximum curvature (angle) between neightbour points. */
  int findPointOfMaxCurvature(int start = -1, int end = -1)
  {
    if (start < 0)
      start = 0;
    if (end < 0)
      end = int(this->cpoints.size()) - 1;

    T maxangle = 0.0;
    int index = -1;
    for (int i = start; i <= end; i++)
    {
      TPoint<T> p0,p1,p2;
      if (findThreePoints(this->cpoints,i,p0,p1,p2,TOLERANCE(T)))
      {
        T angle = ((p1 - p0) < (p2 - p1)) * PCI;

        if (angle > maxangle)
        {
          maxangle = angle;
          index = i;
        }
      }
    }

    return index;
  }

public:

  // parameterisation by numbers or by length
  bool parmsbynumbers = false;

private:
  // at least two non-identical points
  bool ok_ = false;

  // length
  T len = 0.0;

  // point parameters [0..1] parameterised by length
  std::vector<T> parms;
};

}