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

  tsplinecurve.h

  B-spline curve.

*******************************************************************************/

#pragma once

#include "tbasecurve.h"
#include "tmisc.h"
#include "tpoints.h"
#include "tvalues.h"
#include "tpointcurve.h"
#include "tbeziersegment.h"

namespace tcad {

/** Curve of many points with linear approximation between. */
template <class T> class TSplineCurve : public TBaseCurve<T> {
public:
  /** Constructor. */
  TSplineCurve() : TBaseCurve<T>()
  {
  }

  /** Copy constructor. */
  TSplineCurve(const TSplineCurve &other)  
  {
    points = other.points;
    clampedstart = other.clampedstart;
    clampedend = other.clampedend;
    allocate(other.K1,other.M1,other.interpolate);

    this->update();
  }

  /** Assignment operator. */
  TSplineCurve &operator = (const TSplineCurve &other)  
  {
    points = other.points;
    clampedstart = other.clampedstart;
    clampedend = other.clampedend;
    allocate(other.K1,other.M1,other.interpolate);

    this->update();

    return *this;
  }

  /** Constructor for point approximation. */
  TSplineCurve(std::vector<TPoint<T>> &ppoints, int k1, int m1, 
    bool pclampedstart = true, bool pclampedend = true) : TBaseCurve<T>()
  {
    // set parameters
    allocate(k1,m1,false);

    // points are not cpoints 
    points = ppoints;
    clampedstart = pclampedstart;
    clampedend = pclampedend;

    this->update();
  }

  /** Constructor for point interpolation. */
  TSplineCurve(std::vector<TPoint<T>> &ppoints, int m1,
    bool pclampedstart = true, bool pclampedend = true) : TBaseCurve<T>()
  {
    // set parameters
    allocate(int(ppoints.size()) - 1,m1,true);

    // points are not cpoints 
    points = ppoints;
    clampedstart = pclampedstart;
    clampedend = pclampedend;

    this->update();
  }

  /** Constructor for point approximation. */
  TSplineCurve(std::vector<TPoint<T>> &ppoints, int k1, int m1, 
    CurveEndType start, CurveEndType end) : 
  TSplineCurve(ppoints,k1,m1,(start == END_CLAMPED),(end == END_CLAMPED))
  {
  }

  /** Constructor for point interpolation. */
  TSplineCurve(std::vector<TPoint<T>> &ppoints, int m1, 
    CurveEndType start, CurveEndType end) : 
  TSplineCurve(ppoints,m1,(start == END_CLAMPED),(end == END_CLAMPED))
  {
  }

  /** Destructor. */
  virtual ~TSplineCurve() 
  {
    DELETE_CLASS(Uderivative);
    DELETE_CLASS(UUderivative);
  }

  /** Get k-th derivative on on U [0..1]. 0-derivative is position. */
  virtual TPoint<T> derivative(T U, int k)
  {
    assert(!this->cpoints.empty());

    LIMIT(U,0.0,1.0);

    if (k == 0)
    {
      TPoint<T> result;

      std::vector<T> bs(K1 + M1 + 2,0.0);

      splineBasis(M1 + 1,U,Uknots,bs);

      for (int i = 0; i <= K1; i++)
      {
        result = result + this->cpoints[i] * bs[i];
      }

      return result;
    } else if (k == 1)
    {
      if (Uderivative)
      {
        TPoint<T> result = Uderivative->derivative(U,0);
        return result;
      } else
      {
        return TPoint<T>();
      }
    } else if (k == 2)
    {
      if (UUderivative)
      {
        TPoint<T> result = UUderivative->derivative(U,0);
        return result;
      } else
      {
        return TPoint<T>();
      }
    } else
    {
      return TPoint<T>();
    }
  }

  /** Update after any change in control points. */
  virtual void update()
  {
    // number of points here must equal K1 + 1
    std::vector<TPoint<T>> kpoints;
    if (points.size() == K1 + 1)
    {
      kpoints = points;
    } else
    {
      TPointCurve<T> pcurve(points);
      pcurve.redivide(K1 + 1,false,45.0,1.0); 
      kpoints = pcurve.controlPoints();
    }

    assert(kpoints.size() == K1 + 1);

    makeKnots();

    if (interpolate)
    {
      interpolatePoints(kpoints);
    } else
    {
      approximatePoints(kpoints);
    }

    if (clampedstart)
    {
      setClampedStart(points);
    }

    if (clampedend)
    {
      setClampedEnd(points);
    }

    // prepare two derivatives
    if (Uderivative == nullptr)
      Uderivative = new TSplineCurve();
    if (UUderivative == nullptr)
      UUderivative = new TSplineCurve();

    makeUDerivative(*Uderivative);
    Uderivative->makeUDerivative(*UUderivative);
  }

  /** Set clamped direction at the start. */
  void setStartDerivative(TPoint<T> der)
  {
    setClampedStart(der);
  }

  /** Set clamped direction at the end. */
  void setEndDerivative(TPoint<T> der)
  {
    setClampedEnd(der);
  }
  /** Set first derivative at curve starting end from points. It moves control point 1. 
    https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node17.html */
  void setClampedStart(TPoint<T> dir)
  {
    // k 
    int k = M1 + 1;
    T k1 = T(k - 1);
    T tk1 = Uknots[k] - Uknots[1];
    TPoint<T> p1p0 = dir * (tk1 / k1);
    this->cpoints[1] = this->cpoints[0] + p1p0;
  }

  /** Set first derivative at curve starting end from points. It moves control point 1. 
    https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node17.html */
  void setClampedStart(std::vector<TPoint<T>> &points)
  {
    // direction at the start
    TPoint<T> dir = startDirection(points);

    setClampedStart(dir);
  }

  /** Set first derivative at curve end from points. It moves control point cpoints.size() - 2. */
  void setClampedEnd(TPoint<T> dir)
  {
    int n0 = int(this->cpoints.size()) - 1;
    int n1 = n0 - 1;

    // k 
    int k = M1 + 1;
    T k1 = T(k - 1);
    T tk1 = Uknots[K1 + k - 1] - Uknots[K1];
    TPoint<T> pnpn1 = dir * (tk1 / k1); 
    this->cpoints[n1] = this->cpoints[n0] + pnpn1;
  }

  /** Set first derivative at curve end from points. It moves control point cpoints.size() - 2. */
  void setClampedEnd(std::vector<TPoint<T>> &points)
  {
    // direction at the end
    TPoint<T> dir = endDirection(points);

    setClampedEnd(dir);
  }

protected:

  /** Set main parameters. */
  void allocate(int k1, int m1, bool interpolation)
  {
    K1 = k1;
    M1 = m1;
    interpolate = interpolation;
  }

  /** Make uniform knots. */
  void makeKnots()
  {
    tcad::makeKnots(K1,M1,Uknots);
  }

  /** Calculate control points for points to approximate them. Reliable but the 
    curve is approximate. The curve ends always coincide with the specified points
    and first derivatives at ends are exactly defined by directions between pairs
    of original points at ends. */
  bool approximatePoints(std::vector<TPoint<T>> &points)
  {
    this->cpoints.clear();
    this->cpoints.resize(K1 + 1,TPoint<T>());

    // K1 + M1 + 2 = A + 1 - #knots
    std::vector<T> bs(K1 + M1 + 2,0.0);

    for (int i = 0; i <= K1; i++)
    {
      T U = T(i) / T(K1);
      LIMIT(U,0.0,1.0);

      splineBasis(M1 + 1,U,Uknots,bs);

      for (int j = 0; j <= K1; j++)
      {
        this->cpoints[i] = this->cpoints[i] + points[j] * bs[j];
      }
    }

    return true;
  }

  /** Calculate control points for points to interpolate them. The curve passes through 
    all specified points. Requires a solution to a system of equations which is not 
    always successful. In such cases, approximation is used with false returned. */
  bool interpolatePoints(std::vector<TPoint<T>> &points, T parmtolerance = PARM_TOLERANCE, 
    int maxinterpolationpoints = MAX_SPLINEINTERPOLATIONPOINTS)
  {
    if (points.size() > maxinterpolationpoints)
    {
      approximatePoints(points);
      return false;
    }

    //this->cpoints.resize(K1 + 1,TPoint<T>());

    this->cpoints = points;


    std::vector<T> bs(K1 + M1 + 2,0.0);

#if 1
    // calculate bandwidth, just to be sure, looks like it must be (M1 - 1) * 2 + 1
    int halfbandwidth = 0;

    for (int i = 0; i <= K1; i++)
    {
      T U = T(i) / T(K1);
      LIMIT(U,0.0,1.0);

      splineBasis(M1 + 1,U,Uknots,bs); 

      for (int ii = 0; ii <= K1; ii++)
      {
        if (std::abs(bs[ii]) > TOLERANCE(T))
        {
          int d = std::abs(i - ii);
          halfbandwidth = std::max<int>(halfbandwidth,d);
        }
      }
    }

    int bandwidth = halfbandwidth * 2 + 1;

    // banded non-symmetric matrix
    int C = K1 + 1;
    BandedMatrixSimple<T,TPoint<T>> Ab(C,bandwidth);
//!!!    BandedMatrixSimple<T,TPoint<T>> Ab(C,(M1 - 1) * 2 + 1);

    for (int i = 0; i <= K1; i++)
    {
      T U = T(i) / T(K1);
      LIMIT(U,0.0,1.0);

      splineBasis(M1 + 1,U,Uknots,bs); 

      for (int ii = 0; ii <= K1; ii++)
      {
        if (std::abs(bs[ii]) > TOLERANCE(T))
          Ab.element(i,ii) += bs[ii];
      }
    }

    T tolerance = TOLERANCE(T);

    // solve system
    bool ok = Ab.solveSystem(&this->cpoints[0],tolerance);
#else
    // matrix
    std::vector<T> A;
    int C = K1 + 1;
    A.resize(C * C,0.0);

    for (int i = 0; i <= K1; i++)
    {
      T U = T(i) / T(K1);
      LIMIT(U,0.0,1.0);

      splineBasis(M1 + 1,U,Uknots,bs); 

      for (int ii = 0; ii <= K1; ii++)
      {
        A[i * C + ii] += bs[ii]; 
      }
    }

    T tolerance = TOLERANCE(T);

    // solve system
    bool ok = solveSystemVec<T,TPoint<T>>(C,&A[0],&this->cpoints[0],tolerance);
#endif

    if (!ok)
    {
      approximatePoints(points);
      return false;
    }

    // calculate residuals,
    // it is still not clear why solution becomes rubbish with N > ~1000,
    std::vector<T> residuals(K1 + 1,0.0);
    std::vector<TPoint<T>> vresiduals(K1 + 1,TPoint<T>());

    for (int i = 0; i <= K1; i++)
    {
      T U = T(i) / T(K1);
      LIMIT(U,0.0,1.0);

      std::vector<T> equation(K1 + 1,0.0);

      splineBasis(M1 + 1,U,Uknots,bs); 

      for (int ii = 0; ii <= K1; ii++)
      {
        equation[ii] += bs[ii]; 
      }

      TPoint<T> sum;
      for (int ii = 0; ii <= K1; ii++)
      {
        sum += this->cpoints[ii] * equation[ii];
      }

      vresiduals[i] = (points[i] - sum);
      residuals[i] = (!(points[i] - sum));
    }

    // max value of residual
    std::pair<T,T> minmax = calculateMinMax(residuals);

    T len = calculateLength(points);
    T ltolerance = len * parmtolerance;

    // check number 1
    if (minmax.second > ltolerance)
    {
      approximatePoints(points);
      ok = false;
    }

    // check number 2
    T dtolerance = len;
    T diff = difference(points,this->cpoints);
    if (diff > dtolerance)
    {
      approximatePoints(points);
      ok = false;
    }

    return ok;
  }

  /** Make another spline curve to calculate derivatives from the current. */
  void makeUDerivative(TSplineCurve &Uderivative)
  {
                                // create curve 1 order less
    Uderivative.allocate(K1 - 1,M1 - 1,interpolate);
    Uderivative.makeKnots();
                                // prepare control points
    Uderivative.cpoints.resize(Uderivative.K1 + 1,TPoint<T>());
    T m1 = T(M1);
    for (int i = 0; i <= Uderivative.K1; i++)
    {
      T du = Uknots[i + M1 + 1] - Uknots[i + 1];
      Uderivative.cpoints[i] = (this->cpoints[i + 1] - this->cpoints[i]) * (m1 / du); 
    }

    // DO NOT call update() here, cpoints are ready
  }

public:

  // spline K1 : number of intervals
  int K1 = 0;

  // degree of basis functions, if == 2, it is qubic
  int M1 = 0;

  // knots
  std::vector<T> Uknots;

private:
  // interpolate or approximate
  bool interpolate = false; 

  // original points
  std::vector<TPoint<T>> points;

  // clamped (first derivative specified from points) or 
  // natural (second derivative zero) ends
  bool clampedstart = false;
  bool clampedend = false;

  // initialised in update(); they are underinitialised, can be called for 
  // derivative(U,0) but not ,1) or ,2).
  TSplineCurve<T> *Uderivative = nullptr;
  TSplineCurve<T> *UUderivative = nullptr;
};

/** Make a spline curve control points from points. */
template <class T> void makeSplineCurve(std::vector<TPoint<T>> &points, int K, int M, std::vector<TPoint<T>> &cpoints, 
  CurveEndType start = END_FREE, CurveEndType end = END_FREE, bool interpolate = false)
{
  if (interpolate)
  {
    TSplineCurve<T> curve(points,K,M,start,end);
    assert(points.size() == curve.controlPoints().size());
    cpoints = curve.controlPoints();
  } else
  {
    TSplineCurve<T> curve(points,M,start,end);
    assert(points.size() == curve.controlPoints().size());
    cpoints = curve.controlPoints();
  }
}

}