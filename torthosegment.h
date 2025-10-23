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

  torthosegment.h

  Segment on orthogonal Jacobi poly of any type and power, better up to ~20
  (gamma function grows very quickly)

*******************************************************************************/

#pragma once

#include "tbasecurve.h"
#include "tjacobipoly.h"

namespace tcad {

template <class T> class TOrthoSegment : public TBaseCurve<T> {
public:
  /** Constructor. */
  TOrthoSegment() : TBaseCurve<T>()
  {
  }

  /** Constructor. */
  TOrthoSegment(std::vector<TPoint<T>> &points, T palpha, T pbeta, int ppower, 
    int pintegration = OTHER_INTEGRATION, bool parameterisationbynumbers = false) : TBaseCurve<T>()
  {
    alpha = palpha;
    beta = pbeta;
    integration = pintegration;

    power = ppower;
    parmsbynumbers = parameterisationbynumbers;

    this->cpoints = points;

    this->update();
  }

  /** Constructor. */
  TOrthoSegment(std::vector<TPoint<T>> &points, CurveEndType start, CurveEndType end, int ppower, 
    int pintegration = OTHER_INTEGRATION, bool parameterisationbynumbers = false) : 
  TOrthoSegment(points,(end == END_ROUNDED) ? 0.5 : 0.0,(start == END_ROUNDED) ? 0.5 : 0.0,ppower, 
    pintegration,parameterisationbynumbers)
  {
  }

  /** Copy constructor. */
  TOrthoSegment(const TOrthoSegment &other)  
  {
    this->cpoints = other.cpoints;

    alpha = other.alpha;
    beta = other.beta;
    integration = other.integration;

    parmsbynumbers = other.parmsbynumbers;
    power = other.power;

    this->update();
  }

  /** Assignment operator. */
  TOrthoSegment &operator = (const TOrthoSegment &other)  
  {
    this->cpoints = other.cpoints;

    alpha = other.alpha;
    beta = other.beta;
    integration = other.integration;

    parmsbynumbers = other.parmsbynumbers;
    power = other.power;

    this->update();

    return *this;
  }

  /** Destructor. */
  virtual ~TOrthoSegment() {}

  /** Get k-th derivative on on U [0..1]. 0-derivative is poisition. */
  virtual TPoint<T> derivative(T U, int k)
  {
    assert((this->cpoints.size() >= 2) && "Wrong number of control points in LSQ segment");

    if (!_ok)
      return TPoint<T>();

    LIMIT(U,0.0,1.0);

    T u = U * 2.0 - 1.0;

    if (k == 0)
    {
      TPoint<T> point(_fx.getValue(u),_fy.getValue(u),_fz.getValue(u));
      return point;
    } else
    {
      TPoint<T> point(_fx.getDerK(u,k),_fy.getDerK(u,k),_fz.getDerK(u,k));
      return point;
    }
  }

  /** Update after any change in control points. */
  virtual void update()
  {
    // get parameters
    prepareParameters(this->cpoints,parms,true,parmsbynumbers); 

    // split into x,y,z
    std::vector<T> x,y,z;
    splitXYZ(this->cpoints,x,y,z);

    // make fittings
    _fx.alpha = _fy.alpha = _fz.alpha = alpha;
    _fx.beta = _fy.beta = _fz.beta = beta;

    _ok =
      _fx.fit(power,parms,x,integration) && 
      _fy.fit(power,parms,y,integration) &&
      _fz.fit(power,parms,z,integration);

#ifdef _DEBUG
      T ax = _fx.accuracy(parms,x);
      T ay = _fy.accuracy(parms,y);
      T az = _fz.accuracy(parms,z);
#endif
  }

  /** Approximation ok? */
  bool ok()
  {
    return _ok;
  }

  /** Approximation x(parm) */
  TJacobiPoly<T> &fx()
  {
    return _fx;
  }

  /** Approximation y(parm) */
  TJacobiPoly<T> &fy()
  {
    return _fy;
  }

  /** Approximation z(parm) */
  TJacobiPoly<T> &fz()
  {
    return _fz;
  }

private:

  // parameterisation by numbers or by length;
  // for ortho poly parameterisation by numbers is much better
  bool parmsbynumbers = true;

  // Jacobi poly type, Legendre if 0.0,0.0
  T alpha = 0.0;
  T beta = 0.0;

  // integration
  int integration = OTHER_INTEGRATION;

  // poly power
  int power = -1;

private:
  // ok?
  bool _ok = false;

  // parameters used
  std::vector<T> parms;

  // actual fittings
  TJacobiPoly<T> _fx,_fy,_fz;
};

}