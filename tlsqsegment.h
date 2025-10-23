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

  tlsqsegment.h

  LSQ segment of low power up to 3

*******************************************************************************/

#pragma once

#include "tbasecurve.h"
#include "tlsqfitting.h"

namespace tcad {

/** LSQ segment of power up to TLSQFitting::LSFITTING_MAXPOWER - 1. */
template <class T> class TLSQSegment : public TBaseCurve<T> {
public:

  /** Constructor. */
  TLSQSegment() : TBaseCurve<T>()
  {
  }

  /** Constructor. */
  TLSQSegment(std::vector<TPoint<T>> &points, int ppower, bool parametersbynumbers = false) : TBaseCurve<T>()
  {
    power = ppower;
    LIMIT(power,0,TLSQFitting<T>::LSFITTING_MAXPOWER - 1);
    this->cpoints = points;

    parmsbynumbers = parametersbynumbers;

    this->update();
  }

  /** Copy constructor. */
  TLSQSegment(const TLSQSegment &other)  
  {
    this->cpoints = other.cpoints;

    parmsbynumbers = other.parmsbynumbers;
    power = other.power;

    this->update();
  }

  /** Assignment operator. */
  TLSQSegment &operator = (const TLSQSegment &other)  
  {
    this->cpoints = other.cpoints;

    parmsbynumbers = other.parmsbynumbers;
    power = other.power;

    this->update();

    return *this;
  }

  /** Destructor. */
  virtual ~TLSQSegment() {}

  /** Get k-th derivative on on U [0..1]. 0-derivative is poisition. */
  virtual TPoint<T> derivative(T U, int k)
  {
    assert((this->cpoints.size() >= 2) && "Wrong number of control points in LSQ segment");

    if (!_ok)
      return TPoint<T>();

    LIMIT(U,0.0,1.0);

    if (k == 0)
    {
      TPoint<T> point(_fx.GetValue(U),_fy.GetValue(U),_fz.GetValue(U));
      return point;
    } else
    {
      // no derivatives in LSQ
      return TPoint<T>();
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
    _ok =
      _fx.fit(power,int(x.size()),&parms[0],&x[0]) &&
      _fy.fit(power,int(y.size()),&parms[0],&y[0]) &&
      _fz.fit(power,int(z.size()),&parms[0],&z[0]);
  }

  /** Approximation ok? */
  bool ok()
  {
    return _ok;
  }

  /** Approximation x(parm) */
  TLSQFitting<T> &fx()
  {
    return _fx;
  }

  /** Approximation y(parm) */
  TLSQFitting<T> &fy()
  {
    return _fy;
  }

  /** Approximation z(parm) */
  TLSQFitting<T> &fz()
  {
    return _fz;
  }

private:
  // parameterisation by numbers or by length
  bool parmsbynumbers = false;

  // poly power, max 3
  int power = -1;

private:
  // ok?
  bool _ok = false;

  // actual fittings
  TLSQFitting<T> _fx,_fy,_fz;

  // parameters used
  std::vector<T> parms;
};

}