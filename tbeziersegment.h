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

  tbeziersegment.h

  Bezier segment

*******************************************************************************/

#pragma once

#include "tbasecurve.h"
#include "torthosegment.h"
#include "tlsqsegment.h"

namespace tcad {

/** Bezier segment, this->cpoints must contain 4 points. */
template <class T> class TBezierSegment : public TBaseCurve<T> {
public:
  /** Constructor. */
  TBezierSegment() : TBaseCurve<T>()
  {
  }

  /** Copy constructor. */
  TBezierSegment(const TBezierSegment &other)  
  {
    this->cpoints = other.cpoints;

    this->update();
  }

  /** Assignment operator. */
  TBezierSegment &operator = (const TBezierSegment &other)  
  {
    this->cpoints = other.cpoints;

    this->update();

    return *this;
  }

  /** Constructor. */
  TBezierSegment(TPoint<T> p0, TPoint<T> p1, TPoint<T> p2, TPoint<T> p3) : TBaseCurve<T>()
  {
    init(p0,p1,p2,p3);
  }

  /** Constructor. */
  TBezierSegment(TPoint<T> p0, TPoint<T> p3) : TBaseCurve<T>()
  {
    this->cpoints.clear();
    this->cpoints.push_back(p0);
    this->cpoints.push_back(p0 + (p3 - p0) / 3.0);
    this->cpoints.push_back(p0 + (p3 - p0) * (2.0 / 3.0));
    this->cpoints.push_back(p3);
  }

  /** LSQ constructor over points. */
  TBezierSegment(std::vector<TPoint<T>> &points, bool exactleft, bool exactright, 
    bool orthogonalLSQ = true) : TBaseCurve<T>()
  {
    this->cpoints.clear();

    if (points.empty())
      return;

    if (points.size() == 1)
    {
      // just one point in input
      for (int i = 0; i < 4; i++)
        this->cpoints.push_back(points[0]);
    } else if (points.size() == 2)
    {
      // two input points - male a linear segment
      this->cpoints.push_back(points[0]);
      this->cpoints.push_back(points[0] + (points[1] - points[0]) / 3.0);
      this->cpoints.push_back(points[0] + (points[1] - points[0]) * (2.0 / 3.0));
      this->cpoints.push_back(points[1]);
    } else
    {
      // zero length - make degenerated segment
      T length = calculateLength(points);
      if (length < TOLERANCE(T))
      {
        for (int i = 0; i < 4; i++)
          this->cpoints.push_back(points[0]);
      } else
      {
        // this is normal : make LSQ over these points
        if (orthogonalLSQ)
        {
          // if fitting power == 4, new segment will be exact copy of Legendre poly
          TOrthoSegment<T> lsq(points,0.0,0.0,4,GAUSSINT_8);

          if (lsq.ok())
          {
            this->cpoints.resize(4);

            this->cpoints[0] = TPoint<T>(lsq.fx().getValue(-1.0),lsq.fy().getValue(-1.0),lsq.fz().getValue(-1.0));
            this->cpoints[3] = TPoint<T>(lsq.fx().getValue(1.0),lsq.fy().getValue(1.0),lsq.fz().getValue(1.0));

            TPoint<T> der0 = TPoint<T>(lsq.fx().getDerK(-1.0,1),lsq.fy().getDerK(-1.0,1),lsq.fz().getDerK(-1.0,1));
            TPoint<T> der1 = TPoint<T>(lsq.fx().getDerK(1.0,1),lsq.fy().getDerK(1.0,1),lsq.fz().getDerK(1.0,1));

            T len = !(this->cpoints[3] - this->cpoints[0]);

            this->cpoints[1] = this->cpoints[0] + (+der0) * 0.3333333333333333 * len;
            this->cpoints[2] = this->cpoints[3] - (+der1) * 0.3333333333333333 * len;

            // exact starting/ending points
            if (exactleft)
            {
              this->cpoints[0] = points.front();
            }
            if (exactright)
            {
              this->cpoints[3] = points.back();
            }
          } else
          {
            assert(false && "TBezierSegment Jacobi LSQ constructor error.");
          }
        } else
        {
          // standard poly power 3, 4 coefs
          TLSQSegment<T> lsq(points,3);

          if (lsq.ok())
          {
            this->cpoints.resize(4);

            // temporary
            TPoint<T> nodes[4];
            
            nodes[0].X = lsq.fx().coefs[0];
            nodes[1].X = lsq.fx().coefs[1] * 0.3333333333333333 + nodes[0].X;
            nodes[2].X = lsq.fx().coefs[2] * 0.3333333333333333 + 2.0 * nodes[1].X - nodes[0].X;
            nodes[3].X = lsq.fx().coefs[3] + 3.0 * nodes[2].X - 3.0 * nodes[1].X + nodes[0].X;

            nodes[0].Y = lsq.fy().coefs[0];
            nodes[1].Y = lsq.fy().coefs[1] * 0.3333333333333333 + nodes[0].Y;
            nodes[2].Y = lsq.fy().coefs[2] * 0.3333333333333333 + 2.0 * nodes[1].Y - nodes[0].Y;
            nodes[3].Y = lsq.fy().coefs[3] + 3.0 * nodes[2].Y - 3.0 * nodes[1].Y + nodes[0].Y;

            nodes[0].Z = lsq.fz().coefs[0];
            nodes[1].Z = lsq.fz().coefs[1] * 0.3333333333333333 + nodes[0].Z;
            nodes[2].Z = lsq.fz().coefs[2] * 0.3333333333333333 + 2.0 * nodes[1].Z - nodes[0].Z;
            nodes[3].Z = lsq.fz().coefs[3] + 3.0 * nodes[2].Z - 3.0 * nodes[1].Z + nodes[0].Z;

            TPoint<T> offset = points.front() - nodes[0];

            nodes[0] = exactleft ? (points.front()) : (nodes[0] + offset);
            nodes[1] = nodes[1] + offset;
            nodes[2] = nodes[2] + offset;
            nodes[3] = exactright ? (points.back()) : (nodes[3] + offset);

            this->cpoints[0] = nodes[0];
            this->cpoints[1] = nodes[1];
            this->cpoints[2] = nodes[2];
            this->cpoints[3] = nodes[3];
          } else
          {
            assert(false && "TBezierSegment LSQ constructor error.");
          }
        }
      }
    }
  }

  /** LSQ constructor over points. */
  TBezierSegment(std::vector<TPoint<T>> &points, CurveEndType start, CurveEndType end, 
    bool orthogonalLSQ = true) : 
  TBezierSegment(points,(start == END_FIXED),(end == END_FIXED),orthogonalLSQ)
  {
  }

  /** Constructor : cut part from another segment. */
  TBezierSegment(TBezierSegment *copy, T U0, T U1) : TBaseCurve<T>()
  {
    this->cpoints.resize(4);

    LIMIT(U0,0.0,1.0);
    LIMIT(U1,0.0,1.0);

    T U01 = 1.0 - U0;
    T U11 = 1.0 - U1;

    TPoint<T> v[4];
    int i;
    for (i = 0; i < 4; i++) 
      v[i] = copy->controlPoints()[i];

    for (i = 0; i < 4; i++)
    {
      if (i < 3)
      {
        this->cpoints[0].XYZ[i] =
          U01*U01*U01*  v[0].XYZ[i] +
          3.0*U0*U01*U01* v[1].XYZ[i] +
          3.0*U0*U0*U01*  v[2].XYZ[i] +
          U0*U0*U0*     v[3].XYZ[i];

        this->cpoints[1].XYZ[i] =
          U01*U01*U11*           v[0].XYZ[i] +
          U01*(2.0*U0 + U1 - 3.0*U0*U1)* v[1].XYZ[i] +
          U0*(2.0*U1 + U0 - 3.0*U0*U1)*  v[2].XYZ[i] +
          U0*U0*U1*              v[3].XYZ[i];

        this->cpoints[2].XYZ[i] =
          U01*U11*U11*           v[0].XYZ[i] +
          U11*(2.0*U1 + U0 - 3.0*U0*U1)* v[1].XYZ[i] +
          U1*(2.0*U0 + U1 - 3.0*U0*U1)*  v[2].XYZ[i] +
          U1*U1*U0*              v[3].XYZ[i];

        this->cpoints[3].XYZ[i] =
          U11*U11*U11*  v[0].XYZ[i] +
          3.0*U1*U11*U11* v[1].XYZ[i] +
          3.0*U1*U1*U11*  v[2].XYZ[i] +
          U1*U1*U1*     v[3].XYZ[i];
      } else
      {
        this->cpoints[0].XYZ[i] = 1.0;
        this->cpoints[1].XYZ[i] = 1.0;
        this->cpoints[2].XYZ[i] = 1.0;
        this->cpoints[3].XYZ[i] = 1.0;
      }
    }
  }

  /** Constructor : cut part from two points and two directions. len is segment length, 
    dirs are normalized. */
  TBezierSegment(TPoint<T> p0, TPoint<T> p1, TPoint<T> dir0, TPoint<T> dir1, T len) :
    TBaseCurve()
  {
    this->cpoints.resize(4);

    if (len < 0.0)
    {
      len = !(p0 - p1);
    }

    this->cpoints[0] = p0;
    this->cpoints[3] = p1;
    this->cpoints[1] = dir0 * (len * 0.3333333333333333) + p0;
    this->cpoints[2] = p1 - dir1 * (len * 0.3333333333333333);
  }

  /** Create a round segment; dir0, dir1 may not be normalised; they are adjusted and radius returned. */
  TBezierSegment(TPoint<T> p0, TPoint<T> p1, TPoint<T> dir0, TPoint<T> dir1,
    T *R) : TBaseCurve<T>()
  {
    this->cpoints.resize(4);

    // make right derivatives
    dir0 = +dir0;
    dir1 = +dir1;
    T alpha = (dir0 < dir1) * 0.5;
    TPoint<T> dp = p1 - p0;
    T d0 = (!dp) * 0.5;

    T sa = sin(alpha);
    if (std::abs(sa) < TOLERANCE(T))
    {
      *R = 0.0; // this is actually infinity
      this->cpoints[0] = p0;
      this->cpoints[3] = p1;
      TPoint<T> d = p1 - p0;
      this->cpoints[1] = p0 + d * 0.3333333333333333;
      this->cpoints[2] = p1 - d * 0.3333333333333333;
    } else
    {
      *R = d0 / sa;
      T len = 2.0 * (*R) * alpha;
      dir0 = dir0 * len;
      dir1 = dir1 * len;

      this->cpoints[0] = p0;
      this->cpoints[3] = p1;
      this->cpoints[1] = dir0 * 0.3333333333333333 + p0;
      this->cpoints[2] = p1 - dir1 * 0.3333333333333333;
    }
  }

  /** Create ellipse from -90 to 0 (lower == true) or from 0 to 90 with X axis directed along Xdir. */
  TBezierSegment(TPoint<T> centre, T a, T b, bool lower, TPoint<T> Xdir) : TBaseCurve<T>()
  {
    this->cpoints.resize(4);

    T turnangle = Xdir.GetAngleBetween0360(TPoint<T>(1,0,0));

    if (lower)
    {
      this->cpoints[0] = TPoint<T>(0,-b);
      this->cpoints[1] = TPoint<T>(a * 0.55,-b);
      this->cpoints[2] = TPoint<T>(a,-b * 0.55);
      this->cpoints[3] = TPoint<T>(a,0);
    } else
    {
      this->cpoints[0] = TPoint<T>(a,0);
      this->cpoints[1] = TPoint<T>(a,b * 0.55);
      this->cpoints[2] = TPoint<T>(a * 0.55,b);
      this->cpoints[3] = TPoint<T>(0,b);
    }

    for (int i = 0; i < 4; i++)
    {
      Turn(turnangle * CPI,0.0,0.0,this->cpoints[i].X,this->cpoints[i].Y,TOLERANCE(T));
      this->cpoints[i].X += centre.X;
      this->cpoints[i].Y += centre.Y;
    }
  }

  /** Initialise. */
  void init(TPoint<T> p0, TPoint<T> p1, TPoint<T> p2, TPoint<T> p3)
  {
    this->cpoints.clear();
    this->cpoints.push_back(p0);
    this->cpoints.push_back(p1);
    this->cpoints.push_back(p2);
    this->cpoints.push_back(p3);
  }

  /** Destructor. */
  virtual ~TBezierSegment() {}

  /** Get k-th derivative on on U [0..1]. 0-derivative is poisition. */
  virtual TPoint<T> derivative(T U, int k)
  {
    assert((this->cpoints.size() == 4) && "Wrong number of control points in Bezier segment");

    LIMIT(U,0.0,1.0);

    if (k == 0)
    {
      if (U == 0.0)
      {
        return this->cpoints[0];
      } else if (U == 1.0)
      {
        return this->cpoints[3];
      } else
      {
		    TPoint<T> coefs = coefs0(U);
        return multiplyCoefs(coefs);
      }
    } else if (k == 1)
    {
		  TPoint<T> coefs = coefs1(U);
      return multiplyCoefs(coefs);
    } else if (k == 2)
    {
		  TPoint<T> coefs = coefs2(U);
      return multiplyCoefs(coefs);
    } else
    {
      return TPoint<T>();
    }
  }

  /** Update after any change in control points. */
  virtual void update()
  {
  }

private:

  /** Calculate coefs for position vector. */
  TPoint<T> coefs0(T U)
  {
    LIMIT(U,0.0,1.0);
    T U1 = 1.0 - U;

    TPoint<T> coefs;

    coefs.X = U1 * U1 * U1;
    coefs.Y = 3.0 * U * U1 * U1;
    coefs.Z = 3.0 * U * U * U1;
    coefs.W = U * U * U;

    return coefs;
  }

  /** Calculate coefs for first derivative. */
  TPoint<T> coefs1(T U)
  {
    LIMIT(U,0.0,1.0);

    TPoint<T> coefs;

    if (U == 0.0)
    {
      coefs = TPoint<T>(-3.0,3.0,0.0,0.0);
    } else if (U == 1.0)
    {
      coefs = TPoint<T>(0.0,0.0,-3.0,3.0);
    } else
    {
		  T U1 = 1.0 - U; 
		  T U2 = U * U;

		  coefs.X = -3.0 * U1 * U1;
		  coefs.Y = 3.0 - 12.0 * U + 9.0 * U2;
		  coefs.Z = 6.0 * U - 9.0 * U2;
		  coefs.W = 3.0 * U2;
    }

    return coefs;
  }

  /** Calculate coefs for second derivative. */
  TPoint<T> coefs2(T U)
  { 
    LIMIT(U,0.0,1.0);

    TPoint<T> coefs;

	  T U1 = 1.0 - U; 
	  T U2 = U * U;

	  coefs.X = 6.0 * U1;
	  coefs.Y = -12.0 + 18.0 * U;
	  coefs.Z = 6.0 - 18.0 * U;
	  coefs.W = 6.0 * U;

    return coefs;
  }

  /** Multiply nodes by coefficients. */
  TPoint<T> multiplyCoefs(TPoint<T> coefs)
  {
    TPoint<T> r = TPoint<T>(0,0,0,0);

    for (int i = 0; i < 4; i++) 
    {
      r = r + this->cpoints[i] * coefs.XYZW[i];
    }

    return r;
  }

};

}