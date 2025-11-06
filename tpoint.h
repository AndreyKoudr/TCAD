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

  tpoint.h

  Single 3D point with a homogeneous coordinate

  dimensions : 0 (point)

*******************************************************************************/

#pragma once

#include "tbasics.h"

namespace tcad {

/**
  4-component vector template, no SIMD. The fourth component is a homogeneous 
coordinate.
  It does not make much sense to make a template for a N-component (e.g. 3) 
vector as currently there are such things as SIMD enabling to make operations
on all 4 components at once.
 
  Operations
  ----------
  "float" below is 4 or 8 byte float

  float = V[AxisX]  - get vector coordinate
  float = !V        - get vector length                         
  v3 = v1 ^ v2      - cross product                     
  v2 = +v1          - normalisation
  float = v1 * v2   - dot product, W is ignored
  (v1 > v2)         - vectors co-directed?    
  bool(v1==v2)      - vectors equal?                            
  v2 = -v1          - change sign of components                 
  v3 = v1 + v2      - addition
  v1 += v2;         - addition
  v3 = v1 - v2      - subtraction
  v1 -= v2;         - subtraction
  v3 = v1 * float   - multiply by scalar
  v3 = float * v    - multiply by scalar
  v1 *= float       - multiply by scalar
  v3 = v1 / float   - divide by scalar
  v1 /= float       - divide by scalar
*/

                              // axes
enum Axes {AxisX,AxisY,AxisZ,AxisW};                  

template <class T> class TPoint {
public:
                              // data
  union {
    struct { T X,Y,Z,W; };
    T XYZW[4];
    T XYZ[3];
  };
                              // constructors
  TPoint() 
  { 
                              // better to avoid any initialisation here as it takes time
    init();
  }
                              // initialise to zero vector
  void init()
  {
    X = Y = Z = (T) (0.0); 
    W = (T) (0.0);
  }

  TPoint(T x, T y = 0.0, T z = 0.0, T w = 0.0) 
  { 
    X = x; Y = y; Z = z; W = w; 
  }

  TPoint(const TPoint<T>& other) 
  { 
    X = other.X; 
    Y = other.Y; 
    Z = other.Z; 
    W = other.W; 
  }

  TPoint operator=(const TPoint<T>& other) 
  { 
    X = other.X; 
    Y = other.Y; 
    Z = other.Z; 
    W = other.W; 

    return *this;
  }
                              // get coordinate
  T& operator[](Axes A) {
    switch (A) {
    case AxisX : return X ;
    case AxisY : return Y ;
    case AxisZ : return Z ;
    case AxisW : return W ;
    default : assert(false && "Wrong axis value");
    };
  }
                              // get coordinate
  T& operator[](size_t index) {
    assert(index >= 0 && index <= 3);
    switch (index) {
    case 0 : return X ;
    case 1 : return Y ;
    case 2 : return Z ;
    case 3 : return W ;
    default : return X;
    };
  }
                              // operator ! - get vector length
  T operator!() const
  {
    return static_cast<T>(sqrt(X * X + Y * Y + Z * Z));
  }
                              // length
  T length() const
  {
    return static_cast<T>(sqrt(X * X + Y * Y + Z * Z));
  }
                              // length squared
  T length2() const
  {
    return X * X + Y * Y + Z * Z;
  }
                              // operator ^ for two vectors (cross product)
  TPoint operator^(TPoint v2) const
  {
    TPoint v3;
		v3.X = Y*v2.Z-Z*v2.Y;
		v3.Y = Z*v2.X-X*v2.Z;
		v3.Z = X*v2.Y-Y*v2.X;
    return v3;
  }
                              // operator + - normalization
  TPoint operator+() {
    TPoint v3;
    T Len = static_cast<T>(sqrt(X*X+Y*Y+Z*Z));

		if (Len > tolerance_)
		{
			v3.X = X/Len;
			v3.Y = Y/Len;
			v3.Z = Z/Len;
		} else
		{
      v3.X = 0;
      v3.Y = 0;
      v3.Z = 0;
		}
    return v3;
  }
                              // dot product, W is ignored
  T operator*(TPoint v2) const {
    return X * v2.X + Y * v2.Y + Z * v2.Z;
  }
                              // matrix multiply 1 x 4 x 1, including W
  T matrixMult(TPoint v2) const {
    return X * v2.X + Y * v2.Y + Z * v2.Z + W * v2.W;
  }
                              // operator > - two vectors co-directed?
  bool operator>(TPoint v2) {
    return ((*this * v2) >= 0);
  }
                              // operator == - two vectors equal?
  bool operator==(TPoint v2) const;
                              // operator - (change direction)
  TPoint operator-() {
    return (*this) * (-1.0);
  }
                              // addition
  TPoint operator+=(TPoint v2) {
    for (int i = 0; i < 4; i++) XYZW[i] += v2.XYZW[i];
    return *this;
  }

  TPoint inline operator+(TPoint v2) const {
    TPoint v3;
    for (int i = 0; i < 4; i++) v3.XYZW[i] = XYZW[i] + v2.XYZW[i];
    return v3;
  }
			                        // subtraction
  TPoint operator-(TPoint v2) const {
    TPoint v3;
    for (int i = 0; i < 4; i++) v3.XYZW[i] = XYZW[i] - v2.XYZW[i];
    return v3;
  }
                              // subtract
  TPoint operator-=(TPoint v2) {
    for (int i = 0; i < 4; i++) XYZW[i] -= v2.XYZW[i];
    return *this;
  }
                              // operator *scalar
  TPoint operator*(T multiplier) const {
    TPoint v3;
    v3.X = X*multiplier;
    v3.Y = Y*multiplier;
    v3.Z = Z*multiplier;
    v3.W = W*multiplier;
    return v3;
  }
                              // operator *= scalar
  TPoint operator*=(T multiplier) {
    X *= multiplier;
    Y *= multiplier;
    Z *= multiplier;
    W *= multiplier;
    return *this;
  }
                              // operator / scalar
  TPoint operator/(T divider) const {
    TPoint v3;
    assert(std::abs(divider) > tolerance_ && "Division by zero in TPoint");
    v3.X = X/divider;
    v3.Y = Y/divider;
    v3.Z = Z/divider;
    v3.W = W/divider;
    return v3;
  }
                              // operator /= scalar
  TPoint operator/=(T divider) {
    assert(std::abs(divider) > tolerance_ && "Division by zero in TPoint");
    X /= divider;
    Y /= divider;
    Z /= divider;
    W /= divider;
    return *this;
  }
                              // operator /= vector
  TPoint operator/=(TPoint divider) {
    assert((std::abs(divider.X) > tolerance_ && std::abs(divider.Y) > tolerance_ &&
      std::abs(divider.Z) > tolerance_ && std::abs(divider.W) > tolerance_) && "Division by zero in TPoint");
    X /= divider.X;
    Y /= divider.Y;
    Z /= divider.Z;
    W /= divider.W;
    return *this;
  }

  static TPoint zero()
  {
    return TPoint(0,0,0,0);
  }

  /** Get angle between vectors in radians. */
  T operator<(TPoint V2) {
    T L1 = !(*this);
    T L2 = !V2;
    T LL = L1 * L2;
    if (L1 < TOLERANCE(T) || L2 < TOLERANCE(T))
    {
      return (M_PI * T(0.5));
    } else
    {
      T S = (*this * V2) / LL;
			LIMIT(S,T(-1.0),T(+1.0));
      return (acos(S));
    }
  }

	/** Get angle in radians between this vector and v, 
    measured counter-clockwise; plane normal is
    plane with normal to make rotation comparisons. */
	T GetAngleSigned(TPoint planennormal, TPoint v)
	{
		TPoint p = (+(*this)) ^ (+v);
    T plen = !p;
    //assert(plen <= 1.0);
    LIMIT_MAX(plen,1.0);
    T a = fabs(asin(plen));

    if (!(p > planennormal))
      a = -a;

		return a;
	}

  /** To polar coords. X is angle! */
  TPoint toPolar()
  {
	  return TPoint(atan2(Y,X),sqrt(X * X + Y * Y),Z,W);
  }

  /** To Cartesian coords. X is angle! */
  TPoint toCartesian()
  {
	  return TPoint(Y * cos(X),Y * sin(X),Z,W);
  }

  /** Get angle between this vector and v, measured counter-clockwise,
    in the range [0..360]; correct in X-Y ONLY!!! */
	T getAngleBetween0360(TPoint v, bool counterclockwise = true)
	{
		T a = ((*this) < v) * PCI;

		TPoint p = *this ^ v;

    if (counterclockwise)
    {
		  if (p.Z < 0) a = 360.0 - a;
    } else
    {
		  if (p.Z > 0) a = 360.0 - a;
    }
		
		LIMIT(a,T(0.0),T(360.0));

		return a;
	}

  /** tolerance. */
  static T tolerance()
  {
    return tolerance_;
  }

  /** Set tolerance. */
  static void setTolerance(const T tolerance)
  {
    tolerance_ = tolerance;
  }

  /** Set defualt tolerance. */
  static void setDefaultTolerance(const T tolerance)
  {
    tolerance_ = TOLERANCE(T);
  }

private:
  // tolerance
  static T tolerance_;
};

template <typename T> T TPoint<T>::tolerance_ = TOLERANCE(T);

/** operator : scalar * vector. */
template <class T> TPoint<T> operator*(const double scalar, const TPoint<T> v)
{
  return v * static_cast<T>(scalar);
}

/** Get min from components. */
template <typename T> TPoint<T> pointMin(TPoint<T> &first, TPoint<T> &second)
{
  return TPoint<T>(std::min<T>(first.X,second.X),std::min<T>(first.Y,second.Y),std::min<T>(first.Z,second.Z));
}

/** Get max from components. */
template <typename T> TPoint<T> pointMax(TPoint<T> &first, TPoint<T> &second)
{
  return TPoint<T>(std::max<T>(first.X,second.X),std::max<T>(first.Y,second.Y),std::max<T>(first.Z,second.Z));
}



}
