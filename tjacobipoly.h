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

  tjacobipoly.h

  Fitting with Jacobi orthogonal polynomials

*******************************************************************************/

#pragma once

#include "tvalues.h"
#include "tpoints.h"
#include "tgaussint.h"

namespace tcad {

/**
  Fitting points with Jacobi orthogonal polynomials 
  -------------------------------------------------
(https://en.wikipedia.org/wiki/Jacobi_polynomials)

  Compiler
  --------
  VS 2019, can be easily converted into Linux : just std:: and C++ 11/14.

  Parameters
  ----------
  Input : pairs of x-y points, Output : Jacobi poly. The code builds ortho poly 
approximation on points f(x). x must be monotonically increasing.
  Parameter U is within [-1.0 .. +1.0] everywhere.
  Parameters alpha and beta define type of Jacobi polynomial; the case 
alpha = beta = 0.0 corresponds to a Legendre poly. This case is common for
approximation of most curves. But if you need to fit a curve y(x) with infinite 
derivatives at the ends, like a half circle, or an aerofoil surface, 
set alpha = beta = 0.5.
  The trick is that actually f(x) / ((1 - x) ^ alpha) * (1 + x) ^ beta), not 
f(x), is approximated.

  Integration
  -----------
  As the poly is orthogonal, it means that approxomation coefs are calculated by 
integration. There are two ways to integrate : by numerical Gauss (set fit() 
parameter as GAUSSINT_8 for example) or by trapezoid rule (set OTHER_INTEGRATION). 
Keep in mind that Gauss integration is built upon a single poly itself and it 
cannot well integrate curves of multiple bends. Use trapezoid rule for this, but 
it requires many points.

  Poly degree
  -----------
  Keep in mind that a single poly is built upon the whole region and the poly 
degree must not be very high, say, up to 20; otherwise gamma function / factorial
calculation problems will appear (no checks).

  Tests
  -----
Test 1 : horizontal curve, Legendre poly (0.0,0.0) and Gauss integration
Accuracy 0.000000
Test 2 : horizontal curve, Legendre poly (0.0,0.0) and trapezoid integration, few points
Accuracy 0.000000
Test 3 : horizontal curve, Legendre poly (0.0,0.0) and trapezoid integration, many points
Accuracy 0.000000
Test 4 : half of circle, infinite derivative at ends, Legendre poly (0.0,0.0) and Gauss integration
Accuracy 0.384346
Test 5 : half of circle, infinite derivative at ends, use Jacoby poly (0.5,0.5) and Gauss integration
Accuracy 0.006320
Test 6 : fitting cos curve from 0 to 4 Pi. Legendre poly. Trapezoidal integration on many points
Accuracy 0.052173
Test 7 : fitting cos curve from 0 to 4 Pi. Legendre poly. Gauss integration on many points
Accuracy 0.052173
*/


template <class T> class TJacobiPoly {
public:
                            // "alpha" and "beta"
  T alpha = T(0.0);
  T beta = T(0.0);

                            // y min max before scaling
  T ymin = T(0.0);
  T ymax = T(0.0);
                            // x min max before scaling
  T xmin = T(0.0);
  T xmax = T(0.0);

                            // constructors
  TJacobiPoly() = default;
  TJacobiPoly(const T palpha, const T pbeta);
                            // fit arbitrary function y(x), use GAUSSINT_... 
                            // for Gauss integration or OTHER_INTEGRATION for 
                            // trapezoid rule (for many points and complicated
                            // curve shape)
  bool fit(const int degree, std::vector<T> &x, std::vector<T> &y, 
    int integration);
                            // get value on U [-1..+1]!!!
  T getValue(const T U);
                            // get derivative on K on U [-1..+1]
  T getDerK(const T U, const int k);
                            // get accuracy of fitting as max y deflection from
                            // original data
  T accuracy(const std::vector<T> &x, const std::vector<T> &y);

private:
                            // coefficients; real degree is one less than the number of coefs
  std::vector<T> coefs;
                            // get value, degree is one less than the number of coefs,
                            // n being degree,
                            // a - Jacobi poly alpha, b - beta, U [-1..+1]
  static T getPolyValuePrim(int n, const T a, const T b, const T U);
                            // get k-derivative
  static T getPolyDerKPrim(int n, T a, T b, T U, int k);
                            // get value, degree is one less than the number of coefs,
                            // U [-1..+1]
  T getPolyValue(const int degree, const T U);
                            // k-th derivative
  T getDerK(const int degree, const T U, const int k);
                            // get poly values recurrently
  static void getPolyValues(const int n, const T a, const T b, const T U, std::vector<T> &values);
                            // get orthogonality coefs 
  static void getOrthogonalityCoefs(const int n, const T a, const T b, std::vector<T> &values);
};

// Find value inside montonic table
template <class T> static int FindInterval(const std::vector<T> table, const T value)
{
                              // size, signed integer to avoid problems with size() - 2 (last interval)
  int size = static_cast<int>(table.size());
                              // increasing?
  bool increasing = ((table[size - 1] - table[0]) >= 0.0);
                              // table must be monotonic
  #ifdef _DEBUG
    if (increasing)
    {
      for (int i = 0; i < size - 1; i++)
      {
        if (table[i + 1] < table[i])
          assert(false && "Table not monotonic");
      }
    } else
    {
      for (int i = 0; i < size - 1; i++)
      {
        if (table[i + 1] > table[i])
          assert(false && "Table not monotonic");
      }
    }
  #endif
                              // lower and upper, originally out of scope
  int lower = -1;
  int upper = size - 1;
                              // min/max
  T min = increasing ? table[0] : table[size - 1];
  T max = increasing ? table[size - 1] : table[0];
                              // tolerance
  T tolerance = TOLERANCE(T);
                              // obvious outcome
  if (value < min - tolerance)
    return -1;
  if (value > max + tolerance)
    return -1;

  while ((upper - lower) > 1)
  {
                              // middle point
    int middle = (lower + upper) >> 1;
    if ((value >= table[middle]) == increasing)
    {
      lower = middle;
    } else
    {
      upper = middle;
    }
  }
                              // we are going to return lower, correct it just in case
  if (lower < 0) lower = 0;
  if (lower > size - 2) lower = size - 2;
                              // done!
  return lower;
}

// Scale array of points to the range smin,smax and find its actual min/max 
template <class T> void rescale(std::vector<T> &points, const T smin, const T smax, T &min, T &max)
{
  assert(points.size() > 0);
  assert(smax > smin);
                            // get min/max
  min = max = points[0];
  for (size_t i = 1; i < points.size(); i++)
  {
    if (points[i] < min)
      min = points[i];
    if (points[i] > max)
      max = points[i];
  }
                            // scale to smin,smax
  T d = max - min;
  T sd = smax - smin;

  if (std::abs(d) > TOLERANCE(T))
  {
    for (size_t i = 0; i < points.size(); i++)
    {
      points[i] = smin + sd * (points[i] - min) / d;
    }
  } else
  {
    std::fill(points.begin(),points.end(),T(0.0));
  }
}

// Get y value for integration, xx must be from xmin to xmax
template <class T> static T getY(const std::vector<T> &x, const std::vector<T> &y, const T xx, T tolerance)
{
                            // find segment by bisection
  int i = FindInterval<T>(x,xx);
  assert(i >= 0);
                            // find value by linear interpolation within segment
  T dx = x[i + 1] - x[i];
  if (std::abs(dx) < tolerance)
  {
    return y[i];
  } else
  {
    return y[i] + (y[i + 1] - y[i]) * (xx - x[i]) / dx;
  }
}

// Factorial
template <class T> static T factorial(int n)
{
  return tgamma(n + 1);
}

// Gamma-function   
template <class T> static T gamma(T x)
{
  return tgamma(x);
}

// Binomial coefficient (z / n)
template <class T> static T binomialCoef(T z, T n)
{
  assert(n >= 0 && z >= n);

  T res = gamma<T>(z + T(1.0)) / (gamma<T>(z - n + T(1.0)) * gamma<T>(n + T(1.0)));

  return res;
}

template <class T> TJacobiPoly<T>::TJacobiPoly(const T palpha, const T pbeta) : alpha(palpha), beta(pbeta)
{
}

template <class T> T TJacobiPoly<T>::getValue(const T U)
{
  T sum = T(0.0);

  for (int i = 0; i < static_cast<int>(coefs.size()); i++)
  {
    sum += coefs[i] * getPolyValue(i,U);
  }
                              // approximate value is without (1 - x)^a * (1 + x)^b
  sum *= pow(T(1.0) - U,alpha) * pow(T(1.0) + U,beta);
                              // y scaled to 0..1
  sum = ymin + (ymax - ymin) * sum;

  return sum;
}

template <class T> T TJacobiPoly<T>::getDerK(const T U, const int k)
{
  T V2 = 0.0;
  T V2prim = 0.0;

  for (int i = 0; i < static_cast<int>(coefs.size()); i++)
  {
    V2 += coefs[i] * getPolyValue(i,U);
    V2prim += coefs[i] * getDerK(i,U,k);
  }

  // important : at U = -1 and U = 1 with alpha/beta = 0.5
  // pow(1.0 + U,beta - 1.0) and pow(1.0 - U,alpha - 1.0) easily become infinity
  // by 0 ^ -0.5. therefore,
  T u = U;
  LIMIT(u,-0.999999,+0.999999);

  T V0 = pow(1.0 - u,alpha);
  T V0prim = (std::abs(alpha) > TOLERANCE(T)) ? -alpha * pow(1.0 - u,alpha - 1.0) : 0.0;
  T V1 = pow(1.0 + u,beta);
  T V1prim = (std::abs(beta) > TOLERANCE(T)) ? beta * pow(1.0 + u,beta - 1.0) : 0.0;

  T der = V0prim * V1 * V2 + V0 * V1prim * V2 + V0 * V2prim * V1;

  assert(der == der);
                              // y scaled to 0..1
                              // x scaled to -1 .. +1
  der = der * 2.0 * (ymax - ymin) / (xmax - xmin); 

  return der;
}

template <class T> T TJacobiPoly<T>::getPolyValuePrim(int n, T a, T b, T U)
{
//  assert(n >= 0);
  assert(U >= T(-1.0) && U <= T(1.0));

#if 0 // this was for testing
    switch (n) {
    case 0 : return 1.0;
    case 1 : return U;
    case 2 : return 0.5 * (3.0 * U * U - 1.0);
    case 3 : return 0.5 * (5.0 * U * U * U - 3.0 * U);
    case 4 : return 0.125 * (35.0 * U * U * U * U - 30.0 * U * U + 3.0);
    case 5 : return 0.125 * (63.0 * pow(T(U),int(5)) - 70.0 * pow(T(U),int(3)) + 15.0 * U);
    case 6 : return 0.0625 * (231.0 * pow(T(U),int(6)) - 315.0 * pow(T(U),int(4)) + 105.0 * 
      pow(T(U),int(2)) - 5.0);
    default : return 0;
    }
#else
  T sum = T(0.0);

  for (int s = 0; s <= n; s++)
  {
    sum += binomialCoef<T>(n + a, n - s) * binomialCoef<T>(n + b, s) * 
      pow((U - T(1.0)) * T(0.5),s) * pow((U + T(1.0)) * T(0.5),n - s);
  }

  return sum;
#endif
}

template <class T> T TJacobiPoly<T>::getPolyDerKPrim(int n, T a, T b, T U, int k)
{
  assert(n >= 0);
  assert(U >= T(-1.0) && U <= T(1.0));

  T value = gamma<T>(a + b + n + 1 + k) * 
    getPolyValuePrim((n - k),a + k,b + k,U) / 
    (pow(2.0,k) * gamma<T>(a + b + n + 1));

  return value;
}

template <class T> void TJacobiPoly<T>::getOrthogonalityCoefs(const int degree, const T a, const T b, 
  std::vector<T> &values)
{
  assert(degree >= 0);
  values.clear();

  T ab1 = a + b + T(1.0);
  for (int n = 0; n <= degree; n++)
  {
    T fn = factorial<T>(n);
    T na1 = n + a + 1;
    T nb1 = n + b + 1;
    T n2ab1 = n * 2 + ab1;
    T nab1 = n + ab1;
    T v = pow(T(2.0),ab1) * gamma<T>(na1) * gamma<T>(nb1) / (n2ab1 * fn * gamma<T>(nab1));
    values.push_back(v);
  }
}

template <class T> void TJacobiPoly<T>::getPolyValues(const int n, const T a, const T b, 
  const T U, std::vector<T> &values)
{
  assert(n >= 0);
  assert(U >= -1.0 && U <= 1.0);
  values.clear();

#if 0 // this was for testing
  for (int m = 0; m <= n; m++)
  {
    T value = getPolyValuePrim(m,a,b,U);
    values.push_back(value);
  }
#else
  T ab = a + b;
  int njm = n - 1;
  values.push_back(T(1.0));
  values.push_back((a - b + (ab + 2) * U) * 0.5);

  for (int m = 2; m <= n; m++)
  {
    T nab = m + ab;
    T nab2 = m * 2 + ab;
    T nab21 = nab2 - 1;
    T nab22 = nab2 - 2;
    T na1 = m + a - 1;
    T nb1 = m + b - 1;
    T bracket = nab2 * nab22 * U + a * a - b * b;
    T c0 = 2 * m * nab *  nab22;
    T c1 = nab21 * bracket;
    T c2 = 2 * na1 * nb1 * nab2;
    T v = (c1 * values[m - 1] - c2 * values[m - 2]) / c0;

    values.push_back(v);
  }
#endif
}

template <class T> T TJacobiPoly<T>::getPolyValue(const int degree, const T U)
{
  assert(degree >= 0);
  assert(U >= -1.0 && U <= 1.0);

  return getPolyValuePrim(degree,alpha,beta,U);
}

template <class T> T TJacobiPoly<T>::getDerK(const int degree, const T U, const int k)
{
  assert(degree >= 0);
  assert(U >= -1.0 && U <= 1.0);

  return getPolyDerKPrim(degree,alpha,beta,U,k);
}

template <class T> bool comparePoints(const std::pair<T,T> &a, const std::pair<T,T> &b) {
    return a.first < b.first;
}

template <class T> bool TJacobiPoly<T>::fit(const int degree, std::vector<T> &x, std::vector<T> &y, 
  int integration)
{
  assert(x.size() > 1);
  assert(x.size() == y.size());

  if (x.size() > 1)
  {
    // we need to integrate by monotonically increasing x
    std::vector<T> xold = x;
    std::vector<T> yold = y;

    // sort by increasing x
    std::vector<std::pair<double,double>> xy;
    for (int i = 0; i < x.size(); i++)
    {
      xy.push_back(std::pair<double,double>(x[i],y[i]));
    }

    std::sort(xy.begin(),xy.end(),comparePoints<T>);

    for (int i = 0; i < x.size(); i++)
    {
      x[i] = xy[i].first;
      y[i] = xy[i].second;
    }

    std::vector<T> yscaled = y;
    rescale(yscaled,T(0.0),T(1.0),ymin,ymax);

    coefs.clear();
    coefs.resize(degree + 1,T(0.0));

    T tolerance = TOLERANCE(T);

    std::pair<T,T> minmax = calculateMinMax<T>(x);

    xmin = minmax.first;
    xmax = minmax.second;
    assert(xmax > xmin);
    T dx = xmax - xmin;
                              // orthogonality coefs
    std::vector<T> ocoefs;
    getOrthogonalityCoefs(degree,alpha,beta,ocoefs);

                              // integrate by Gauss
    if (integration == GAUSSINT_1 ||
      integration == GAUSSINT_2 ||
      integration == GAUSSINT_4 ||
      integration == GAUSSINT_8 ||
      integration == GAUSSINT_20)
    {
                          // if y-x function is complex,
                          // do not use Gauss integration here, as
                          // only a few points (max 8 in out case) are 
                          // taken into account
      for (int k = 0; k < GaussInt[integration].numpoints; k++)
      {
        T e = static_cast<T>(GaussInt[integration].knots[k]);
        T w = static_cast<T>(GaussInt[integration].weights[k]);

        T xx = xmin + dx * (e + T(1.0)) * T(0.5);
        T ymiddle = getY(x,yscaled,xx,tolerance);

        std::vector<T> values;
        getPolyValues(degree,alpha,beta,e,values);

        for (int j = 0; j <= degree; j++)
        {
          coefs[j] += ymiddle * values[j] * w / ocoefs[j];
        }
      }
    } else
    {
                          // compute coefficients by trapezoidal rule
      for (int j = 1; j < static_cast<int>(x.size()); j++)
      {
        assert(x[j] >= x[j - 1]);

        T xx = (x[j] + x[j - 1]) * T(0.5);
        T u = (xx - xmin) / dx;
        T u1 = (x[j - 1] - xmin) / dx;
        T u2 = (x[j] - xmin) / dx;
        T e = u + u - 1;

        std::vector<T> values;
        getPolyValues(degree,alpha,beta,e,values);

        for (int i = 0; i <= degree; i++)
        {
          coefs[i] += (yscaled[j] + yscaled[j - 1]) * T(0.5) * values[i] * (u2 - u1) * 
            T(2.0) /* -1..+1*/ / ocoefs[i];
        }
      }
    }

    //printf("coefs\n");
    //for (int i = 0; i < coefs.size(); i++)
    //{
    //  printf("i %d %f\n",i,coefs[i]);
    //}

    // restore old values
    x = xold;
    y = yold;

    return true;
  } else
  {
    coefs.clear();
    return false;
  }
}

template <class T> T TJacobiPoly<T>::accuracy(const std::vector<T> &x, const std::vector<T> &y)
{
  assert(x.size() > 1);
  assert(x.size() == y.size());
  assert(coefs.size() > 0);

  if (coefs.size() > 1)
  {
    T accuracy = T(0.0);

    // xmin/xmax must be set
    assert(xmax > xmin);
    T dx = xmax - xmin;

#ifdef _DEBUG
    std::vector<T> diffs;
    std::vector<T> polys;
#endif

    for (int j = 0; j < static_cast<int>(x.size()); j++)
    {
      T u = (x[j] - xmin) / dx;
                              // -1..+1
      u = u + u - 1;

      T poly = getValue(u);
      T der = getDerK(u,1);
      T diff = std::abs(poly - y[j]);
      T dxdy = (j > 0) ? ((y[j] - y[j - 1]) / (x[j] - x[j - 1])) : 0.0;

//      printf("U %f x %f y %f value %f dxdy %f der %f diff %f\n",u,x[j],y[j],poly,dxdy,der,diff);

#ifdef _DEBUG
      polys.push_back(poly);
      diffs.push_back(diff);
#endif

      if (diff > accuracy)
        accuracy = diff;
    }

    return accuracy;
  } else
  {
    return T(0.0);
  }
}

}