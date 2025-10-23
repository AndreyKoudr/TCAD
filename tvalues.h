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

  tvalues.h

  Operations on values as std::vector<T> like min/max etc.

*******************************************************************************/

#pragma once

#include "tbasics.h"

#include <vector>

namespace tcad {

/** Get min/max. */
template <class T> std::pair<T,T> calculateMinMax(const std::vector<T> &values, std::pair<int,int> *iminmax = nullptr)
{
  T min = std::numeric_limits<T>::max();
  T max = -std::numeric_limits<T>::max();
  if (iminmax)
  {
    iminmax->first = iminmax->second = -1;
  }

  for (int i = 0; i < values.size(); i++)
  {
    if (values[i] < min)
    {
      min = values[i];
      if (iminmax)
        iminmax->first = i;
    }

    if (values[i] > max)
    {
      max = values[i];
      if (iminmax)
        iminmax->second = i;
    }
  }

  std::pair<T,T> minmax(min,max);
  return minmax;
}

/** Get mean value. */
template <class T> T getMean(std::vector<T> &x)
{
  if (x.empty())
    return 0.0;

  T sum = 0.0;
  for (int i = 0; i < int(x.size()); i++)
  {
    sum += x[i];
  }

  sum /= T(x.size());

  return sum;
}

/** Get variance. */
template <class T> T getVariance(std::vector<T> &x)
{
  T mean = getMean(x);

  T sum = 0.0;
  for (int i = 0; i < x.size(); i++)
  {
    T t = x[i] - mean;
    sum += t * t;
  }

  return sum / T(x.size());
}

/** Get co-variance. */
template <class T> T getCovariance(std::vector<T> &x, std::vector<T> &y)
{
  T xmean = getMean(x);
  T ymean = getMean(y);

  T sum = 0.0;
  for (int i = 0; i < x.size(); i++)
  {
    sum += (x[i] - xmean) * (y[i] - ymean);
  }

  return sum / T(x.size());
}

}
