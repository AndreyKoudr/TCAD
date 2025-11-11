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

  tbeziervolume.h

  Volume with Bezier basis. Parameterisation is uniform.

  dimensions : 3 (U,V,W parameters)

*******************************************************************************/

#pragma once

#include "tbasevolume.h"

namespace tcad {

template <class T> class TBezierVolume : public TBaseVolume<T> {
public:

  //===== Construction =========================================================

  /** Constructor. */
  TBezierVolume() : TBaseVolume<T>() {}

  /** Constrictor for box. numsegs... - number of Bezier segments in each direction. */
  TBezierVolume(TPoint<T> min, TPoint<T> max, int numsegsU, int numsegsV, int numsegsW) : TBaseVolume<T>() 
  {
    initBox(min,max,numsegsU,numsegsV,numsegsW);
  }

  /** Constructor. */
  TBezierVolume(const TBezierVolume &other) : TBaseVolume<T>()
  {
    this->K1 = other.K1;
    this->K2 = other.K2;
    this->K3 = other.K3;
    this->cpoints = other.cpoints;
    this->Uknots = other.Uknots;
    this->Vknots = other.Vknots;
    this->Wknots = other.Wknots;

    update();
  }

  /** Assignment operator. */
  TBezierVolume &operator = (const TBezierVolume &other)  
  {
    this->K1 = other.K1;
    this->K2 = other.K2;
    this->K3 = other.K3;
    this->cpoints = other.cpoints;
    this->Uknots = other.Uknots;
    this->Vknots = other.Vknots;
    this->Wknots = other.Wknots;

    update();

    return *this;
  }

  /** Initialize from box. */
  void initBox(TPoint<T> min, TPoint<T> max, int numsegsU = 1, int numsegsV = 1, int numsegsW = 1) 
  {
    assert(numsegsU > 0);
    assert(numsegsV > 0);
    assert(numsegsW > 0);

    LIMIT_MIN(numsegsU,1);
    LIMIT_MIN(numsegsV,1);
    LIMIT_MIN(numsegsW,1);

    this->cpoints.clear();

    int pnumU = numsegsU * 3 + 1;
    int pnumV = numsegsV * 3 + 1;
    int pnumW = numsegsW * 3 + 1;

    this->K1 = pnumU - 1;
    this->K2 = pnumV - 1;
    this->K3 = pnumW - 1;

    TPoint<T> d = max - min;
    T DZ = d.Z / T(pnumW - 1);
    T DY = d.Y / T(pnumV - 1);
    T DX = d.X / T(pnumU - 1);

    for (int i = 0; i < pnumW; i++)
    {
      T Z = min.Z + T(i) * DZ;
      for (int j = 0; j < pnumV; j++)
      {
        T Y = min.Y + T(j) * DY;
        for (int k = 0; k < pnumU; k++)
        {
          T X = min.X + T(k) * DX;
          this->cpoints.push_back(TPoint<T>(X,Y,Z));
        }
      }
    }

    // make knots between segments
    makeKnots();
  }

  /** Destructor. */
  virtual ~TBezierVolume() {}

  //===== Abstract =============================================================

  /** Get k-th derivative on parameter U,V [0..1]. 0-derivative is position(U,V) (see below). */
  virtual TPoint<T> derivative(T U, T V, T W, Parameter onparameter, int k)
  {
    assert(this->cpoints.size() > 0);

    LIMIT(U,0.0,1.0);
    LIMIT(V,0.0,1.0);
    LIMIT(W,0.0,1.0);

    int numU = this->K1 + 1;
    int numV = this->K2 + 1;
    int numW = this->K3 + 1;

    // basis functions at point U,V,W
    std::vector<T> funcsU,funcsV,funcsW;

    // only 0-th derivative done yet
    if (k == 0)
    {
      getUVWfuncs(U,V,W,0,0,0,funcsU,funcsV,funcsW);
    } else if (k == 1)
    {
      if (onparameter == PARAMETER_U)
      {
        getUVWfuncs(U,V,W,1,0,0,funcsU,funcsV,funcsW);
      } else if (onparameter == PARAMETER_V)
      {
        getUVWfuncs(U,V,W,0,1,0,funcsU,funcsV,funcsW);
      } else if (onparameter == PARAMETER_W)
      {
        getUVWfuncs(U,V,W,0,0,1,funcsU,funcsV,funcsW);
      }
    } else if (k == 2) //!!! this stuff is untested
    {
      if (onparameter == PARAMETER_UU)
      {
        getUVWfuncs(U,V,W,2,0,0,funcsU,funcsV,funcsW);
      } else if (onparameter == PARAMETER_UV)
      {
        getUVWfuncs(U,V,W,1,1,0,funcsU,funcsV,funcsW);
      } else if (onparameter == PARAMETER_UW)
      {
        getUVWfuncs(U,V,W,1,0,1,funcsU,funcsV,funcsW);

      } else if (onparameter == PARAMETER_VU)
      {
        getUVWfuncs(U,V,W,1,1,0,funcsU,funcsV,funcsW);
      } else if (onparameter == PARAMETER_VV)
      {
        getUVWfuncs(U,V,W,0,2,0,funcsU,funcsV,funcsW);
      } else if (onparameter == PARAMETER_VW)
      {
        getUVWfuncs(U,V,W,0,1,1,funcsU,funcsV,funcsW);

      } else if (onparameter == PARAMETER_WU)
      {
        getUVWfuncs(U,V,W,1,0,1,funcsU,funcsV,funcsW);
      } else if (onparameter == PARAMETER_WV)
      {
        getUVWfuncs(U,V,W,0,1,1,funcsU,funcsV,funcsW);
      } else if (onparameter == PARAMETER_WW)
      {
        getUVWfuncs(U,V,W,0,0,2,funcsU,funcsV,funcsW);
      }
    }

    TPoint<T> result;

    if (!funcsU.empty() && !funcsV.empty() && !funcsW.empty())
    {
      // temp
      int numVU = numV * numU;
      for (int i = 0; i < numW; i++)
      {
        for (int j = 0; j < numV; j++)
        {
          for (int k = 0; k < numU; k++)
          {
            result = result + this->cpoints[i * numVU + j * numU + k] * (funcsW[i] * funcsV[j] * funcsU[k]);
          }
        }
      }
    }

    return result;
  }

  /** Update after any change in control points. */
  virtual void update()
  {
  }

  /** C1 smooth, Bezier style. */
  void smooth(bool Usmooth, bool Vsmooth, bool Wsmooth)
  {
    int numU = this->K1 + 1;
    int numV = this->K2 + 1;
    int numW = this->K3 + 1;

    if (Usmooth)
    {
      for (int j = 0; j < numV; j++)
      {
        for (int k = 0; k < numW; k++)
        {
          std::vector<int> points;
          this->getULinePoints(j,k,points);

          for (int l = 3; l < points.size() - 3; l += 3)
          {
            straightenThreePoints(this->cpoints[points[l - 1]],this->cpoints[points[l]], 
              this->cpoints[points[l + 1]]);
          }
        }
      }
    }

    if (Vsmooth)
    {
      for (int j = 0; j < numU; j++)
      {
        for (int k = 0; k < numW; k++)
        {
          std::vector<int> points;
          this->getVLinePoints(j,k,points);

          for (int l = 3; l < points.size() - 3; l += 3)
          {
            straightenThreePoints(this->cpoints[points[l - 1]],this->cpoints[points[l]], 
              this->cpoints[points[l + 1]]);
          }
        }
      }
    }

    if (Wsmooth)
    {
      for (int j = 0; j < numU; j++)
      {
        for (int k = 0; k < numV; k++)
        {
          std::vector<int> points;
          this->getWLinePoints(j,k,points);

          for (int l = 3; l < points.size() - 3; l += 3)
          {
            straightenThreePoints(this->cpoints[points[l - 1]],this->cpoints[points[l]], 
              this->cpoints[points[l + 1]]);
          }
        }
      }
    }
  }

  
public:

///*
//            ---------------------
//         /                      /|
//         ----------------------  |
//        |                      | |
//  V (j) |                      | |
//  ^     |                      | |
//  |     |                      |/
//         ----------------------
// /   -> U (i)
//W (k)
//
//*/
//
//  std::vector<T> Uknots;     
//  std::vector<T> Vknots;   
//  std::vector<T> Wknots;    
//
//  // number of columns minus 1
//  int K1 = 0;
//
//  // number of rows munus 1
//  int K2 = 0;
//
//  // number of levels minus 1
//  int K3 = 0;

protected:

  //// (K1 + 1) * (K2 + 1) * (K3 + 1) control points, call update() after every change
  //std::vector<TPoint<T>> cpoints;

  /** Make uniformly ditributed knots between every 4 Bezier segments. */
  static void makeUniformKnots(int numPoints, std::vector<T> &knots)
  {
    knots.clear();

    int numSegs = (numPoints - 1) / 3;
    T du = 1.0 / T(numSegs);
    for (int i = 0; i <= numSegs; i++)
    {
      knots.push_back(du * T(i));
    }
  }

  /** Make all knots. */
  void makeKnots()
  {
    // make uniformly spaced knots
    makeUniformKnots(this->K1 + 1,this->Uknots);
    makeUniformKnots(this->K2 + 1,this->Vknots);
    makeUniformKnots(this->K3 + 1,this->Wknots);
  }

private:

  /** Create basis for derivatives. */
  void getUVWfuncs(T U, T V, T W, int derU, int derV, int derW, 
    std::vector<T> &funcsU, std::vector<T> &funcsV, std::vector<T> &funcsW)
  {
    int segmentU,segmentV,segmentW;
    if (derU == 0)
    {
      BezierBasis(U,this->Uknots,funcsU,segmentU);
    } else if (derU == 1)
    {
      BezierBasisDer1(U,this->Uknots,funcsU,segmentU);
    } else if (derU == 2)
    {
      BezierBasisDer2(U,this->Uknots,funcsU,segmentU);
    }

    if (derV == 0)
    {
      BezierBasis(V,this->Vknots,funcsV,segmentV);
    } else if (derV == 1)
    {
      BezierBasisDer1(V,this->Vknots,funcsV,segmentV);
    } else if (derV == 2)
    {
      BezierBasisDer2(V,this->Vknots,funcsV,segmentV);
    }

    if (derW == 0)
    {
      BezierBasis(W,this->Wknots,funcsW,segmentW);
    } else if (derW == 1)
    {
      BezierBasisDer1(W,this->Wknots,funcsW,segmentW);
    } else if (derW == 2)
    {
      BezierBasisDer2(W,this->Wknots,funcsW,segmentW);
    }
  }

};

}
