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

  tsplinesurface.h

  B-spline surface

  dimensions : 2 (U,V parameters)

*******************************************************************************/

#pragma once

#include "tbasesurface.h"
#include "tsplinecurve.h"

namespace tcad {

template <class T> class TSplineSurface : public TBaseSurface<T> {
public:

  //===== Construction =========================================================

  /** Constructor. */
  TSplineSurface() : TBaseSurface<T>()
  {
  }

  /** Constructor. */
  TSplineSurface(const TSplineSurface &other) : TBaseSurface<T>()
  {
    this->K1 = other.K1;
    this->K2 = other.K2;
    this->M1 = other.M1;
    this->M2 = other.M2;
    this->cpoints = other.cpoints;

    interpolate = other.interpolate;
    clampedstartU = other.clampedstartU;
    clampedendU = other.clampedendU;
    clampedstartV = other.clampedstartV;
    clampedendV = other.clampedendV;

    points = other.points;

    update();
  }

  /** Assignment operator. */
  TSplineSurface &operator = (const TSplineSurface &other)  
  {
    this->K1 = other.K1;
    this->K2 = other.K2;
    this->M1 = other.M1;
    this->M2 = other.M2;
    this->cpoints = other.cpoints;

    interpolate = other.interpolate;
    clampedstartU = other.clampedstartU;
    clampedendU = other.clampedendU;
    clampedstartV = other.clampedstartV;
    clampedendV = other.clampedendV;

    points = other.points;

    update();

    return *this;
  }

  /** Constructor. Every points[i] is a row of points (U-changing). */
  TSplineSurface(std::vector<std::vector<TPoint<T>>> &ppoints, 
    int k1, int m1, int k2, int m2, 
    CurveEndType startU, CurveEndType endU,
    CurveEndType startV, CurveEndType endV) : TBaseSurface<T>()
  {
    this->points = ppoints;

    this->cpoints.clear();

    allocate(k1,m1,k2,m2,false);

    clampedstartU = (startU == END_CLAMPED);
    clampedendU = (endU == END_CLAMPED);
    clampedstartV = (startV == END_CLAMPED);
    clampedendV = (endV == END_CLAMPED);

    interpolate = false;

    update();
  }

  /** Constructor. Every points[i] is a row of points (U-changing fastest). */
  TSplineSurface(std::vector<std::vector<TPoint<T>>> &ppoints, 
    int m1, int m2, 
    CurveEndType startU, CurveEndType endU,
    CurveEndType startV, CurveEndType endV) : TBaseSurface<T>()
  {
    int k1 = int(ppoints[0].size()) - 1;
    int k2 = int(ppoints.size()) - 1;

    this->points = ppoints;

    this->cpoints.clear();

    allocate(k1,m1,k2,m2,false);

    clampedstartU = (startU == END_CLAMPED);
    clampedendU = (endU == END_CLAMPED);
    clampedstartV = (startV == END_CLAMPED);
    clampedendV = (endV == END_CLAMPED);

    interpolate = true;

    update();
  }

  /** Constructor from any other surface. Set refine... for start/end of parametric
    directions to 0.5 to handle rounded edges. */
  TSplineSurface(TBaseSurface<T> &surface, 
    int k1, int m1, int k2, int m2, 
    CurveEndType startU, CurveEndType endU,
    CurveEndType startV, CurveEndType endV,
    T refinestartU = 1.0, T refineendU = 1.0, 
    T refinestartV = 1.0, T refineendV = 1.0) : TBaseSurface<T>()
  {
    std::vector<TPoint<T>> temp;
    int tempK1 = 0;
    int tempK2 = 0;
    surface.createPoints(temp,nullptr,&tempK1,&tempK2, 
       k1 + 1,k2 + 1,refinestartU,refineendU,refinestartV,refineendV);

    // make 2D points
    points1Dto2D(temp,tempK1,tempK2,this->points);

    this->cpoints.clear();

    allocate(k1,m1,k2,m2,false);

    clampedstartU = (startU == END_CLAMPED);
    clampedendU = (endU == END_CLAMPED);
    clampedstartV = (startV == END_CLAMPED);
    clampedendV = (endV == END_CLAMPED);

    interpolate = false;

    update();
  }

  /** Destructor. */
  virtual ~TSplineSurface() 
  {
    DELETE_CLASS(Uderivative);
    DELETE_CLASS(UUderivative);
    DELETE_CLASS(UVderivative);
    DELETE_CLASS(Vderivative);
    DELETE_CLASS(VVderivative);
    DELETE_CLASS(VUderivative);
  }

  //===== Abstract =============================================================

  /** Get k-th derivative on parameter U,V [0..1]. 0-derivative is position(U,V) (see below). 
    No UV derivatives. */
  virtual TPoint<T> derivative(T U, T V, Parameter onparameter, int k)
  {
    LIMIT(U,0.0,1.0);
    LIMIT(V,0.0,1.0);

    if (k == 0)
    {
      TPoint<T> result;

      std::vector<T> bs(this->K1 + this->M1 + 2,0.0);
      std::vector<T> bt(this->K2 + this->M2 + 2,0.0);

      splineBasis(this->M1 + 1,U,Uknots,bs);
      splineBasis(this->M2 + 1,V,Vknots,bt);

      for (int i = 0; i <= this->K1; i++)
      {
        for (int j = 0; j <= this->K2; j++)
        {
          result = result + this->cpoints[j * (this->K1 + 1) + i] * (bs[i] * bt[j]);
        }
      }

      return result;
    } else if (k == 1)
    {
      if (onparameter == PARAMETER_U)
      {
        return Uderivative->position(U,V);
      } else if (onparameter == PARAMETER_V)
      {
        return Vderivative->position(U,V);
      } else
      {
        return TPoint<T>();
      }
    } else if (k == 2)
    {
      if (onparameter == PARAMETER_UU)
      {
        return UUderivative->position(U,V);
      } else if (onparameter == PARAMETER_UV)
      {
        return UVderivative->position(U,V);
      } else if (onparameter == PARAMETER_VV)
      {
        return VVderivative->position(U,V);
      } else if (onparameter == PARAMETER_VU)
      {
        return VUderivative->position(U,V);
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
    // number of points here must equal (K1 + 1) x (K2 + 1)
    std::vector<std::vector<TPoint<T>>> kpoints;

    if (points.size() == this->K2 + 1 && points[0].size() == this->K1 + 1)
    {
      kpoints = points;
    } else
    // it is a problem : we need new points (K1 + 1) x (K2 + 1)
    {
      kpoints.resize((this->K2 + 1),std::vector<TPoint<T>>((this->K1 + 1),TPoint<T>()));

      std::vector<TPoint<T>> temp;
      int tempK1 = 0;
      int tempK2 = int(points.size()) - 1;

      // loop by rows
      for (int i = 0; i < int(points.size()); i++)
      {
        TSplineCurve<T> crow(points[i],this->K1,this->M1,clampedstartU,clampedendU);
        temp.insert(temp.end(),crow.controlPoints().begin(),crow.controlPoints().end());
        if (i == 0)
        {
          tempK1 = int(crow.controlPoints().size()) - 1;
        }
      }

      // tempK2 is wrong here; we must redivide columns

      // now we need to redivide temp columns in V direction
      for (int i = 0; i <= tempK1; i++)
      {
        std::vector<TPoint<T>> col;
        tcad::getColumn(temp,tempK1,tempK2,i,col);
        TSplineCurve<T> ccol(col,this->K2,this->M2,clampedstartV,clampedendV);

        // set this column in kpoints
        for (int k = 0; k <= this->K2; k++)
        {
          kpoints[k][i] = ccol.controlPoints()[k];
        }
      }
    }

    // knots
    makeKnots();

    assert(kpoints.size() == this->K2 + 1 && kpoints[0].size() == this->K1 + 1);

    if (interpolate)
    {
      // clamping is done inside, do not do it twice; if interpolatePoints() fails
      // it means that the points are approximated
      if (!interpolatePoints(kpoints))
      {
        if (clampedstartU || clampedendU)
        {
          setUClamping(points,clampedstartU,clampedendU);
        }
  
        if (clampedstartV || clampedendV)
        {
          setVClamping(points,clampedstartV,clampedendV);
        }
      }
    } else
    {
      approximatePoints(kpoints);
  
      if (clampedstartU || clampedendU)
      {
        setUClamping(points,clampedstartU,clampedendU);
      }
  
      if (clampedstartV || clampedendV)
      {
        setVClamping(points,clampedstartV,clampedendV);
      }
    }

    // prepare two derivatives
    if (Uderivative == nullptr)
      Uderivative = new TSplineSurface();
    if (UUderivative == nullptr)
      UUderivative = new TSplineSurface();
    if (UVderivative == nullptr)
      UVderivative = new TSplineSurface();
    if (Vderivative == nullptr)
      Vderivative = new TSplineSurface();
    if (VVderivative == nullptr)
      VVderivative = new TSplineSurface();
    if (VUderivative == nullptr)
      VUderivative = new TSplineSurface();

    makeDerivatives(*Uderivative,*Vderivative);

    Uderivative->makeDerivatives(*UUderivative,*UVderivative);
    Vderivative->makeDerivatives(*VUderivative,*VVderivative);
  }

  /** Make transform. */
  virtual void makeTransform(TTransform<T> *transform)
  {
    for (auto &line : points)
    {
      for (auto &p : line)
      {
        p = transform->applyTransform(p);
      }
    }

    // call virtual update()
    this->update();
  }

  /** Set U clamping at start / end. points are the original points from where
    the derivatives are taken. */
  void setUClamping(std::vector<std::vector<TPoint<T>>> &points, bool start, bool end)
  {
    // get end derivatives from the original points
    std::vector<TPoint<T>> temp;
    getDerivativesU0(points,temp);
    TPointCurve<T> U0(temp);
    getDerivativesU1(points,temp);
    TPointCurve<T> U1(temp);

    if (!U0.ok() || !U1.ok())
      return;

    // set rows
    for (int i = 0; i <= this->K2; i++)
    {
      T V = T(i) / T(this->K2);
      TPoint<T> der0 = U0.derivative(V,0);
      TPoint<T> der1 = -U1.derivative(V,0);

      std::vector<TPoint<T>> temp;
      this->getRow(i,temp);

      // approximate them along U
      TSplineCurve<T> curve(temp,this->M1,start,end); //!!! exact

      if (start)
        curve.setClampedStart(der0);
      if (end)
        curve.setClampedEnd(der1);

      // copy curve control points into surface control points
      this->setRow(i,curve.controlPoints());
    }
  }

  /** Set V clamping at start / end. points are the original points from where
    the derivatives are taken. */
  void setVClamping(std::vector<std::vector<TPoint<T>>> &points, bool start, bool end)
  {
    std::vector<TPoint<T>> temp;
    getDerivativesV0(points,temp);
    TPointCurve<T> V0(temp);
    getDerivativesV1(points,temp);
    TPointCurve<T> V1(temp);

    if (!V0.ok() || !V1.ok())
      return;

    for (int i = 0; i <= this->K1; i++)
    {
      T U = T(i) / T(this->K1);
      TPoint<T> der0 = V0.derivative(U,0);
      TPoint<T> der1 = -V1.derivative(U,0);

      std::vector<TPoint<T>> temp;
      this->getColumn(i,temp);

      // approximate them along U
      TSplineCurve<T> curve(temp,this->M2,start,end);

      if (start)
        curve.setClampedStart(der0);
      if (end)
        curve.setClampedEnd(der1);

      // copy curve control points into surface control points
      this->setColumn(i,curve.controlPoints());
    }
  }

  /** Calculate control points for points to approximate them. Reliable but the 
    curve is approximate. The curve ends always coincide with the specified points
    and first derivatives at ends are exactly defined by directions between pairs
    of original points at ends. */
  void approximatePoints(std::vector<std::vector<TPoint<T>>> &points)
  {
    this->cpoints.clear();
    this->cpoints.resize((this->K1 + 1) * (this->K2 + 1),TPoint<T>());

    // K1 + M1 + 2 = A + 1 - #knots
    std::vector<T> bs(this->K1 + this->M1 + 2,0.0);
    std::vector<T> bt(this->K2 + this->M2 + 2,0.0);

    for (int i = 0; i <= this->K1; i++)
    {
      T U = T(i) / T(this->K1);
      LIMIT(U,0.0,1.0);

      splineBasis(this->M1 + 1,U,Uknots,bs);

      for (int j = 0; j <= this->K2; j++)
      {
        T V = T(j) / T(this->K2);
        LIMIT(V,0.0,1.0);

        splineBasis(this->M2 + 1,V,Vknots,bt);

        this->cpoints[j * (this->K1 + 1) + i] = TPoint<T>();

        for (int ii = 0; ii <= this->K1; ii++)
        {
          for (int jj = 0; jj <= this->K2; jj++)
          {
            this->cpoints[j * (this->K1 + 1) + i] = this->cpoints[j * (this->K1 + 1) + i] +
              points[jj][ii] * (bs[ii] * bt[jj]);
          }
        }
      }
    }
  }

  /** Calculate control points for points to interpolate them. The surface passes through 
    all specified points.  
    https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/INT-APP/SURF-INT-global.html. */
  bool interpolatePoints(std::vector<std::vector<TPoint<T>>> &points, T parmtolerance = PARM_TOLERANCE, 
    int maxinterpolationpoints = MAX_SPLINEINTERPOLATIONPOINTS)
  {
#if 0 // this works fine, but control points after a correct solution are rubbish at start/end
    if (points.size() > maxinterpolationpoints)
    {
      approximatePoints(points);
      return false;
    }

    points2Dto1D(points,this->cpoints);

    std::vector<T> bs(this->K1 + this->M1 + 2,0.0);
    std::vector<T> bt(this->K2 + this->M2 + 2,0.0);

#if 1
    // banded non-symmetric matrix
    int C = (this->K1 + 1) * (this->K2 + 1);

    // calculate bandwidth, loook like it must be std::max<int>(this->K1 + 1,this->K2 + 1) * 2 + 2
    int halfbandwidth = 0;

    for (int j = 0; j <= this->K2; j++)
    {
      T V = T(j) / T(this->K2);
      LIMIT(V,0.0,1.0);

      splineBasis(this->M2 + 1,V,this->Vknots,bt); 

      for (int i = 0; i <= this->K1; i++)
      {
        T U = T(i) / T(this->K1);
        LIMIT(U,0.0,1.0);

        splineBasis(M1 + 1,U,this->Uknots,bs); 

        int index = j * (this->K1 + 1) + i;

        for (int jj = 0; jj <= this->K2; jj++)
        {
          for (int ii = 0; ii <= this->K1; ii++)
          {

            T product = bs[ii] * bt[jj];

            if (std::abs(product) > TOLERANCE(T))
            {
              int aindex = index * C + jj * (this->K1 + 1) + ii;

              int ai = aindex / C;
              int aj = aindex % C;
              int d = std::abs(ai - aj);

              halfbandwidth = std::max<int>(halfbandwidth,d);
            }
          }
        }
      }
    }

    int bandwidth = halfbandwidth * 2 + 1;
    BandedMatrixSimple<T,TPoint<T>> Ab(C,bandwidth);

    for (int j = 0; j <= this->K2; j++)
    {
      T V = T(j) / T(this->K2);
      LIMIT(V,0.0,1.0);

      splineBasis(this->M2 + 1,V,this->Vknots,bt); 

      for (int i = 0; i <= this->K1; i++)
      {
        T U = T(i) / T(this->K1);
        LIMIT(U,0.0,1.0);

        splineBasis(M1 + 1,U,this->Uknots,bs); 

        int index = j * (this->K1 + 1) + i;

        for (int jj = 0; jj <= this->K2; jj++)
        {
          for (int ii = 0; ii <= this->K1; ii++)
          {
            int aindex = index * C + jj * (this->K1 + 1) + ii;

            int ai = aindex / C;
            int aj = aindex % C;

            T product = bs[ii] * bt[jj];

            if (std::abs(product) > TOLERANCE(T))
              Ab.element(ai,aj) += product;
          }
        }
      }
    }

#ifdef _DEBUG
    // matrix
    std::vector<T> A;
    A.resize(C * C,0.0);

    for (int j = 0; j <= this->K2; j++)
    {
      T V = T(j) / T(this->K2);
      LIMIT(V,0.0,1.0);

      splineBasis(this->M2 + 1,V,this->Vknots,bt); 

      for (int i = 0; i <= this->K1; i++)
      {
        T U = T(i) / T(this->K1);
        LIMIT(U,0.0,1.0);

        splineBasis(M1 + 1,U,this->Uknots,bs); 

        int index = j * (this->K1 + 1) + i;

        for (int jj = 0; jj <= this->K2; jj++)
        {
          for (int ii = 0; ii <= this->K1; ii++)
          {
            int aindex = index * C + jj * (this->K1 + 1) + ii;

            T product = bs[ii] * bt[jj];

            if (std::abs(product) > TOLERANCE(T))
              A[aindex] += product;
          }
        }
      }
    }

    // compare
    for (int i = 0; i < C; i++)
    {
      int i1 = i - halfbandwidth;
      LIMIT_MIN(i1,0);
      int i2 = i + halfbandwidth; 
      LIMIT_MAX(i2,C - 1);

      for (int j = i1; j <= i2; j++)
      {
        T e0 = Ab.element(i,j);
        T e1 = A[C * i + j];
        assert(e0 == e1);

        Ab.element(i,j) = e1;
      }
    }
#endif

    T tolerance = TOLERANCE(T);

    // solve system
    bool ok = Ab.solveSystem(&this->cpoints[0],tolerance);

#else

    // matrix
    std::vector<T> A;
    int C = (this->K1 + 1) * (this->K2 + 1);
    A.resize(C * C,0.0);

    for (int j = 0; j <= this->K2; j++)
    {
      T V = T(j) / T(this->K2);
      LIMIT(V,0.0,1.0);

      splineBasis(this->M2 + 1,V,this->Vknots,bt); 

      for (int i = 0; i <= this->K1; i++)
      {
        T U = T(i) / T(this->K1);
        LIMIT(U,0.0,1.0);

        splineBasis(M1 + 1,U,this->Uknots,bs); 

        int index = j * (this->K1 + 1) + i;

        for (int jj = 0; jj <= this->K2; jj++)
        {
          for (int ii = 0; ii <= this->K1; ii++)
          {
            int aindex = index * C + jj * (this->K1 + 1) + ii;

            A[aindex] += bs[ii] * bt[jj];
          }
        }
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

#if 0
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

#endif

    return true;

#else

    this->cpoints.clear();
    this->cpoints.resize((this->K1 + 1) * (this->K2 + 1),TPoint<T>());

//    T usize = Usize(points);
//    T vsize = Vsize(points);

    // interpolate along shortest size

    if (clampedstartU || clampedendU)
    {
      // interpolate rows from initial column control points
      for (int i = 0; i <= this->K2; i++)
      {
        std::vector<TPoint<T>> temp;
        temp = points[i];

        // approximate them along U
        TSplineCurve<T> curve(temp,this->M1,clampedstartU,clampedendU);

        // copy curve control temp into surface control temp
        this->setRow(i,curve.controlPoints());
      }

      // interpolate column control points 
      for (int i = 0; i <= this->K1; i++)
      {
        // get column 
        std::vector<TPoint<T>> temp;
        for (int k = 0; k <= this->K2; k++)
        {
          temp.push_back(points[k][i]); 
        }

        // make curve approximation along V, we need control points from it
        TSplineCurve<T> curve(temp,this->M2,clampedstartV,clampedendV);

        // copy curve control points into surface control points
        this->setColumn(i,curve.controlPoints());
      }
    } else
    {
      // interpolate column control points 
      for (int i = 0; i <= this->K1; i++)
      {
        // get column 
        std::vector<TPoint<T>> temp;
        for (int k = 0; k <= this->K2; k++)
        {
          temp.push_back(points[k][i]); 
        }

        // make curve approximation along V, we need control points from it
        TSplineCurve<T> curve(temp,this->M2,clampedstartV,clampedendV);

        // copy curve control points into surface control points
        this->setColumn(i,curve.controlPoints());
      }

      // interpolate rows from initial column control points
      for (int i = 0; i <= this->K2; i++)
      {
        std::vector<TPoint<T>> temp;
        temp = points[i];

        // approximate them along U
        TSplineCurve<T> curve(temp,this->M1,clampedstartU,clampedendU);

        // copy curve control temp into surface control temp
        this->setRow(i,curve.controlPoints());
      }
    }

    return true;

#endif
  }

  /** Make knots in both directions. */
  void makeKnots()
  {
    tcad::makeKnots(this->K1,this->M1,Uknots);
    tcad::makeKnots(this->K2,this->M2,Vknots);
  }

  /** Make derivatives from the current surface. */
  void makeDerivatives(TSplineSurface &Uderivative, TSplineSurface &Vderivative)
  {
    // create surfaces 1 order less
    Uderivative.allocate(this->K1 - 1,this->M1 - 1,this->K2 - 1,this->M2 - 1,interpolate);

    Vderivative.allocate(this->K1 - 1,this->M1 - 1,this->K2 - 1,this->M2 - 1,interpolate);

    Uderivative.cpoints.resize(this->K1 * this->K2);
    Vderivative.cpoints.resize(this->K1 * this->K2);

    // prepare control points
    T m1 = T(this->M1);
    T m2 = T(this->M2);
    for (int i = 0; i <= Uderivative.K1; i++)
    {
      T du = Uknots[i + this->M1 + 1] - Uknots[i + 1];

      if (abs(du) < PARM_TOLERANCE)
      {
        du = PARM_TOLERANCE;
      }

      std::vector<TPoint<T>> pv;
      this->getColumn(i,pv);

      for (int j = 0; j <= Uderivative.K2; j++)
      {
        T dv = Vknots[j + this->M2 + 1] - Vknots[j + 1];

        if (abs(dv) < PARM_TOLERANCE)
        {
          dv = PARM_TOLERANCE;
        }

        std::vector<TPoint<T>> pu;
        this->getRow(j,pu);

        Uderivative.cpoints[j * (Uderivative.K1 + 1) + i] = (pu[i + 1] - pu[i]) * (m1 / du);
        Vderivative.cpoints[j * (Vderivative.K1 + 1) + i] = (pv[j + 1] - pv[j]) * (m2 / dv);
      }
    }

    // drop two end points in derivatives
    // drop two end points in derivatives
    removeFirstLast(Uknots,Uderivative.Uknots);
    removeFirstLast(Vknots,Uderivative.Vknots);

    Vderivative.Uknots = Uderivative.Uknots; 
    Vderivative.Vknots = Uderivative.Vknots; 

    // DO NOT call update() here, cpoints are ready
  }

  /** Reverse U. Normal is changed to opposite. */
  virtual void reverseU()
  {
    // reverse rows of control poins
    for (int i = 0; i <= this->K2; i++)
    {
      std::vector<TPoint<T>> row;
      this->getRow(i,row);
      std::reverse(row.begin(),row.end());
      this->setRow(i,row);
    }

    for (int i = 0; i < int(points.size()); i++)
    {
      std::vector<TPoint<T>> row;
      getRow(points,i,row);
      std::reverse(row.begin(),row.end());
      setRow(points,i,row);
    }

    SWAP(bool,clampedstartU,clampedendU);
  }

  /** Reverse V. Normal is changed to opposite. */
  virtual void reverseV()
  {
    // reverse columns of control poins
    for (int i = 0; i <= this->K1; i++)
    {
      std::vector<TPoint<T>> col;
      this->getColumn(i,col);
      std::reverse(col.begin(),col.end());
      this->setColumn(i,col);
    }

    for (int i = 0; i < int(points[0].size()); i++)
    {
      std::vector<TPoint<T>> col;
      getColumn(points,i,col);
      std::reverse(col.begin(),col.end());
      setColumn(points,i,col);
    }

    SWAP(bool,clampedstartV,clampedendV);
  }
  
  /** Clamped start U? */
  bool clampedStartU()
  {
    return clampedstartU;
  }
  
  /** Clamped end U? */
  bool clampedEndU()
  {
    return clampedendU;
  }
  
  /** Clamped start V? */
  bool clampedStartV()
  {
    return clampedstartV;
  }
  
  /** Clamped end V? */
  bool clampedEndV()
  {
    return clampedendV;
  }

protected:

  /** Set main parameters. */
  void allocate(int k1, int m1, int k2, int m2, bool interpolation)
  {
    this->K1 = k1;
    this->M1 = m1;
    this->K2 = k2;
    this->M2 = m2;
    interpolate = interpolation;
  }

public: //!!!!!!!
  // degree in U direction
  int M1 = 0;
  // degree in V direction
  int M2 = 0;

  // U knots
  std::vector<T> Uknots;
  // V knots
  std::vector<T> Vknots;

  // original points
  std::vector<std::vector<TPoint<T>>> points;

private:

  // every row/columns makes a Bezier curve of K1/K2 points;
  // every curve contains (K1/K2 + 1) / 4 Bezier segments

  //// number of columns minus 1
  //int K1 = 0;
  //// number of rows minus 1
  //int K2 = 0;

  //// (K1 + 1) * (K2 + 1) control points, call update() after every change
  //std::vector<TPoint<T>> cpoints;

  // interpolate or approximate
  bool interpolate = false; 

  // clamped (first derivative specified from points) or 
  // natural (second derivative zero) ends
  bool clampedstartU = false;
  bool clampedendU = false;
  bool clampedstartV = false;
  bool clampedendV = false;

  // initialised in update(); they are underinitialised, can be called for 
  // derivative(U,V,0) but not ,1) or ,2).
  TSplineSurface<T> *Uderivative = nullptr;
  TSplineSurface<T> *UUderivative = nullptr;
  TSplineSurface<T> *UVderivative = nullptr;
  TSplineSurface<T> *Vderivative = nullptr;
  TSplineSurface<T> *VVderivative = nullptr;
  TSplineSurface<T> *VUderivative = nullptr;
};

/** Name surfaces with name + number. */
template <class T> void nameSurfaces(std::vector<TSplineSurface<T> *> &surfaces,
  std::string prefix)
{
  for (int i = 0; i < int(surfaces.size()); i++)
  {
    surfaces[i]->name = prefix + " " + to_string(i);
  }
}

/** Close outer boundary. */
template <class T> void closeOuterBoundary(std::vector<TSplineSurface<T> *> &surfaces,
  // surface  // loop 0   // 4 pieces // piece contents
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV)
{
  for (int i = 0; i < int(surfaces.size()); i++)
  {
    std::vector<std::vector<TPoint<T>>> loop;
    surfaces[i]->closeOuterBoundaryLoop(loop);

    boundariesUV.push_back(std::vector<std::vector<std::vector<tcad::TPoint<T>>>>());
    boundariesUV.back().push_back(loop);
  }
}

/** Prepare empty outer boundary. */
template <class T> void prepareOuterBoundary(std::vector<TSplineSurface<T> *> &surfaces,
  // surface  // loop 0   // 4 pieces // piece contents
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV)
{
  for (int i = 0; i < int(surfaces.size()); i++)
  {
    boundariesUV.push_back(std::vector<std::vector<std::vector<tcad::TPoint<T>>>>());
  }
}

}
