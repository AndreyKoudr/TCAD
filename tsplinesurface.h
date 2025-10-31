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

  Spline surface

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

  /** Constructor. Every points[i] is a row of points (U-changing). */
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
        return Uderivative->derivative(U,V,PARAMETER_ANY,0);
      } else if (onparameter == PARAMETER_V)
      {
        return Vderivative->derivative(U,V,PARAMETER_ANY,0);
      } else
      {
        return TPoint<T>();
      }
    } else if (k == 2)
    {
      if (onparameter == PARAMETER_U)
      {
        return UUderivative->derivative(U,V,PARAMETER_ANY,0);
      } else if (onparameter == PARAMETER_V)
      {
        return VVderivative->derivative(U,V,PARAMETER_ANY,0);
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
      // clamping is done inside, do not do it twice
      interpolatePoints(kpoints);
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

  /** Set U clamping at start / end. points are the original points from where
    the dirivatives are taken. */
  void setUClamping(std::vector<std::vector<TPoint<T>>> &points, bool start, bool end)
  {
    // get end derivatives from the original points
    std::vector<TPoint<T>> temp;
    getDerivativesU0(points,temp);
    TPointCurve<T> U0(temp);
    getDerivativesU1(points,temp);
    TPointCurve<T> U1(temp);

    // set rows
    for (int i = 0; i <= this->K2; i++)
    {
      T V = T(i) / T(this->K2);
      TPoint<T> der0 = U0.derivative(V,0);
      TPoint<T> der1 = -U1.derivative(V,0);

      std::vector<TPoint<T>> temp;
      this->getRow(i,temp);

      // approximate them along U
      TSplineCurve<T> curve(temp,this->M1,start,end);

      if (start)
        curve.setClampedStart(der0);
      if (end)
        curve.setClampedEnd(der1);

      // copy curve control points into surface control points
      this->setRow(i,curve.controlPoints());
    }
  }

  /** Set V clamping at start / end. points are the original points from where
    the dirivatives are taken. */
  void setVClamping(std::vector<std::vector<TPoint<T>>> &points, bool start, bool end)
  {
    std::vector<TPoint<T>> temp;
    getDerivativesV0(points,temp);
    TPointCurve<T> V0(temp);
    getDerivativesV1(points,temp);
    TPointCurve<T> V1(temp);

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
  void interpolatePoints(std::vector<std::vector<TPoint<T>>> &points)
  {
    bool ok = true;

    this->cpoints.clear();
    this->cpoints.resize((this->K1 + 1) * (this->K2 + 1),TPoint<T>());

    // K1 + M1 + 2 = A + 1 - #knots
    std::vector<T> bs(this->K1 + this->M1 + 2,0.0);
    std::vector<T> bt(this->K2 + this->M2 + 2,0.0);

    // step 1 : interpolate rows from initial column control points
    for (int i = 0; i <= this->K2; i++)
    {
      std::vector<TPoint<T>> temp;
      temp = points[i];

      // approximate them along U
      TSplineCurve<T> curve(temp,this->M1,clampedstartU,clampedendU);

      // copy curve control temp into surface control temp
      this->setRow(i,curve.controlPoints());
    }

    // step 2 : interpolate column control points 
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
    Uderivative.makeKnots();

    Vderivative.allocate(this->K1 - 1,this->M1 - 1,this->K2 - 1,this->M2 - 1,interpolate);
    Vderivative.makeKnots();

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

    // correct knots
    for (int i = 0; i <= Uderivative.K1; i++)
    {
      Uderivative.Uknots[i] = Uknots[i + 1];
    }

    for (int i = 0; i <= Vderivative.K2; i++)
    {
      Vderivative.Vknots[i] = Vknots[i + 1];
    }

    // DO NOT call update() here, cpoints are ready
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

private:

  // every row/columns makes a Bezier curve of K1/K2 points;
  // every curve contains (K1/K2 + 1) / 4 Bezier segments

  //// number of columns minus 1
  //int K1 = 0;
  //// number of rows munus 1
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

  // original points
  std::vector<std::vector<TPoint<T>>> points;

  // initialised in update(); they are underinitialised, can be called for 
  // derivative(U,0) but not ,1) or ,2).
  TSplineSurface<T> *Uderivative = nullptr;
  TSplineSurface<T> *UUderivative = nullptr;
  TSplineSurface<T> *UVderivative = nullptr;
  TSplineSurface<T> *Vderivative = nullptr;
  TSplineSurface<T> *VVderivative = nullptr;
  TSplineSurface<T> *VUderivative = nullptr;
};

}
