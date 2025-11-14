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

  tsplinevolume.h

  Volume on B-spline basis. 

  dimensions : 3 (U,V,W parameters)

*******************************************************************************/

#pragma once

#include "tbasevolume.h"

namespace tcad {

template <class T> class TSplineVolume : public TBaseVolume<T> {
public:

  //===== Construction =========================================================

  /** Constructor. */
  TSplineVolume() : TBaseVolume<T>() {}

  /** Constrictor for box. numpoints... - number of points in each direction. */
  TSplineVolume(TPoint<T> min, TPoint<T> max, int m1 = SPLINE_DEGREE, int m2 = SPLINE_DEGREE, int m3 = SPLINE_DEGREE,
    int numpointsU = MANY_POINTS3D, int numpointsV = MANY_POINTS3D, int numpointsW = MANY_POINTS3D) : TBaseVolume<T>() 
  {
    initBox(min,max,m1,m2,m3,numpointsU,numpointsV,numpointsW);

    update();
  }

  /** Constructor from any TBaseVolume. */
  TSplineVolume(TBaseVolume<T> &other, int k1, int k2, int k3,
    int m1 = SPLINE_DEGREE, int m2 = SPLINE_DEGREE, int m3 = SPLINE_DEGREE) : TBaseVolume<T>() 
  {
    // set main parameters
    allocate(k1,m1,k2,m2,k3,m3);

    // make knots
    makeKnots();

    this->cpoints.clear();
    this->cpoints.resize((this->K1 + 1) * (this->K2 + 1) * (this->K3 + 1),TPoint<T>());

    // K1 + M1 + 2 = A + 1 - #knots
    std::vector<T> bs(this->K1 + this->M1 + 2,0.0);
    std::vector<T> bt(this->K2 + this->M2 + 2,0.0);
    std::vector<T> bx(this->K3 + this->M3 + 2,0.0);

    for (int i = 0; i <= this->K1; i++)
    {
      T U = T(i) / T(this->K1);
      LIMIT(U,0.0,1.0);

      splineBasis(this->M1 + 1,U,this->Uknots,bs);

      for (int j = 0; j <= this->K2; j++)
      {
        T V = T(j) / T(this->K2);
        LIMIT(V,0.0,1.0);

        splineBasis(this->M2 + 1,V,this->Vknots,bt);

        for (int k = 0; k <= this->K3; k++)
        {
          T W = T(k) / T(this->K3);
          LIMIT(W,0.0,1.0);

          splineBasis(this->M3 + 1,W,this->Wknots,bx);

          int index = this->getIndex(i,j,k);
          this->cpoints[index] = TPoint<T>();

          TPoint<T> pos = other.position(U,V,W);

          for (int ii = 0; ii <= this->K1; ii++)
          {
            for (int jj = 0; jj <= this->K2; jj++)
            {
              for (int kk = 0; kk <= this->K3; kk++)
              {
                this->cpoints[index] = this->cpoints[index] + pos * (bs[ii] * bt[jj] * bx[kk]);
              }
            }
          }
        }
      }
    }

    update();
  }

  /** Constructor. */
  TSplineVolume(const TSplineVolume &other) : TBaseVolume<T>()
  {
    this->K1 = other.K1;
    this->K2 = other.K2;
    this->K3 = other.K3;
    this->M1 = other.M1;
    this->M2 = other.M2;
    this->M3 = other.M3;
    this->cpoints = other.cpoints;
    this->Uknots = other.Uknots;
    this->Vknots = other.Vknots;
    this->Wknots = other.Wknots;

    update();
  }

  /** Assignment operator. */
  TSplineVolume &operator = (const TSplineVolume &other)  
  {
    this->K1 = other.K1;
    this->K2 = other.K2;
    this->K3 = other.K3;
    this->M1 = other.M1;
    this->M2 = other.M2;
    this->M3 = other.M3;
    this->cpoints = other.cpoints;
    this->Uknots = other.Uknots;
    this->Vknots = other.Vknots;
    this->Wknots = other.Wknots;

    update();

    return *this;
  }

  /** Initialize from box. */
  void initBox(TPoint<T> min, TPoint<T> max, int m1 = SPLINE_DEGREE, int m2 = SPLINE_DEGREE, int m3 = SPLINE_DEGREE, 
    int numpointsU = MANY_POINTS3D, int numpointsV = MANY_POINTS3D, int numpointsW = MANY_POINTS3D) 
  {
    this->cpoints.clear();

    allocate(numpointsU - 1,m1,numpointsV - 1,m2,numpointsW - 1,m3);

    TPoint<T> d = max - min;
    T DZ = d.Z / T(numpointsW - 1);
    T DY = d.Y / T(numpointsV - 1);
    T DX = d.X / T(numpointsU - 1);

    for (int i = 0; i < numpointsW; i++)
    {
      T Z = min.Z + T(i) * DZ;
      for (int j = 0; j < numpointsV; j++)
      {
        T Y = min.Y + T(j) * DY;
        for (int k = 0; k < numpointsU; k++)
        {
          T X = min.X + T(k) * DX;
          this->cpoints.push_back(TPoint<T>(X,Y,Z));
        }
      }
    }

    // make knots between segments
    makeKnots();

    // make spline control points from simply regular mesh
    makeSplineControlPoints(true);
  }

  /** Destructor. */
  virtual ~TSplineVolume() 
  {
    DELETE_CLASS(Uderivative);
    DELETE_CLASS(Vderivative);
    DELETE_CLASS(Wderivative);

    DELETE_CLASS(UUderivative);
    DELETE_CLASS(UVderivative);
    DELETE_CLASS(UWderivative);

    DELETE_CLASS(VUderivative);
    DELETE_CLASS(VVderivative);
    DELETE_CLASS(VWderivative);

    DELETE_CLASS(WUderivative);
    DELETE_CLASS(WVderivative);
    DELETE_CLASS(WWderivative);
  }

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

    TPoint<T> result;

    // only 0-th derivative done yet
    if (k == 0)
    {
      std::vector<T> bs(this->K1 + this->M1 + 2,0.0);
      std::vector<T> bt(this->K2 + this->M2 + 2,0.0);
      std::vector<T> bx(this->K3 + this->M3 + 2,0.0);

      splineBasis(this->M1 + 1,U,this->Uknots,bs);
      splineBasis(this->M2 + 1,V,this->Vknots,bt);
      splineBasis(this->M3 + 1,W,this->Wknots,bx);

      for (int i = 0; i <= this->K1; i++)
      {
        for (int j = 0; j <= this->K2; j++)
        {
          for (int k = 0; k <= this->K3; k++)
          {
            int index = this->getIndex(i,j,k);
            result = result + this->cpoints[index] * (bs[i] * bt[j] * bx[k]);
          }
        }
      }
    } else if (k == 1)
    {
      if (onparameter == PARAMETER_U)
      {
        return Uderivative->position(U,V,W);
      } else if (onparameter == PARAMETER_V)
      {
        return Vderivative->position(U,V,W);
      } else if (onparameter == PARAMETER_W)
      {
        return Wderivative->position(U,V,W);
      }
    } else if (k == 2) 
    {
      if (onparameter == PARAMETER_UU)
      {
        return UUderivative->position(U,V,W);
      } else if (onparameter == PARAMETER_UV)
      {
        return UVderivative->position(U,V,W);
      } else if (onparameter == PARAMETER_UW)
      {
        return UWderivative->position(U,V,W);

      } else if (onparameter == PARAMETER_VU)
      {
        return VUderivative->position(U,V,W);
      } else if (onparameter == PARAMETER_VV)
      {
        return VVderivative->position(U,V,W);
      } else if (onparameter == PARAMETER_VW)
      {
        return VWderivative->position(U,V,W);

      } else if (onparameter == PARAMETER_WU)
      {
        return WUderivative->position(U,V,W);
      } else if (onparameter == PARAMETER_WV)
      {
        return WVderivative->position(U,V,W);
      } else if (onparameter == PARAMETER_WW)
      {
        return WWderivative->position(U,V,W);
      }
    }

    return result;
  }

  /** Update after any change in control points. */
  virtual void update()
  {
    // prepare two derivatives
    if (Uderivative == nullptr)
      Uderivative = new TSplineVolume();
    if (Vderivative == nullptr)
      Vderivative = new TSplineVolume();
    if (Wderivative == nullptr)
      Wderivative = new TSplineVolume();

    if (UUderivative == nullptr)
      UUderivative = new TSplineVolume();
    if (UVderivative == nullptr)
      UVderivative = new TSplineVolume();
    if (UWderivative == nullptr)
      UWderivative = new TSplineVolume();

    if (VUderivative == nullptr)
      VUderivative = new TSplineVolume();
    if (VVderivative == nullptr)
      VVderivative = new TSplineVolume();
    if (VWderivative == nullptr)
      VWderivative = new TSplineVolume();

    if (WUderivative == nullptr)
      WUderivative = new TSplineVolume();
    if (WVderivative == nullptr)
      WVderivative = new TSplineVolume();
    if (WWderivative == nullptr)
      WWderivative = new TSplineVolume();

    makeDerivatives(*Uderivative,*Vderivative,*Wderivative);

    Uderivative->makeDerivatives(*UUderivative,*UVderivative,*UWderivative);
    Vderivative->makeDerivatives(*VUderivative,*VVderivative,*VWderivative);
    Wderivative->makeDerivatives(*WUderivative,*WVderivative,*WWderivative);
  }

  /** Make U,V,W derivatives from the current volume. */
  void makeDerivatives(TSplineVolume &Uderivative, TSplineVolume &Vderivative, TSplineVolume &Wderivative)
  {
    // create surfaces 1 order less
    Uderivative.allocate(this->K1 - 1,this->M1 - 1,this->K2 - 1,this->M2 - 1,this->K3 - 1,this->M3 - 1);
    Vderivative.allocate(this->K1 - 1,this->M1 - 1,this->K2 - 1,this->M2 - 1,this->K3 - 1,this->M3 - 1);
    Wderivative.allocate(this->K1 - 1,this->M1 - 1,this->K2 - 1,this->M2 - 1,this->K3 - 1,this->M3 - 1);

    Uderivative.cpoints.resize(this->K1 * this->K2 * this->K3);
    Vderivative.cpoints.resize(this->K1 * this->K2 * this->K3);
    Wderivative.cpoints.resize(this->K1 * this->K2 * this->K3);

    // prepare control points
    T m1 = T(this->M1);
    T m2 = T(this->M2);
    T m3 = T(this->M3);
    for (int i = 0; i <= Uderivative.K1; i++)
    {
      T du = this->Uknots[i + this->M1 + 1] - this->Uknots[i + 1];

      if (abs(du) < PARM_TOLERANCE)
      {
        du = PARM_TOLERANCE;
      }

      for (int j = 0; j <= Uderivative.K2; j++)
      {
        T dv = this->Vknots[j + this->M2 + 1] - this->Vknots[j + 1];

        if (abs(dv) < PARM_TOLERANCE)
        {
          dv = PARM_TOLERANCE;
        }

        for (int k = 0; k <= Uderivative.K3; k++)
        {
          T dw = this->Wknots[k + this->M3 + 1] - this->Wknots[k + 1];

          if (abs(dw) < PARM_TOLERANCE)
          {
            dw = PARM_TOLERANCE;
          }

          std::vector<TPoint<T>> row;
          this->getRow(j,k,row);

          std::vector<TPoint<T>> col;
          this->getColumn(i,k,col);

          std::vector<TPoint<T>> lay;
          this->getLayer(i,j,lay);

          int index = Uderivative.getIndex(i,j,k);

          Uderivative.cpoints[index] = (row[i + 1] - row[i]) * (m1 / du);
          Vderivative.cpoints[index] = (col[j + 1] - col[j]) * (m2 / dv);
          Wderivative.cpoints[index] = (lay[k + 1] - lay[k]) * (m3 / dw);
        }
      }
    }

    // drop two end points in derivatives
    removeFirstLast(this->Uknots,Uderivative.Uknots);
    removeFirstLast(this->Vknots,Uderivative.Vknots);
    removeFirstLast(this->Wknots,Uderivative.Wknots);

    Vderivative.Uknots = Wderivative.Uknots = Uderivative.Uknots; 
    Vderivative.Vknots = Wderivative.Vknots = Uderivative.Vknots; 
    Vderivative.Wknots = Wderivative.Wknots = Uderivative.Wknots; 

    // DO NOT call update() here, cpoints are ready
  }

  /** Improve positions of spline control points to achieve uniform parameterisation. */
  void makeSplineControlPoints(bool uniformparms, 
    CurveEndType Ustart = END_FREE, CurveEndType Uend = END_FREE,
    CurveEndType Vstart = END_FREE, CurveEndType Vend = END_FREE,
    CurveEndType Wstart = END_FREE, CurveEndType Wend = END_FREE)
  {
    for (int j = 0; j <= this->K2; j++)
    {
      for (int k = 0; k <= this->K3; k++)
      {
        std::vector<TPoint<T>> row;
        this->getRow(j,k,row);
        makeSplineCurve(row,this->K1,this->M1,row,Ustart,Uend,!uniformparms);
        this->setRow(j,k,row);
      }
    }

    for (int i = 0; i <= this->K1; i++)
    {
      for (int k = 0; k <= this->K3; k++)
      {
        std::vector<TPoint<T>> col;
        this->getColumn(i,k,col);
        makeSplineCurve(col,this->K2,this->M2,col,Vstart,Vend,!uniformparms);
        this->setColumn(i,k,col);
      }
    }

    for (int i = 0; i <= this->K1; i++)
    {
      for (int j = 0; j <= this->K2; j++)
      {
        std::vector<TPoint<T>> lay;
        this->getLayer(i,j,lay);
        makeSplineCurve(lay,this->K3,this->M3,lay,Wstart,Wend,!uniformparms);
        this->setLayer(i,j,lay);
      }
    }
  }

/**
  Faces (all normals point outside) :
              7-------------------6 
             /|        /         /|    
            / |      (5)->      / |     
           /  |        |       /  |          
          4---------<-(3)-----5   |    
          |   |               |   |     
          |(0)|     ^         |(1)|    
          | | |     |         ||/ |
          |/  3----(2)->------|---2 
          |  /                |  / 
  WZk     ^ /VYj       (4)->  | /
  layers  |/ columns   /      |/
          0-->----------------1
          UXi rows

  cpoints are ALL control points, (K1 + 1) * (K2 + 1) * (K3 + 1)
  Knots are supposed to be standard set, uniform with no multiplicity.
*/
  TSplineSurface<T> *getFace(int faceno, CurveEndType startU = END_FREE, CurveEndType endU = END_FREE,
    CurveEndType startV = END_FREE, CurveEndType endV = END_FREE)
  {
    std::vector<std::vector<TPoint<T>>> points;
    tcad::getFace<T>(this->cpoints,faceno,this->K1,this->K2,this->K3,points);

    // take spline degrees from volume
    int m1 = 0;
    int m2 = 0;
    if (faceno == 0)
    {
      m1 = this->M2; 
      m2 = this->M3;
    } else if (faceno == 1)
    {
      m1 = this->M2; 
      m2 = this->M3;
    } else if (faceno == 2)
    {
      m1 = this->M1; 
      m2 = this->M3;
    } else if (faceno == 3)
    {
      m1 = this->M1; 
      m2 = this->M3;
    } else if (faceno == 4)
    {
      m1 = this->M1; 
      m2 = this->M2;
    } else if (faceno == 5)
    {
      m1 = this->M1; 
      m2 = this->M2;
    }

    TSplineSurface<T> *face = new TSplineSurface<T>(points,m1,m2,startU,endU,startV,endV);

    return face;
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

  // spline degree
  int M1 = 0;
  int M2 = 0;
  int M3 = 0;

protected:

  /** Set main parameters. */
  void allocate(int k1, int m1, int k2, int m2, int k3, int m3)
  {
    this->K1 = k1;
    this->M1 = m1;
    this->K2 = k2;
    this->M2 = m2;
    this->K3 = k3;
    this->M3 = m3;
  }

  //// (K1 + 1) * (K2 + 1) * (K3 + 1) control points, call update() after every change
  //std::vector<TPoint<T>> cpoints;

  /** Make all knots. */
  void makeKnots()
  {
    // make uniformly spaced knots
    tcad::makeKnots(this->K1,this->M1,this->Uknots);
    tcad::makeKnots(this->K2,this->M2,this->Vknots);
    tcad::makeKnots(this->K3,this->M3,this->Wknots);
  }

private:
  // initialised in update(); they are underinitialised, can be called for 
  // derivative(U,W,V,0) but not ,1) or ,2).
  TSplineVolume<T> *Uderivative = nullptr;
  TSplineVolume<T> *Vderivative = nullptr;
  TSplineVolume<T> *Wderivative = nullptr;

  TSplineVolume<T> *UUderivative = nullptr;
  TSplineVolume<T> *UVderivative = nullptr;
  TSplineVolume<T> *UWderivative = nullptr;

  TSplineVolume<T> *VUderivative = nullptr;
  TSplineVolume<T> *VVderivative = nullptr;
  TSplineVolume<T> *VWderivative = nullptr;

  TSplineVolume<T> *WUderivative = nullptr;
  TSplineVolume<T> *WVderivative = nullptr;
  TSplineVolume<T> *WWderivative = nullptr;
};

}
