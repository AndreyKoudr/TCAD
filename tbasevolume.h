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

  tbasevolume.h

  Basic abstract class for volumes

  dimensions : 3 (U,V,W parameters)

*******************************************************************************/

#pragma once

#include "tbasics.h"
#include "tpoint.h"
#include "ttransform.h"

namespace tcad {

template <class T> class TBaseVolume {
public:

  //===== Construction =========================================================

  /** Constructor. */
  TBaseVolume() {}

  /** Destructor. */
  virtual ~TBaseVolume() {}

  //===== Abstract =============================================================

  /** Get k-th derivative on parameter U,V [0..1]. 0-derivative is position(U,V) (see below). */
  virtual TPoint<T> derivative(T U, T V, T W, Parameter onparameter, int k) = 0;

  /** Update after any change in control points. */
  virtual void update() = 0;

  //===== Operations ===========================================================

  /** Same as 0-th derivative*/
  virtual TPoint<T> position(T U, T V, T W)
  {
    return this->derivative(U,V,W,PARAMETER_ANY,0);
  }

  /** Get control points. */
  virtual std::vector<TPoint<T>> &controlPoints()
  {
    return cpoints;
  }

  /** Make transform. */
  virtual void makeTransform(TTransform<T> *transform)
  {
    for (auto &p : cpoints)
    {
      p = transform->applyTransform(p);
    }

    // call virtual update()
    this->update();
  }

  /** Are equal? */
  virtual bool equal(TBaseVolume &other, T tolerance)
  {
    T diff = difference(cpoints,other.cpoints);
    return (diff >= 0.0 && diff < tolerance);
  }

  /** Location index for a control point. */
  int getIndex(int i, int j, int k)
  {
    return tcad::getIndex<T>(K1,K2,K3,i,j,k);
  }

  /** Get row of control points (numbers). */
  void getULinePoints(int j, int k, std::vector<int> &points)
  {
    tcad::getULinePoints<T>(K1,K2,K3,j,k,points);
  }

  /** Get column of control points (numbers). */
  void getVLinePoints(int i, int k, std::vector<int> &points)
  {
    tcad::getVLinePoints<T>(K1,K2,K3,i,k,points);
  }

  /** Get layer of control points (numbers). */
  void getWLinePoints(int i, int j, std::vector<int> &points)
  {
    tcad::getWLinePoints<T>(K1,K2,K3,i,j,points);
  }

  /** Get row of control points. */
  void getRow(int j, int k, std::vector<TPoint<T>> &points)
  {
    tcad::getRow<T>(this->cpoints,K1,K2,K3,j,k,points);
  }

  /** Set row of control points. */
  void setRow(int j, int k, std::vector<TPoint<T>> &points)
  {
    tcad::setRow<T>(this->cpoints,K1,K2,K3,j,k,points);
  }

  /** Get column of control points. */
  void getColumn(int i, int k, std::vector<TPoint<T>> &points)
  {
    tcad::getColumn<T>(this->cpoints,K1,K2,K3,i,k,points);
  }

  /** Set column of control points. */
  void setColumn(int i, int k, std::vector<TPoint<T>> &points)
  {
    tcad::setColumn<T>(this->cpoints,K1,K2,K3,i,k,points);
  }

  /** Get layer of control points. */
  void getLayer(int i, int j, std::vector<TPoint<T>> &points)
  {
    tcad::getLayer<T>(this->cpoints,K1,K2,K3,i,j,points);
  }

  /** Set layer of control points. */
  void setLayer(int i, int j, std::vector<TPoint<T>> &points)
  {
    tcad::setLayer<T>(this->cpoints,K1,K2,K3,i,j,points);
  }

  /** Generate a uniform set of actual points inside the volume. 

    A set of according parameter values is generated as well.
    Set refine... to 0.5 at a corresponding end to refine,
    1.0 has no effect.

    UVWpoints, if not null, contain U,V,W parameter values in X,Y
  */
  virtual void createPoints(std::vector<TPoint<T>> &points, 
    std::vector<TPoint<T>> *UVWpoints,
    int *k1 = nullptr, int *k2 = nullptr, int *k3 = nullptr,
    int numpointsU = MANY_POINTS3D, int numpointsV = MANY_POINTS3D, int numpointsW = MANY_POINTS3D,
    T refinestartU = 1.0, T refineendU = 1.0, 
    T refinestartV = 1.0, T refineendV = 1.0,
    T refinestartW = 1.0, T refineendW = 1.0)
  {
    points.clear();

    for (int i = 0; i < numpointsW; i++)
    {
      T W = T(i) / T(numpointsW - 1);

      W = refineParameter(W,refinestartW,refineendW);

      for (int j = 0; j < numpointsV; j++)
      {
        T V = T(j) / T(numpointsV - 1);

        V = refineParameter(V,refinestartV,refineendV);

        for (int k = 0; k < numpointsU; k++)
        {
          T U = T(k) / T(numpointsU - 1);

          U = refineParameter(U,refinestartU,refineendU);

          TPoint<T> p = this->position(U,V,W);
          points.push_back(p);

          if (UVWpoints)
            UVWpoints->push_back(TPoint<T>(U,V,W));
        }
      }
    }

    if (k1)
      *k1 = numpointsU - 1;
    if (k2)
      *k2 = numpointsV - 1;
    if (k3)
      *k3 = numpointsW - 1;
  }

  /** Find values of parameters U,V,W for a point inside the volume. */
  virtual TPoint<T> findUVWforPoint(std::vector<TPoint<T>> &points,
    std::vector<TPoint<T>> &UVWpoints, int k1, int k2, int k3, TPoint<T> p) 
  {
    // temp
    TPoint<T> proj;
    int seg = 0;
    T u = 0.0;

    // look for closest point on rows
    T U = -1.0;
    T mindist = std::numeric_limits<T>::max();
    for (int j = 0; j < k2 + 1; j++)
    {
      for (int k = 0; k < k3 + 1; k++)
      {
        std::vector<TPoint<T>> row;
        tcad::getRow(points,k1,k2,k3,j,k,row);

        if (projectPointOnPoints(row,p,proj,&seg,&u))
        {
          T dist = !(p - proj);
          if (dist < mindist)
          {
            TPoint<T> p0 = row[seg];
            TPoint<T> p1 = row[seg + 1];
            int index0 = tcad::getIndex<T>(k1,k2,k3,seg,j,k);
            int index1 = tcad::getIndex<T>(k1,k2,k3,seg + 1,j,k);
            U = UVWpoints[index0].X + (UVWpoints[index1].X - UVWpoints[index0].X) * u;
            mindist = dist;
          }
        }
      }
    }

    // look for closest point on columns
    T V = -1.0;
    mindist = std::numeric_limits<T>::max();
    for (int i = 0; i < k1 + 1; i++)
    {
      for (int k = 0; k < k3 + 1; k++)
      {
        std::vector<TPoint<T>> col;
        tcad::getColumn(points,k1,k2,k3,i,k,col);

        if (projectPointOnPoints(col,p,proj,&seg,&u))
        {
          T dist = !(p - proj);
          if (dist < mindist)
          {
            TPoint<T> p0 = col[seg];
            TPoint<T> p1 = col[seg + 1];
            int index0 = tcad::getIndex<T>(k1,k2,k3,i,seg,k);
            int index1 = tcad::getIndex<T>(k1,k2,k3,i,seg + 1,k);
            V = UVWpoints[index0].Y + (UVWpoints[index1].Y - UVWpoints[index0].Y) * u;
            mindist = dist;
          }
        }
      }
    }

    // look for closest point on columns
    T W = -1.0;
    mindist = std::numeric_limits<T>::max();
    for (int i = 0; i < k1 + 1; i++)
    {
      for (int j = 0; j < k2 + 1; j++)
      {
        std::vector<TPoint<T>> lay;
        tcad::getLayer(points,k1,k2,k3,i,j,lay);

        if (projectPointOnPoints(lay,p,proj,&seg,&u))
        {
          T dist = !(p - proj);
          if (dist < mindist)
          {
            TPoint<T> p0 = lay[seg];
            TPoint<T> p1 = lay[seg + 1];
            int index0 = tcad::getIndex<T>(k1,k2,k3,i,j,seg);
            int index1 = tcad::getIndex<T>(k1,k2,k3,i,j,seg + 1);
            W = UVWpoints[index0].Z + (UVWpoints[index1].Z - UVWpoints[index0].Z) * u;
            mindist = dist;
          }
        }
      }
    }

    return TPoint<T>(U,V,W);
  }

  /** Calculate min/max. */
  virtual bool calculateMinMax(TPoint<T> *min, TPoint<T> *max, TPoint<T> *imin = nullptr, 
    TPoint<T> *imax = nullptr, 
    int numpointsU = MANY_POINTS3D, int numpointsV = MANY_POINTS3D, int numpointsW = MANY_POINTS3D)
  {
    std::vector<TPoint<T>> points;
    createPoints(points,nullptr,nullptr,nullptr,nullptr,numpointsU,numpointsV,numpointsW);

    return tcad::calculateMinMax(points,min,max,imin,imax);
  }

public:

/*
            ---------------------
         /                      /|
         ----------------------  |
        |                      | |
  V (j) |                      | |
  ^     |                      | |
  |     |                      |/
         ----------------------
 /   -> U (i)
W (k)

*/

  std::vector<T> knotsU;     
  std::vector<T> knotsV;   
  std::vector<T> knotsW;    

  // number of columns minus 1
  int K1 = 0;

  // number of rows munus 1
  int K2 = 0;

  // number of levels minus 1
  int K3 = 0;

protected:

  // (K1 + 1) * (K2 + 1) * (K3 + 1) control points, call update() after every change
  std::vector<TPoint<T>> cpoints;

};

}
