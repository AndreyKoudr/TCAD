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

  FFD.h

  Free-form deformation based on TBezierVolume. The algorithm is taken from 
  "S.-M.Hu et al. Direct manipulation of FFD: efficient explicit solutions and 
  decomposible multiple point constraints"

*******************************************************************************/

#pragma once

#include "tbeziervolume.h"
#include "tbasecurve.h"
#include "tbasesurface.h"
#include "tsystems.h"

namespace tcad {

template <class T> class FFD : public TBezierVolume<T> {
public:

  //===== Construction =========================================================

  /** Constructor. */
  FFD() : TBezierVolume<T>() {}

  /** Constrictor for FFD box. numsegs... - number of Bezier segments in each direction. 
    If numsegs... = 1, it means that lattice is a single qubic in this direction (i.e.
    very smooth). */
  FFD(TPoint<T> min, TPoint<T> max, int numsegsU = 1, int numsegsV = 1, int numsegsW = 1) : 
    TBezierVolume<T>(min,max,numsegsU,numsegsV,numsegsW) 
  {
  }

  /** Constructor for a list of surfaces. Surfaces are FFD-distorted right in the 
    constructor. oldpositions points are moved to newpositions to define FFD
    distortion. */
  FFD(std::vector<TBaseSurface<T> *> &surfaces, std::vector<TPoint<T>> &oldpositions, 
    std::vector<TPoint<T>> &newpositions, int numsegsU = 1, int numsegsV = 1, int numsegsW = 1,
    int numpointsU = 2, int numpointsV = 2, int numpointsW = 2,
    T boxextension = 1.01) : TBezierVolume<T>() 
  {
    // get min/max for all surfaces
    TPoint<T> min,max;
    if (calculateMinMax<T>(surfaces,min,max))
    {
      // extend box to include move points
      for (TPoint<T> p : oldpositions)
      {
        min = pointMin<T>(min,p);
        max = pointMax<T>(max,p);
      }

      // extend box to make sure all control points are inside the box
      extendMinMax<T>(min,max,boxextension);

      // init lattice
      this->initBox(min,max,numsegsU,numsegsV,numsegsW);

      // uvw are parametric positions for lattice points to be moved,
      // disp are XYZ displacements
      std::vector<TPoint<T>> uvw,disp;

      assert(oldpositions.size() == newpositions.size());

      // create fine mesh of points
      std::vector<TPoint<T>> points,UVWpoints;
      int k1 = 0;
      int k2 = 0;
      int k3 = 0;
      this->createPoints(points,&UVWpoints,&k1,&k2,&k3,numpointsU,numpointsV,numpointsW); 

      for (int i = 0; i < int(oldpositions.size()); i++)
      {
        TPoint<T> parm = this->findUVWforPoint(points,UVWpoints,k1,k2,k3,oldpositions[i]);

        uvw.push_back(parm);
        disp.push_back(newpositions[i] - oldpositions[i]);
      }

      // distort and smooth lattice, get accuracy
      T acc = distortLattice(uvw,disp); 

      // apply to all control points
      for (int i = 0; i < int(surfaces.size()); i++)
      {
        // if spline, it keeps original points inside
        TSplineSurface<T> *spline = dynamic_cast<TSplineSurface<T> *>(surfaces[i]);

        if (spline)
        {
          for (int j = 0; j < int(spline->points.size()); j++)
          {
            for (int k = 0; k < int(spline->points[j].size()); k++)
            {
              TPoint<T> parm = this->findUVWforPoint(points,UVWpoints,k1,k2,k3,spline->points[j][k]);
              spline->points[j][k] = this->position(parm.X,parm.Y,parm.Z);
            }
          }

          // call update() as original points were transformed
          spline->update();
        } else
        {
          for (int j = 0; j < int(surfaces[i]->controlPoints().size()); j++)
          {
            TPoint<T> parm = this->findUVWforPoint(points,UVWpoints,k1,k2,k3,surfaces[i]->controlPoints()[j]);
            surfaces[i]->controlPoints()[j] = this->position(parm.X,parm.Y,parm.Z);

            // do not call update()
          }
        }
      }
    }
  }

  /** Constructor. */
  FFD(const FFD &other) : TBezierVolume<T>(other)
  {
  }

  /** Destructor. */
  virtual ~FFD() {}

  /** Distort lattice (i.e. volume) : move points at parametric positions uvw by 
    displacements in XYZ. Returns max distance from target points to estimate accuracy. */
  T distortLattice(std::vector<TPoint<T>> &uvw, std::vector<TPoint<T>> &disp, 
    bool Usmooth = true, bool Vsmooth = true, bool Wsmooth = true)
  {
    if (disp.size() == 0)
    {
      return 0.0;
    }

    // assertions
    assert(uvw.size() > 0);
    assert(disp.size() > 0);
    assert(disp.size() == uvw.size());

    int numU = this->K1 + 1;
    int numV = this->K2 + 1;
    int numW = this->K3 + 1;

    // temp
    int numVU = numV * numU;

    // keep target points
    std::vector<TPoint<T>> targetPoints;
    for (int i = 0; i < uvw.size(); i++)
    {
      targetPoints.push_back(this->position(uvw[i].X,uvw[i].Y,uvw[i].Z) + disp[i]);
    }

    // the algorithm is taken from "S.-M.Hu et al. Direct
    // manipulation of FFD: efficient explicit solutions and decomposible multiple
    // point constraints"

    // To test matrices with high precision
    #define MATRIX_FLOAT T

    // matrix R
    std::vector<MATRIX_FLOAT> R(numW * numV * numU * uvw.size());
    // right-hand side
    std::vector<TPoint<T>> D(disp.size());
    for (int i = 0; i < disp.size(); i++)
    {
      D[i].X = static_cast<MATRIX_FLOAT>(disp[i].X);
      D[i].Y = static_cast<MATRIX_FLOAT>(disp[i].Y);
      D[i].Z = static_cast<MATRIX_FLOAT>(disp[i].Z);
    }

    int count = 0;

    for (int i = 0; i < numW; i++)
    {
      for (int j = 0; j < numV; j++)
      {
        for (int k = 0; k < numU; k++)
        {
          for (int q = 0; q < static_cast<int>(uvw.size()); q++)
          {
            // compute basis functions at source point
            int segmentU,segmentV,segmentW;
            std::vector<T> funcsU;
            BezierBasis(uvw[q].X,this->Uknots,funcsU,segmentU);
            std::vector<T> funcsV;
            BezierBasis(uvw[q].Y,this->Vknots,funcsV,segmentV);
            std::vector<T> funcsW;
            BezierBasis(uvw[q].Z,this->Wknots,funcsW,segmentW);

            T f = funcsW[i] * funcsV[j] * funcsU[k];
            R[count++] = static_cast<MATRIX_FLOAT>(f);
          }
        }
      }
    }

    // Solve system
    bool res = solveATASystem<MATRIX_FLOAT,int>(numW * numV * numU,int(uvw.size()),&R[0],&D[0],
      static_cast<MATRIX_FLOAT>(1.0e-20),false);
    assert(res);

    // Move point if only system has been solved
    if (!res) {
    } else
    {
      count = 0;
      for (int i = 0; i < numW; i++)
      {
        for (int j = 0; j < numV; j++)
        {
          for (int k = 0; k < numU; k++)
          {
            TPoint<T> sum(0,0,0);

            for (int q = 0; q < static_cast<int>(uvw.size()); q++)
            {
              // compute basis functions at source point
              int segmentU,segmentV,segmentW;
              std::vector<T> funcsU;
              BezierBasis(uvw[q].X,this->Uknots,funcsU,segmentU);
              std::vector<T> funcsV;
              BezierBasis(uvw[q].Y,this->Vknots,funcsV,segmentV);
              std::vector<T> funcsW;
              BezierBasis(uvw[q].Z,this->Wknots,funcsW,segmentW);

              T f = funcsW[i] * funcsV[j] * funcsU[k];

              sum.X += f * static_cast<T>(D[q].X);
              sum.Y += f * static_cast<T>(D[q].Y);
              sum.Z += f * static_cast<T>(D[q].Z);
            }

            int index = i * numVU + j * numU + k;
            this->cpoints[index] += sum;
          }
        }
      }
    }

    #undef MATRIX_FLOAT

    // smooth lattice
    this->smooth(Usmooth,Vsmooth,Wsmooth);

    // calculate max distance from target points
    T maxDist = 0;
    for (int i = 0; i < uvw.size(); i++)
    {
      TPoint<T> oldPoint = targetPoints[i];
      TPoint<T> newPoint = this->position(uvw[i].X,uvw[i].Y,uvw[i].Z);
      T dist = !(newPoint - oldPoint);
      if (dist > maxDist)
      {
        maxDist = dist;
      }
    }

    return maxDist;
  }
};

}
