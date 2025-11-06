/*
BSD 2-Clause License

Copyright (c) 2020, Andrey Kudryavtsev (andrewkoudr@hotmail.com)
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

#pragma once
                              
#include "tpoint.h"

#include <limits>
#include <algorithm>
#include <array>
#include <assert.h>

using namespace std;
using namespace tcad;

/** Very important class to define integer positions of cells and nodes within octree. */
#define IPosition TPoint<LINT>

/** 
  Background parameters for a balanced octree.

  Background is a set of largest cells obtained by uniform division of the whole cuboid region.
  Background cells are cells of level 0. Every cell of level 0 can be subdivided (refined)
to make a hierarchy of cells inside a background cell. 

  A search of every 3D point inside octree is two-step :
  (1) which background cell? 
  (2) starting from this background cell, search within a hierarchy
Both operations are fast. This is the way an octree provides a kind of indexing for 3D points.

  T is float type (float or double) for real coordinates.

  It is used to speedup triangle/triangle intersections by spacial partitioning.
*/

template <typename T> class OBackground {
public:

  OBackground() = delete;

  /** Contructor */
  OBackground(const TPoint<T> &min, const TPoint<T> &max, 
    LINT Im, LINT Jm, LINT Km, const int maxOctreeLevel);

  /** Return 8 boundary nodes */
  std::array<TPoint<T>,8> boundaryNodes();

  /** Convert back cell I,J,K (not normal integer!) position into back cell number */
  LINT positionToBackCellIndex(const IPosition &position);

  /** Convert back cell number into I,J,K (not normal integer!) position */
  IPosition backCellIndexToPosition(const LINT index);

  /** Find cell number, -1 in failure. */
  LINT findCell(TPoint<T> position);

  /** Number of cells. */
  LINT numBackgroundCells();

  /** Num background cells in I,J,K directions. */
  IPosition IJKnumCells = IPosition(0,0,0);

  /** Size of whole region in I,J,K directions in integer units. */
  IPosition IJKintSizes = IPosition(0,0,0);

  /** Box X,Y,Z min coordinate. */
  TPoint<T> boxMin = TPoint<T>(0.0,0.0,0.0);
  /** Box X,Y,Z max coordinate. */
  TPoint<T> boxMax = TPoint<T>(0.0,0.0,0.0);

  /** Max octree level. */
  int maxLevel = 0;

  /** Real sizes of integer unit in 3 directions. */
  TPoint<T> intUnitSizes = TPoint<T>(0.0,0.0,0.0);

  /** Cell size in I,J,K directions. */
  TPoint<T> cellSize = TPoint<T>(0.0,0.0,0.0);

  /** Max box size. */
  T scale = 0.0;
};

template <typename T> OBackground<T>::OBackground(const TPoint<T> &min, const TPoint<T> &max, 
  LINT Im, LINT Jm, LINT Km, const int maxOctreeLevel)
{
  LIMIT_MIN(Im,1);
  LIMIT_MIN(Jm,1);
  LIMIT_MIN(Km,1);

  //assert(Im > 0);
  //assert(Jm > 0);
  //assert(Km > 0);
  assert(max.XYZ[0] > min.XYZ[0]);
  assert(max.XYZ[1] > min.XYZ[1]);
  assert(max.XYZ[2] > min.XYZ[2]);

  IJKnumCells = IPosition(Im,Jm,Km);
  maxLevel = maxOctreeLevel;

  LINT l = 1i64 << maxLevel;
  IJKintSizes = IJKnumCells * l;

  assert(IJKintSizes.XYZ[0] > 0);
  assert(IJKintSizes.XYZ[1] > 0);
  assert(IJKintSizes.XYZ[2] > 0);

  TPoint<T> d = max - min;

  cellSize.XYZ[0] = d.XYZ[0] / static_cast<T>(Im);
  cellSize.XYZ[1] = d.XYZ[1] / static_cast<T>(Jm);
  cellSize.XYZ[2] = d.XYZ[2] / static_cast<T>(Km);

#ifndef NDEBUG
  T tolerance = std::numeric_limits<T>::epsilon() * static_cast<T>(10.0);

  assert(cellSize.XYZ[0] > tolerance);
  assert(cellSize.XYZ[1] > tolerance);
  assert(cellSize.XYZ[2] > tolerance);
#endif
  
  boxMin = min;
  boxMax = max;

  intUnitSizes.XYZ[0] = cellSize.XYZ[0] * pow(0.5,maxLevel);
  intUnitSizes.XYZ[1] = cellSize.XYZ[1] * pow(0.5,maxLevel);
  intUnitSizes.XYZ[2] = cellSize.XYZ[2] * pow(0.5,maxLevel);

  scale = std::max<T>(d.XYZ[0],std::max<T>(d.XYZ[1],d.XYZ[2]));
}

template <typename T> std::array<TPoint<T>,8> OBackground<T>::boundaryNodes()
{
  std::array<TPoint<T>,8> nodes = {
    TPoint<T>(boxMin.XYZ[0],boxMin.XYZ[1],boxMin.XYZ[2]),
    TPoint<T>(boxMax.XYZ[0],boxMin.XYZ[1],boxMin.XYZ[2]),
    TPoint<T>(boxMax.XYZ[0],boxMax.XYZ[1],boxMin.XYZ[2]),
    TPoint<T>(boxMin.XYZ[0],boxMax.XYZ[1],boxMin.XYZ[2]),
    TPoint<T>(boxMin.XYZ[0],boxMin.XYZ[1],boxMax.XYZ[2]),
    TPoint<T>(boxMax.XYZ[0],boxMin.XYZ[1],boxMax.XYZ[2]),
    TPoint<T>(boxMax.XYZ[0],boxMax.XYZ[1],boxMax.XYZ[2]),
    TPoint<T>(boxMin.XYZ[0],boxMax.XYZ[1],boxMax.XYZ[2])
  };

  return nodes;
}

template <typename T> LINT OBackground<T>::positionToBackCellIndex(const IPosition &position)
{
  assert(position.XYZ[0] >= 0 && position.XYZ[0] < IJKnumCells.XYZ[0]);
  assert(position.XYZ[1] >= 0 && position.XYZ[1] < IJKnumCells.XYZ[1]);
  assert(position.XYZ[2] >= 0 && position.XYZ[2] < IJKnumCells.XYZ[2]);

  LINT index = static_cast<LINT>(position.XYZ[0]) * 
      (static_cast<LINT>(IJKnumCells.XYZ[1]) * static_cast<LINT>(IJKnumCells.XYZ[2])) +
    static_cast<LINT>(position.XYZ[1]) * 
      static_cast<LINT>(IJKnumCells.XYZ[2]) +
    static_cast<LINT>(position.XYZ[2]);

  return index;
}

template <typename T> IPosition OBackground<T>::backCellIndexToPosition(const LINT index)
{
  assert(index < IJKnumCells.XYZ[0] * IJKnumCells.XYZ[1] * IJKnumCells.XYZ[2]);

  LINT rest = index;
  LINT x = index / (IJKnumCells.XYZ[1] * IJKnumCells.XYZ[2]);
  rest -= x * (IJKnumCells.XYZ[1] * IJKnumCells.XYZ[2]);
  LINT y = rest / IJKnumCells.XYZ[2];
  rest -= y * IJKnumCells.XYZ[2];
  LINT z = rest;

  return IPosition(x,y,z);
}

template <typename T> LINT OBackground<T>::findCell(TPoint<T> position)
{
  // get background cell number
  IPosition backCell;

  for (LINT j = 0; j < 3; j++)
  {
    backCell.XYZ[j] = static_cast<LINT>((position.XYZ[j] - boxMin[j]) / cellSize[j]);
    LIMIT(backCell.XYZ[j],0,IJKnumCells.XYZ[j] - 1);
  }

#ifdef DEBUG_FIND
  cout << "backCell " << backCell.XYZ[0] << " " << backCell.XYZ[1] << " " << backCell.XYZ[2] << endl;
#endif

  // outside
  if (backCell.XYZ[0] < 0 || backCell.XYZ[0] >= IJKnumCells.XYZ[0])
    return -1;
  if (backCell.XYZ[1] < 0 || backCell.XYZ[1] >= IJKnumCells.XYZ[1])
    return -1;
  if (backCell.XYZ[2] < 0 || backCell.XYZ[2] >= IJKnumCells.XYZ[2])
    return -1;

  // get starting background cell
  LINT cellNo = positionToBackCellIndex(backCell);

  return cellNo;
}

template <typename T> LINT OBackground<T>::numBackgroundCells()
{
  return IJKnumCells.X * IJKnumCells.Y * IJKnumCells.Z;
}