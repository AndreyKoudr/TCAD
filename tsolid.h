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

  tsolid.h

  Construct edges and vertices for a solid model from spline surfaces and 
  boundary loops

*******************************************************************************/

#pragma once

#include "tpoints.h"
#include "tpointcurve.h"
#include "tsplinesurface.h"
#include "tedge.h"
#include <assert.h>

namespace tcad {

/** Find all bad edges in the map which do not have exactly two neightbour faces. 
  Returns number of bad edges. */
template <class T> int findBadEdges(std::map<std::array<LINT,3>,std::vector<std::array<LINT,4>>,TEdgeCompare<T>> &edgemap,
  std::vector<std::array<LINT,3>> &badedges)
{
  badedges.clear();

  for (auto &e : edgemap)
  {
    if (e.second.size() != 2)
    {
      badedges.push_back(e.first);
    }
  }

  return int(badedges.size());
}

/** Get piece of boundary curve in XYZ from parametric boundariesUV. */
template <class T> void getBoundaryPartXYZ(std::vector<tcad::TSplineSurface<T> *> &surfaces, 
//surface     loop        part        points
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV,
  int surface, int loop, int part, 
  std::vector<TPoint<T>> &points)
{
  points.clear();

  // this curve should be XYZ, not UV
  for (auto &UV : boundariesUV[surface][loop][part])
  {
    points.push_back(surfaces[surface]->position(UV.X,UV.Y));
  }
}

/** Create solid model as edges; each edge has two end nodes (from vertices) 
  pluse one middle node (from middlevertices)
  and two surfaces from left and right; the first has right edge direction 
  (surface to the left), second has opposite.
*/
template <class T> bool createSolidEdges(std::vector<tcad::TSplineSurface<T> *> &surfaces, 
//surface     loop        part        points
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV,
  std::vector<TPoint<T>> &vertices, std::vector<TPoint<T>> &middlevertices, 
  std::map<std::array<LINT,3>,std::vector<std::array<LINT,4>>,TEdgeCompare<T>> &edgemap,
  T tolerance, std::vector<std::vector<TPoint<T>>> *pbadedges = nullptr)
{
  // step 1 : make vertices and middlevertices
  vertices.clear();
  middlevertices.clear();
  edgemap.clear();
  if (pbadedges)
    pbadedges->clear();

  assert(surfaces.size() == boundariesUV.size());

  // original list of edges with duplicate nodes
  std::vector<std::array<LINT,3>> edges;
  // these are indices to identify where the edges are from (surface-loop-boundarypart 
  // location); the last index defines if the edge is reversed (1) or not
  std::vector<std::array<LINT,4>> edgelocations;

  for (int i = 0; i < int(surfaces.size()); i++)
  {
    for (int j = 0; j < int(boundariesUV[i].size()); j++)
    {
      for (int k = 0; k < int(boundariesUV[i][j].size()); k++)
      {
        // this curve should be XYZ, not UV
        std::vector<TPoint<T>> points;
        getBoundaryPartXYZ<T>(surfaces,boundariesUV,i,j,k,points);

        TPointCurve<T> curve(points);

        if (curve.length() > tolerance)
        {
          TPoint<T> s = curve.start();
          TPoint<T> e = curve.end();
          TPoint<T> m = curve.middle();

          vertices.push_back(s);
          vertices.push_back(e);
          middlevertices.push_back(m);

          edges.push_back(std::array<LINT,3>{(LINT) vertices.size() - 2,(LINT) vertices.size() - 1,(LINT) middlevertices.size() - 1});
          edgelocations.push_back(std::array<LINT,4>{i,j,k,0});
        }
      } // boundary part (curve)
    } // loops
  } // surfaces

//!!!!!!! #ifdef _DEBUG
#if 1
  outputDebugString(std::string("BEFORE tolerance ") + to_string(tolerance) + 
    std::string(" vertices ") + to_string(int(vertices.size())) + 
    std::string(" middlevertices ") + to_string(int(middlevertices.size())));
#endif

  // step 2 : exclude duplicate vertices
  std::vector<LINT> vreplacement,mreplacement;
  removeDuplicates(vertices,true,tolerance,&vreplacement);
  removeDuplicates(middlevertices,true,tolerance,&mreplacement);

//!!!!!!! #ifdef _DEBUG
#if 1
  outputDebugString(std::string("AFTER tolerance ") + to_string(tolerance) + 
    std::string(" vertices ") + to_string(int(vertices.size())) + 
    std::string(" middlevertices ") + to_string(int(middlevertices.size())));
#endif

  // step 3 : renumber indices in edges
  for (auto &e : edges)
  {
    e[0] = vreplacement[e[0]];
    e[1] = vreplacement[e[1]];
    e[2] = mreplacement[e[2]];
  }

  // step 4 : make a map : list of boundary pieces for every edge
  edgemap.clear();

  for (int i = 0; i < int(edges.size()); i++)
  {
    // new edge
    if (edgemap.find(edges[i]) == edgemap.end())
    {
      edgemap.insert(std::pair<std::array<LINT,3>,std::vector<std::array<LINT,4>>>(edges[i],{edgelocations[i]})); 
    } else
    // same edge, add new boundary part
    {
      edgemap[edges[i]].push_back(edgelocations[i]);
    }
  }

  // set "reversed" for every edge location
  for (auto &e : edgemap)
  {
    // reversed?
    TPoint<T> v0 = vertices[e.first[0]];
    TPoint<T> v1 = vertices[e.first[1]];

    for (int i = 0; i < int(e.second.size()); i++)
    {
      std::vector<TPoint<T>> points;
      getBoundaryPartXYZ<T>(surfaces,boundariesUV,
        int(e.second[i][0]),int(e.second[i][1]),int(e.second[i][2]),points);

      TPoint<T> p0 = points.front();
      TPoint<T> p1 = points.back();

      T d = !(v0 - p0) + !(v1 - p1);
      T dr = !(v0 - p1) + !(v1 - p0);

      // reversed
      e.second[i][3] = (LINT) (dr < d);
    }
  }

  // list of bad edges which do not have two face neighbours
  std::vector<std::array<LINT,3>> badedges;

  int n = findBadEdges(edgemap,badedges);

  if (n && pbadedges != nullptr)
  {
    for (int k = 0; k < int(badedges.size()); k++)
    {
      std::vector<TPoint<T>> points;
      points.push_back(vertices[badedges[k][0]]);
      points.push_back(middlevertices[badedges[k][2]]);
      points.push_back(vertices[badedges[k][1]]);
      
      pbadedges->push_back(points);
    }
  }

//!!!!!!! #ifdef _DEBUG
#if 1
  outputDebugString(std::string("num  bad edges ") + to_string(n));
#endif

  return (n == 0);
}

/** Create solid model as edges; each edge has two end nodes (from vertices) 
  pluse one middle node (from middlevertices)
  and two surfaces from left and right; the first has right edge direction 
  (surface to the left), second has opposite.
  Every edge contains
  [0] - first vertex
  [1] - second vertex
  [2] - middle vertex

  [3] - first surface
  [4] - first surface loop
  [5] - first surface loop piece (4 pieces for plain untrimmed face)
  [6] - 1 if reversed (compared to [0]->[1] vertex direction)

  [7] - second surface
  [8] - second surface loop
  [9] - second surface loop piece (4 pieces for plain untrimmed face)
  [10] - 1 if reversed (compared to [0]->[1] vertex direction)
*/
template <class T> bool createSolidEdges(std::vector<tcad::TSplineSurface<T> *> &surfaces, 
//surface     loop        part        points
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV,
  std::vector<TPoint<T>> &vertices, std::vector<TPoint<T>> &middlevertices, 
  std::vector<std::array<LINT,11>> &edges,
  T tolerance, std::vector<std::vector<TPoint<T>>> *pbadedges = nullptr, int attempts = 10)
{
  edges.clear();

  std::map<std::array<LINT,3>,std::vector<std::array<LINT,4>>,TEdgeCompare<T>> edgemap;

  bool ok = false;
  T atolerance = tolerance;

  for (int i = 0; i < attempts; i++)
  {
    ok = createSolidEdges<T>(surfaces,boundariesUV,vertices,middlevertices,edgemap,atolerance,pbadedges);
    if (ok)
      break;
    atolerance *= 2.0;
  }

  if (!ok)
  {
    return false;
  }

  // only two faces per edge are allowed here
  for (auto &e : edgemap)
  {
    std::array<LINT,11> a = {e.first[0],e.first[1],e.first[2],
      e.second[0][0],e.second[0][1],e.second[0][2],e.second[0][3],
      e.second[1][0],e.second[1][1],e.second[1][2],e.second[1][3]};
    edges.push_back(a);
  }

  return true;
}

/** Get a list of edges for surface and its loop. */
template <class T> int findLoopEdges(std::vector<std::array<LINT,11>> &edges, int surface, 
  int loop, int loopsize, std::vector<std::pair<int,bool>> &iedges)
{
  for (int i = 0; i < loopsize; i++)
  {
    for (int j = 0; j < edges.size(); j++)
    {
      if (edges[j][3] == surface && edges[j][4] == loop && edges[j][5] == i)
      {
        iedges.push_back(std::pair<int,bool>(j,bool(edges[j][6])));
      }
      if (edges[j][7] == surface && edges[j][8] == loop && edges[j][9] == i)
      {
        iedges.push_back(std::pair<int,bool>(j,bool(edges[j][10])));
      }
    }
  }

  return int(iedges.size());
}

}