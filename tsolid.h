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
#include "toperations.h"
#include <assert.h>

//#define DEBUG_SOLID
//#ifdef NDEBUG
//  #undef DEBUG_SOLID
//#endif

namespace tcad {

/** Find all bad edges in the map which do not have exactly two neightbour faces. 
  Returns number of bad edges. */
template <class T> int findBadEdges(std::vector<std::vector<std::array<LINT,3>>> &edgepairs,
  std::vector<std::array<LINT,3>> &badedges)
{
  badedges.clear();

  for (auto &e : edgepairs)
  {
    if (e.size() != 2)
    {
      badedges.push_back(e[0]);
    }
  }

  return int(badedges.size());
}

/** Remove degenerated pieces of boundary. */
template <class T> int removeDegeneratedBoundaryPieces(std::vector<tcad::TSplineSurface<T> *> &surfaces, 
//surface     loop        part        points
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV,
  T tolerance)
{
  assert(surfaces.size() == boundariesUV.size());

  int count = 0;
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

//outputDebugString("len " + to_string(curve.length()) + " tol " + to_string(tolerance));
        if (curve.length() < tolerance) 
        {
          // remove this degenerated piece of boundary, only ONE! piece of boundary can 
          // be degenerated //!!!!!! - not sure
          boundariesUV[i][j].erase(boundariesUV[i][j].begin() + k);
          count++;
        }
      }
    } 
  }

//outputDebugString("count " + to_string(count)); 

  return count;
}

/** Find boundary piece by control sum. */
template <class T> bool findBoundaryPiece(
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV,
  unsigned int sum,  std::array<LINT,3> &location)
{
  for (int i = 0; i < int(boundariesUV.size()); i++)
  {
    for (int j = 0; j < int(boundariesUV[i].size()); j++)
    {
      for (int k = 0; k < int(boundariesUV[i][j].size()); k++)
      {
        unsigned int s = setChecksum(boundariesUV[i][j][k],false);
        if (s == sum)
        {
          location = {i,j,k};
          return true;
        }
      }
    }
  }

  return false;
}

/** Find edge by location. */
template <class T> bool findEdgeByLocation(
  std::vector<std::array<LINT,3>> &edges,
  std::vector<std::array<LINT,4>> &edgelocations,
  std::array<LINT,3> &location,
  std::array<LINT,3> &edge)
{
  for (int i = 0; i < int(edgelocations.size()); i++)
  {
    if (edgelocations[i][0] == location[0] &&
      edgelocations[i][1] == location[1] &&
      edgelocations[i][2] == location[2])
    {
      edge = edges[i];
      return true;
    }
  }

  return false;
}

/** Find location by control sum. */
template <class T> bool findLocation(
  std::vector<std::vector<std::vector<unsigned int>>> &controlsums,
  unsigned int sum, std::array<LINT,3> &location)
{
  for (int i = 0; i < int(controlsums.size()); i++)
  {
    for (int j = 0; j < int(controlsums[i].size()); j++)
    {
      for (int k = 0; k < int(controlsums[i][j].size()); k++)
      {
        if (controlsums[i][j][k] == sum)
        {
          location = {i,j,k};
          return true;
        }
      }
    }
  }

  return false;
}

/** Find piece closest to middle point. */
template <class T> bool findClosestPiece(
  std::vector<std::vector<std::vector<TPoint<T>>>> &middles,
  std::vector<std::vector<std::vector<bool>>> &busy,
  TPoint<T> middle, std::array<LINT,3> midlocation, 
  std::array<LINT,3> &location, T bigtolerance)
{
  T mindist = std::numeric_limits<T>::max();
  location = {-1,-1,-1}; 
  for (int i = 0; i < int(middles.size()); i++)
  {
    // do not touch same surfaces
    if (midlocation[0] == i)
      continue;

    for (int j = 0; j < int(middles[i].size()); j++)
    {
      for (int k = 0; k < int(middles[i][j].size()); k++)
      {
        if (busy[i][j][k])
          continue;

        T dist = !(middles[i][j][k] - middle);
        if (dist < mindist)
        {
          mindist = dist;
          location = {i,j,k};
        }
      }
    }
  }

#ifdef DEBUG_SOLID
  if (mindist >= bigtolerance)
  {
    outputDebugString(std::string("failure! ") +
      " mindist " + to_string(mindist) +  " bigtolerance " + to_string(bigtolerance));
  }
#endif

  return (mindist < bigtolerance);
}

/** Find first not busy. */
template <class T> bool findFirstNotBusy(
  std::vector<std::vector<std::vector<bool>>> &busy,
  std::array<LINT,3> &location)
{
  for (int i = 0; i < int(busy.size()); i++)
  {
    for (int j = 0; j < int(busy[i].size()); j++)
    {
      for (int k = 0; k < int(busy[i][j].size()); k++)
      {
        if (!busy[i][j][k])
        {
          location = {i,j,k};
          return true;
        }
      }
    }
  }

  return false;
}

/** Find locations (surface number + loop number + loop piece number) for one edge, 
  normally 2 for manifold, maybe one in case of failure.
  piecesXYZ - boundary pieces in XYZ coordinates.
  middleXYZ - middle points of boundary pieces in XYZ coordinates. */
template <class T> void makeEdgePairs(std::vector<tcad::TSplineSurface<T> *> &surfaces, 
//surface     loop        part        points
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV,
  std::vector<std::vector<std::array<LINT,3>>> &edgepairs,
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &piecesXYZ,
  std::vector<std::vector<std::vector<TPoint<T>>>> &middlesXYZ,
  T closestcoef = 0.1)
{
  edgepairs.clear();
  piecesXYZ.clear();
  middlesXYZ.clear();

  assert(surfaces.size() == boundariesUV.size());

  // busy
  std::vector<std::vector<std::vector<bool>>> busy;

  for (int i = 0; i < int(boundariesUV.size()); i++)
  {
    piecesXYZ.push_back(std::vector<std::vector<std::vector<tcad::TPoint<T>>>>());
    middlesXYZ.push_back(std::vector<std::vector<tcad::TPoint<T>>>());
    busy.push_back(std::vector<std::vector<bool>>());

    for (int j = 0; j < int(boundariesUV[i].size()); j++)
    {
      piecesXYZ.back().push_back(std::vector<std::vector<tcad::TPoint<T>>>());
      middlesXYZ.back().push_back(std::vector<tcad::TPoint<T>>());
      busy.back().push_back(std::vector<bool>());

      for (int k = 0; k < int(boundariesUV[i][j].size()); k++)
      {
        // this curve should be XYZ, not UV
        std::vector<TPoint<T>> points;
        getBoundaryPartXYZ<T>(surfaces,boundariesUV,i,j,k,points);
        piecesXYZ.back().back().push_back(points);

        TPointCurve<T> curve(points);
        TPoint<T> middle = curve.position(0.5);
        middlesXYZ.back().back().push_back(middle);

        busy.back().back().push_back(false);
      }
    }
  }

  // making boundary piece pairs, using "the closest" approach
  do {
    std::array<LINT,3> start;
    if (!findFirstNotBusy<T>(busy,start))
      break;

    TPoint<T> middle = middlesXYZ[start[0]][start[1]][start[2]];
    std::vector<TPoint<T>> pieceXYZ = piecesXYZ[start[0]][start[1]][start[2]];
    std::vector<TPoint<T>> pieceUV = boundariesUV[start[0]][start[1]][start[2]];

    T len = calculateLength(pieceXYZ);
    T bigtolerance = len * closestcoef;

    // mark as busy, find a second piece
    busy[start[0]][start[1]][start[2]] = true;
    edgepairs.push_back(std::vector<std::array<LINT,3>>());
    edgepairs.back().push_back(start);

    std::array<LINT,3> location;
    if (findClosestPiece(middlesXYZ,busy,middle,start,location,bigtolerance))
    {
      busy[location[0]][location[1]][location[2]] = true;
      // second location, good
      edgepairs.back().push_back(location);
    } else
    {
      // single location, not good

#ifdef DEBUG_SOLID
      std::vector<TPoint<T>> piece = piecesXYZ[location[0]][location[1]][location[2]];
      T dist0 = !(piece.front() - pieceXYZ.front());
      T dist1 = !(piece.back() - pieceXYZ.back());
      T dist2 = !(piece.front() - pieceXYZ.back());
      T dist3 = !(piece.back() - pieceXYZ.front());
      outputDebugString(std::string("failure! ") +
        " dist 0 " + to_string(dist0) +
        " dist 1 " + to_string(dist1) +
        " dist 2 " + to_string(dist2) +
        " dist 3 " + to_string(dist3)); 
#endif
    }
  } while (1);
}

/** Divide parametric piece of boundary by location and parameter U. */
template <class T> void divideBoundary(std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV,
  std::array<LINT,3> loc, std::vector<std::vector<TPoint<T>>> &newpointsUV)
{
  boundariesUV[loc[0]][loc[1]].erase(boundariesUV[loc[0]][loc[1]].begin() + loc[2]);
  boundariesUV[loc[0]][loc[1]].insert(boundariesUV[loc[0]][loc[1]].begin() + loc[2],newpointsUV.begin(),newpointsUV.end());
}

/** Find an overlapping edge for a bad edge and divide either bad edge an overlapped edge. */
template <class T> int fixBadEdge(std::array<LINT,3> badedge, 
  std::vector<tcad::TSplineSurface<T> *> &surfaces, 
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV,
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &piecesXYZ,
  T tolerance, T parmtolerance = PARM_TOLERANCE)
{
  int fixed = 0;

  std::vector<bool> busy;

  // find overlapping edges 
  std::vector<std::array<LINT,4>> list;
  int count = findOverlapping(badedge,piecesXYZ,list,tolerance * 100.0,parmtolerance); //!!!!!!!

  for (int i = 0; i < count; i++)
  {
    std::vector<tcad::TPoint<T>> badXYZ = piecesXYZ[badedge[0]][badedge[1]][badedge[2]];
    std::vector<tcad::TPoint<T>> otherXYZ = piecesXYZ[list[i][0]][list[i][1]][list[i][2]];
    std::vector<tcad::TPoint<T>> badUV = boundariesUV[badedge[0]][badedge[1]][badedge[2]];
    std::vector<tcad::TPoint<T>> otherUV = boundariesUV[list[i][0]][list[i][1]][list[i][2]];

    if (list[i][3] == 1) // bad edge is smaller, divide other edge
    {
      std::array<LINT,3> loc({list[i][0],list[i][1],list[i][2]});

      std::vector<std::vector<TPoint<T>>> newpoints,newpointsUV;

      // try both ends
      if (divide<T>(otherXYZ,otherUV,badXYZ.front(),newpoints,newpointsUV,tolerance,parmtolerance))
      {
        divideBoundary(boundariesUV,loc,newpointsUV);
        fixed++;
      } else if (divide<T>(otherXYZ,otherUV,badXYZ.back(),newpoints,newpointsUV,tolerance,parmtolerance))
      {
        divideBoundary(boundariesUV,loc,newpointsUV);
        fixed++;
      }
    } else if (list[i][3] == 2) // bad edge is larger, divide bad edge
    {
      std::array<LINT,3> loc = badedge;

      std::vector<std::vector<TPoint<T>>> newpoints,newpointsUV;

      // try both ends
      if (divide<T>(badXYZ,badUV,otherXYZ.front(),newpoints,newpointsUV,tolerance,parmtolerance))
      {
        divideBoundary(boundariesUV,loc,newpointsUV);
        fixed++;
      } else if (divide<T>(badXYZ,badUV,otherXYZ.back(),newpoints,newpointsUV,tolerance,parmtolerance))
      {
        divideBoundary(boundariesUV,loc,newpointsUV);
        fixed++;
      }
    }
  }

  return fixed;
}

/** Fix edge pairs. */
template <class T> int fixEdgePairs(std::vector<tcad::TSplineSurface<T> *> &surfaces, 
//surface     loop        part        points
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV,
  std::vector<std::vector<std::array<LINT,3>>> &edgepairs,
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &piecesXYZ,
  T tolerance, T parmtolerance, T maxedgecoef = 0.001)
{
  int fixed = 0;

  // these are locations as surface+loop+piece
  for (int i = 0; i < int(edgepairs.size()); i++)
  {
    if (edgepairs[i].size() == 1)
    {
      std::array<LINT,3> loc = edgepairs[i][0];
      int count = fixBadEdge(loc,surfaces,boundariesUV,piecesXYZ,tolerance,parmtolerance);
      fixed += count;
    }
  }

  return fixed;
}

/** Create solid model as edges; each edge has two end nodes (from vertices). 
  and two surfaces from left and right; the first has right edge direction 
  (surface to the left), second has opposite.

  Every edge contains
  [0] - first vertex
  [1] - second vertex
  [2] - middle vertex //!!! NOT USED

  [3] - first surface
  [4] - first surface loop
  [5] - first surface loop piece (4 pieces for plain untrimmed face)
  [6] - 1 if reversed (compared to [0]->[1] vertex direction)

  [7] - second surface
  [8] - second surface loop
  [9] - second surface loop piece (4 pieces for plain untrimmed face)
  [10] - 1 if reversed (compared to [0]->[1] vertex direction)
*/
template <class T> bool createSolidEdgesPrim(std::vector<tcad::TSplineSurface<T> *> &surfaces, 
//surface     loop        part        points
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV,
  std::vector<TPoint<T>> &vertices, 
  std::vector<std::array<LINT,11>> &edges,
  T tolerance, T parmtolerance, std::vector<std::vector<TPoint<T>>> *pbadedges = nullptr,
  T closestcoef = 0.1, bool makefix = true, int attempts = 50)
{
  // step 1 : make vertices and middlevertices
  vertices.clear();
  edges.clear();
  if (pbadedges)
    pbadedges->clear();

  assert(surfaces.size() == boundariesUV.size());

  // two locations for one edge, maybe one in case of failure
  std::vector<std::vector<std::array<LINT,3>>> edgepairs;
  // boundary pieces in XYZ coordinates
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> piecesXYZ;
  // middle points of boundary pieces in XYZ coordinates
  std::vector<std::vector<std::vector<TPoint<T>>>> middlesXYZ;

  for (int a = 0; a < attempts; a++)
  {
    // prepare edge pairs
    makeEdgePairs(surfaces,boundariesUV,edgepairs,piecesXYZ,middlesXYZ,closestcoef);

    if (!makefix)
      break;

    // try to fix edge pairs (changes boundariesUV)
    if (!fixEdgePairs(surfaces,boundariesUV,edgepairs,piecesXYZ,tolerance,parmtolerance))
    {
      break;
    }
  }

  // fill edges array
  for (int i = 0; i < int(edgepairs.size()); i++)
  {
    if (edgepairs[i].size() == 2)
    {
      std::array<LINT,3> loc0 = edgepairs[i][0];
      std::array<LINT,3> loc1 = edgepairs[i][1];

      std::vector<TPoint<T>> piece0 = piecesXYZ[loc0[0]][loc0[1]][loc0[2]];
      std::vector<TPoint<T>> piece1 = piecesXYZ[loc1[0]][loc1[1]][loc1[2]];

      TPoint<T> s0 = piece0.front();
      TPoint<T> e0 = piece0.back();
      TPoint<T> s1 = piece1.front();
      TPoint<T> e1 = piece1.back();

      // s0/e0 go as edge direction
      vertices.push_back(s0);
      vertices.push_back(e0);

      LINT reversed0 = 0;

      // vertex direction
      TPoint<T> v0 = s0;
      TPoint<T> v1 = e0;
      // second edge direction
      TPoint<T> p0 = s1;
      TPoint<T> p1 = e1;

      T d = !(v0 - p0) + !(v1 - p1);
      T dr = !(v0 - p1) + !(v1 - p0);

      // reversed, must be true
      LINT reversed1 = (LINT) (dr < d);

      edges.push_back(std::array<LINT,11>{(LINT) vertices.size() - 2,(LINT) vertices.size() - 1,-1,
        loc0[0],loc0[1],loc0[2],reversed0,loc1[0],loc1[1],loc1[2],reversed1});
    }
  } 

  // step 2 : exclude duplicate vertices
  std::vector<LINT> vreplacement;
  removeDuplicates(vertices,true,tolerance,&vreplacement);

  // step 3 : renumber indices in edges
  for (auto &e : edges)
  {
    e[0] = vreplacement[e[0]];
    e[1] = vreplacement[e[1]];
  }

  // list of bad edges which do not have two face neighbours
  std::vector<std::array<LINT,3>> badedges;

  int n = findBadEdges<T>(edgepairs,badedges);

  if (n && pbadedges != nullptr)
  {
    for (int k = 0; k < int(badedges.size()); k++)
    {
      std::vector<TPoint<T>> points = piecesXYZ[badedges[k][0]][badedges[k][1]][badedges[k][2]];
      pbadedges->push_back(points);
    }
  }

#ifdef DEBUG_SOLID
  outputDebugString(std::string("num bad edges ") + to_string(n));
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
  [2] - middle vertex //!!! NOT USED

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
  std::vector<TPoint<T>> &vertices,  
  std::vector<std::array<LINT,11>> &edges,
  T tolerance, T parmtolerance, std::vector<std::vector<TPoint<T>>> *pbadedges = nullptr)
{
  edges.clear();

  bool ok = createSolidEdgesPrim<T>(surfaces,boundariesUV,vertices,edges,tolerance,parmtolerance,pbadedges);

  return ok;
}

/** Get a list of edges for surface and its loop. */
template <class T> int findLoopEdges(std::vector<std::array<LINT,11>> &edges, int surface, 
  int loop, int loopsize, std::vector<std::pair<int,bool>> &iedges)
{
  for (int i = 0; i < loopsize; i++)
  {
    for (int j = 0; j < int(edges.size()); j++)
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

/** Check loop. */
template <class T> bool loopOK(std::vector<tcad::TSplineSurface<T> *> &surfaces, 
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV,
  std::vector<TPoint<T>> &vertices, 
  std::vector<std::array<LINT,11>> &edges, int nsurface, 
  int nloop, int loopsize, std::vector<std::pair<int,bool>> &iedges, T tolerance,
  std::vector<std::vector<TPoint<T>>> &dloop0, std::vector<std::vector<TPoint<T>>> &dloop1)
{
  // remember edges involved
  std::vector<int> nedges;

  // make iedges for this nsurface and nloop
  for (int i = 0; i < loopsize; i++)
  {
    for (int j = 0; j < int(edges.size()); j++)
    {
      if (edges[j][3] == nsurface && edges[j][4] == nloop && edges[j][5] == i)
      {
        iedges.push_back(std::pair<int,bool>(j,bool(edges[j][6])));
        nedges.push_back(j);
      }
      if (edges[j][7] == nsurface && edges[j][8] == nloop && edges[j][9] == i)
      {
        iedges.push_back(std::pair<int,bool>(j,bool(edges[j][10])));
        nedges.push_back(j);
      }
    }
  }

  // edge maybe degenerated
  //assert(loopsize == iedges.size());
  //assert(loopsize == nedges.size());
  //if (loopsize != iedges.size())
  //  return false;

  bool ok = true;

  // take surface and loop, ther are constant here
  tcad::TSplineSurface<T> *surface = surfaces[nsurface];
  std::vector<std::vector<tcad::TPoint<T>>> &loop = boundariesUV[nsurface][nloop];

#ifdef DEBUG_SOLID
  outputDebugString(string(""));
#endif

  dloop0.clear();
  dloop1.clear();

  // these are edges i.e. pieces of boundary
  for (int i = 0; i < int(iedges.size()); i++)
  {
    int nedge = iedges[i].first;
    bool reversed = iedges[i].second;

    // ending/middle vertices from edge
    TPoint<T> v0 = vertices[edges[nedge][0]];
    TPoint<T> v1 = vertices[edges[nedge][1]];

    // piece of boundary in UV
    std::vector<tcad::TPoint<T>> &pieceUV = loop[i];

    // into XYZ
    std::vector<tcad::TPoint<T>> points;
    surface->UVIntoPoints(pieceUV,points);

    TPoint<T> p0 = points.front();
    TPoint<T> p1 = points.back();

    T d = !(v0 - p0) + !(v1 - p1);
    T dr = !(v0 - p1) + !(v1 - p0);

    bool reversed01 = (dr < d);

    assert(reversed == reversed01);

    if (reversed)
    {
      std::vector<TPoint<T>> line;
      line.push_back(v1);
      line.push_back(v0);
      dloop0.push_back(line);
    } else
    {
      std::vector<TPoint<T>> line;
      line.push_back(v0);
      line.push_back(v1);
      dloop0.push_back(line);
    }

    dloop1.push_back(points);

#ifdef DEBUG_SOLID
    outputDebugString(std::string("rev ") + to_string(int(reversed)) + std::string(" rev01 ") + to_string(int(reversed01)));
#endif
  }

  return ok;
}

}