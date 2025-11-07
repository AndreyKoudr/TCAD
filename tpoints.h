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

  tpoints.h

  Operations on 3D points as std::vector<TPoint<T>> like min/max etc.

*******************************************************************************/

#pragma once

#include "tbasics.h"
#include "tpoint.h"
#include "tmisc.h"
#include "tplane.h"
#include "ttransform.h"

#include <vector>
#include <map>
#include <algorithm>

namespace tcad {

/** Operators on lists of points. */

template <class T> std::vector<TPoint<T>> operator + (const std::vector<TPoint<T>> &points0, 
  const std::vector<TPoint<T>> &points1)
{
  assert(points0.size() == points1.size());

  std::vector<TPoint<T>> points;
  for (int i = 0; i < int(points0.size()); i++)
  {
    points.push_back(points0[i] + points1[i]);
  }

  return points;
}

template <class T> std::vector<TPoint<T>> operator - (const std::vector<TPoint<T>> &points0, 
  const std::vector<TPoint<T>> &points1)
{
  assert(points0.size() == points1.size());

  std::vector<TPoint<T>> points;
  for (int i = 0; i < int(points0.size()); i++)
  {
    points.push_back(points0[i] - points1[i]);
  }

  return points;
}

template <class T> std::vector<TPoint<T>> operator * (const std::vector<TPoint<T>> &points0, 
  const std::vector<T> &coefs)
{
  assert(points0.size() == coefs.size());

  std::vector<TPoint<T>> points;
  for (int i = 0; i < int(points0.size()); i++)
  {
    points.push_back(points0[i] * coefs[i]);
  }

  return points;
}

/** Normalisation. */
template <class T> std::vector<TPoint<T>> operator + (std::vector<TPoint<T>> &points0)
{
  std::vector<TPoint<T>> points;
  for (int i = 0; i < int(points0.size()); i++)
  {
    points.push_back((+points0[i]));
  }

  return points;
}

/** Get start direction (first derivative on U) for point list. */
template <class T> TPoint<T> startDirection(std::vector<TPoint<T>> &points)
{
  TPoint<T> d = points[1] - points[0];
  T len = !d;
  T curvelen = calculateLength(points);
  T DU = len / curvelen;
  TPoint<T> dir = d / DU;
  return dir;
}

/** Get end direction (first derivative on U) for point list, directed from end "inside". */
template <class T> TPoint<T> endDirection(std::vector<TPoint<T>> &points)
{
  TPoint<T> d = points[points.size() - 2] - points[points.size() - 1];
  T len = !d;
  T curvelen = calculateLength(points);
  T DU = len / curvelen;
  TPoint<T> dir = d / DU;
  return dir;
}

/** Calculate min/max among a list of points; imin, imax contain corresponding indices as reals. */
template <class T> bool calculateMinMax(std::vector<TPoint<T>> &points, TPoint<T> *min, TPoint<T> *max,
  TPoint<T> *imin = nullptr, TPoint<T> *imax = nullptr)
{
  if (points.empty())
    return false;

  *min = *max = points[0];
  TPoint<T> point;

  if (imin)
  {
    *imin = TPoint<T>(0,0,0);
  }

  if (imax)
  {
    *imax = TPoint<T>(0,0,0);
  }

  for (int i = 1; i < int(points.size()); i++)
  {
    point = points[i];
    if (point.X > max->X) 
    {
      max->X = point.X;
      if (imax)
        imax->X = i;
    }
    if (point.Y > max->Y) 
    {
      max->Y = point.Y;
      if (imax)
        imax->Y = i;
    }
    if (point.Z > max->Z) 
    {
      max->Z = point.Z;
      if (imax)
        imax->Z = i;
    }

    if (point.X < min->X) 
    {
      min->X = point.X;
      if (imin)
        imin->X = i;
    }
    if (point.Y < min->Y) 
    {
      min->Y = point.Y;
      if (imin)
        imin->Y = i;
    }
    if (point.Z < min->Z) 
    {
      min->Z = point.Z;
      if (imin)
        imin->Z = i;
    }
  }

  return true;
}

/** Calculate length. */
template <class T> T calculateLength(std::vector<TPoint<T>> &points)
{
  T len = 0.0;

  if (points.size() > 0)
  {
    for (int i = 1; i < points.size(); i++)
    {
      len += !(points[i] - points[i - 1]);
    }
  }

  return len;
}

/** Calculate max difference. */
template <class T> T difference(std::vector<TPoint<T>> &points, std::vector<TPoint<T>> &other)
{
  if (points.size() != other.size())
    return -1.0;

  T maxdiff = 0.0;
  if (points.size() > 0)
  {
    for (int i = 0; i < points.size(); i++)
    {
      T diff = !(other[i] - points[i]);
      maxdiff = std::max<T>(diff,maxdiff);
    }
  }

  return maxdiff;
}

/** Get parameters [0..1]. Note : it is better to use parameterisation by numbers
  everywhere, especially in TPointCurve and TOrthoSegment! */
template <class T> bool prepareParameters(std::vector<TPoint<T>> &points, 
  std::vector<T> &parms, bool normalise = true, bool bynumbers = false, T tolerance = TOLERANCE(T))
{
  parms.clear();

  if (points.size() < 2)
    return false;
  if (points.size() == 2 && !(points[1] - points[0]) < tolerance)
  {
    return false;
  }

  T len = calculateLength(points);
  if (len < tolerance)
  {
    for (int i = 0; i < points.size(); i++)
    {
      T U = T(i) / T(points.size() - 1);
      parms.push_back(U);
    }
  } else
  {
    if (bynumbers)
    {
      for (int i = 0; i < points.size(); i++)
      {
        T U = T(i) / T(points.size() - 1);
        parms.push_back(U);
      }

      // these parms are already scaled
    } else
    {
      T len = 0.0;
      parms.push_back(len);
      for (int i = 1; i < points.size(); i++)
      {
        len += !(points[i] - points[i - 1]);
        parms.push_back(len);
      }
                          
      // scale parms from 0 to 1
      if (normalise)
      {    
        for (auto &p : parms)
        {
          p /= parms.back();
        }
      }
    }
  }

  return true;
}

/** Split points into x,y,z arrays. */
template <class T> void splitXYZ(std::vector<TPoint<T>> &points, std::vector<T> &x, std::vector<T> &y, std::vector<T> &z)
{
  x.clear();
  y.clear();
  z.clear();

  for (auto &p : points)
  {
    x.push_back(p.X);  
    y.push_back(p.Y);
    z.push_back(p.Z);
  }
}

/** Combine points x,y,z into points. */
template <class T> void combineXYZ(std::vector<T> &x, std::vector<T> &y, std::vector<T> &z, std::vector<TPoint<T>> &points)
{
  points.clear();

  for (int i = 0; i < int(x.size()); i++)
  {
    points.push_back(TPoint<T>(x[i],y[i],z[i]));
  }
}

/** Find closest vector. */
template <class T> int findClosest(std::vector<TPoint<T>> &points, TPoint<T> value, T *mindist = nullptr)
{
  T min = std::numeric_limits<T>::max();
  int index = -1;

  for (int i = 0; i < points.size(); i++)
  {
    T dist = !(points[i] - value);
    if (dist < min)
    {
      min = dist;
      index = i;
    }
  }

  if (mindist)
    *mindist = min;

  return index;
}

/** Find closest vector among not busy. */
template <class T> int findClosestNotBusy(std::vector<TPoint<T>> &points, TPoint<T> value,
  std::vector<bool> &busy, T *mindist = nullptr)
{
  T min = std::numeric_limits<T>::max();
  int index = -1;

  for (int i = 0; i < points.size(); i++)
  {
    if (busy[i])
      continue;

    T dist = !(points[i] - value);
    if (dist < min)
    {
      min = dist;
      index = i;
    }
  }

  if (index < 0)
    min = 0.0;

  if (mindist)
    *mindist = min;

  return index;
}

/** Find closest vector among not busy. */
template <class T> int findClosestNotBusyExcept(std::vector<TPoint<T>> &points, TPoint<T> value,
  std::vector<bool> &busy, int except, T *mindist = nullptr)
{
  T min = std::numeric_limits<T>::max();
  int index = -1;

  for (int i = 0; i < points.size(); i++)
  {
    if (busy[i] || i == except)
      continue;

    T dist = !(points[i] - value);
    if (dist < min)
    {
      min = dist;
      index = i;
    }
  }

  if (index < 0)
    min = 0.0;

  if (mindist)
    *mindist = min;

  return index;
}

/** Calculate number of intersections of segment p0->p1 with edges in XY plane. */
template <class T> int numEdgeIntersectionsXY(std::vector<std::pair<TPoint<T>,TPoint<T>>> &edges, 
  TPoint<T> p0, TPoint<T> p1, T tolerance)
{
  int count = 0;

  for (int i = 0; i < edges.size(); i++)
  {
    TPoint<T> v0 = edges[i].first;
    TPoint<T> v1 = edges[i].second;
    
    T t1 = 0.0;
    T t2 = 0.0;
    T Xi = p1.X;
    T Yi = p1.Y;

    if (intersectSegmentsXY(p0.X,p0.Y,p1.X,p1.Y,v0.X,v0.Y,v1.X,v1.Y,&t1,&t2,&Xi,&Yi) &&
      (t1 >= 0.0) && (t1 <= 1.0) && (t2 >= 0.0) && (t2 <= 1.0))
    {
      TPoint<T> intr(Xi,Yi);
      bool startend = (
        (!(p0 - intr) < tolerance) || 
        (!(p1 - intr) < tolerance) || 
        (!(v0 - intr) < tolerance) || 
        (!(v1 - intr) < tolerance));

      if (!startend) 
        count++;
    }
  }

  return count;
}

/** Calculate number of intersections of segment p0->p1 with edges. */
template <class T> int numEdgeIntersections(std::vector<std::pair<TPoint<T>,TPoint<T>>> &edges, 
  TPoint<T> p0, TPoint<T> p1, T tolerance)
{
  int count = 0;

  for (int i = 0; i < edges.size(); i++)
  {
    TPoint<T> v0 = edges[i].first;
    TPoint<T> v1 = edges[i].second;
    
    T t1 = 0.0;
    T t2 = 0.0;
    T dist = 0.0;
    TPoint<T> ip,iv;
    if (intersectSegments(p0,p1,v0,v1,t1,t2,dist,&ip,&iv) &&
      (t1 >= 0.0) && (t1 <= 1.0) && (t2 >= 0.0) && (t2 <= 1.0) && dist < tolerance) //!!!
    {
      bool startend = (
        (!(p0 - ip) < tolerance) || 
        (!(p1 - ip) < tolerance) || 
        (!(v0 - iv) < tolerance) || 
        (!(v1 - iv) < tolerance));

      if (!startend) 
        count++;
    }
  }

  return count;
}

/** Find intersections between two curves represented as points. Every point MUST contain 
  their U parameter in W as per createPoints() in tbasecurve. UV contains intersection U 
  parameters for both curves in X and Y. */
template <class T> int findIntersections(std::vector<TPoint<T>> &points0, std::vector<TPoint<T>> &points1, 
  std::vector<TPoint<T>> &UV, T tolerance, T parmtolerance = PARM_TOLERANCE)
{
  UV.clear();

  for (int i = 0; i < points0.size() - 1; i++)
  {
    TPoint<T> p0 = points0[i];
    TPoint<T> p1 = points0[i + 1];

    for (int j = 0; j < points1.size() - 1; j++)
    {
      TPoint<T> v0 = points1[j];
      TPoint<T> v1 = points1[j + 1];
    
      T t1 = 0.0;
      T t2 = 0.0;
      T dist = 0.0;
      TPoint<T> ip,iv;
      if (intersectSegments(p0,p1,v0,v1,t1,t2,dist,&ip,&iv) &&
        (t1 >= 0.0 - parmtolerance) && (t1 <= 1.0 + parmtolerance) && 
        (t2 >= 0.0 - parmtolerance) && (t2 <= 1.0 + parmtolerance) && dist < tolerance) //!!!
      {
        LIMIT(t1,0.0,1.0);
        LIMIT(t2,0.0,1.0);

        T U1 = p0.W + (p1.W - p0.W) * t1;
        T U2 = v0.W + (v1.W - v0.W) * t2;
        UV.push_back(TPoint<T>(U1,U2));
      }
    }
  }

  // remove all duplicates
  removeDuplicates(UV,true,parmtolerance * 4.0);

  return int(UV.size());
}

/** Find edge closest to edges. Set maxedge to -1.0 or/and alledges to empty to avoid corresponding checks. */
template <class T> bool findClosestNotBusy(std::vector<std::pair<TPoint<T>,TPoint<T>>> &edges, 
  std::vector<std::pair<TPoint<T>,TPoint<T>>> &alledges,
  TPoint<T> value, std::vector<bool> &busy, int &index, bool &reversed, T maxedge, T tolerance, 
  T *mindist = nullptr)
{
  T min = std::numeric_limits<T>::max();
  index = -1;
  bool ok = false;

  for (int i = 0; i < edges.size(); i++)
  {
    if (busy[i])
      continue;

    T dist0 = !(edges[i].first - value);
    T dist1 = !(edges[i].second - value);
    if (dist0 < min)
    {
      min = dist0;
      reversed = false;
      index = i;
    }
    if (dist1 < min)
    {
      min = dist1;
      reversed = true;
      index = i;
    }
  }

  if (index >= 0)
  {
    ok = true;

    if (min > maxedge)
    {
      ok = false;
    } else
    {
      if (!alledges.empty())
      {
        int count = numEdgeIntersections(alledges,value,reversed ? edges[index].second : 
          edges[index].first,tolerance);

        if (count)
        {
          ok = false;
        }
      }
    }
  }

  if (!ok)
    min = 0.0;

  if (mindist)
    *mindist = min;

  return ok;
}

/** Make up curve(s) from unordered pieces, like hanging edges on triangles boundary
  or intersection curve pieces. */
template <class T> bool makeUpCurves(std::vector<std::pair<TPoint<T>,TPoint<T>>> &edges,
  std::vector<std::pair<TPoint<T>,TPoint<T>>> &alledges, T maxedge, 
  std::vector<std::vector<TPoint<T>>> &lines, T tolerance, bool bothways = true)
{
  if (edges.empty())
    return false;

  std::vector<bool> busy(edges.size(),false);

  // we need to start at an edge close to the middle as in case of
  // degenerated edges (separate points) we cannot calculate maxedge,
  // it is set as a large fraction of the whole size
  TPoint<T> middle;
  for (int i = 0; i < int(edges.size()); i++)
  {
    TPoint<T> m = (edges[i].first + edges[i].second) * 0.5;
    middle += m;
  }
  middle /= T(edges.size());

  bool reversed = false;

  // this may return negative is dist > maxedge, we do not need it here
  int startedge = -1;
  findClosestNotBusy(edges,alledges,middle,busy,startedge,reversed,maxedge,tolerance);

  std::vector<TPoint<T>> line;

  line.push_back(edges[startedge].first);
  line.push_back(edges[startedge].second);
  busy[startedge] = true;

  while (!allBusy(busy))
  {
    // this fails only if all edges are busy; if distance > maxedge, a negative index is returned
    bool reversed = false;

    int index = -1;
    bool ok0 = findClosestNotBusy(edges,alledges,line.back(),busy,index,reversed,maxedge,tolerance);

    if (!ok0)
    {
      if (bothways)
      {
        bool reversed1 = false;
        int index1 = -1;
        bool ok1 = findClosestNotBusy(edges,alledges,line.front(),busy,index1,reversed1,maxedge,tolerance);
        if (!ok1)
        {
          // restart
          lines.push_back(line);
          line.clear();

          startedge = findFirstNotBusy(busy);
          if (startedge < 0)
            break;

          line.push_back(edges[startedge].first);
          line.push_back(edges[startedge].second);
          busy[startedge] = true;
        } else
        {
          if (reversed1)
          {
            line.insert(line.begin(),edges[index1].first);
            busy[index1] = true;
          } else
          {
            line.insert(line.begin(),edges[index1].second);
            busy[index1] = true;
          }
        }
      } else
      {
        if (index >= 0)
          busy[index] = true;
      }
    } else
    {
      if (reversed)
      {
        line.push_back(edges[index].first);
        busy[index] = true;
      } else
      {
        line.push_back(edges[index].second);
        busy[index] = true;
      }
    }
  }

  if (!line.empty())
    lines.push_back(line);

  return true;
}

/** Make up curves from pieces. In case of NOT degenerateedges,
  pieces contain many pieces ot two points generated by intersections.
  If degenerateedges, pieces.size() is 1 which contains many isolated
  points each of which will be duplicated into degenarted edges.
  maxedge ratio is a guessed max edge size to the model size,
  it is unknown in case edges are points (degenerated edges). */
template <class T> bool curvesFromPieces(std::vector<std::vector<TPoint<T>>> &pieces,
  std::vector<std::vector<TPoint<T>>> &lines, T tolerance, 
  bool degenerateedges = false, T maxedgeratio = 0.1)
{
  // make edges
  std::vector<std::pair<TPoint<T>,TPoint<T>>> edges;

  T maxedge = 0.0;
  for (int i = 0; i < pieces.size(); i++)
  {
    if (degenerateedges)
    {
      // make degenerated edges
      for (int j = 0; j < int(pieces[i].size()); j++)
      {
        TPoint<T> p0 = pieces[i][j];
        edges.push_back(std::pair<TPoint<T>,TPoint<T>>(p0,p0));
      }

      // we cannot get a reliable egde length here, so set it as a large 
      // part of the whole length and start making line from a node
      // close to the middle
      TPoint<T> min,max;
      if (calculateMinMax(pieces[0],&min,&max))
      {
        maxedge = !(max - min) * maxedgeratio;
      }
    } else
    {
      for (int j = 0; j < int(pieces[i].size()) - 1; j++)
      {
        TPoint<T> p0 = pieces[i][j];
        TPoint<T> p1 = pieces[i][j + 1];
        T dist = !(p1 - p0);
        maxedge = std::max<T>(maxedge,dist);

        edges.push_back(std::pair<TPoint<T>,TPoint<T>>(p0,p1));
      }
    }
  }

  std::vector<std::pair<TPoint<T>,TPoint<T>>> alledges;
  bool ok = false;

  if (degenerateedges)
  {
    ok = makeUpCurves(edges,alledges,maxedge + tolerance,lines,tolerance,true); 
  } else
  {
    ok = makeUpCurves(edges,alledges,tolerance,lines,tolerance,true); 
  }

  return ok;
}

/** Make up single curve from pieces. In case of NOT degenerateedges,
  pieces contain many pieces ot two points generated by intersections.
  If degenerateedges, pieces.size() is 1 which contains many isolated
  points each of which will be duplicated into degenerated edges. 
  maxedge ratio is a guessed max edge size to the model size,
  it is unknown in case edges are points (degenerated edges). */
template <class T> bool curveFromPieces(std::vector<std::vector<TPoint<T>>> &pieces,
  std::vector<TPoint<T>> &line, T tolerance, bool degenerateedges = false, T maxedgeratio = 0.1)
{
  std::vector<std::vector<TPoint<T>>> lines;
  if (curvesFromPieces(pieces,lines,tolerance,degenerateedges,maxedgeratio) && lines.size() == 1)
  {
    line = lines[0];
    return true;
  } else
  {
    return false;
  }
}

/** Closed? */
template <class T> bool closed(const std::vector<TPoint<T>> &points, T tolerance)
{
  if (points.empty())
    return false;

  T d = !(points.front() - points.back());

  return (d < tolerance);
}

/** Get three neighbours around point index. */
template <class T> bool findHeighbours(const std::vector<TPoint<T>> &points, int index, 
  int &prev, int &next, T tolerance)
{
  if (points.size() < 3)
    return false;

  bool cl = closed(points,tolerance);

  if (cl)
  {
    prev = index - 1;
    if (prev < 0)
      prev += int(points.size());
    next = index + 1;
    if (next >= points.size())
      next -= int(points.size());
  } else
  {
    prev = index - 1;
    if (prev < 0)
      return false;
    next = index + 1;
    if (next >= int(points.size()))
      return false;
  }

  return true;
}

/** Get three neighbours around point index. */
template <class T> bool findThreePoints(const std::vector<TPoint<T>> &points, int index, 
  TPoint<T> &p0, TPoint<T> &p1, TPoint<T> &p2, T tolerance)
{
  int prev = -1;
  int next = -1;

  if (findHeighbours(points,index,prev,next,tolerance))
  {
    p0 = points[prev];
    p1 = points[index];
    p2 = points[next];

    return true;
  } else
  {
    return false;
  }
}

/** Get angle in degrees around point index. Returns -1.0 in case of failure. */
template <class T> T calcThreeAngle(const std::vector<TPoint<T>> &points, int index, T tolerance)
{
  if (points.size() < 3)
    return false;

  TPoint<T> p0,p1,p2;

  if (findThreePoints(points,index,p0,p1,p2,tolerance))
  {
    T angle = ((p1 - p0) < (p2 - p1)) * PCI;
    return angle;
  } else
  {
    return -1.0;
  }
}

/** Find sharp corners. First, try to find ratio if angle at sharp corner to angles 
  of its closest neighbour, typically 2.0-40.0-3.0 (max ~6.0) or so. 
  If ratio is > 5.0, it is a sharp corner. Otherwise use just sharp angle. */
template <class T> bool findSharpCorners(const std::vector<TPoint<T>> &points,
  std::vector<int> &indices, T tolerance, T sharpangledeg = 45.0, T sharpangleratio = 0.0)
{ 
  indices.clear();

  if (points.empty())
    return false;

  for (int i = 1; i < points.size() - 1; i++)
  {
    bool sharp = false;

    T a0 = calcThreeAngle(points,i - 1,tolerance);
    T a1 = calcThreeAngle(points,i,tolerance);
    T a2 = calcThreeAngle(points,i + 1,tolerance);

    T ratio = 0.0;

    if ((sharpangleratio > 0.0) && (a0 >= 0.0) && (a1 >= 0.0) && (a2 >= 0.0))
    {
      T a02 = std::max<T>(a0,a2);
      if (a02 < tolerance)
      {
        if (a1 < tolerance)
        {
          // both are zeroes
        } else
        {
          // sharp edge between flat parts
          if (a1 > sharpangledeg)
          {
            indices.push_back(i);
            sharp = true;
          }
        }
      } else 
      {
        if (a1 < tolerance)
        {
          // flat
        } else
        {
          ratio = a1 / a02;
          if (ratio > sharpangleratio)
          {
            indices.push_back(i);
            sharp = true;
          }
        }
      }
    } else
    {
      if (a1 > sharpangledeg)
      {
        indices.push_back(i);
        sharp = true;
      }
    }
  }

  return !indices.empty();
}

/** Find sharp corners, redivide all curves. */
template <class T> void redividePoints(std::vector<std::vector<TPoint<T>>> &points,
  std::vector<std::vector<TPoint<T>>> &newpoints,
  T tolerance, T sharpangledeg = 45.0, T sharpangleratio = 0.0)
{ 
  newpoints.clear();

  for (int i = 0; i < int(points.size()); i++)
  {
    std::vector<int> indices; 
    if (findSharpCorners(points[i],indices,tolerance,sharpangledeg,sharpangleratio))
    {
      std::vector<std::pair<int,int>> division;
      divideByIndices(indices,int(points[i].size()),division);
      
      for (int j = 0; j < int(division.size()); j++)
      {
        std::vector<TPoint<T>> part(points[i].begin() + division[j].first,
          points[i].begin() + division[j].second + 1);

        removeDuplicates(part,false,tolerance);

        if (!part.empty())
          newpoints.push_back(part);
      }
    } else
    {
      newpoints.push_back(points[i]);
    }
  }
}

template <class T> bool segmentLenMinMax(std::vector<TPoint<T>> &line, T &min, T &max, 
  int *imin = nullptr, int *imax = nullptr)
{
  if (line.size() < 2)
    return false;

  min = std::numeric_limits<T>::max();
  max = 0.0;
  int index0 = -1;
  int index1 = -1;

  for (int i = 0; i < int(line.size() - 1); i++)
  {
    TPoint<T> p0 = line[i];
    TPoint<T> p1 = line[i + 1];
    T len = !(p1 - p0);

    if (len < min)
    {
      min = len;
      index0 = i;
    }
    if (len > max)
    {
      max = len;
      index1 = i;
    }
  }

  if (imin)
    *imin = index0;
  if (imax)
    *imax = index1;

  return true;
} 

/** Find projection of point on points segments. */
template <class T> bool projectPointOnPoints(std::vector<TPoint<T>> &points, TPoint<T> p, TPoint<T> &proj, 
  int *seg = nullptr, T *Useg = nullptr, T parmtolerance = PARM_TOLERANCE)
{
  if (seg)
    *seg = -1;
  if (Useg)
    *Useg = 0.0;
  T minDist = std::numeric_limits<T>::max();

  bool found = false;
  for (int i = 0; i < points.size() - 1; i++)
  {
    T t = 0;
    TPoint<T> intr;
    if (projectPointOnSegment(p,points[i],points[i + 1],&intr,&t,parmtolerance))
    {
      T dist = !(intr - p);
      if (dist < minDist)
      {
        minDist = dist;
        if (seg)
          *seg = i;
        if (Useg)
        {
          *Useg = t;
          LIMIT(*Useg,0.0,1.0);
        }
        proj = intr;
        found = true;
      }
    }
  }

  return found;
}

/** Straight line?. It is assumed the line points are ordered and its first/last points are
  line ends. */
template <class T> bool straightLine(std::vector<TPoint<T>> &line, T toleranceCoef = 0.001)
{
  if (line.size() < 3)
    return true;

  T len = calculateLength(line);
  T tolerance = len * toleranceCoef;

  TPoint<T> V0 = line.front();
  TPoint<T> V1 = line.back();

  for (int i = 1; i < line.size() - 1; i++)
  {
    TPoint<T> intr;
    T t;
    TPoint<T> V = line[i];
    if (projectPointOnSegment(V,V0,V1,&intr,&t,0.0))
    {
      T dist = !(V - intr);
      if (dist > tolerance)
        return false;
    } else
    {
      return false;
    }
  }

  return true;
}

/** Get mass centre for points. */
template <class T> bool calculateCentre(std::vector<TPoint<T>> &points, TPoint<T> &centre)
{
  centre = TPoint<T>(0.0,0.0,0.0);

  if (points.size() < 1)
    return false;

  for (auto c : points)
  {
    centre += c;
  }

  centre /= T(points.size());

  return true;
}

/** Close points. */
template <class T> bool closePoints(std::vector<TPoint<T>> &points, T tolerance)
{
  T dist = !(points.back() - points.front());
  if (dist > tolerance)
  {
    points.push_back(points.front());
    return true;
  } else
  {
    return false;
  }
}

/** Find all intersections by plane, returns list of parameters along point line.
  //!!! The points MUST have a parameter value in W. */
template <class T> int intersectByPlane(std::vector<TPoint<T>> &points, TPlane<T> &plane, 
  std::vector<T> &Upoints, T tolerance, T parmtolerance = PARM_TOLERANCE)
{
  Upoints.clear();

  if (points.size() < 2)
    return 0;

  for (int i = 0; i < points.size() - 1; i++)
  {
    TPoint<T> p0 = points[i];
    TPoint<T> p1 = points[i + 1];
    TPoint<T> intr;
    T U = 0.5;
    bool ok = plane.segmentIntersect(p0,p1,&intr,&U,tolerance);
    if (ok && U > 0.0 - parmtolerance && U < 1.0 + parmtolerance)
    {
      LIMIT(U,0.0,1.0);
  
      T Up = p0.W + (p1.W - p0.W) * U;
      Upoints.push_back(Up);
    }
  }

  return int(Upoints.size());
}

// comparator for equality
template <class T> bool comparePoint(TPoint<T> p1, TPoint<T> p2)
{
  if (p1.X < p2.X)
  {
    return true;
  } else if (p1.X > p2.X)
  {
    return false;
  } else
  {
    if (p1.Y < p2.Y)
    {
      return true;
    } else if (p1.Y > p2.Y)
    {
      return false;
    } else
    {
      if (p1.Z < p2.Z)
      {
        return true;
      } else if (p1.Z > p2.Z)
      {
        return false;
      } else
      {
        return false;
      }
    }
  }
}

/** Exclude duplicates between points. 

  sortcoords = false for exclusion of node duplicates among close point neighbours (use it for
    lines of points)

  sortcoords = true for exclusion and renumbering all nodes like in triangles with duplicate 
    nodes : nodes are sorted by XYZ, their order is not restored

  replacement is an array to establish corresponence between old and new node numbers after 
  node exclusions :
  int newnode = replacement[oldnode];
*/
template <class T> bool removeDuplicates(std::vector<TPoint<T>> &points, bool sortcoords, 
  T tolerance, std::vector<LINT> *replacement = nullptr)
{
                              // list is empty
  if (points.size() < 2) 
    return true;
                              // just one vector
  if (points.size() == 1)
  {
    if (replacement)
      replacement->push_back(0);
    return true;
  }
                              // put index number in W as original point indices
  for (size_t i = 0; i < points.size(); i++)
  {
    points[i].W = static_cast<T>(i);
  }
                              // sort points by X,Y,Z
  if (sortcoords)
    std::sort(points.begin(),points.end(),comparePoint<T>);

                              // replacement indices for node numbering
  std::vector<LINT> rep;
  rep.resize(points.size(),-1);

#if 1 // for very big arrays : temporarily replace vectors by map, 
// vectors erase is too slow for e.g. 240000 elements

  std::map<int,TPoint<T>> list;
  for (int i = 0; i < int(points.size()); i++)
  {
    list.insert(std::pair<int,TPoint<T>>(i,points[i]));
  }

  // now move down and remove duplicates which are neighbours to each other
  for (int i = int(points.size()) - 1; i >= 0; i--)
  {
    // len is number of same points closed to each other
    int len = 0;
    for (int j = i; j >= 0; j--)
    {
      // len cannot be 0, it includes the leftmost node as well
      T dist = !(points[i] - points[j]);
      if (dist <= tolerance)
      {
        len++;
      } else
      {
        break;
      }
    }

    // this number is leftmost remaining, all the rest to the right (len - 1) 
    // points are deleted
    int i1 = i - len + 1;
    int count = len;
    for (int j = i; j >= 0; j--)
    {
      rep[ROUND(points[j].W)] = i1;
      if (--count == 0)
        break;
    }

    // delete len - 1 nodes to the right from i1;
    // we are deleting list - correct replacement indices 
    if (len > 1)
    {
      for (int j = len - 1; j >= 1; j--)
      {
        // removing node n
        int n = i1 + j;
        list.erase(n);

        // correct replacements
        for (auto &k : rep)
        {
          if (k >= n)
            k--;
        }
      }
    }

    // skip len - 1 nodes
    i -= (len - 1);
  }

  points.clear();
  for (auto p : list)
  {
    points.push_back(p.second);
  }

#else

  // now move down and remove duplicates which are neighbours to each other
  for (int i = int(points.size()) - 1; i >= 0; i--)
  {
    // len is number of same points closed to each other
    int len = 0;
    for (int j = i; j >= 0; j--)
    {
      // len cannot be 0, it includes the leftmost node as well
      T dist = !(points[i] - points[j]);
      if (dist <= tolerance)
      {
        len++;
      } else
      {
        break;
      }
    }

    // this number is leftmost remaining, all the rest to the right (len - 1) 
    // points are deleted
    int i1 = i - len + 1;
    int count = len;
    for (int j = i; j >= 0; j--)
    {
      rep[ROUND(points[j].W)] = i1;
      if (--count == 0)
        break;
    }

    // delete len - 1 nodes to the right from i1;
    // we are deleting points - correct replacement indices 
    if (len > 1)
    {
      for (int j = len - 1; j >= 1; j--)
      {
        // removing node n
        int n = i1 + j;
        points.erase(points.begin() + n);

        // correct replacements
        for (auto &k : rep)
        {
          if (k >= n)
            k--;
        }
      }
    }

    // skip len - 1 nodes
    i -= (len - 1);
  }

#endif

  // all replacements must be filled up
#ifdef _DEBUG
  for (auto r : rep)
  {
    assert(r >= 0);
  }
#endif

  if (replacement)
    *replacement = rep;

  return true;
}

/** Make transform. */
template <class T> void makeTransform(std::vector<TPoint<T>> &points, TTransform<T> *transform)
{
  for (auto &p : points)
  {
    p = transform->applyTransform(p);
  }
}

/** Make circle in XY plane. */
template <class T> void makeCircleXY(int numpoints, TPoint<T> centre, T R, std::vector<TPoint<T>> &points,
  T adegfrom = 0.0, T adegto = 360.0)
{
  int numsegments = numpoints - 1;
  T da = (adegfrom - adegto) / T(numsegments);

  for (int i = 0; i <= numsegments; i++)
  {
    T a = T(i) * da;
    TPoint<T> p(R * cos(a * CPI),R * sin(a * CPI));
    points.push_back(centre + p);
  }
}

/** Make ellipse in XY plane. */
template <class T> void makeEllipseXY(int numpoints, TPoint<T> centre, T a, T b, std::vector<TPoint<T>> &points,
  T adegfrom = 0.0, T adegto = 360.0)
{
  makeCircleXY(numpoints,TPoint<T>(0,0,0),a,points,adegfrom,adegto);

  TTransform<T> t0;
  t0.Resize(TPoint<T>(1.0,b / a,1.0));
  makeTransform(points,&t0);

  TTransform<T> t1;
  t1.Translate(centre);
  makeTransform(points,&t1);
}

}
