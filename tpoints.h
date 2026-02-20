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
#include "strings.h"

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

template <class T> std::vector<TPoint<T>> operator * (const std::vector<TPoint<T>> &points0, const T coef)
{
  std::vector<TPoint<T>> points;
  for (int i = 0; i < int(points0.size()); i++)
  {
    points.push_back(points0[i] * coef);
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
  TPoint<T> dir = (DU > PARM_TOLERANCE) ? (d / DU) : TPoint<T>();
  return dir;
}

/** Get end direction (first derivative on U) for point list, directed from end "inside". */
template <class T> TPoint<T> endDirection(std::vector<TPoint<T>> &points)
{
  TPoint<T> d = points[points.size() - 2] - points[points.size() - 1];
  T len = !d;
  T curvelen = calculateLength(points);
  T DU = len / curvelen;
  TPoint<T> dir = (DU > PARM_TOLERANCE) ? (d / DU) : TPoint<T>();
  return dir;
}

/** Get normalised direction at the middle. */
template <class T> TPoint<T> midDirection(std::vector<TPoint<T>> &points)
{
  assert(points.size() > 1);

  int i0 = int(points.size() / 2) - 1;
  LIMIT_MIN(i0,0);
  int i1 = i0 + 1;
  LIMIT_MAX(i1,int(points.size() - 1));

  // suppose no duplicates
  TPoint<T> dir = +(points[i1] - points[i0]);
  return dir;
}

/** Get point at the middle. */
template <class T> TPoint<T> midPoint(std::vector<TPoint<T>> &points)
{
  assert(points.size() > 0);

  int i = int(points.size() / 2);

  return points[i];
}

/** Get next increasing element of an array. */
template <class T> int next(int size, int index)
{
  if (++index >= int(size))
  {
    index = 0;
  }
  return index;
}

/** Get previuos element of an array. */
template <class T> int prev(int size, int index)
{
  if (--index < 0)
  {
    index = size - 1;
  }
  return index;
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

/** Calculate length TILL segment, i.e. length is 0 if seg = 0, the whole line length
  if seg = points.size() - 1. */
template <class T> T calculateLength(std::vector<TPoint<T>> &points, int seg)
{
  T len = 0.0;

  if (points.size() > 0)
  {
    for (int i = 0; i < seg; i++)
    {
      len += !(points[i + 1] - points[i]);
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

//!!!!!!! very important
  //if (points.size() == 2 && !(points[1] - points[0]) < tolerance)
  //{
  //  return false;
  //}

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

/** Get uniform parameters [0..1]. */
template <class T> void prepareUniformParameters(int numpoints, 
  std::vector<T> &parms, T refinestartU = 1.0, T refineendU = 1.0)
{
  parms.clear();

  for (int i = 0; i < numpoints; i++)
  {
    T U = T(i) / T(numpoints - 1);

    U = refineParameter(U,refinestartU,refineendU);

    parms.push_back(U);
  }
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

/** Is point inside boundary? Ray casting. */
template <class T> bool insideBoundary(std::vector<std::vector<TPoint<T>>> &boundary, TPoint<T> point, 
  T parmtolerance = PARM_TOLERANCE)
{
  std::vector<TPoint<T>> intrs;

  for (int i = 0; i < int(boundary.size()); i++)
  {
    for (int j = 0; j < int(boundary[i].size()) - 1; j++)
    {
      TPoint<T> p0 = boundary[i][j];
      TPoint<T> p1 = boundary[i][j + 1];

      // right on the boundary?
      T t = 0;
      TPoint<T> intr;
      if (projectPointOnSegment(point,p0,p1,&intr,&t,parmtolerance))
      {
        T d = !(point - intr);
        if (d < parmtolerance)
          return true;
      }

      // make ray casting
      T t1 = 0.0;
      T t2 = 0.0;
      T Xi = p1.X;
      T Yi = p1.Y;

//      if (intersectSegmentsXY(p0.X,p0.Y,p1.X,p1.Y,point.X,point.Y,point.X + 10.0,point.Y + 10.0,
      if (intersectSegmentsXY(p0.X,p0.Y,p1.X,p1.Y,point.X,point.Y,point.X,point.Y + 10.0, // important to cast ray along V
        &t1,&t2,&Xi,&Yi) && (t1 >= 0.0) && (t1 <= 1.0) && (t2 >= 0.0)) // t2 >= 0 - single ray direction
      {
        intrs.push_back(TPoint<T>(Xi,Yi));
      }

      // old code :
      //T t1,t2,Ui,Vi;
      //bool res = s->IntersectSegmentUV(0,U,V,U,V + 10.0,&t1,&t2,&Ui,&Vi,0.0); 

      //if (res)
      //{
      //  T dist = GetDistFromPointToLine(TPoint<T>(U,V),TPoint<T>(s->tu[0][0],s->tv[0][0]),
      //    +(TPoint<T>(s->tu[0][1],s->tv[0][1]) - TPoint<T>(s->tu[0][0],s->tv[0][0])));
      //  if (dist > minparmdist)
      //  { 
      //    intrs.push_back(TPoint<T>(Ui,Vi));
      //  }
      //}
    }
  }

  // exclude equal intersections
  if (intrs.size() == 0)
  {
    return false;
  } else if (intrs.size() == 1)
  {
    return true;
  } else
  {
    bool found = false;

    do {
      found = false;
      for (int i = 0; i < intrs.size() - 1; i++)
      {
        for (int j = i + 1; j < intrs.size(); j++)
        {
          T dist = !(intrs[i] - intrs[j]);
          if (dist < parmtolerance)
          {
            intrs.erase(intrs.begin() + j);
            found = true;
            break;
          }
        }
        if (found)
          break;
      }
    } while (found);

    return (intrs.size() % 2 == 1);
  }
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

  if (points0.empty() || points1.empty())
    return 0;

  for (int i = 0; i < points0.size() - 1; i++)
  {
    TPoint<T> p0 = points0[i];
    TPoint<T> p1 = points0[i + 1];
    T d0 = !(p1 - p0);

    for (int j = 0; j < points1.size() - 1; j++)
    {
      TPoint<T> v0 = points1[j];
      TPoint<T> v1 = points1[j + 1];
      T d1 = !(v1 - v0);
    
      T t1 = 0.0;
      T t2 = 0.0;
      T dist = 0.0;
      TPoint<T> ip,iv;
      if (intersectSegments(p0,p1,v0,v1,t1,t2,dist,&ip,&iv))
      {
#if 0
        T parmtolerance0 = (d0 < TOLERANCE(T)) ? parmtolerance : parmtolerance / d0; //!!!!!!!
        T parmtolerance1 = (d1 < TOLERANCE(T)) ? parmtolerance : parmtolerance / d1;
#else
        T parmtolerance0 = parmtolerance;
        T parmtolerance1 = parmtolerance;
#endif

        if ((t1 >= 0.0 - parmtolerance0) && (t1 <= 1.0 + parmtolerance0) && 
          (t2 >= 0.0 - parmtolerance1) && (t2 <= 1.0 + parmtolerance1))
        { 
          if (dist < tolerance) //!!!
          {
            LIMIT(t1,0.0,1.0);
            LIMIT(t2,0.0,1.0);

            T U1 = p0.W + (p1.W - p0.W) * t1;
            T U2 = v0.W + (v1.W - v0.W) * t2;
            UV.push_back(TPoint<T>(U1,U2));
          }
        }
      }
    }
  }

  // remove all duplicates
  removeDuplicates(UV,true,parmtolerance * 4.0); //!!!!!!!

  return int(UV.size());
}

/** Find edge closest to edges. Set maxedge to -1.0 to avoid corresponding checks. */
template <class T> bool findClosestNotBusy(std::vector<std::pair<TPoint<T>,TPoint<T>>> &edges, 
  TPoint<T> value, std::vector<bool> &busy, int &index, bool &reversed, T maxedge, 
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
    }
  }

  if (!ok)
    min = 0.0;

  if (mindist)
    *mindist = min;

  return ok;
}

/** Edges are directed toward intersection line, so they are defined only by
  coordinates of the first node. */
template <class T> void removeEdgeDuplicates(std::vector<std::pair<TPoint<T>,TPoint<T>>> &edges,
  T tolerance)
{
  // take only first points, they uniquely define a directed edge
  std::vector<TPoint<T>> points;
  for (auto &e : edges)
  {
    points.push_back(e.first);
  }

  // exclude duplicates, take replacement
  std::vector<LINT> replacement;
  removeDuplicates(points,true,tolerance,&replacement);

  std::vector<std::pair<TPoint<T>,TPoint<T>>> newedges = edges;

  for (int i = 0; i < int(edges.size()); i++)
  {
    newedges[replacement[i]] = edges[i];
  }

  // cut out tail
  newedges.resize(points.size());
  edges = newedges;
}

/** Make up curve(s) from unordered pieces, like hanging edges on triangles boundary
  or intersection of curve pieces. maxedge and tolerance are usually big, like 0.01 
  of model size or 0.1 of max edge, this value to separate impossible cases while 
  the connection algorithm being "find the closest", not "take the first within tolerance". */
template <class T> bool makeUpCurves(std::vector<std::pair<TPoint<T>,TPoint<T>>> &pedges,
  T maxedge, std::vector<std::vector<TPoint<T>>> &lines, T tolerance, bool bothways = true, 
  bool removeedgeduplicates = false, bool checkdirection = false)
{
  std::vector<std::pair<TPoint<T>,TPoint<T>>> edges = pedges;

  if (edges.empty())
    return false;

  if (removeedgeduplicates)
    removeEdgeDuplicates(edges,tolerance);

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
  findClosestNotBusy(edges,middle,busy,startedge,reversed,maxedge);

  std::vector<TPoint<T>> line;

  line.push_back(edges[startedge].first);
  line.push_back(edges[startedge].second);
  busy[startedge] = true;

  while (!allBusy(busy))
  {
#if 1
    // this fails only if all edges are busy; if distance > maxedge, a negative index is returned
    bool reversed0 = false;
    int index0 = -1;
    T mindist0 = 0.0;
    bool ok0 = findClosestNotBusy(edges,line.back(),busy,index0,reversed0,maxedge,&mindist0);

    bool reversed1 = false;
    int index1 = -1;
    T mindist1 = 0.0;
    bool ok1 = findClosestNotBusy(edges,line.front(),busy,index1,reversed1,maxedge,&mindist1);

    if (ok0 && ok1)
    {
      if (mindist0 < mindist1)
      {
        ok1 = false;
      } else
      {
        ok0 = false;
      }
    }

    if (!ok0)
    {
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
      if (reversed0)
      {
        line.push_back(edges[index0].first);
        busy[index0] = true;
      } else
      {
        line.push_back(edges[index0].second);
        busy[index0] = true;
      }
    }

#else
    // this fails only if all edges are busy; if distance > maxedge, a negative index is returned
    bool reversed = false;

    int index = -1;
    T mindist0 = 0.0;
    bool ok0 = findClosestNotBusy(edges,line.back(),busy,index,reversed,maxedge,&mindist0);

    if (!ok0 || (checkdirection && reversed))
    {
      if (bothways)
      {
        bool reversed1 = false;
        int index1 = -1;
        bool ok1 = findClosestNotBusy(edges,line.front(),busy,index1,reversed1,maxedge);
        if (!ok1 || (checkdirection && reversed1))
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
#endif
  }

  if (!line.empty())
    lines.push_back(line);

  return true;
}

/** Make up curves from pieces. In case of NOT degenerateedges,
  pieces contain many pieces ot two points generated by intersections.
  If degenerateedges, pieces.size() is 1 which contains many isolated
  points each of which will be duplicated into degenerated edges.
  maxedge ratio is a guessed max edge size to the model size,
  it is unknown in case edges are points (degenerated edges). */
template <class T> bool curvesFromPieces(std::vector<std::vector<TPoint<T>>> &pieces,
  std::vector<std::vector<TPoint<T>>> &lines, T tolerance, 
  bool degenerateedges = false, T maxedgeratio = MAXEDGE_RATIO, bool removeedgeduplicates = false)
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

  bool ok = makeUpCurves(edges,maxedge,lines,tolerance,true,removeedgeduplicates);

  return ok;
}

/** Make up single curve from pieces. In case of NOT degenerateedges,
  pieces contain many pieces ot two points generated by intersections.
  If degenerateedges, pieces.size() is 1 which contains many isolated
  points each of which will be duplicated into degenerated edges. 
  maxedge ratio is a guessed max edge size to the model size,
  it is unknown in case edges are points (degenerated edges). */
template <class T> bool curveFromPieces(std::vector<std::vector<TPoint<T>>> &pieces,
  std::vector<TPoint<T>> &line, T tolerance, bool degenerateedges = false, T maxedgeratio = MAXEDGE_RATIO)
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

/** Exclude duplicate pieces and order points. */
template <class T> bool orderPoints(std::vector<TPoint<T>> &points,
  std::vector<std::vector<TPoint<T>>> &lines, T tolerance, T maxedgeratio = MAXEDGE_RATIO)
{
  // step 1 : exclude duplicates
  removeDuplicates(points,true,tolerance);

  // step 2 : order points
  std::vector<std::vector<TPoint<T>>> pieces(1);
  for (auto &p : points)
  {
    pieces[0].push_back({p});
  }

  bool ok = curvesFromPieces(pieces,lines,tolerance,true,maxedgeratio);

  return ok;
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
template <class T> TPoint<T> calculateCentre(std::vector<TPoint<T>> &points)
{
  assert(!points.empty());

  TPoint<T> centre = TPoint<T>(0.0,0.0,0.0);

  for (auto c : points)
  {
    centre += c;
  }

  centre /= T(points.size());

  return centre;
}

/** Get normal for curve orientation direction. https://en.wikipedia.org/wiki/Curve_orientation */
template <class T> TPoint<T> calculateNormalXY(std::vector<TPoint<T>> &points)
{
  assert(points.size() > 2);

  // find point with min X and min Y
  int index = -1;
  T xmin = std::numeric_limits<T>::max();
  T ymin = std::numeric_limits<T>::max();
  for (int i = 0; i < int(points.size()); i++)
  {
    if (points[i].X < xmin)
    {
      xmin = points[i].X;
      ymin = points[i].Y;
      index = i;
    } else if (points[i].X == xmin && points[i].Y < ymin)
    {
      xmin = points[i].X;
      ymin = points[i].Y;
      index = i;
    }
  }

  int i0 = prev<T>(int(points.size()),index);
  int i1 = next<T>(int(points.size()),index);

  TPoint<T> d0,d1;
  for (int a = 0; a < 5; a++)
  {
    d0 = points[index] - points[i0];
    if ((!d0) > TOLERANCE(T))
      break;
    i0 = prev<T>(int(points.size()),i0);
  }
  for (int a = 0; a < 5; a++)
  {
    d1 = points[i1] - points[index];
    if ((!d1) > TOLERANCE(T))
      break;
    i1 = next<T>(int(points.size()),i1);
  }

  assert((!d0) > TOLERANCE(T));
  assert((!d1) > TOLERANCE(T));

  TPoint<T> normal = +(d0 ^ d1);

  return normal;
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

// important
static double ptolerance = 0.0;

template <class T> bool comparePoint(TPoint<T> p1, TPoint<T> p2)
{
  if (p1.X < p2.X - T(ptolerance))
  {
    return true;
  } else if (p1.X > p2.X + T(ptolerance))
  {
    return false;
  } else
  {
    if (p1.Y < p2.Y - T(ptolerance))
    {
      return true;
    } else if (p1.Y > p2.Y + T(ptolerance))
    {
      return false;
    } else
    {
      if (p1.Z < p2.Z - T(ptolerance))
      {
        return true;
      } else if (p1.Z > p2.Z + T(ptolerance))
      {
        return false;
      } else
      {
        return false;
      }
    }
  }
}

/** Sort two long integers. */
static int TwoIntSort(const void *E0, const void *E1)
{
  LINT *i0 = (LINT *) E0;
  LINT *i1 = (LINT *) E1;

  if (*i0 > *i1) 
  {
    return +1; 
  } else if (*i0 < *i1)
  {
    return -1; 
  } else 
  {
    return 0;
  }
}

/** Exclude duplicates from 3D array (fast and reliable) */
template <class T> static bool removeDupNodes(std::vector<TPoint<T>> &points, 
  std::vector<LINT> &replacement, T tolerance)
{
                              // list is empty
  if (points.size() < 1) 
    return true;
                              // just one vector
  if (points.size() == 1)
  {
    replacement.push_back(0);
    return true;
  }
                              // grid step
  T gridstep = tolerance * T(2.0);
                              // max possible coordinate value
  T maxcoord = 0.0;

  TPoint<T> *v = &points[0];
  for (size_t i = 0; i < points.size(); i++) 
  {
    if (fabs(v->X) > maxcoord) maxcoord = fabs(v->X);
    if (fabs(v->Y) > maxcoord) maxcoord = fabs(v->Y);
    if (fabs(v->Z) > maxcoord) maxcoord = fabs(v->Z);
    v++;
  }
                              // get max integer value we need
  T numcellsneeded = maxcoord / gridstep;
  T numcellsavailable = static_cast<T>(std::numeric_limits<LINT>::max());
  if (numcellsavailable < numcellsneeded)
  {
                              // increase tolerance
    tolerance = maxcoord / numcellsavailable;
    gridstep = tolerance * T(2.0);
  }
                              // array will hold integer coordinate + vertex
                              // number and will be sorted by coordinate
  size_t numvectors = points.size();
  LINT knsize = sizeof(LINT) * (LINT) numvectors * 2;
  LINT *cnumber = (LINT *) malloc(knsize);
  if (cnumber == nullptr) return false;
                              // zero out
  memset(cnumber,0,knsize);
  replacement.resize(points.size(),0);
                              // loop over all 8 possible positions of the 
                              // 2*tolerance cube
  static T incs[8][3] =
  {
    -0.5,-0.5,-0.5,
    +0.5,-0.5,-0.5,
    +0.5,+0.5,-0.5,
    -0.5,+0.5,-0.5,
    -0.5,-0.5,+0.5,
    +0.5,-0.5,+0.5,
    +0.5,+0.5,+0.5,
    -0.5,+0.5,+0.5
  };

  for (int pos = 0; pos < 8; pos++) 
  {
    T dx = tolerance * incs[pos][0];
    T dy = tolerance * incs[pos][1];
    T dz = tolerance * incs[pos][2];
                              // set up sort array
    LINT numtosort = 0;
    for (size_t i = 0; i < numvectors; i++) 
    {
      if (replacement[i] == 0) 
      {
        cnumber[numtosort * 2] = static_cast<LINT>((points[i].X + dx) / gridstep);
        cnumber[numtosort * 2 + 1] = (LINT) i;
        numtosort++;
      }
    }
                              // sort then match x
    qsort(cnumber,numtosort,sizeof(LINT) * 2,TwoIntSort);

    LINT xmin = cnumber[0];
    for (LINT i = 0; i < numtosort; ) 
    {
      LINT j = i;
      while (j < numtosort && (cnumber[j * 2] == xmin)) 
      {
                              // replace with y
        cnumber[j * 2] = static_cast<LINT>((points[cnumber[j * 2 + 1]].Y + dy) / gridstep);
        j++;
      }
                              // sort then match y
      if ((j - i) > 1) 
      {
        qsort(cnumber + i * 2,j - i,sizeof(LINT) * 2,TwoIntSort);

        LINT ymin = cnumber[2 * i];
        for (LINT k = i; k < j; ) 
        {
          LINT l = k;
          while (l < j && (cnumber[l * 2] == ymin)) 
          {
                              // replace with z
            cnumber[l * 2] = static_cast<LINT>((points[cnumber[l * 2 + 1]].Z + dz) / gridstep);
            l++;
          }
                              // sort then match z
          if ((l - k) > 1) 
          {
            qsort(cnumber + k * 2,l - k,sizeof(LINT) * 2,TwoIntSort);

            LINT zmin = cnumber[k * 2];
            for (LINT m = k; m < l; ) 
            {
              LINT p = cnumber[m * 2 + 1];
              LINT n = m;
              while (n < l && (cnumber[n * 2] == zmin)) 
              {
                if (n > m) 
                {
                              // point at first node in matches; do not move if 
                              // face indices not specified
                  LINT index = cnumber[n * 2 + 1];
                  replacement[index] = -1 - p;
                }
                              // increment
                n++;
              }
                              // mark node as paired
              if ((n - m) > 1) replacement[p] = 1;

              m = n; 
              if (m < l) zmin = cnumber[m * 2];
            }
          }
          k = l; 
          if (k < j) ymin = cnumber[k * 2];
        }
      }
      i = j; 
      if (i < numtosort) xmin = cnumber[i * 2];
    }
  } // 8 positions
                              // eliminate nodes from the points array and 
                              // complete replacement
  LINT k = 0;

  for (size_t i = 0; i < numvectors; i++) 
  {
    if (replacement[i] >= 0) 
    {
      replacement[i] = k;

      points[k] = points[i];
      k++;
    }
  }

  for(int i = 0; i < numvectors; i++) 
  {
    if (replacement[i] < 0) 
    {
      replacement[i] = replacement[-1 - replacement[i]];
    }
  }

  points.resize(k);

  free(cnumber);

  return true;
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
  // use a more reliable way
  if (sortcoords)
  {
    std::vector<LINT> rep;
    bool ok = removeDupNodes(points,rep,tolerance);
    if (replacement)
      *replacement = rep;

    return ok;
  }
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
  {
//!!!!!!    ptolerance = TOLERANCE(T);
    ptolerance = double(tolerance);
    std::sort(points.begin(),points.end(),comparePoint<T>);
  }
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

/** Make transform. */
template <class T> void makeTransform(std::vector<std::vector<TPoint<T>>> &points, TTransform<T> *transform)
{
  for (auto &list : points)
  {
    for (auto &p : list)
    {
      p = transform->applyTransform(p);
    }
  }
}

/** Make circle in XY plane. */
template <class T> void makeCircleXY(int numpoints, TPoint<T> centre, T R, std::vector<TPoint<T>> &points,
  T adegfrom = 0.0, T adegto = 360.0)
{
  int numsegments = numpoints - 1;
  T da = (adegto - adegfrom) / T(numsegments);

  for (int i = 0; i <= numsegments; i++)
  {
    T a = adegfrom + T(i) * da;
    TPoint<T> p(R * cos(a * CPI),R * sin(a * CPI));
    points.push_back(centre + p);
  }
}

/** Make ellipse in XY plane. */
template <class T> void makeEllipseXY(int numpoints, TPoint<T> centre, T a, T b, std::vector<TPoint<T>> &points,
  T adegfrom = 0.0, T adegto = 360.0)
{
  makeCircleXY(numpoints,TPoint<T>(0,0,0),a,points,adegfrom,adegto);

  if (a > TOLERANCE(T))
  {
    TTransform<T> t0;
    t0.Resize(TPoint<T>(1.0,b / a,1.0));
    makeTransform(points,&t0);
  }

  TTransform<T> t1;
  t1.Translate(centre);
  makeTransform(points,&t1);
}

/** Make ellipse in XY plane. */
template <class T> void makeEllipseXY(int numpoints, int numparts, TPoint<T> centre, T a, T b, 
  std::vector<std::vector<TPoint<T>>> &points, T adegfrom = 0.0, T adegto = 360.0)
{
  assert(numparts >= 1);
  LIMIT_MIN(numparts,1);

  T da = (adegto - adegfrom) / T(numparts);

  for (int i = 0; i < numparts; i++)
  {
    T a0 = T(i) * da;
    T a1 = T(i + 1) * da;

    std::vector<TPoint<T>> part;
    makeEllipseXY(numpoints,centre,a,b,part,a0,a1);

    points.push_back(part);
  }
}

/** Cut out piece of points from 0 to 1. Last possible seg1 must be points.size() - 2. */
template <class T> void cutOut(std::vector<TPoint<T>> &points, int seg0, T U0, int seg1, T U1,
  std::vector<TPoint<T>> &newpoints)
{
  if (seg1 == points.size() - 1)
  {
    seg1 = int(points.size()) - 2;
    U1 = 1.0;
  }

  for (int i = seg0; i <= seg1; i++)
  {
    TPoint<T> p0 = points[i];
    TPoint<T> p1 = points[i + 1];
    TPoint<T> dp = p1 - p0;
    TPoint<T> p = p0;

    if (seg0 == seg1)
    {
      if (i == seg0)
      {
        p = p0 + dp * U0;
        newpoints.push_back(p);
      }

      if (i == seg1)
      {
        p = p0 + dp * U1;
        newpoints.push_back(p);
      }
    } else
    {
      if (i == seg0)
      {
        p = p0 + dp * U0;
      }

      if (i == seg1)
      {
        p = p0 + dp * U1;
      }

      newpoints.push_back(p);
    }
  }
}

/** Extend curve length by moving start or end by dlen. */
template <class T> void extendByLength(std::vector<TPoint<T>> &points, bool extendstart, T dlen)
{
  assert(points.size() > 1);

  if (dlen > 0.0)
  {
    T len = calculateLength(points);

    if (extendstart)
    {
      TPoint<T> p0 = points[0];
      TPoint<T> p1 = points[1];
      TPoint<T> dir = p0 - p1;
      T d = !dir;
      T dnew = d + dlen;
      dir = +dir;
      points[0] = p1 + dir * dnew;
    } else
    {
      TPoint<T> p0 = points[points.size() - 2];
      TPoint<T> p1 = points[points.size() - 1];
      TPoint<T> dir = p1 - p0;
      T d = !dir;
      T dnew = d + dlen;
      dir = +dir;
      points[points.size() - 1] = p0 + dir * dnew;
    }
  }
}

/** Extend curve by moving start or end by coef, which must be > 1.0. */
template <class T> void extend(std::vector<TPoint<T>> &points, bool extendstart, T coef)
{
  assert(points.size() > 1);
  assert(coef >= 1.0);
  LIMIT_MIN(coef,1.0);

  if (coef > 1.0)
  {
    T len = calculateLength(points);
    T lennew = len * coef;
    T dlen = lennew - len;

    extendByLength(points,extendstart,dlen);
  }
}

/** Get direction at node. */
template <class T> TPoint<T> direction(std::vector<TPoint<T>> &points, int index, bool normalise = true)
{
  int i0 = index - 1;
  int i1 = index + 1;
  if (i0 < 0)
  {
    i0++;
    i1++;
  }
  if (i1 > int(points.size()) - 1)
  {
    i1--;
    i0--;
  }

  if (i0 >= 0 && i0 < int(points.size()) &&
    i1 >= 0 && i1 < int(points.size()) && i1 > i0)
  {
    TPoint<T> p0 = points[i0];
    TPoint<T> p1 = points[i1];
    TPoint<T> d = normalise ? (+(p1 - p0)) : (p1 - p0);
    return d;
  } else
  {
    assert(false);
    return TPoint<T>();
  }
}

/** Get plane on three points. */
template <class T> bool getPlanePrim(std::vector<TPoint<T>> &points, 
  int i0, int i, int i1, TPlane<T> &plane)
{
  if (points.size() < 3)
    return false;

  assert(i > i0 && i1 > i);

  if (i0 >= 0 && i0 < int(points.size()) &&
    i1 >= 0 && i1 < int(points.size()) && i1 > i0)
  {
    bool ok = true;
    plane = TPlane<T>(points[i0],points[i],points[i1],ok);
    return ok;
  } else
  {
    assert(false);
    return false;
  }
}

/** Get plane on three neighbour points. */
template <class T> bool getPlane(std::vector<TPoint<T>> &points, int index, TPlane<T> &plane)
{
  if (points.size() < 3)
    return false;

  int i = index;
  int i0 = index - 1;
  int i1 = index + 1;
  if (i0 < 0)
  {
    i0++;
    i++;
    i1++;
  }
  if (i1 > int(points.size()) - 1)
  {
    i1--;
    i--;
    i0--;
  }

  return getPlanePrim(points,i0,i,i1,plane);
}

/** Get plane on three points around index which make a plane. */
template <class T> bool getClosestPlane(std::vector<TPoint<T>> &points, int index, TPlane<T> &plane)
{
  if (points.size() < 3)
    return false;

  int i = index;
  int i0 = index - 1;
  int i1 = index + 1;
  if (i0 < 0)
  {
    i0++;
    i++;
    i1++;
  }
  if (i1 > int(points.size()) - 1)
  {
    i1--;
    i--;
    i0--;
  }

  if (getPlanePrim(points,i0,i,i1,plane))
  {
    return true;
  }

  bool done = false;
  do {
    i0--;
    LIMIT_MIN(i0,0);
    i1++;
    LIMIT_MAX(i1,int(points.size()) - 1);

    if (getPlanePrim(points,i0,i,i1,plane))
    {
      done = true;
    } else
    {
      if (i0 <= 0 && i1 >= int(points.size()) - 1)
        break;
    }
  } while (!done);

  return done;
}

/** Inflate/deflate curve in normal direction by dist. Plane on curve points
  is calculated for all points (impossible for a straight curve). */
template <class T> bool inflatePoints(std::vector<TPoint<T>> &points, T dist)
{
  assert(points.size() > 2);

  // go to along the curve and move all points "to the right"
  std::vector<TPoint<T>> newpoints;
  for (int i = 0; i < int(points.size()); i++)
  {
    // get plane on neighbour points
    TPlane<T> plane;
    if (!getClosestPlane(points,i,plane))
      return false;

    // get normalised direction along the curve
    TPoint<T> dir = direction(points,i);

    // move direction, both normal and dir are normalised
    TPoint<T> movedir = dir ^ plane.normal;

    // new position
    TPoint<T> p = points[i] + movedir * dist;

    newpoints.push_back(p);
  }

  points = newpoints;
  return true;
}

/** Make a straight line. */
template <class T> void makeStraightLine(TPoint<T> p0, TPoint<T> p1, int numpoints,
  std::vector<TPoint<T>> &points)
{
  points.clear();

  TPoint<T> d = p1 - p0;
  for (int i = 0; i < numpoints; i++)
  {
    T U = T(i) / T (numpoints - 1);
    TPoint<T> p = p0 + d * U;
    points.push_back(p);
  }
}

/** Divide into two halves. */
template <class T> void divideHalf(std::vector<TPoint<T>> &points,
  std::vector<TPoint<T>> &left, std::vector<TPoint<T>> &right)
{
  left.clear();
  right.clear();

  left.insert(left.end(),points.begin(),points.begin() + points.size() / 2 + 1);
  right.insert(right.end(),points.begin() + points.size() / 2,points.end());
}

/** Calculate Adler checksum. setInFront - set this value in W od the first element. */
template <class T> unsigned int setChecksum(std::vector<TPoint<T>> &points, bool setInFront)
{
  unsigned int sum = 0;

  for (auto &p : points)
  {
    sum += p.checksum();
  }

  if (setInFront)
    points.front().W = sum;

  return sum;
}

/** Two starting points coincident? */
template <class T> bool startCoincident(std::vector<TPoint<T>> &points, T parmtolerance = PARM_TOLERANCE)
{
  if (points.size() < 2)
    return false;

  TPoint<T> p0 = points[0];
  TPoint<T> p1 = points[1];
  T d = !(p1 - p0);
  return (d < parmtolerance);
}

/** Two ending points coincident? */
template <class T> bool endCoincident(std::vector<TPoint<T>> &points, T parmtolerance = PARM_TOLERANCE)
{
  if (points.size() < 2)
    return false;

  TPoint<T> p0 = points[points.size() - 1];
  TPoint<T> p1 = points[points.size() - 2];
  T d = !(p1 - p0);
  return (d < parmtolerance);
}

/** Correct coincident (it happens after improvement) parametric values at ends by moving 
  the second point to the middle between neighbours. */
template <class T> void correctEndingParms(std::vector<TPoint<T>> &points,
  bool correctstart, bool correctend)
{
  if (points.size() < 3)
    return;

  if (correctstart)
  {
    TPoint<T> p0 = points[0];
    TPoint<T> p1 = points[1];
    points[1] = (points[0] + points[2]) * 0.5;
  }

  if (correctend)
  {
    TPoint<T> p0 = points[points.size() - 1];
    TPoint<T> p1 = points[points.size() - 2];
    points[points.size() - 2] = (points[points.size() - 1] + points[points.size() - 3]) * 0.5;
  }
}

/** Set parm value as middle between neighbours. */
template <class T> void correctParm(std::vector<TPoint<T>> &points, int index)
{
  if (index > 0 && index < int(points.size()) - 1)
  {
    points[index] = (points[index - 1] + points[index + 1]) * 0.5;
  }
}

/** Distance between points. */
template <class T> inline T distanceBwPoints(std::vector<TPoint<T>> &points, int i0, int i1)
{
  TPoint<T> p0 = points[i0];
  TPoint<T> p1 = points[i1];
  T d = !(p1 - p0);
  return d;
}

/** Correct coincident (it happens after improvement) parametric values by moving 
  the second point to the middle between neighbours. All changes are applied to
  both equal-size parametric arrays at once. */
template <class T> bool correctParmsParallel(std::vector<TPoint<T>> &points0, std::vector<TPoint<T>> &points1,
  T parmtolerance)
{
  assert(points0.size() == points1.size());
  if (points0.size() != points1.size())
    return false;

  for (int i = 1; i < points0.size(); i++)
  {
    int i0 = i - 1;
    int i1 = i;

    T d0 = distanceBwPoints(points0,i0,i1);
    T d1 = distanceBwPoints(points1,i0,i1);

    if (d0 < parmtolerance || d1 < parmtolerance)
    {
      if (i < int(points0.size()) / 2)
      {
        correctParm<T>(points0,i);
        correctParm<T>(points1,i);
      } else
      {
        correctParm<T>(points0,i - 1);
        correctParm<T>(points1,i - 1);
      }
    }
  }

  return true;
}

/** Calculate number of duplicates. */
template <class T> int numDuplicates(std::vector<TPoint<T>> &points, T tolerance = PARM_TOLERANCE)
{
  int count = 0;
  for (int i = 0; i < int(points.size()) - 1; i++)
  {
    T d = !(points[i + 1] - points[i]);

    if (d < tolerance)
      count++;
  }

  return count;
}

/** Calculate number of triplicates. */
template <class T> int numTriplicates(std::vector<TPoint<T>> &points, T tolerance = PARM_TOLERANCE)
{
  int count = 0;
  for (int i = 1; i < int(points.size()) - 1; i++)
  {
    int i0 = i - 1;
    int i1 = i;
    int i2 = i + 1;

    T d0 = !(points[i1] - points[i0]);
    T d1 = !(points[i2] - points[i1]);

    if (d0 < tolerance && d1 < tolerance)
      count++;
  }

  return count;
}

/** Remove triplicates. */
template <class T> int removeTriplicates(std::vector<TPoint<T>> &points, T tolerance = PARM_TOLERANCE)
{
  if (points.size() < 3)
    return 0;

  int count = 0;
  for (int i = int(points.size()) - 2; i >= 1; i--)
  {
    int i0 = i - 1;
    int i1 = i;
    int i2 = i + 1;

    T d0 = !(points[i1] - points[i0]);
    T d1 = !(points[i2] - points[i1]);

    if (d0 < tolerance && d1 < tolerance)
    {
      points.erase(points.begin() + i2);
      count++;
    }
  }

  return count;
}

/** i -> i + 1 are duplicate points. */
template <class T> bool isDuplicate(std::vector<TPoint<T>> &points, int i, T tolerance = PARM_TOLERANCE)
{
  int i1 = i + 1;

  if (i < 0)
    return false;

  if (i1 >= points.size())
    return false;

  T dist = !(points[i1] - points[i]);
  return (dist < tolerance);
}


/** Reverse two layers of points. */
template <class T> void reverse(std::vector<std::vector<TPoint<T>>> &points)
{
  for (int i = 0; i < int(points.size()); i++)
  {
    std::reverse(points[i].begin(),points[i].end());
  }
  std::reverse(points.begin(),points.end());
}

#if 1

/** Divide points into two by a divider. */
template <class T> bool divide(std::vector<TPoint<T>> &points, TPoint<T> &divider,
  std::vector<std::vector<TPoint<T>>> &newpoints, T &U, 
  T tolerance, T parmtolerance = PARM_TOLERANCE)
{
  int index = findClosest<T>(points,divider);

  // do not take front/backs
  if (index < 0 || index == 0 || index == points.size() - 1)
    return false;

  TPoint<T> proj; 
  int seg = -1;
  T Useg = 0.5;
  if (projectPointOnPoints(points,divider,proj,&seg,&Useg,parmtolerance))
  {
    // total length
    T totallen = calculateLength(points);

    // current length to projection point
    T seglen = !(proj - points[seg]);
    T len = calculateLength(points,seg);

    U = (len + seglen) / totallen;
    LIMIT(U,0.0,1.0);

    std::vector<TPoint<T>> part0(points.begin(),points.begin() + seg + 1);
    part0.push_back(proj);
    removeDuplicates(part0,false,tolerance);

    std::vector<TPoint<T>> part1(points.begin() + seg,points.end());
    part1.insert(part1.begin(),proj);
    removeDuplicates(part1,false,tolerance);

    newpoints.clear();
    newpoints.push_back(part0);
    newpoints.push_back(part1);

    return true;
  } else
  {
    return false;
  }
}

/** Divide points into two by a divider. Together with corresponding UV array. */
template <class T> bool divide(std::vector<TPoint<T>> &points, std::vector<TPoint<T>> &pointsUV, TPoint<T> &divider, 
  std::vector<std::vector<TPoint<T>>> &newpoints, std::vector<std::vector<TPoint<T>>> &newpointsUV,
  T tolerance, T parmtolerance = PARM_TOLERANCE)
{
  int index = findClosest<T>(points,divider);

  // do not take front/backs
  if (index < 0 || index == 0 || index == points.size() - 1)
    return false;

  if (points.size() != pointsUV.size())
    return false;

  TPoint<T> proj; 
  int seg = -1;
  T Useg = 0.5;
  if (projectPointOnPoints(points,divider,proj,&seg,&Useg,parmtolerance))
  {
    std::vector<TPoint<T>> part0(points.begin(),points.begin() + seg + 1);
    part0.push_back(proj);
    removeDuplicates(part0,false,tolerance);

    std::vector<TPoint<T>> part1(points.begin() + seg,points.end());
    part1.insert(part1.begin(),proj);
    removeDuplicates(part1,false,tolerance);

    newpoints.clear();
    newpoints.push_back(part0);
    newpoints.push_back(part1);

    TPoint<T> projUV = pointsUV[seg] + (pointsUV[seg + 1] - pointsUV[seg]) * Useg;

    std::vector<TPoint<T>> part0UV(pointsUV.begin(),pointsUV.begin() + seg + 1);
    part0UV.push_back(projUV);
    removeDuplicates(part0UV,false,tolerance);

    std::vector<TPoint<T>> part1UV(pointsUV.begin() + seg,pointsUV.end());
    part1UV.insert(part1UV.begin(),projUV);
    removeDuplicates(part1UV,false,tolerance);

    newpointsUV.clear();
    newpointsUV.push_back(part0UV);
    newpointsUV.push_back(part1UV);

    return true;
  } else
  {
    return false;
  }
}

#else

/** Divide points by into two by a divider. This division is INACCURATE!
  as the closest division points are sought. */
template <class T> bool divide(std::vector<TPoint<T>> &points, TPoint<T> &divider,
  std::vector<std::vector<TPoint<T>>> &newpoints)
{
  int index = findClosest<T>(points,divider);

  // do not take front/backs
  if (index < 0 || index == 0 || index == points.size() - 1)
    return false;

  newpoints.push_back(std::vector<TPoint<T>>(points.begin(),points.begin() + index + 1));
  newpoints.push_back(std::vector<TPoint<T>>(points.begin() + index,points.end()));

  return true;
}

/** Divide points by a divider. This division is INACCURATE!
  as the closest division points are sought. */
template <class T> void divide(std::vector<std::vector<TPoint<T>>> &points, TPoint<T> &divider)
{
  for (int i = int(points.size()) - 1; i >= 0; i--)
  {
    std::vector<std::vector<TPoint<T>>> newpoints;
    if (divide(points[i],divider,newpoints))
    {
      points.erase(points.begin() + i);
      points.insert(points.begin() + i,newpoints.begin(),newpoints.end());
    }
  }
}

#endif

/** Find longest piece. */
template <class T> int longestPiece(std::vector<std::vector<TPoint<T>>> &points, T *pmaxlen = nullptr)
{
  int longest = -1;

  T maxlen = 0.0;
  for (int i = 0; i < int(points.size()); i++)
  {
    T len = calculateLength(points[i]);
    if (len > maxlen)
    {
      maxlen = std::max(len,maxlen);
      longest = i;
    }
  }

  if (pmaxlen)
    *pmaxlen = maxlen;

  return longest;
}

/** Find longest piece longest by number of points. */
template <class T> int longestSizePiece(std::vector<std::vector<TPoint<T>>> &points)
{
  int longest = -1;

  int maxlen = 0;
  for (int i = 0; i < int(points.size()); i++)
  {
    int len = int(points[i].size());
    if (len > maxlen)
    {
      maxlen = std::max(len,maxlen);
      longest = i;
    }
  }

  return longest;
}

/** Get loop normal. Is it a hole or space around the hole? */
template <class T> TPoint<T> getLoopNormal(std::vector<std::vector<TPoint<T>>> &points, 
  T parmtolerance = PARM_TOLERANCE)
{
  std::vector<TPoint<T>> wholeloop;
  for (int i = 0; i < int(points.size()); i++)
  {
    wholeloop.insert(wholeloop.end(),points[i].begin(),points[i].end());
  }

  removeDuplicates<T>(wholeloop,false,parmtolerance);

  TPoint<T> loopnormal = calculateNormalXY(wholeloop);

  return loopnormal;
}

/** Find a pair of two closest values. n^2. */
template <class T> bool findTwoClosest(std::vector<TPoint<T>> &points,
  int &index0, int &index1)
{
  if (points.size() <= 2)
    return false;

  T mindist = std::numeric_limits<T>::max();
  for (int i = 0; i < int(points.size()); i++)
  {
    for (int j = i + 1; j < int(points.size()); j++)
    {
      T dist = !(points[i] - points[j]);
      if (dist < mindist)
      {
        mindist = dist;
        index0 = i;
        index1 = j;
      }
    }
  }

  return true;
}

/** Exclude pairs of closest numbers till the list reaches numpoints. Used in
  intersectTriangleByTriangle() to get exactly two intersection points. n^2. */
template <class T> void excludeClosestNumbers(std::vector<TPoint<T>> &points,
  T tolerance, int numpoints = 2)
{
  assert(numpoints >= 2);

  while (points.size() > numpoints)
  {
    int index0 = -1;
    int index1 = -1;
    if (!findTwoClosest(points,index0,index1))
      break;

    // we need to exlude a "less rounded" value
    int count0 = digitRoundCount(points[index0],tolerance);
    int count1 = digitRoundCount(points[index1],tolerance);

    if (count0 > count1)
    {
      points.erase(points.begin() + index1);
    } else
    {
      points.erase(points.begin() + index0);
    }
  }
}

/** Decimate points. */
template <class T> void decimatePoints(std::vector<TPoint<T>> &points, int numpoints)
{
  if (numpoints >= int(points.size()))
    return;

  std::vector<TPoint<T>> newpoints;

#if 1
  std::vector<int> ranges;
  stretchBresenhamPoints<T>(0,0,numpoints - 1,int(points.size()) - 1,ranges);

  for (int i = 0; i < int(ranges.size()); i++)
  {
    newpoints.push_back(points[ranges[i]]);
  }

  assert(!(points.front() - newpoints.front()) < TOLERANCE(T));
  assert(!(points.back() - newpoints.back()) < TOLERANCE(T));

#else
  std::vector<int> ranges;
  getRanges(int(points.size()) - 1,numpoints - 1,ranges);

  for (int i = 0; i < int(ranges.size()); i++)
  {
    newpoints.push_back(points[ranges[i]]);
  }

  assert(!(points.front() - newpoints.front()) < TOLERANCE(T));
  assert(!(points.back() - newpoints.back()) < TOLERANCE(T));
#endif

  points = newpoints;
}


}
