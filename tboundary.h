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

  tboundary.h

  Surface boundary

*******************************************************************************/

#pragma once

#include "tmisc.h"
#include "tpoints.h"

namespace tcad {

/** Surface corner UV values. 
  parm - boundary parameter 0.0..4.0

      parm = 3.0        side 2      
            3------------------------------2 parm = 2.0
            |                              |
            |                              |
            |                              |
   side 3   |                              | side 1
            ^ V                            |
            |    U                         |
 parm = 4.0 0---->-------------------------1
       parm = 0.0       side 0        parm = 1.0

*/

template <class T> const std::array<TPoint<T>,4> cornerUV =
{
  TPoint<T>(0.0,0.0,0.0),
  TPoint<T>(1.0,0.0,0.0),
  TPoint<T>(1.0,1.0,0.0),
  TPoint<T>(0.0,1.0,0.0)
};

/** Next boundary corner or next side. */
template <class T> int nextCorner(int i)
{
  int i1 = i + 1;
  if (i1 > 3)
    i1 = 0;

  return i1;
}

/** Next boundary corner or next side. */
template <class T> int nextSide(int i)
{
  return nextCorner<T>(i);
}

/** Return correct boundary side direction. */
template <class T> TPoint<T> boundarySideDir(int i)
{
  assert (i >= 0 && i <= 3);

  int i1 = nextCorner<T>(i);

  TPoint<T> dir = cornerUV<T>[i1] - cornerUV<T>[i];
  return dir;
}

/** Is it a boundary point? */
template <class T> bool boundaryPoint(TPoint<T> p, T parmtolerance = PARM_TOLERANCE)
{
  return (
    std::abs(p.X - 0.0) < parmtolerance ||
    std::abs(p.X - 1.0) < parmtolerance ||
    std::abs(p.Y - 0.0) < parmtolerance ||
    std::abs(p.Y - 1.0) < parmtolerance);
}

/** All points are boundary? */
template <class T> bool boundaryPoints(std::vector<TPoint<T>> &points, T parmtolerance = PARM_TOLERANCE)
{
  for (auto &p : points)
  {
    if (!boundaryPoint<T>(p,parmtolerance))
      return false;
  }

  return true;
}

/** All points are boundary? */
template <class T> bool boundaryPoints(std::vector<std::vector<TPoint<T>>> &points, T parmtolerance = PARM_TOLERANCE)
{
  for (auto &p : points)
  {
    if (!boundaryPoints<T>(p,parmtolerance))
      return false;
  }

  return true;
}

/** Boundary piece number? */
template <class T> int boundarySide(TPoint<T> p, T parmtolerance = PARM_TOLERANCE)
{
  if (std::abs(p.X - 0.0) < parmtolerance)
  {
    // U == 0
    return 3;
  } else if (std::abs(p.X - 1.0) < parmtolerance)
  {
    // U == 1
    return 1;
  } else if (std::abs(p.Y - 0.0) < parmtolerance)
  {
    // V == 0
    return 0;
  } else if (std::abs(p.Y - 1.0) < parmtolerance)
  {
    // V == 1
    return 2;
  } else
  {
    return -1;
  }
}

/** Boundary side ends. */
template <class T> void boundarySideEnds(int side, TPoint<T> &p0, TPoint<T> &p1)
{
  int i = side;

  int i1 = nextCorner<T>(i);

  p0 = cornerUV<T>[i];
  p1 = cornerUV<T>[i1];
}

/** Is it a boundary line? To check bad intersection lines along boundary. */
template <class T> bool boundaryLine(std::vector<TPoint<T>> &boundary, T parmtolerance = PARM_TOLERANCE)
{
  for (auto &p : boundary)
  {
    if (!boundaryPoint(p))
      return false;
  }

  return true;
}

/** Piece of boundary line on side number. */
template <class T> int boundaryLine(std::vector<TPoint<T>> &boundary, T parmtolerance = PARM_TOLERANCE, 
  bool *reversed = nullptr)
{
  assert(boundary.size() > 1); //!!!
  if (boundary.size() <= 1)
    return -1;

  int index = -1;

  for (int i = 1; i < int(boundary.size() - 1); i++)
  {
    int in = boundarySide(boundary[i],parmtolerance);

    if (in >= 0)
    {
      if (index < 0)
      {
        index = in;
      } else
      {
        if (in != index)
          return -1;
      }
    } else
    {
      return -1;
    }
  }

  if (index >= 0 && reversed)
  {
    TPoint<T> mdir = midDirection<T>(boundary);
    TPoint<T> sdir = boundarySideDir<T>(index);

    *reversed = !(mdir > sdir);
  }

  return index;
}

/** If full boundary line returns side number else -1. */
template <class T> int fullBoundaryLine(std::vector<TPoint<T>> &boundary, T parmtolerance = PARM_TOLERANCE, 
  bool *reversed = nullptr)
{
  int index = boundaryLine(boundary,parmtolerance,reversed);

  if (index < 0)
    return index;

  TPoint<T> p0,p1;
  boundarySideEnds(index,p0,p1);

  if (!(boundary.front() - p0) < parmtolerance && !(boundary.back() - p1) < parmtolerance)
  {
    if (reversed)
      *reversed = false;
    return index;
  } else if (!(boundary.front() - p1) < parmtolerance && !(boundary.back() - p0) < parmtolerance)
  {
    if (reversed)
      *reversed = true;
    return index;
  } else
  {
    return -1;
  }
}

/** Get (outer) boundary lines numbers [0..3] which are NOT in the loop. */
template <class T> void getUsedBoundaryPieces(std::vector<std::vector<TPoint<T>>> &loop, 
  std::array<int,4> &used, T parmtolerance = PARM_TOLERANCE)
{
  used.fill(false);

  for (int i = 0; i < int(loop.size()); i++)
  {
    bool reversed = false;
    int index = boundaryLine(loop[i],parmtolerance,&reversed);
    if (index >= 0)
    {
      used[index] = true;
    }
  }
}

/** Set accurate values for boundary point. */
template <class T> void correctBoundaryPoint(TPoint<T> &p, T parmtolerance = PARM_TOLERANCE)
{
  if (std::abs(p.X - 0.0) < parmtolerance)
    p.X = 0.0;
  if (std::abs(p.X - 1.0) < parmtolerance)
    p.X = 1.0;
  if (std::abs(p.Y - 0.0) < parmtolerance)
    p.Y = 0.0;
  if (std::abs(p.Y - 1.0) < parmtolerance)
    p.Y = 1.0;
}

/** Set accurate values for boundary point. */
template <class T> void correctBoundaryPoints(std::vector<std::vector<TPoint<T>>> &loop, 
  T parmtolerance = PARM_TOLERANCE)
{
  // correct points on the boundary
  for (int i = 0; i < int(loop.size()); i++)
  {
    if (boundaryPoint(loop[i].front(),parmtolerance))
    {
      correctBoundaryPoint(loop[i].front(),parmtolerance);
    }
    if (boundaryPoint(loop[i].back(),parmtolerance))
    {
      correctBoundaryPoint(loop[i].back(),parmtolerance);
    }
  }
}

/** Boundary parameters are measured from node 0 (U = V = 0.0)
  counter-clockwise, max value being 4.0. Returns -1.0 in case of failure. */
template <class T> T UVToBoundaryParm(T U, T V, T parmtolerance = PARM_TOLERANCE)
{
  assert(boundaryPoint(TPoint<T>(U,V)));

  if (std::abs(V - 0.0) < parmtolerance)
  {
    return U;
  } else if (std::abs(U - 1.0) < parmtolerance)
  {
    return 1.0 + V;
  } else if (std::abs(V - 1.0) < parmtolerance)
  {
    return 2.0 + (1.0 - U);
  } else if (std::abs(U - 0.0) < parmtolerance)
  {
    return 3.0 + (1.0 - V);
  } else
  {
    assert(false);
    return -1.0;
  }
}

/** Convert boundary parameter [0.0..4.0] to U,V. */
template <class T> TPoint<T> boundaryParmToUV(T parm)
{
  T U = 0.0;
  T V = 0.0;

  if (parm >= 0.0 && parm < 1.0)
  {
    U = parm;
    V = 0.0;
  } else if (parm >= 1.0 && parm < 2.0)
  {
    U = 1.0;
    V = parm - 1.0;
  } else if (parm >= 2.0 && parm < 3.0)
  {
    U = 1.0 - (parm - 2.0);
    V = 1.0;
  } else if (parm >= 3.0 && parm < 4.0)
  {
    U = 0.0;
    V = 1.0 - (parm - 3.0);
  } else if (parm >= 4.0)
  {
    U = 0.0;
    V = 0.0;
  } else
  {
    assert(false);
  }
  assert(U >= 0 && U <= 1.0);
  assert(V >= 0 && V <= 1.0);
  LIMIT(U,0.0,1.0);
  LIMIT(V,0.0,1.0);

  assert(boundaryPoint(TPoint<T>(U,V)));

  return TPoint<T>(U,V);
}

/** Get parameter on the opposite side. parm being [0..4] */
template <class T> T oppositeParm(T parm)
{
  int i = (int) parm;
  T r = parm - T(i);
  T oparm = T(i) + 2.0;
  oparm += (1.0 - r);
  if (oparm >= 4.0) oparm -= 4.0;
  return oparm;
}

/** Get next parameter value when closing the intersection curve. */
template <class T> bool nextBoundaryParm(T parm, T endparm, T &nextparm,
  T parmtolerance = PARM_TOLERANCE)
{
  // ended?
  T dist = std::abs(endparm - parm);
  if (dist < parmtolerance)
    return false;

  int side = int(parm);
  T r = parm - T(side);

  if (std::abs(r - 0.0) < parmtolerance || std::abs(r - 1.0) < parmtolerance)
  {
    nextparm = parm + 1.0;
  } else
  {
    nextparm = parm + (1.0 - r);
  }

  if (endparm >= parm && endparm <= nextparm)
  {
    nextparm = endparm;
    return true;
  }

  if (nextparm >= 4.0)
    nextparm -= 4.0;

  if (endparm >= parm && endparm <= nextparm)
    nextparm = endparm;

  return true;
}


/** Go round from p0 (U,V in X,Y) to p1 counter-clockwise. Both p0 and p1 must
  be on the boundary. Points in the result list are corner points passed
  when going round from p0 to p1. p0 is first (right), p1 is second! */
template <class T> bool goRoundBoundary(TPoint<T> p0, TPoint<T> p1, std::vector<TPoint<T>> &points, 
  T parmtolerance = PARM_TOLERANCE)
{
  points.clear();

  assert(boundaryPoint(p0,parmtolerance));
  assert(boundaryPoint(p1,parmtolerance));

  if (!boundaryPoint(p0,parmtolerance) || !boundaryPoint(p1,parmtolerance))
    return false;

  // p is current point
  TPoint<T> p = p0;

  // list start point
  points.push_back(p);

  int side = boundarySide(p,parmtolerance);

  bool found = false;
  for (int i = 0; i < 5; i++)
  {
    // get side ends
    TPoint<T> v0,v1;
    boundarySideEnds(side,v0,v1);

    // to the right only
    v0 = p;

    // project ending point p1 on segment v0 -> v1 (a boundary rectangle side)
    TPoint<T> intr;
    T t = 0.5;
    if (projectPointOnSegment(p1,v0,v1,&intr,&t,parmtolerance))
    {
      T dist = !(p1 - intr);
      // all done, found
      if (dist < parmtolerance)
      {
        points.push_back(p1);
        found = true;
        break;
      }
    }

    side = nextSide(side);
    p = v1;
    points.push_back(p);
  }

  // duplicates
  removeDuplicates(points,false,parmtolerance);

  assert(found);

  return found && (points.size() >= 2);
}


/** Prepare points before intersection. Divide pieces by duplicate points. */
template <class T> void cutIntoPoints(std::vector<std::vector<TPoint<T>>> &cut, std::vector<TPoint<T>> &points)
{
  points.clear();

  int count = 0;
  for (int i = 0; i < int(cut.size()); i++)
  {
    for (int j = 0; j < int(cut[i].size()); j++) // points may be duplicate, they mark new part
    {
      TPoint<T> p = cut[i][j];
      p.W = T(count++);
      points.push_back(p);
    }
  }
}

/** Divide points into pieces by duplicates. */
template <class T> void divideByDuplicates(std::vector<TPoint<T>> &points,
  std::vector<std::vector<TPoint<T>>> &loop, T tolerance)
{
  int start = 0; 
  for (int i = 0; i < int(points.size()); i++)
  {
    TPoint<T> p = points[i];

    if (i < int(points.size()) - 1)
    {
      TPoint<T> p1 = points[i + 1];
      T d = !(p1 - p);

      if (d < tolerance)
      {
        std::vector<TPoint<T>> piece(points.begin() + start,points.begin() + i + 1);
        if (!piece.empty())
        {
          loop.push_back(piece);
        }

        start = i + 1;
      }
    } else
    {
      std::vector<TPoint<T>> piece(points.begin() + start,points.end());
      if (!piece.empty())
      {
        loop.push_back(piece);
      }
    }
  } 
}

/** Convert cut and loop into points with duplicate points to mark boundaries between peieces.
  Intersect loop by cut. Return number of intersections. */
template <class T> int intersectLoopByCut(
  std::vector<std::vector<TPoint<T>>> &loop, std::vector<std::vector<TPoint<T>>> &cut, 
  std::vector<TPoint<T>> &looppoints, std::vector<TPoint<T>> &cutpoints, 
  std::vector<TPoint<T>> &UV, T parmtolerance = PARM_TOLERANCE)
{
  UV.clear();

  // make a single continuous curve from outer loop with marking sharp corners :
  // they mark every START of new curve piece
  cutIntoPoints(loop,looppoints);

  // set W for cut
  cutIntoPoints(cut,cutpoints);

  // find intersections, UV contain segment number + fraction of intersections
  // cut curve specifies a correct direction of the loop, it must go first here
  int numintrs = findIntersections(cutpoints,looppoints,UV,parmtolerance); 

  // cleanup duplicate points
  if (UV.size() > 1)
  {
    std::vector<TPoint<T>> UVclean;

    for (int k = 0; k < int(UV.size()); k++)
    {
      bool found = false;
      for (int l = k + 1; l < int(UV.size()); l++)
      {
        T d0 = UV[l].X - UV[k].X;
        T d1 = UV[l].Y - UV[k].Y;
        T dist0 = std::abs(d0);
        T dist1 = std::abs(d1);

        found = (
          dist0 < parmtolerance / T(cutpoints.size()) ||
          dist1 < parmtolerance / T(looppoints.size())
        );

        if (found)
          break;
      }
      if (!found)
        UVclean.push_back(UV[k]);
    }

    UV = UVclean;

    numintrs = int(UV.size());
  }

  return numintrs;
}

/** Get segments and segment UVs from intersection parms for TWO points. */
template <class T> void getSegment2UVs(std::vector<TPoint<T>> &UV,
  int &Useg0, T &U0, int &Useg1, T &U1, int &Vseg0, T &V0, int &Vseg1, T &V1)
{
  assert(UV.size() == 2);

  if (UV[1].X < UV[0].X)
  {
    SWAP(TPoint<T>,UV[0],UV[1]);
  }

  // all points
  Useg0 = int(UV[0].X);
  U0 = UV[0].X - T(Useg0);
  Useg1 = int(UV[1].X);
  U1 = UV[1].X - T(Useg1);

  // now go all cut points
  Vseg0 = int(UV[0].Y);
  V0 = UV[0].Y - T(Vseg0);
  Vseg1 = int(UV[1].Y);
  V1 = UV[1].Y - T(Vseg1);
}

/** Cut out newpoints from points by going round the contour from newpoints end to 
  newpoints start. */
template <class T> bool cutOutGoingRound(std::vector<TPoint<T>> &points, int seg0, T U0,
  TPoint<T> startUV, TPoint<T> endUV, std::vector<TPoint<T>> &newpoints, int &n3, T parmtolerance = PARM_TOLERANCE)
{
  bool found = false;

  // starting point must be duplicate
  newpoints.push_back(startUV);

  // go from seg0 + 1 to the end, look for endUV inside intervals
  //!!! BUG for (int i = seg0 + 1; i < points.size() - 1; i++) 

  // check FIRST same segment, this is not quite right, as the small segment may be from
  // already trimmed outer boundary
  if (UVToBoundaryParm(endUV.X,endUV.Y,parmtolerance) > UVToBoundaryParm(startUV.X,startUV.Y,parmtolerance))
  {
    TPoint<T> p0 = points[seg0];
    TPoint<T> p1 = points[seg0 + 1];
    TPoint<T> dp = p1 - p0;
    T len = !dp;

    if (len > parmtolerance)
    {
      // find end point between p0 and p1
      TPoint<T> intr;
      T t = 0.5;
      if (projectPointOnSegment(endUV,p0,p1,&intr,&t,parmtolerance))
      {
        T dist = !(endUV - intr);

        if (dist < parmtolerance)
        {
          newpoints.push_back(endUV);
          return true;
        }
      }
    }
  }

  // go from seg0 to the end, look for endUV inside intervals
  for (int i = seg0 + 1; i < points.size() - 1; i++)
  {
    TPoint<T> p0 = points[i];
    TPoint<T> p1 = points[i + 1];
    TPoint<T> dp = p1 - p0;
    T len = !dp;

    if (len > parmtolerance)
    {
      // find end point between p0 and p1
      TPoint<T> intr;
      T t = 0.5;
      if (projectPointOnSegment(endUV,p0,p1,&intr,&t,parmtolerance))
      {
        T dist = !(endUV - intr);

        if (dist < parmtolerance)
        {
          // non-trivial bug
          if (isDuplicate(points,i - 1,parmtolerance))
            newpoints.push_back(points[i - 1]);

          newpoints.push_back(endUV);
          found = true;
          break;
        } else
        {
          newpoints.push_back(p0);
        }
      } else
      {
        newpoints.push_back(p0);
      }
    } else
    {
      newpoints.push_back(p0);
    }
  }

  // another part
  if (!found)
  {
    // add a duplicate to the end
    newpoints.push_back(points.back());
    newpoints.push_back(points.back());

    // test starting part of points from 0 to seg0
    for (int i = 0; i <= seg0; i++)
    {
      TPoint<T> p0 = points[i];
      TPoint<T> p1 = points[i + 1];
      TPoint<T> dp = p1 - p0;
      T len = !dp;

      if (len > parmtolerance)
      {
        // find end point between p0 and p1
        TPoint<T> intr;
        T t = 0.5;
        if (projectPointOnSegment(endUV,p0,p1,&intr,&t,parmtolerance))
        {
          T dist = !(endUV - intr);

          if (dist < parmtolerance)
          {
            // non-trivial bug
            if (isDuplicate(points,i - 1,parmtolerance))
              newpoints.push_back(points[i - 1]);

            newpoints.push_back(endUV);
            found = true;
            break;
          } else
          {
            newpoints.push_back(p0);
          }
        } else
        {
          newpoints.push_back(p0);
        }
      } else
      {
        newpoints.push_back(p0);
      }
    }
  }

  removeTriplicates(newpoints,parmtolerance);

  n3 = numTriplicates(newpoints);

#if _DEBUG
  int n2 = numDuplicates(newpoints);

  //assert(n3 == 0);
#endif

  return found;
}

/** Cut out newpoints from points by going round the contour from newpoints end to 
  newpoints start. */
template <class T> bool cutOutGoingRound(std::vector<TPoint<T>> &points, int seg0, T U0,
  std::vector<TPoint<T>> &newpoints, int &n3, T parmtolerance = PARM_TOLERANCE)
{
  TPoint<T> startUV = newpoints.back();
  TPoint<T> endUV = newpoints.front();

  return cutOutGoingRound(points,seg0,U0,startUV,endUV,newpoints,n3,parmtolerance);
}

/** Combine new loop points for 2 intersection points. They have duplicates bewteen pieces. */
template <class T> bool makeNewLoopPoints(std::vector<TPoint<T>> &cutpoints,
  std::vector<TPoint<T>> &looppoints, std::vector<TPoint<T>> &UV, 
  std::vector<TPoint<T>> &newpoints, bool insertcut, T parmtolerance = PARM_TOLERANCE)
{
  newpoints.clear();

  // first intersection point is from cutpoints start, cutpoints direction from
  // first to second intersection point defines cutting direction : void is to
  // the right, surface is to the left

  int Useg0,Useg1,Vseg0,Vseg1;
  T U0,U1,V0,V1;
  getSegment2UVs<T>(UV,Useg0,U0,Useg1,U1,Vseg0,V0,Vseg1,V1);

  // just take all cut points void is to the right, surface is to the left
  cutOut(cutpoints,Useg0,U0,Useg1,U1,newpoints);

  TPoint<T> startUV = newpoints.back();
  TPoint<T> endUV = newpoints.front();

  if (!insertcut)
    newpoints.clear();

  // we need to proceed from newpoints last point to newpoints first
  // point; we always leave surface to the left, so there are TWO CASES here : 
  // point 1 has "greater" UV position when going along the boundary or "less"
  int n3 = 0;
  bool ok = cutOutGoingRound<T>(looppoints,Vseg1,V1,startUV,endUV,newpoints,n3,parmtolerance);
  if (n3 > 0)
  {
    ok = false;
  }

  // it must be always ok
  assert(ok);

  return ok;
}

/** Find a cut piece (except busy) close to p. p is "tail", cut front is "head" to be
  connected to tail. */
template <class T> int findClosestFront(std::vector<std::vector<TPoint<T>>> &cut,
  std::vector<bool> &busy, TPoint<T> p, T tolerance)
{
  // find closest piece to the current loop end among not busy
  int closest = -1;

  T mindist = std::numeric_limits<T>::max();
  for (int i = 0; i < int(cut.size()); i++)
  {
    if (busy[i])
      continue;

    T dist = !(p - cut[i].front());
    if (dist < tolerance && dist < mindist)
    {
      mindist = dist;
      closest = i;
    }
  }

  return closest;
}

/** Find a cut piece (except busy) by one end close to p. */
template <class T> int findClosest(std::vector<std::vector<TPoint<T>>> &cut,
  std::vector<bool> &busy, TPoint<T> p, bool &front, T tolerance)
{
  // find closest piece to the current loop end among not busy
  int closest = -1;
  front = false;

  T mindist = std::numeric_limits<T>::max();
  for (int i = 0; i < int(cut.size()); i++)
  {
    if (busy[i])
      continue;

    T dist = !(p - cut[i].front());
    if (dist < tolerance && dist < mindist)
    {
      front = true;
      mindist = dist;
      closest = i;
    }

    dist = !(p - cut[i].back());
    if (dist < tolerance && dist < mindist)
    {
      front = false;
      mindist = dist;
      closest = i;
    }
  }

  return closest;
}

/** Boundary parameter from 0.0 to 4.0. */
template <class T> struct TBoundary {

  /** Parameter from 0.0 to 4.0, counterclockwise around normal. */
  T parm = 0.0;

  /** Default constructor. */
  TBoundary()
  {
  }

  /** Constructor. */
  TBoundary(T pparm)
  {
    parm = pparm;
  }

  /** Constructor from U,V (U,V must be on boundary). */
  TBoundary(T U, T V, T parmtolerance = PARM_TOLERANCE)
  {
    assert(boundaryPoint(TPoint<T>(U,V)));

    parm = UVToBoundaryParm(U,V,parmtolerance);
  }

  /** Copy constructor. */
  TBoundary(const TBoundary &copy)
  {
    parm = copy.parm;
  }

  /** Assignment. */
  TBoundary operator = (const TBoundary &copy)
  {
    parm = copy.parm;
    return *this;
  }

  /** Get U,V in return X,Y. */
  TPoint<T> UV()
  {
    return boundaryParmToUV(parm);
  }

  ///** Boundary direction is measured counter-clockwise passing face with
  //  surface to the left. */
  //operator > (TBoundary first) const
  //{
  //  TPoint<T> p0 = first.UV();
  //  TPoint<T> p1 = UV();
  //  TPoint<T> d = p1 - p0;
  //  TPoint<T> normal();
  //}
};

}
