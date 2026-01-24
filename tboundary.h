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

  p0 = cornerUV[i];
  p1 = cornerUV[i1];
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
template <class T> int boundaryLine(std::vector<TPoint<T>> &boundary, T parmtolerance = PARM_TOLERANCE, bool *reversed = nullptr)
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
