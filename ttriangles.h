/**
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
#include "tmisc.h"
#include "tpoints.h"
#include "ttransform.h"
#include "tplane.h"
#include "tjacobipoly.h"
#include "tdata.h"
// used to speedup triangle/triangle intersections by spacial partitioning
#include "obackground.h"

#include "strings.h"

#include <assert.h>
#include <vector>
#include <array>
#include <set>
#include <unordered_set>
#include <map>
#include <string>
#include <algorithm>
#include <fstream>

#include <thread>
#include <mutex>

#define NOMINMAX
#include "windows.h"

// output all debug info
//#define DEBUG_TRIS

#ifdef NDEBUG
  #undef DEBUG_TRIS
#endif

// slower but more reliable
//!!!!!!! #define GLOBAL_INTRCHECK

namespace tcad {

/**
  Class TTriangles - a collection of triangles
  --------------------------------------------

  Simple templated class to keep a surface as a set of 3D triangles represented by node coordinates <I>coords</I>
and connectivity array <I>corners</I> (3 coordinate indices for every triangle (face)). So, the
three coordinates for an i-th triangle can be extracted as

  coord0 = coords[corners[i * 3]];
  coord1 = coords[corners[i * 3 + 1]];
  coord2 = coords[corners[i * 3 + 2]];

  There is a massive infrastructure to handle face topology like getCornerTris() (find all triangles around
corner), getEdgeTris() (find all triangles around an edge) etc. There are two functions to define
if a triangulation is manifold() or solid(). buildConnectivityArray() excludes node duplicates and renumbers
corners accordingly. Tolerance there is very important. intersectByPlane() cuts triangles and reconstructs
the intersection line. getBoundary() makes a boundary as an ordered closed set of points.
  
  Triangles can be loaded and saved from/into STL files. No checks about node numeration :
it is supposed that face normal is defined by counter-clockwise numeration of nodes when
looking from the normal sharp end (see faceNormal()).

  When loading, duplicate nodes are excluded with a tolerance and face corners renumbered.
The tolerance specified in loadSTL() is important, it must be large enough - it is
used in exclusion of duplicate coordinates. If tolerance too large, degenerated triangles 
may appear; if too small - duplicate nodes may not be all excluded.

  STL files in the binary form may contain multiple parts by simply glueing multiple STL files
together; there is a provision in the code for that; not well tested though. Binary STL
files contain floats in only 4-byte format.

  In text STL files the code selects proper number of digits in text representation for
templated float and double classes.

*/

/** To compare edges. */
bool edgeComp(const std::pair<LINT,LINT> &a, const std::pair<LINT,LINT> &b);

/** Edges equal? */
bool edgesEqual(const std::pair<LINT,LINT> &a, const std::pair<LINT,LINT> &b, bool *reversed = nullptr);

/** Comparison function for this edges. */
class EdgeCompare {
public:
  bool operator()(const std::pair<LINT,LINT> &a, const std::pair<LINT,LINT> &b) const
  {
    std::pair<LINT,LINT> aa = a;
    std::pair<LINT,LINT> bb = b;

    if (aa.second < aa.first)
    {
      LINT temp = aa.first; aa.first = aa.second; aa.second = temp;
    }

    if (bb.second < bb.first)
    {
      LINT temp = bb.first; bb.first = bb.second; bb.second = temp;
    }

    if (aa.first < bb.first)
    {
      return true;
    } else if (aa.first > bb.first)
    {
      return false;
    } else
    {
      if (aa.second < bb.second)
      {
        return true;
      } else if (aa.second > bb.second)
      {
        return false;
      } else
      {
        return false;
      }
    }
  }
};

/** Does edge contains this tri? */
bool edgeContainsTri(std::pair<LINT,LINT> edge, LINT tri, std::map<std::pair<LINT,LINT>,std::vector<LINT>,EdgeCompare> &edgeTris,
  LINT &othertri);

/** Find intersection by tri number among not busy edges. */
int findIntersectionByTri(LINT tri, std::vector<std::pair<LINT,LINT>> &iedges, std::vector<bool> &ibusy, 
  std::map<std::pair<LINT,LINT>,std::vector<LINT>,EdgeCompare> &edgeTris, LINT &othertri);

/** Find free hanging edge among not busy. */
bool findFreeEdge(std::map<std::pair<LINT,LINT>,std::vector<LINT>,EdgeCompare> &edgeTris,
  std::set<std::pair<LINT,LINT>,EdgeCompare> &busy, std::pair<LINT,LINT> &startedge);

/** Cut triangle by triangle. tri and other being corner coordinates. */
template <class T> int intersectTriangleByTriangle(std::array<TPoint<T>,3> &tri,
  std::array<TPoint<T>,3> &other, std::vector<TPoint<T>> &intrs, T tolerance, T parmtolerance = PARM_TOLERANCE)
{
  std::pair<TPoint<T>,TPoint<T>> edges[3];
  std::pair<TPoint<T>,TPoint<T>> oedges[3];
  edges[0] = std::pair<TPoint<T>,TPoint<T>>(tri[0],tri[1]);
  edges[1] = std::pair<TPoint<T>,TPoint<T>>(tri[1],tri[2]);
  edges[2] = std::pair<TPoint<T>,TPoint<T>>(tri[2],tri[0]);
  oedges[0] = std::pair<TPoint<T>,TPoint<T>>(other[0],other[1]);
  oedges[1] = std::pair<TPoint<T>,TPoint<T>>(other[1],other[2]);
  oedges[2] = std::pair<TPoint<T>,TPoint<T>>(other[2],other[0]);

//!!!!!!!
//bool stop = false;
//TPoint<T> err0(-29.799847932004926,0.14142135623730948,0.14142135623730764);
//TPoint<T> err1(-29.802149779470824,0.14142135623730953,0.14142135623730767);
//T d0 = !(tri[1] - err0);
//T d1 = !(tri[2] - err0);
//T d2 = !(tri[1] - err1);
//T d3 = !(tri[2] - err1);
//if (d0 < 0.0000001 || d1 < 0.0000001 || d2 < 0.0000001 || d3 < 0.0000001)
//{
//  int gsgsg = 0;
//  stop = true;
//}

  intrs.clear();

  // tri edges with other triangle
  for (int k = 0; k < 3; k++)
  {
    T U = 0.0;
    TPoint<T> intersection;
    // tolerance not used inside
    if (segTriIntersect(edges[k].first,edges[k].second,other,U,intersection,tolerance,parmtolerance))
    {
      intrs.push_back(intersection);
    }
  }

  // other edges with tri
  for (int k = 0; k < 3; k++)
  {
    T U = 0.0;
    TPoint<T> intersection;
    // tolerance not used inside
    if (segTriIntersect(oedges[k].first,oedges[k].second,tri,U,intersection,tolerance,parmtolerance))
    {
      intrs.push_back(intersection);
    }
  }

  removeDuplicates(intrs,true,tolerance); // true is correct here

  //!!! important : set correct orientation of intersection curve : part to the right is
  // thrown off
  if (intrs.size() == 2)
  {
    TPoint<T> normal = (tri[1] - tri[0]) ^ (tri[2] - tri[1]);
    TPoint<T> onormal = (other[1] - other[0]) ^ (other[2] - other[1]);
    TPoint<T> d = intrs[1] - intrs[0];
    TPoint<T> cross = d ^ normal;

    if (cross > onormal)
    {
      std::reverse(intrs.begin(),intrs.end());
    }
  }

  return int(intrs.size());
}

template <class T> class TTriangles {
public:
                              
  // node coordinates
  std::vector<TPoint<T>> coords;
                              
  // corner indices into vectors, 3 per face (triangle)
  std::vector<LINT> corners;
                              
  // replacement after exclusion of node duplicates
  std::vector<LINT> replacement;

  // auxiliary data on UV surface coordinates of every triangle
  // node if created from a surface 
  std::vector<std::array<TPoint<T>,3>> UVcorners;
      
  /** Constructor. */
  TTriangles() = default;

  /** Constructor. */
  TTriangles(const TTriangles &copy) 
  {
    coords = copy.coords;
    corners = copy.corners;
    replacement = copy.replacement;
  }

  /** Assignment. */
  TTriangles &operator=(const TTriangles &copy)
  {
    coords = copy.coords;
    corners = copy.corners;
    replacement = copy.replacement;

    return *this;
  }
                              
  /** Destructor. */
  ~TTriangles() = default;

  /** Number of triangles (faces) */
  inline LINT numFaces() const;

  /** Clear all, save memory. */
  void clear();

  /** Add flat triangle of three corners like that from STL file; call coords.reserve()
    before adding.
    NB to avoid checks (much faster), set tolerance to zero.
    angledegtolerance is set to throw off sliver tris. */
  bool addTri(TPoint<T> v0, TPoint<T> v1, TPoint<T> v2, T tolerance, T angletolerancedeg = T(0.0), 
    T sliverdeg = T(0.0), bool revert = false);

  /** Build connectivity array (with unique corners) by exclusion of
    duplicate coords and renumbering node indices (corners).
    a "false" return means that tolerance was too high but was fixed
    by removing tris with duplicate nodes. */
  bool buildConnectivityArray(const T tolerance, std::vector<LINT> *preplacement = nullptr, bool fixsliver = true);

  /** Node duplicates excluded? */
  bool duplicatesExcluded()
  {
    return (coords.size() < corners.size());
  }

  /** Get min/max of node coordinates. */
  std::pair<TPoint<T>,TPoint<T>> minmax();

  /** Get size. */
  T maxSize();

  /** Get triangles around every corner. */
  void getCornerTris(std::vector<std::vector<LINT>> &cornerTris);

  /** Find neighbour tris in specified number of layers 
    around tris already in triNumbers list (must not be empty). */
  void getTriNeighbours(const int numLayersAround, 
    const std::vector<std::vector<LINT>> &cornerTris, std::set<LINT> &triNumbers);

  /** Get edges and triangles on sides. */
  void getEdgeTris(std::map<std::pair<LINT,LINT>, std::vector<LINT>,EdgeCompare> &edgeTris);

  /** Get unique corners around a node */
  void getCornersAround(std::vector<std::set<LINT>> &cornersAround);

  /** Face pierced by segment. */
  bool segIntersect(const TPoint<T> &point0, const TPoint<T> &point1, const int tri,
    T &U, TPoint<T> &intersection, const T tolerance, const T parmtolerance = PARM_TOLERANCE);

  /** Get tri normal. */
  TPoint<T> faceNormal(LINT faceNo) const;

  /** Get tri centre. */
  TPoint<T> faceCentre(LINT faceNo, T *R = nullptr) const;

  /** Get tri area. */
  T faceArea(LINT faceNo) const;

  /** Area of all tris. */
  T area() const;
                               
  /** Centre of all tris. */
  TPoint<T> centre() const;

  /** Get tri size as max of sides. */
  T faceMaxSize(LINT faceNo) const;

  /** Get tri size as min of sides. */
  T faceMinSize(LINT faceNo) const;

  /** Get tri corners. */
  std::array<TPoint<T>,3> threeCorners(LINT faceNo) const;

  /** Get face min/max. */
  void faceMinMax(LINT faceNo, TPoint<T> &min, TPoint<T> &max) const;

  /** Make face box. */
  void faceBox(LINT faceNo, std::array<TPoint<T>,8> &box) const;

  /** Make all face boxes. */
  void makeFaceBoxes(std::vector<std::array<TPoint<T>,8>> &boxes) const;

  /** Apply Laplace smooth to nodes; coef = 0 - old nodes;
    0.5 - mean; 1.0 - pure Laplace (dangerous, may produce degenarate
    tris). */
  void smoothLaplace(const T coef);

  /** Load tris from STL file, exclude duplicate nodes with
    tolerance and renumber corners, throw off sliver tris. */
  bool loadSTL(const std::string filename, std::string &partname, bool &binary, T &tolerance, 
    T angledegtolerance = T(0.0), T sliverdegtolerance = T(0.0));

  /** Save as STL file. */
  bool saveSTL(const std::string filename, const std::string partname, bool binary) const;

  /** Save triangles into OBJ file. */
  bool saveOBJ(const std::string filename, const std::string partname) const;

  /** Save a face as STL file. */
  bool saveFaceSTL(int faceNo, const std::string filename, const std::string partname, bool binary) const;

  /** Save box as STL. */
  bool saveBoxSTL(std::array<TPoint<T>,8> &box, const std::string filename, 
    const std::string partname, bool binary) const;

  /** Generate points from triangle surfaces. pointsPerTri - 
    number of points per traingle (approximate). Points may be duplicate. */
  bool generatePoints(int pointsPerTri, bool excludeDuplicates, T tolerance, 
    std::vector<TPoint<T>> &points);

  /** Generate points from triangle surfaces. pointsPerTri - 
    number of points per traingle (approximate). Points may be duplicate. */
  bool generatePointsByRefine(T fraction, bool excludeDuplicates, T tolerance, 
    std::vector<TPoint<T>> &points);

  /** Get free boundary edges with triangle multiplicity 1. */
  void getBoundaryEdges(T tolerance, std::vector<std::pair<LINT,LINT>> &edges);

  /** Edges to unique nodes. */
  void edgesToNodes(std::vector<std::pair<LINT,LINT>> &edges, std::set<LINT> &nodes);

  /** Manifold? - no edges with more than 2 neighbour triangles, just 1 or 2? */
  bool manifold(T tolerance, std::vector<std::pair<LINT,LINT>> &badedges,
    std::vector<std::pair<LINT,LINT>> *boundaryedges = nullptr);

  /** Solid? - all edges with 2 neighbour triangles. */
  bool solid(T tolerance, std::vector<std::pair<LINT,LINT>> &badedges,
    std::vector<std::pair<LINT,LINT>> *boundaryedges = nullptr);

  /** Remove edges from coords. */
  int removeBadEdges(std::vector<std::pair<LINT,LINT>> &badedges,
    std::vector<std::pair<LINT,LINT>> &boundaryedges);

  /** Safe tolerance calculated om min edge size to be used on removeDuplicates().
    It MUST be called BEFORE removeDuplicates() and buildConnectivityArray(). */
  T minEdge();

  /** Max edge length. */
  T maxEdge();

  /** Make NACA0012 airfoil X[-0.5..0.5], Z[-span * 0.5, +span * 0.5]
    numX - #intervals along chord X, numZ - along span. */
  bool makeNACA0012(LINT numX, LINT numZ, T span, int fitdegree = 4);

  /** Make triangles from a (twisted) quad between 4 corners. */
  void makeQuad(TPoint<T> corners[4], LINT numU, LINT numV);

  /** Transform all triangles. */
  void makeTransform(TTransform<T> &transform);

  /** Get closest point on triangle. */
  enum {
    CLOSEST_EDGE01, // closest point is on edge 01
    CLOSEST_EDGE12, // ..12
    CLOSEST_EDGE20, // ..20
    CLOSEST_INSIDE  // inside triangle
  };

  TPoint<T> closestPoint(LINT faceNo, const TPoint<T> &point, int &closestpos) const;

  /** Point is inside triangle in XY plane. */
  bool pointInsideFaceXY(LINT faceNo, const TPoint<T> &point, T tolerance) const;

  /** Point is a corner or inside a triangle in XY plane. */
  bool pointInsideXY(const TPoint<T> &point, T tolerance) const;

  /** Add arrow to display velocity. */
  void addArrow(TPoint<T> p, TPoint<T> velocity);

  /** Get signed distance to face plane. */
  T signedDist(LINT faceNo, const TPoint<T> &point) const;

  /** Cut by plane into intersection line. boundary is UV boundary curve ofr surface. */
  bool intersectByPlane(TPlane<T> &plane, std::vector<std::vector<TPoint<T>>> &lines, T tolerance,
    T parmtolerance = PARM_TOLERANCE,
    std::vector<std::vector<TPoint<T>>> *boundary = nullptr);

  /** Intersect with other tris. boundary0,1 contain in U,V parameters in X,Y for 
    first and second surfaces for trimming. Set BOTH boundaries to null or not null 
    at the same time. */
  bool intersect(TTriangles<T> &other, std::vector<std::vector<TPoint<T>>> &lines, 
    T parmtolerance = TOLERANCE(T),
    std::vector<std::vector<TPoint<T>>> *boundary0 = nullptr,
    std::vector<std::vector<TPoint<T>>> *boundary1 = nullptr,
    int numthreads = NUM_THREADS); 

  /** Intersect with other tris. boundary0,1 contain in U,V parameters in X,Y for 
    first and second surfaces for trimming. Set BOTH boundaries to null or not null 
    at the same time. boxes are ready for every triangle. */
  bool intersect(TTriangles<T> &other, 
    std::vector<std::array<TPoint<T>,8>> &boxes,
    std::vector<std::array<TPoint<T>,8>> &oboxes,
    std::vector<std::vector<TPoint<T>>> &lines, T parmtolerance,
    std::vector<std::vector<TPoint<T>>> *boundary0,
    std::vector<std::vector<TPoint<T>>> *boundary1,
    int numthreads);

  /** Get pieces (normally one) of boundary as ordered boundary nodes. maxlinelength
    is max line length to prevent crashes (not clear). */
  bool getBoundaryNodes(std::vector<std::vector<LINT>> &ilines, T tolerance, LINT maxlinelength = 10000);

  /** Convert lines of indices into that of XYZ coords. */
  void indexLinesToNodes(std::vector<std::vector<LINT>> &ilines,
    std::vector<std::vector<TPoint<T>>> &lines);

  /** Get rough boundary points, duplicates NOT excluded. */
  bool getBoundaryPoints(std::vector<TPoint<T>> &boundary, T tolerance);

  /** Get edge min/max. */
  bool getEdgeMinMax(T &min, T &max);

  /** Get edge min/max for XYZ components. */
  bool getEdgeMinMax(TPoint<T> &dmin, TPoint<T> &dmax);

  /** Connect edges into a closed curve. */
  bool edgesIntoCurves(std::vector<std::pair<TPoint<T>,TPoint<T>>> &edges,
    std::vector<std::vector<TPoint<T>>> &line, T tolerance);

  /** Get boundary from free edges, maybe not a single close line. */
  bool getBoundary(std::vector<std::vector<TPoint<T>>> &boundary, T tolerance);

  /** Delete face, no checks. */
  void deleteFace(LINT faceNo);

  /** Delete triangles with duplicate nodes (good to call after buildConnectivity()). */
  bool checkSliver(bool fix = false);

  /** Centre of mass as average of face centres. */
  TPoint<T> getCentre() const;

  /** Get corner numbers. */
  inline void getFaceCorners(LINT faceNo, LINT &i0, LINT &i1, LINT &i2);

  /** Face edges. */
  inline void getFaceEdges(LINT faceNo, std::pair<LINT,LINT> edges[3]);

  /** Face edges as coordinates. */
  void getFaceEdges(LINT faceNo, std::pair<TPoint<T>,TPoint<T>> pedges[3]);

  /** Cut face. */
  int cutFaceByPlane(TPlane<T> &plane, LINT faceNo, std::vector<TPoint<T>> &intrs, T tolerance);

  /** Get all face centres with radius at W. R coef is R expansion coefficient. */
  void getCentresAndRadii(std::vector<TPoint<T>> &centres, T Rcoef = 2.1);

private:

  // Subdive a sliver tri, 0 is the sharpest node. Conformity of mesh is destroyed!
  void addSliver(TPoint<T> v0, TPoint<T> v1, TPoint<T> v2, LINT numdivs, bool revert);
                              // save as STL file
  bool saveSTL(FILE *fp, std::string partname, bool binary) const;
                              // save triangles into OBJ file
  bool saveOBJ(FILE *fp, const std::string partname) const;

  /** Find two boundary nodes around boundary node. */
  bool getBoundaryNodesAround(LINT node, std::vector<std::vector<LINT>> &cornerTris, 
    std::map<std::pair<LINT,LINT>,std::vector<LINT>,EdgeCompare> &edgeTris, 
    LINT &node0, LINT &node1);
};

template<class T> LINT TTriangles<T>::numFaces() const
{
  assert(corners.size() % 3 == 0);
  return (LINT)(corners.size()) / 3;
}

template<class T> void TTriangles<T>::clear()
{
  std::vector<TPoint<T>> ocoords;
  std::vector<LINT> ocorners;
  std::vector<std::array<TPoint<T>,3>> ouvcorners;

  coords.swap(ocoords);
  corners.swap(ocorners);
  UVcorners.swap(ouvcorners);

  replacement.clear();
}

template<class T> void TTriangles<T>::addSliver(TPoint<T> v0, TPoint<T> v1, TPoint<T> v2, LINT numdivs, bool revert)
{
  // 0 is the shapest node
 
  // two long sides
  TPoint<T> v01 = v1 - v0;
  TPoint<T> v02 = v2 - v0;

  // make starting tri
  addTri(v0,v0 + v01 / T(numdivs),v0 + v02 / T(numdivs),0,0,0,revert);

  for (LINT i = 1; i < numdivs; i++)
  {
    TPoint<T> c0 = v0 + v01 * (T(i) / T(numdivs));
    TPoint<T> c1 = v0 + v01 * (T(i + 1) / T(numdivs));
    TPoint<T> c2 = v0 + v02 * (T(i + 1) / T(numdivs));
    TPoint<T> c3 = v0 + v02 * (T(i) / T(numdivs));

    T d02 = !(c2 - c0);
    T d13 = !(c3 - c1);

    if (d02 < d13)
    {
      addTri(c0,c1,c2,0,0,0,revert);
      addTri(c0,c2,c3,0,0,0,revert);
    } else
    {
      addTri(c0,c1,c3,0,0,0,revert);
      addTri(c3,c1,c2,0,0,0,revert);
    }
  }
}

template<class T> bool TTriangles<T>::addTri(TPoint<T> v0, TPoint<T> v1, TPoint<T> v2, T tolerance,
  T angletolerancedeg, T sliverdeg, bool revert)
{
  TPoint<T> v01 = v1 - v0;
  TPoint<T> v12 = v2 - v1;
  TPoint<T> v20 = v0 - v2;
                           
  T len01 = !v01;
  T len12 = !v12;
  T len20 = !v20;

  TPoint<T> d01 = +v01;
  TPoint<T> d12 = +v12;
  TPoint<T> d20 = +v20;

  T angle0 = asin(!(d01 ^ d20)) * PCI;
  T angle1 = asin(!(d12 ^ d01)) * PCI;
  T angle2 = asin(!(d20 ^ d12)) * PCI;

  T anglemin = std::min<T>(angle0,std::min<T>(angle1,angle2));
                              // avoid checks
  if (tolerance > T(0.0) || angletolerancedeg > T(0.0))
  {
                              // check degeneration
    if (len01 <= tolerance)
      return false;
    if (len12 <= tolerance)
      return false;
    if (len20 <= tolerance)
      return false;

    T area = !(v01 ^ v12) * T(0.5);
    if (area < tolerance * tolerance)
      return false;

    T angle = (!(v01 ^ v12)) / ((!v01 * !v12)) * PCI;
    if (angle < angletolerancedeg)
      return false;
  }

  if (sliverdeg > T(0.0) && anglemin <= sliverdeg)
  {
    // subdivide sliver
    // 1-2 is min side
    if (len12 <= len01 && len12 <= len20)
    {
      addSliver(v0,v1,v2,10,revert);
      return true;
    // 2-0 is min side
    } else if (len20 <= len01 && len20 <= len12)
    {
      addSliver(v1,v2,v0,10,revert);
      return true;
    } else // 0-1 is min side
    {
      addSliver(v2,v0,v1,10,revert);
      return true;
    }
  }

                              // connectivity array (with coord duplicates)
  if (revert)
  {
    corners.push_back(coords.size());
    corners.push_back(coords.size() + 2);
    corners.push_back(coords.size() + 1);
  } else
  {
    corners.push_back(coords.size());
    corners.push_back(coords.size() + 1);
    corners.push_back(coords.size() + 2);
  }
                            // coordinates contain all three nodes
  coords.push_back(v0);
  coords.push_back(v1);
  coords.push_back(v2);

  assert(corners.size() % 3 == 0);

  return true;
}

template<class T> bool TTriangles<T>::saveSTL(FILE *fp, std::string partname, bool binary) const
{
  if (binary)
  {
    size_t error = 0;

    char header[80] = {0};
                            // remove "solid" if any
    std::string name = upCase(partname);
    if (name.substr(0,5) == "SOLID")
    {
                            // spoil first character by replacing by space
      partname[0] = ' ';
    }
                            // limit name length by 80 characters
    if (partname.length() > 80)
      partname.erase(80,std::string::npos);
                            // move string to header
    assert(partname.length() <= 80);
    memmove(header,partname.c_str(),partname.length());
                            // save header
    error = fwrite(header,sizeof(header),1,fp);
    if (error == 0) return false;
                            // save number of triangles
    int numtriangles = static_cast<int>(numFaces());
    assert(sizeof(numtriangles) == 4);
    error = fwrite(&numtriangles,sizeof(numtriangles),1,fp);
    if (error == 0) return false;
                            // save triangles (vectors must be 4-byte floats)
    float v[3] = {0};
    for (int i = 0; i < numFaces(); i++)
    {
      int i3 = i * 3;
                            // save normal
      TPoint<T> normal = faceNormal(i);
      v[0] = static_cast<float>(normal.X);
      v[1] = static_cast<float>(normal.Y);
      v[2] = static_cast<float>(normal.Z);
      error = fwrite(v,sizeof(v),1,fp);
      if (error == 0) return false;
                            // save three triangle nodes
      TPoint<T> co[3];
      co[0] = coords[corners[i3]];
      co[1] = coords[corners[i3 + 1]];
      co[2] = coords[corners[i3 + 2]];

      for (int j = 0; j < 3; j++)
      {
        v[0] = static_cast<float>(co[j].X);
        v[1] = static_cast<float>(co[j].Y);
        v[2] = static_cast<float>(co[j].Z);
        error = fwrite(v,sizeof(v),1,fp);
        if (error == 0) return false;
      }
                            // save attribute
      short int attr = 0;
      assert(sizeof(attr) == 2);
      error = fwrite(&attr,sizeof(attr),1,fp);
      if (error == 0) return false;
    }

    return true;
  } else
  {
    int error = 0;

    error = fprintf(fp,"solid %s\n",partname.c_str());
    if (error < 0) return false;

    //assert(coords.size() % 3 == 0);

                            // format string depends only on current float type
    int numdigits = std::numeric_limits<T>::digits10;
    std::string digitsstr = std::to_string(numdigits);
    std::string formatstr = std::string("vertex ") + std::string("%.") + digitsstr + "e %." + digitsstr + "e %." + digitsstr + "e\n";
    std::string normalstr = std::string("facet normal ") + std::string("%.") + digitsstr + "e %." + digitsstr + "e %." + digitsstr + "e\n";

    for (int i = 0; i < numFaces(); i++)
    {
      int i3 = i * 3;

      TPoint<T> normal = faceNormal(i);

      error = fprintf(fp,normalstr.c_str(),normal.X,normal.Y,normal.Z);
      if (error < 0) return false;
      error = fprintf(fp,"outer loop\n");
      if (error < 0) return false;
                            // save three triangle nodes
      TPoint<T> co[3];
      co[0] = coords[corners[i3]];
      co[1] = coords[corners[i3 + 1]];
      co[2] = coords[corners[i3 + 2]];

      error = fprintf(fp,formatstr.c_str(),co[0].X,co[0].Y,co[0].Z);
      if (error < 0) return false;
      error = fprintf(fp,formatstr.c_str(),co[1].X,co[1].Y,co[1].Z);
      if (error < 0) return false;
      error = fprintf(fp,formatstr.c_str(),co[2].X,co[2].Y,co[2].Z);
      if (error < 0) return false;

      error = fprintf(fp,"endloop\n");
      if (error < 0) return false;
      error = fprintf(fp,"endfacet\n");
      if (error < 0) return false;
    }

    error = fprintf(fp,"endsolid %s\n",partname.c_str());
    if (error < 0) return false;

    return true;
  }
}

template<class T> bool TTriangles<T>::saveSTL(const std::string filename, const std::string partname, 
  bool binary) const
{
  FILE *fp = nullptr;
  if (fopen_s(&fp,filename.c_str(),"wb") != 0)
    return false;

  bool res = saveSTL(fp,partname,binary);

  fclose(fp);

  return res;
}

template<class T> bool TTriangles<T>::saveFaceSTL(int faceNo, const std::string filename, const std::string partname, 
  bool binary) const
{
  TTriangles<T> tris;
  
  std::array<TPoint<T>,3> corners = threeCorners(faceNo);
  tris.addTri(corners[0],corners[1],corners[2],0.0);

  return tris.saveSTL(filename,partname,binary);
}

template<class T> bool TTriangles<T>::saveBoxSTL(std::array<TPoint<T>,8> &box, const std::string filename, 
  const std::string partname, bool binary) const
{
  TTriangles<T> tris;
  
  tris.addTri(box[0],box[1],box[4],0.0);
  tris.addTri(box[4],box[1],box[5],0.0);

  tris.addTri(box[1],box[2],box[5],0.0);
  tris.addTri(box[5],box[2],box[6],0.0);

  tris.addTri(box[2],box[3],box[6],0.0);
  tris.addTri(box[6],box[3],box[7],0.0);

  tris.addTri(box[3],box[0],box[4],0.0);
  tris.addTri(box[4],box[7],box[3],0.0);

  tris.addTri(box[1],box[0],box[3],0.0);
  tris.addTri(box[3],box[2],box[1],0.0);

  tris.addTri(box[4],box[5],box[6],0.0);
  tris.addTri(box[6],box[7],box[4],0.0);

  return tris.saveSTL(filename,partname,binary);
}

template<class T> bool TTriangles<T>::saveOBJ(FILE *fp, const std::string partname) const
{
                            // we save only corners
  assert(coords.size() > 0);

  int error = 0;
                            // header
  error = fprintf(fp,"# %s\n",partname.c_str());
  if (error < 0) return false;
                            // save coordinates
  error = fprintf(fp,"\n# coordinates\n\n");
  if (error < 0) return false;
                            // format string depends only on current float type
  int numdigits = std::numeric_limits<T>::digits10;
  std::string digitsstr = std::to_string(numdigits);
  std::string formatstr = std::string("v ") + std::string("%.") + digitsstr + "e %." + digitsstr + "e %." + digitsstr + "e\n";

  for (int i = 0; i < coords.size(); i++)
  {
    error = fprintf(fp,formatstr.c_str(),coords[i].X,coords[i].Y,coords[i].Z);
    if (error < 0) return false;
  }
                            // save coordinates
  error = fprintf(fp,"\n# faces\n\n");
  if (error < 0) return false;

  for (int i = 0; i < numFaces(); i++)
  {
    int i3 = i * 3;
    int i0 = corners[i3];
    int i1 = corners[i3 + 1];
    int i2 = corners[i3 + 2];

    error = fprintf(fp,"f %zd %zd %zd\n",i0 + 1,i1 + 1,i2 + 1);
    if (error < 0) return false;
  }

  return true;
}

template<class T> bool TTriangles<T>::saveOBJ(const std::string filename, const std::string partname) const
{
  FILE *fp = nullptr;
  if (fopen_s(&fp,filename.c_str(),"wb") != 0)
    return false;

  bool res = saveOBJ(fp,partname);

  fclose(fp);

  return res;
}

template<class T> bool TTriangles<T>::loadSTL(const std::string filename, std::string &partname, 
  bool &binary, T &tolerance, T angledegtolerance, T sliverdegtolerance)
{
  FILE *fp = nullptr;
  if (fopen_s(&fp,filename.c_str(),"rb") != 0)
    return false;
                            // result
  bool OK = false;
                            // read first bytes to find out if the file is binary or text
  binary = false;

  char header[200];
  size_t bytesread = fread(header,1,sizeof(header),fp);
  OK = (bytesread > 0);
  if (!OK) 
  {
    fclose(fp);
    return false;
  }
                            // text file will contain "solid" and a number of LFs in 200 bytes
                            // truncate header and use std::string to find "solid"
  header[199] = 0;
  auto sheader = std::string(header);
                          // "solid" starting from zero byte
  bool solidfound = (sheader.find("solid",0,5) == 0) ||
    (sheader.find("SOLID",0,5) == 0);
                          // count number of LFs
  int LFcount = 0;
  for (int i = 0; i < bytesread; i++)
  {
    if (header[i] == 10) LFcount++;
  }
                          // 100% the file is ASCII
  binary = !(solidfound && LFcount > 0);
                          // file is binary just seek its beginning
  if (binary)
  {
    partname = std::string(header);
                          // go to start of file
    OK = (fseek(fp,0,SEEK_SET) == 0);
                          // read non-standard file as many files glued together
    while (!feof(fp))
    {
                          // all int-s are 32-bit in STL file
      assert(sizeof(int) == 4);
                          // skip header
      OK = OK && (fseek(fp,80,SEEK_CUR) == 0);
                          // read number of triangles
      int numtriangles = 0;
      OK = OK && (fread(&numtriangles,1,sizeof(numtriangles),fp) == sizeof(numtriangles));
      if (!OK) 
      {
        fclose(fp);

        if (tolerance == 0.0)
        {
          std::pair<TPoint<T>,TPoint<T>> mm = minmax();
          TPoint<T> d = mm.second - mm.first;
          tolerance = !d * 0.000001;
        }
                            // OK, if some triangles already loaded
        bool ok = buildConnectivityArray(tolerance);
        return ok;
      }
      OK = (numtriangles > 0);
      if (!OK) 
      {
        fclose(fp);

        if (tolerance == 0.0)
        {
          std::pair<TPoint<T>,TPoint<T>> mm = minmax();
          TPoint<T> d = mm.second - mm.first;
          tolerance = !d * 0.000001;
        }
                            // OK, if some triangles already loaded
        bool ok = buildConnectivityArray(tolerance);
        return ok;
      }
                          // read all triangles
      assert(sizeof(float) == 4);
      float xyz[3];

      for (int i = 0; i < numtriangles; i++)
      {
        //TPoint<T> normal;

        OK = (fread(xyz,1,sizeof(xyz),fp) == sizeof(xyz));
        if (!OK) 
        {
          fclose(fp);
          return false;
        }

        //normal.X = static_cast<T>(xyz[0]);
        //normal.Y = static_cast<T>(xyz[1]);
        //normal.Z = static_cast<T>(xyz[2]);
        //normal = +normal;
                            // three vertices
        TPoint<T> v[3];

        for (int j = 0; j < 3; j++)
        {
          OK = (fread(xyz,1,sizeof(xyz),fp) == sizeof(xyz));
          if (!OK) 
          {
            fclose(fp);
            return false;
          }

          v[j].X = static_cast<T>(xyz[0]);
          v[j].Y = static_cast<T>(xyz[1]);
          v[j].Z = static_cast<T>(xyz[2]);
        }

        addTri(v[0],v[1],v[2],tolerance,angledegtolerance,sliverdegtolerance);
                            // 2 bytes of unused data
        short int temp = 0;
        OK = (fread(&temp,1,2,fp) == 2);
        if (!OK) 
        {
          fclose(fp);
          return false;
        }
      }
    }

    fclose(fp);

    if (tolerance == 0.0)
    {
      std::pair<TPoint<T>,TPoint<T>> mm = minmax();
      TPoint<T> d = mm.second - mm.first;
      tolerance = !d * 0.000001;
    }

    // this thing spoils connectivity
    if (sliverdegtolerance > T(0.0))
    {
      return false;
    } else
    {
                            // remove duplicate nodes and renumber corners
      bool ok = buildConnectivityArray(tolerance);
      return ok;
    }
  } else
                            // ASCII
  {
                            // close binary file and open as text
    fclose(fp);

    std::ifstream file;
    std::ios_base::openmode oMode(std::ios::in|std::ios::binary);
    file.open(filename.c_str(), oMode);

    if (!file.is_open())
    {
      return false;
    }
                            // line by line...

                            // loop theoretically may contain more than 3 nodes but
                            // it does not seem possible as binary version does not allow
                            // it at all
    int nodecount = 0;
    TPoint<T> normal;
    TPoint<T> v[3];

    while (!file.eof())
    {
      std::string line;
      std::getline(file,line);
                            // last line might well be without ending (CR)/LF
      if (file.fail()) {
        break;
      }
                            // trim line
      line = trim(line," \n\r\t");
                            // skip empty
      if (line.length() == 0)
        continue;
                            // original line
      std::string original = line;
                            // uppercase
      line = upCase(line);
                            // lines starting with SOLID
      if (line.substr(0,5) == "SOLID")
      {
        if (original.length() > 6)
        {
                            // store part name
          partname = original.substr(6);
          partname = trim(partname," \n\r\t");
        }
        continue;
      }
                            // lines starting with ENDSOLID
      if (line.substr(0,8) == "ENDSOLID")
      {
        continue;
      }

      if (line == "ENDFACET")
      {
      } else if (line.substr(0,5) == "FACET")
      {
                            // normal for this facet, this line is expected to be like
                            // "facet normal -0.7074 0.0 0.7074"
        int pos1[100],pos2[100];
        int numwords = parseWords(line,' ',pos1,pos2,100);
                            // skip this line
        if (numwords != 5)
        {
          normal = TPoint<T>();
                            // error, issue warning, do not return
          continue;
        }
                            // do not parse normal

      } else if (line == "OUTER LOOP")
      {
        nodecount = 0;

      } else if (line == "ENDLOOP")
      {
        if (nodecount >= 3)
        {
                            // add triangle
          addTri(v[0],v[1],v[2],tolerance,angledegtolerance,sliverdegtolerance);
        } else
        {
                            // Wrong count of nodes
                            // error, issue warning, do not return
        }

                            // "vertex -1.0 0.5 0.33"
      } else if (line.substr(0,6) == "VERTEX")
      {
        int pos1[100],pos2[100];
        int numwords = parseWords(line,' ',pos1,pos2,100);
                            // "vertex" line contains less than 3 coordinates, take only
                            // first available coordinates, others left zeroes
        if (nodecount < 3)
        {
          if (numwords < 4)
          {
            if (numwords == 2)
            {
              v[nodecount].X = static_cast<T>(atof(line.substr(pos1[1],(pos2[1] - pos1[1] + 1)).c_str()));
              v[nodecount].Y = T(0.0);
              v[nodecount].Z = T(0.0);
            } if (numwords == 3)
            {
              v[nodecount].X = static_cast<T>(atof(line.substr(pos1[1],(pos2[1] - pos1[1] + 1)).c_str()));
              v[nodecount].Y = static_cast<T>(atof(line.substr(pos1[2],(pos2[2] - pos1[2] + 1)).c_str()));
              v[nodecount].Z = T(0.0);
            }
                              // error, issue warning, do not return
          } else
          {
            v[nodecount].X = static_cast<T>(atof(line.substr(pos1[1],(pos2[1] - pos1[1] + 1)).c_str()));
            v[nodecount].Y = static_cast<T>(atof(line.substr(pos1[2],(pos2[2] - pos1[2] + 1)).c_str()));
            v[nodecount].Z = static_cast<T>(atof(line.substr(pos1[3],(pos2[3] - pos1[3] + 1)).c_str()));
          }

          nodecount++;
        }
      } else

                            // unidentified line
      {
                            // error, issue warning, do not return
      }
    } // while loop on all lines

    file.close();

    if (tolerance == 0.0)
    {
      std::pair<TPoint<T>,TPoint<T>> mm = minmax();
      TPoint<T> d = mm.second - mm.first;
      tolerance = !d * PARM_TOLERANCE;
    }

    // this thing spoils connectivity
    if (sliverdegtolerance > T(0.0))
    {
      return false;
    } else
    {

                            // remove duplicate nodes and renumber corners
      bool ok = buildConnectivityArray(tolerance);
      return ok;
    }
  } // ASCII

  return false;
}

template <class T> bool TTriangles<T>::buildConnectivityArray(const T tolerance, std::vector<LINT> *preplacement, bool fixsliver)
{
  if (numFaces() == 0)
    return false;

  assert(coords.size() > 0);
  assert(corners.size() > 0);

#ifdef DEBUG_TRIS
  outputDebugString(std::string("coords size before ") + std::to_string(int(coords.size())));
#endif
                            // exclude duplicate nodes
  if (!removeDuplicates(coords,true,tolerance,&replacement))    
    return false;

#ifdef DEBUG_TRIS
  outputDebugString(std::string("coords size after ") + std::to_string(int(coords.size())));
#endif

                            // renumber all indices
  for (int i = 0; i < corners.size(); i++)
  {
    corners[i] = replacement[corners[i]];
  }

  //// test for correct node numbering : coincident nodes indicate that
  //// tolerance is TOO high
  //for (int i = 0; i < corners.size(); i += 3)
  //{
  //  if (corners[i] == corners[i + 1] || corners[i] == corners[i + 2] || corners[i + 1] == corners[i + 2])
  //    return false;
  //}

  // check and fix instead
  bool check = checkSliver(fixsliver);

#ifdef DEBUG_TRIS
  outputDebugString("NODES");
  for (int i = 0; i < coords.size(); i++)
  {
    outputDebugString(
      std::to_string(i) + std::string(" ") +
      std::to_string(coords[i].X) + std::string(" ") +
      std::to_string(coords[i].Y) + std::string(" ") +
      std::to_string(coords[i].Z)
    );
  }
  outputDebugString("CORNERS");
  for (int i = 0; i < corners.size(); i += 3)
  {
    outputDebugString(
      std::string("face ") + std::to_string(i / 3) + std::string(" corner ") + std::to_string(i) + std::string(" ") +
      std::to_string(int(corners[i])) + std::string(" ") +
      std::to_string(int(corners[i + 1])) + std::string(" ") +
      std::to_string(int(corners[i + 2]))
    );
  }
#endif

  if (preplacement)
    *preplacement = replacement;

  // not necessarily bad
  return check;
}

template <class T> std::pair<TPoint<T>,TPoint<T>> TTriangles<T>::minmax()
{
  assert(coords.size() > 0);
  TPoint<T> min = coords[0];
  TPoint<T> max = coords[0];

  for (int i = 1; i < coords.size(); i++)
  {
    for (int j = 0; j < 3; j++)
    {
      min.XYZ[j] = std::min(min.XYZ[j],coords[i].XYZ[j]);
      max.XYZ[j] = std::max(max.XYZ[j],coords[i].XYZ[j]);
    }
  }

  return std::pair<TPoint<T>,TPoint<T>>(min,max);
}

template <class T> T TTriangles<T>::maxSize()
{
  std::pair<TPoint<T>,TPoint<T>> mm = minmax();
  TPoint<T> d = mm.second - mm.first;
  return !d;
}

template <class T> void TTriangles<T>::getCornerTris(std::vector<std::vector<LINT>> &cornerTris)
{
  assert(coords.size() > 0);
  cornerTris.clear();
  cornerTris.resize(coords.size(),std::vector<LINT>());

  for (int i = 0; i < corners.size(); i += 3)
  {
    int trinumber = i / 3;

    for (int j = 0; j < 3; j++)
    {
      cornerTris[corners[i + j]].push_back(trinumber);
    }
  }
}

template <class T> void TTriangles<T>::getCornersAround(std::vector<std::set<LINT>> &cornersAround)
{
  assert(coords.size() > 0);
  cornersAround.clear();
  cornersAround.resize(coords.size(),std::set<LINT>());

  for (int i = 0; i < numFaces(); i++)
  {
    LINT i0,i1,i2;
    getFaceCorners(i,i0,i1,i2);

    cornersAround[i0].insert(i1);
    cornersAround[i0].insert(i2);
    cornersAround[i1].insert(i0);
    cornersAround[i1].insert(i2);
    cornersAround[i2].insert(i0);
    cornersAround[i2].insert(i1);
  }
}

template <class T> void TTriangles<T>::getEdgeTris(std::map<std::pair<LINT,LINT>,
  std::vector<LINT>,EdgeCompare> &edgeTris)
{
  assert(coords.size() > 0);
  edgeTris.clear();

  for (int i = 0; i < corners.size(); i += 3)
  {
    LINT trinumber = i / 3;

    std::pair<LINT,LINT> e[3] = {
      std::pair<LINT,LINT>(corners[i],corners[i + 1]),
      std::pair<LINT,LINT>(corners[i + 1],corners[i + 2]),
      std::pair<LINT,LINT>(corners[i + 2],corners[i])
    };

    for (int j = 0; j < 3; j++)
    {
      // swap numbers
      if (e[j].second < e[j].first)
      {
        LINT temp = e[j].second; e[j].second = e[j].first; e[j].first = temp;
      }

      auto it = edgeTris.find(e[j]);
      if (it == edgeTris.end())
      {
        edgeTris.insert(std::pair<std::pair<LINT,LINT>,std::vector<LINT>>(e[j],
          std::vector<LINT>(std::initializer_list<LINT>{trinumber})));
      } else
      {
        it->second.push_back(trinumber);
      }
    }
  }
}

template <class T> void TTriangles<T>::smoothLaplace(const T coef)
{
  std::vector<std::vector<LINT>> cornerTris;
  getCornerTris(cornerTris);

  std::vector<TPoint<T>> newCoords(coords.size());

  for (int i = 0; i < coords.size(); i++)
  {
    std::set<LINT> neiNodes;
    for (auto tri : cornerTris[i])
    {
      neiNodes.insert(corners[tri * 3]);
      neiNodes.insert(corners[tri * 3 + 1]);
      neiNodes.insert(corners[tri * 3 + 2]);
    }

    TPoint<T> sum(0,0,0);
    int count = 0;
    for (auto n : neiNodes)
    {
      if (n != i)
      {
        sum += coords[n];
        count++;
      }
    }

    newCoords[i] = coords[i] + ((sum /= static_cast<T>(count)) - coords[i]) * coef;

#if 0
#ifdef _DEBUG
    T diff = !(coords[i] - newCoords[i]);
    cout << diff << endl;
#endif
#endif
  }

  coords = newCoords;
}

template <class T> TPoint<T> TTriangles<T>::faceNormal(LINT faceNo) const
{
  LINT i0 = corners[faceNo * 3];
  LINT i1 = corners[faceNo * 3 + 1];
  LINT i2 = corners[faceNo * 3 + 2];
  TPoint<T> n = +((coords[i1] - coords[i0]) ^ (coords[i2] - coords[i1]));
  return n;
}

template <class T> TPoint<T> TTriangles<T>::faceCentre(LINT faceNo, T *R) const
{
  LINT i0 = corners[faceNo * 3];
  LINT i1 = corners[faceNo * 3 + 1];
  LINT i2 = corners[faceNo * 3 + 2];
  TPoint<T> c = (coords[i0] + coords[i1] + coords[i2]) * 0.333333333;
  if (R)
    *R = !(coords[i0] - c);
  return c;
}

template <class T> T TTriangles<T>::faceMaxSize(LINT faceNo) const
{
  LINT i0 = corners[faceNo * 3];
  LINT i1 = corners[faceNo * 3 + 1];
  LINT i2 = corners[faceNo * 3 + 2];

  T size = std::max<T>(!(coords[i0] - coords[i2]),std::max<T>(!(coords[i0] - coords[i1]),!(coords[i1] - coords[i2])));
  return size;
}

template <class T> T TTriangles<T>::faceMinSize(LINT faceNo) const
{
  LINT i0 = corners[faceNo * 3];
  LINT i1 = corners[faceNo * 3 + 1];
  LINT i2 = corners[faceNo * 3 + 2];

  T size = std::min<T>(!(coords[i0] - coords[i2]),std::min<T>(!(coords[i0] - coords[i1]),!(coords[i1] - coords[i2])));
  return size;
}

template <class T> std::array<TPoint<T>,3> TTriangles<T>::threeCorners(LINT faceNo) const
{
  LINT i0 = corners[faceNo * 3];
  LINT i1 = corners[faceNo * 3 + 1];
  LINT i2 = corners[faceNo * 3 + 2];
  return std::array<TPoint<T>,3>({coords[i0],coords[i1],coords[i2]});
}

template <class T> void TTriangles<T>::faceMinMax(LINT faceNo, TPoint<T> &min, TPoint<T> &max) const
{
  std::array<TPoint<T>,3> corners = threeCorners(faceNo);

  min = pointMin<T>(corners[0],pointMin<T>(corners[1],corners[2]));
  max = pointMax<T>(corners[0],pointMax<T>(corners[1],corners[2]));
}

template <class T> void TTriangles<T>::faceBox(LINT faceNo, std::array<TPoint<T>,8> &box) const
{
  TPoint<T> min,max;
  faceMinMax(faceNo,min,max);

  makeBox<T>(min,max,box);
}

template <class T> void TTriangles<T>::makeFaceBoxes(std::vector<std::array<TPoint<T>,8>> &boxes) const
{
  boxes.clear();

  for (int i = 0; i < numFaces(); i++)
  {
    std::array<TPoint<T>,8> box;
    faceBox(i,box);

    boxes.push_back(box);
  }
}

template <class T> bool TTriangles<T>::segIntersect(const TPoint<T> &point0, const TPoint<T> &point1, 
  const int tri, T &U, TPoint<T> &intersection, const T tolerance, const T parmtolerance)
{
  std::array<TPoint<T>,3> c;
  c[0] = coords[corners[tri * 3 + 0]];
  c[1] = coords[corners[tri * 3 + 1]];
  c[2] = coords[corners[tri * 3 + 2]];

  return segTriIntersect(point0,point1,c,U,intersection,tolerance,parmtolerance);
}

template <class T> void TTriangles<T>::getTriNeighbours(const int numLayersAround, 
  const std::vector<std::vector<LINT>> &cornerTris, std::set<LINT> &triNumbers)
{
  // trinumbers should not be empty (must contain at least one starting tri number)
  assert(triNumbers.size() > 0);

  for (int i = 0; i < numLayersAround; i++)
  {
    std::set<LINT> newtrinumbers;
    for (auto tri : triNumbers)
    {
      for (int j = 0; j < 3; j++)
      {
        LINT corner = corners[tri * 3 + j];
        for (auto it : cornerTris[corner])
        {
          newtrinumbers.insert(it);
        }
      }
    }
    triNumbers.insert(newtrinumbers.begin(),newtrinumbers.end());
  }
}

template <class T> bool TTriangles<T>::generatePoints(int pointsPerTri, bool excludeDuplicates, 
  T tolerance, std::vector<TPoint<T>> &points)
{
  int numSteps = (int) sqrt(double(pointsPerTri));
  LIMIT_MIN(numSteps,1);
                            // got through all tris
  for (int i = 0; i < numFaces(); i++)
  {
     int i3 = i * 3;
                            // three triangle nodes
    TPoint<T> co[3];
    co[0] = coords[corners[i3]];
    co[1] = coords[corners[i3 + 1]];
    co[2] = coords[corners[i3 + 2]];

    T d = 1.0 / T(numSteps);
    for (int j = 0; j < numSteps + 1; j++)
    {
      T u = T(j) * d;
      for (int k = 0; k < numSteps + 1; k++)
      {
        T v = T(k) * d;
        T w = 1.0 - u - v;

        if (w >= 0.0)
        {
          TPoint<T> p = co[0] * u + co[1] * v + co[2] * w;
          points.push_back(p);
        }
      }
    }
  }

  if (excludeDuplicates)
  {
    if (tolerance == 0.0)
    {
      std::pair<TPoint<T>,TPoint<T>> mm = minmax();
      TPoint<T> d = mm.second - mm.first;
      tolerance = !d * PARM_TOLERANCE;
    }
                            // remove duplicate nodes and renumber corners
    removeDuplicates(coords,true,tolerance);
  }

  return true;
}


template <class T> void refineTriangle(const TPoint<T> corner0, 
  const TPoint<T> corner1, const TPoint<T> corner2, 
  const T fraction, std::vector<TPoint<T>> &points)
{
  // for every corner, two other corners
  static int sides[3][2] = {{1,2},{2,0},{0,1}};

  // corners
  TPoint<T> c[3] = {corner0, corner1, corner2};

  // num divisions of all sides
  int divs[3];
  // side length in front of every corner
  T len[3];
  for (int i = 0; i < 3; i++)
  {
    len[i] = !(c[sides[i][1]] - c[sides[i][0]]);
    divs[i] = int(len[i] / fraction);
    LIMIT_MIN(divs[i],1);
  }

  // the shortest side?
  int start = 0; 
  T startLen = len[0];
  for (int i = 1; i < 3; i++)
  {
    if (len[i] < startLen)
    {
      startLen = len[i];
      start = i;
    }
  }

  // now go along two longest side in front of the shortest
  int ds = std::max<int>(divs[sides[start][0]],divs[sides[start][1]]);

  TPoint<T> vstart = c[start];
  TPoint<T> vend0 = c[sides[start][0]];
  TPoint<T> vend1 = c[sides[start][1]];

  points.push_back(vstart);

  for (int i = 1; i <= ds; i++)
  {
    T U = T(i) / T(ds);

    TPoint<T> v0 = vstart + (vend0 - vstart) * U;
    TPoint<T> v1 = vstart + (vend1 - vstart) * U;

    // cross direction
    T clen = !(v1 - v0);
    int cdivs = int(clen / fraction);
    LIMIT_MIN(cdivs,1);

    for (int j = 0; j <= cdivs; j++)
    {
      T V = T(j) / T(cdivs);
      TPoint<T> v = v0 + (v1 - v0) * V;
      points.push_back(v);
    }
  }
}

template <class T> bool TTriangles<T>::generatePointsByRefine(T fraction, bool excludeDuplicates, 
  T tolerance, std::vector<TPoint<T>> &points)
{
                            // got through all tris
  for (int i = 0; i < numFaces(); i++)
  {
     int i3 = i * 3;
                            // three triangle nodes
    TPoint<T> co[3];
    co[0] = coords[corners[i3]];
    co[1] = coords[corners[i3 + 1]];
    co[2] = coords[corners[i3 + 2]];

    refineTriangle(co[0],co[1],co[2],fraction,points);
  }

  if (excludeDuplicates)
  {
    if (tolerance == 0.0)
    {
      std::pair<TPoint<T>,TPoint<T>> mm = minmax();
      TPoint<T> d = mm.second - mm.first;
      tolerance = !d * PARM_TOLERANCE;
    }
                            // remove duplicate nodes and renumber corners
    removeDuplicates(coords,true,tolerance);
  }

  return true;
}

template <class T> void TTriangles<T>::getBoundaryEdges(T tolerance, std::vector<std::pair<LINT,LINT>> &edges)
{
  edges.clear();

  // make unique node numbering
  if (!duplicatesExcluded())
  {
    buildConnectivityArray(tolerance);
  }

  std::map<std::pair<LINT,LINT>, std::vector<LINT>,EdgeCompare> edgeTris;
  getEdgeTris(edgeTris);

  for (auto e : edgeTris)
  {
    if (e.second.size() == 1)
    {
      edges.push_back(e.first);
    }
  }
}

template <class T> void TTriangles<T>::edgesToNodes(std::vector<std::pair<LINT,LINT>> &edges, std::set<LINT> &nodes)
{
  nodes.clear();

  for (auto &e : edges)
  {
    nodes.insert(e.first);
    nodes.insert(e.second);
  }
}

template <class T> bool TTriangles<T>::manifold(T tolerance, std::vector<std::pair<LINT,LINT>> &badedges,
  std::vector<std::pair<LINT,LINT>> *boundaryedges)
{
  bool ok = true;

  badedges.clear();

  // make unique node numbering
  if (!duplicatesExcluded())
  {
    buildConnectivityArray(tolerance);
  }

  if (ok)
  {
    std::map<std::pair<LINT,LINT>, std::vector<LINT>,EdgeCompare> edgeTris;
    getEdgeTris(edgeTris);

    for (auto e : edgeTris)
    {
      if (e.second.size() > 2)
      {
        badedges.push_back(e.first);
        ok = false;
      }
    }

    if (boundaryedges)
      getBoundaryEdges(tolerance,*boundaryedges);
  }

  return ok;
}

template <class T> bool TTriangles<T>::solid(T tolerance, std::vector<std::pair<LINT,LINT>> &badedges,
  std::vector<std::pair<LINT,LINT>> *boundaryedges)
{
  bool ok = true;

  badedges.clear();

  // make unique node numbering
  if (!duplicatesExcluded())
  {
    buildConnectivityArray(tolerance);
  }

  if (ok)
  {
    std::map<std::pair<LINT,LINT>, std::vector<LINT>,EdgeCompare> edgeTris;
    getEdgeTris(edgeTris);

    for (auto e : edgeTris)
    {
      if (e.second.size() != 2)
      {
        badedges.push_back(e.first);
        ok = false;
      }
    }

    if (boundaryedges)
      getBoundaryEdges(tolerance,*boundaryedges);
  }

  return ok;
}

template <class T> int TTriangles<T>::removeBadEdges(std::vector<std::pair<LINT,LINT>> &badedges,
  std::vector<std::pair<LINT,LINT>> &boundaryedges)
{
  int count = 0;

  std::set<LINT> excludenodes;
  edgesToNodes(badedges,excludenodes);

  std::set<LINT> boundarynodes;
  edgesToNodes(boundaryedges,boundarynodes);

  for (auto riter = excludenodes.rbegin(); riter != excludenodes.rend(); ++riter) 
  {
    // do not delete boundary nodes
    if (boundarynodes.find(*riter) == boundarynodes.end())
    {
      coords.erase(coords.begin() + (*riter));
      count++;
    }
  }

  return count;
}

template <class T> bool TTriangles<T>::makeNACA0012(LINT numX, LINT numZ, T span, int fitdegree)
{
  // rescale to -0.5, +0.5]
  std::vector<T> x,y;
  for (auto p : NACA0012xy<T>)
  {
    x.push_back(p.first - 0.5);
    y.push_back(p.second);
  }
 
  TJacobiPoly<T> f(0.0,0.5);

  bool res = f.fit(fitdegree,x,y,GAUSSINT_20);

  T acc = f.accuracy(x,y);

  if (res)
  {
    T DU = 2.0 / numX;
    T Zmin = -span * 0.5;
    T Zmax = +span * 0.5;
    T DZ = span / T(numZ);
    T Z0 = Zmin;

    for (int i = 0; i < numZ; i++)
    {
      T Z1 = Z0 + DZ;

      T U0 = -1.0;
      for (int j = 0; j < numX; j++)
      {
        T U1 = U0 + DU;
        LIMIT_MAX(U1,+1.0);

        T Y0 = f.getValue(U0);
        T Y1 = f.getValue(U1);

        if (j == 0)
        {
          Y0 = 0.0;
        }
        if (j == numX - 1)
        {
          Y1 = 0.0;
        }

        T X0 = U0 * 0.5;
        T X1 = U1 * 0.5;

        // upper surface
        addTri(TPoint<T>(X0,Y0,Z0),TPoint<T>(X1,Y1,Z0),TPoint<T>(X1,Y1,Z1),0.0,0.0,0.0,true);
        addTri(TPoint<T>(X0,Y0,Z0),TPoint<T>(X1,Y1,Z1),TPoint<T>(X0,Y0,Z1),0.0,0.0,0.0,true);

        // lower surface
        addTri(TPoint<T>(X1,-Y1,Z1),TPoint<T>(X1,-Y1,Z0),TPoint<T>(X0,-Y0,Z0),0.0,0.0,0.0,true);
        addTri(TPoint<T>(X0,-Y0,Z1),TPoint<T>(X1,-Y1,Z1),TPoint<T>(X0,-Y0,Z0),0.0,0.0,0.0,true);

        if (i == 0)
        {
          if (j == 0)
          {
            addTri(TPoint<T>(X0,-Y0,Z0),TPoint<T>(X1,-Y1,Z0),TPoint<T>(X1,Y1,Z0),0.0,0.0,0.0,true);
          } else if (j == numX - 1)
          {
            addTri(TPoint<T>(X0,-Y0,Z0),TPoint<T>(X1,-Y1,Z0),TPoint<T>(X0,Y0,Z0),0.0,0.0,0.0,true);
          } else
          {
            addTri(TPoint<T>(X0,-Y0,Z0),TPoint<T>(X1,-Y1,Z0),TPoint<T>(X1,Y1,Z0),0.0,0.0,0.0,true);
            addTri(TPoint<T>(X0,-Y0,Z0),TPoint<T>(X1,Y1,Z0),TPoint<T>(X0,Y0,Z0),0.0,0.0,0.0,true);
          }
        }

        if (i == numZ - 1)
        {
          if (j == 0)
          {
            addTri(TPoint<T>(X1,Y1,Z1),TPoint<T>(X1,-Y1,Z1),TPoint<T>(X0,-Y0,Z1),0.0,0.0,0.0,true);
          } else if (j == numX - 1)
          {
            addTri(TPoint<T>(X0,Y0,Z1),TPoint<T>(X1,-Y1,Z1),TPoint<T>(X0,-Y0,Z1),0.0,0.0,0.0,true);
          } else
          {
            addTri(TPoint<T>(X1,Y1,Z1),TPoint<T>(X1,-Y1,Z1),TPoint<T>(X0,-Y0,Z1),0.0,0.0,0.0,true);
            addTri(TPoint<T>(X0,Y0,Z1),TPoint<T>(X1,Y1,Z1),TPoint<T>(X0,-Y0,Z1),0.0,0.0,0.0,true);
          }
        }

        U0 += DU;
      }

      Z0 += DZ;
    }
  }

  return res;
}

template <class T> T TTriangles<T>::minEdge()
{
  T minedge = std::numeric_limits<T>::max();
  for (int i = 0; i < coords.size(); i += 3)
  {
    T e0 = !(coords[i + 1] - coords[i]);
    T e1 = !(coords[i + 2] - coords[i + 1]);
    T e2 = !(coords[i] - coords[i + 2]);

    minedge = std::min<T>(minedge,std::min<T>(e0,std::min<T>(e1,e2)));
  }

  return minEdge;
}

template <class T> void TTriangles<T>::makeQuad(TPoint<T> corners[4], LINT numU, LINT numV)
{
  for (int i = 0; i < numU; i++)
  {
    T U0 = T(i) / T(numU) * 2.0 - 1.0;
    T U1 = T(i + 1) / T(numU) * 2.0 - 1.0;
    for (int j = 0; j < numV; j++)
    {
      T V0 = T(j) / T(numV) * 2.0 - 1.0;
      T V1 = T(j + 1) / T(numV) * 2.0 - 1.0;

      TPoint<T> coords[4];
      coords[0] = rectCoord(corners,U0,V0);
      coords[1] = rectCoord(corners,U1,V0);
      coords[2] = rectCoord(corners,U1,V1);
      coords[3] = rectCoord(corners,U0,V1);

      T d0 = !(coords[0] - coords[2]);
      T d1 = !(coords[1] - coords[3]);

      if (d0 < d1)
      {
        addTri(coords[0],coords[1],coords[2],0.0);
        addTri(coords[0],coords[2],coords[3],0.0);
      } else
      {
        addTri(coords[0],coords[1],coords[3],0.0);
        addTri(coords[1],coords[2],coords[3],0.0);
      }
    }
  }
}

template <class T> void TTriangles<T>::makeTransform(TTransform<T> &transform)
{
  for (LINT i = 0; i < int(coords.size()); i++)
  {
    coords[i] = transform.applyTransform(coords[i]);
  }
}

template <class T> TPoint<T> TTriangles<T>::closestPoint(LINT faceNo, 
  const TPoint<T> &point, int &closestpos) const
{
  // corners
  std::array<TPoint<T>,3> corners = threeCorners(faceNo);

  // unit normal
  TPoint<T> normal = faceNormal(faceNo);

  // projection on plane
  TPoint<T> p = point - normal * (normal * (point - corners[0]));

  // edges
  TPoint<T> e01 = corners[1] - corners[0];
  TPoint<T> e12 = corners[2] - corners[1];
  TPoint<T> e20 = corners[0] - corners[2];

  // project on edge 01
  T E0 = ((corners[1] - p) ^ e12) * normal;
  if (E0 < 0.0) 
  {
    closestpos = CLOSEST_EDGE01;

    TPoint<T> e12n = +e12;
    p = corners[2] - e12n * (e12n * (corners[2] - p));

    if ((corners[1] - p) * e12 > 0.0) 
    {
      return corners[1];
    } else if ((corners[2] - p) * e12 < 0.0) 
    {
      return corners[2];
    } else 
    {
      return p;
    }
  }

  // project on edge 12
  T E1 = ((corners[2] - p) ^ e20) * normal;
  if (E1 < 0.0) 
  {
    closestpos = CLOSEST_EDGE12;

    TPoint<T> e20n = +e20;
    p = corners[0] - e20n * (e20n * (corners[0] - p));

    if ((corners[2] - p) * e20 > 0.0) 
    {
      return corners[2];
    } else if ((corners[0] - p) * e20 < 0.0) 
    {
      return corners[0];
    } else 
    {
      return p;
    }
  }

  // project on edge 20
  T E2 = ((corners[0] - p) ^ e01) * normal;
  if (E2 < 0.0) 
  {
    closestpos = CLOSEST_EDGE20;

    TPoint<T> e01n = +e01;
    p = corners[1] - e01n * (e01n * (corners[1] - p));

    if ((corners[0] - p) * e01 > 0.0) 
    {
      return corners[0];
    } else if ((corners[1] - p) * e01 < 0.0) 
    {
      return corners[1];
    } else 
    {
      return p;
    }
  }

  // inside triangle
  closestpos = CLOSEST_INSIDE;

  return p;
}

template <class T> bool TTriangles<T>::pointInsideFaceXY(LINT faceNo, const TPoint<T> &point, T tolerance) const
{
  std::array<TPoint<T>,3> c = threeCorners(faceNo);

  // signed distances
  T Z0 = ((+(c[1] - c[0])) ^ (point - c[0])).Z;
  T Z1 = ((+(c[2] - c[1])) ^ (point - c[1])).Z;
  T Z2 = ((+(c[0] - c[2])) ^ (point - c[2])).Z;

  return (sign(Z0) == sign(Z1)) && (sign(Z1) == sign(Z2));
}

template <class T> bool TTriangles<T>::pointInsideXY(const TPoint<T> &point, T tolerance) const
{
  // if point is among corners, it is inside
  for (int i = 0; i < coords.size(); i++)
  {
    T dist = !(coords[i] - point);
    if (dist < tolerance)
      return true;
  }

  for (int i = 0; i < numFaces(); i++)
  {
    if (pointInsideFaceXY(i,point,tolerance))
      return true;
  }

  return false;
}

template <class T> void TTriangles<T>::addArrow(TPoint<T> p, TPoint<T> velocity)
{
  // create two tris

  TPoint<T> pend = p + velocity;
  TPoint<T> dir = +velocity;
  T len = !velocity;
  T taillen = len * 0.05;
  TPoint<T> pdir0 = +(dir ^ TPoint<T>(0,0,1));
  TPoint<T> p0 = p + pdir0 * taillen;
  TPoint<T> p1 = p - pdir0 * taillen;
  TPoint<T> pdir1 = +(dir ^ pdir0);
  TPoint<T> p2 = p + pdir1 * taillen;
  TPoint<T> p3 = p - pdir1 * taillen;

  addTri(p0,p1,pend,0.0);
  addTri(p2,p3,pend,0.0);
}

template <class T> T TTriangles<T>::signedDist(LINT faceNo, const TPoint<T> &point) const
{
  TPoint<T> normal = faceNormal(faceNo);
  TPoint<T> centre = faceCentre(faceNo);
  TPlane<T> plane(normal,centre);

  // signed
  T dist = plane.distance(point);
  return dist;
}

template <class T> bool TTriangles<T>::intersectByPlane(TPlane<T> &plane, 
  std::vector<std::vector<TPoint<T>>> &lines, T tolerance, T parmtolerance,
  std::vector<std::vector<TPoint<T>>> *boundary)
{
  // return value
  bool ok = false;

  // no lines
  lines.clear();

  // intersection pieces each of two points
  std::vector<std::vector<TPoint<T>>> pieces;

  // parametric values for intersection points, X,Y - U,V 
  std::vector<TPoint<T>> UVintrs;

  bool makeboundaries = (boundary && !UVcorners.empty());

  if (makeboundaries)
  {
    for (LINT i = 0; i < numFaces(); i++)
    {
      std::array<TPoint<T>,3> corners = threeCorners(i);

      // UVcorners must be filled in TBaseSurface::createTriangles()
      std::array<TPoint<T>,3> uvcorners = UVcorners[i];

      std::vector<TPoint<T>> intrs;
      if (cutFaceByPlane(plane,i,intrs,tolerance) == 2)
      {
        TPoint<T> coord0 = barycentricCoord(corners,intrs[0]);
        TPoint<T> coord1 = barycentricCoord(corners,intrs[1]);

        // these will contain U,V in X,Y for every intersection point
        TPoint<T> uv0,uv1;
        for (int k = 0; k < 3; k++)
        {
          uv0 += uvcorners[k] * coord0.XYZW[k];
          uv1 += uvcorners[k] * coord1.XYZW[k];
        }

        // keep U,V for every intersection
        UVintrs.push_back(TPoint<T>(uv0.X,uv0.Y));
        intrs[0].W = T(UVintrs.size() - 1);
        UVintrs.push_back(TPoint<T>(uv1.X,uv1.Y));
        intrs[1].W = T(UVintrs.size() - 1);

        pieces.push_back(intrs);
      }
    }
  } else
  {
    for (LINT i = 0; i < numFaces(); i++)
    {
      std::vector<TPoint<T>> intrs;
      if (cutFaceByPlane(plane,i,intrs,tolerance) == 2)
      {
        pieces.push_back(intrs);
      }
    }
  }

  // order a curve from two-point pieces
  ok = curvesFromPieces(pieces,lines,tolerance,false);

  if (makeboundaries)
  {
    for (int i = 0; i < int(lines.size()); i++)
    {
      boundary->push_back(std::vector<TPoint<T>>());
      for (int j = 0; j < int(lines[i].size()); j++)
      {
        boundary->back().push_back(TPoint<T>(UVintrs[ROUND(lines[i][j].W)].X,UVintrs[ROUND(lines[i][j].W)].Y));
      }
    }
  }

  return ok;
}

#ifdef GLOBAL_INTRCHECK
  extern std::mutex imutex;
  extern std::mutex emutex;
#endif

struct pair_hash {
  inline std::size_t operator()(const std::pair<LINT,LINT> &v) const {
    return v.first * 31 + v.second;
  }
};

/** Intersect tris which are in one cell cellindex. */
template <class T> void intersectCells(TTriangles<T> *tris, TTriangles<T> *other, LINT cellindex,
  std::vector<std::set<LINT>> *celltris, std::vector<std::set<LINT>> *ocelltris,
  std::vector<TPoint<T>> *centres, std::vector<TPoint<T>> *ocentres,
  // output :
  std::vector<std::pair<LINT,LINT>> *intersected, std::vector<TPoint<T>> *UVintrs,
  std::vector<std::pair<TPoint<T>,TPoint<T>>> *edges,
  T tolerance, T parmtolerance, std::unordered_set<std::pair<LINT,LINT>,pair_hash> *sintersected = nullptr)
{
  for (LINT index0 : (*celltris)[cellindex])
  {
    TPoint<T> c = (*centres)[index0];

    std::array<TPoint<T>,3> corners = tris->threeCorners(index0);

    // UVcorners must be filled in TBaseSurface::createTriangles()
    std::array<TPoint<T>,3> uvcorners = tris->UVcorners[index0];

    for (LINT index1 : (*ocelltris)[cellindex])
    {
      bool done = false;

#ifdef GLOBAL_INTRCHECK
      if (sintersected)
      {
        imutex.lock();

        done = (sintersected && 
          (sintersected->find(std::pair<LINT,LINT>(index0,index1)) != sintersected->end() ||
          sintersected->find(std::pair<LINT,LINT>(index1,index0)) != sintersected->end()));

        imutex.unlock();
      }
#endif

      if (done)
        continue;

      // check by sizes
      TPoint<T> oc = (*ocentres)[index1];
      T d = !(oc - c);
      if (d > oc.W + c.W)
        continue;

      std::array<TPoint<T>,3> ocorners = other->threeCorners(index1);
      // UVcorners must be filled in TBaseSurface::createTriangles()
      std::array<TPoint<T>,3> ouvcorners = other->UVcorners[index1];

      std::vector<TPoint<T>> intrs;
      if (intersectTriangleByTriangle(corners,ocorners,intrs,tolerance,parmtolerance))
      {

//!!!!!!!
//bool print = false;
//if (print)
//{
//  tris->saveFaceSTL(int(index0),"index0.stl","TCAD",false);
//  other->saveFaceSTL(int(index1),"index1.stl","TCAD",false);
//}
        if (intrs.size() == 2)
        {

// outputDebugString(to_string(int(index0)) + " " + to_string(int(index1)));

          TPoint<T> coord0 = barycentricCoord(corners,intrs[0]);
          TPoint<T> coord1 = barycentricCoord(corners,intrs[1]);
          TPoint<T> ocoord0 = barycentricCoord(ocorners,intrs[0]);
          TPoint<T> ocoord1 = barycentricCoord(ocorners,intrs[1]);

          // these will contain U,V in X,Y for every intersection point
          TPoint<T> uv0,uv1,ouv0,ouv1;
          for (int k = 0; k < 3; k++)
          {
            uv0 += uvcorners[k] * coord0.XYZW[k];
            uv1 += uvcorners[k] * coord1.XYZW[k];
            ouv0 += ouvcorners[k] * ocoord0.XYZW[k];
            ouv1 += ouvcorners[k] * ocoord1.XYZW[k];
          }

#ifdef GLOBAL_INTRCHECK
          emutex.lock();
#endif

          // keep U,V for every intersection on both surfaces
          UVintrs->push_back(TPoint<T>(uv0.X,uv0.Y,ouv0.X,ouv0.Y));
          intrs[0].W = T(UVintrs->size() - 1);
          UVintrs->push_back(TPoint<T>(uv1.X,uv1.Y,ouv1.X,ouv1.Y));
          intrs[1].W = T(UVintrs->size() - 1);

          edges->push_back(std::pair<TPoint<T>,TPoint<T>>(intrs[0],intrs[1]));

#ifdef GLOBAL_INTRCHECK
          emutex.unlock();
#endif

          if (intersected)
            intersected->push_back(std::pair<LINT,LINT>(index0,index1));
        }
      }

#ifdef GLOBAL_INTRCHECK
      if (sintersected)
      {
        imutex.lock();
        sintersected->insert(std::pair<LINT,LINT>(index0,index1));
        imutex.unlock();
      }
#endif
    }
  }
}

/** Intersect range of cells. */
template <class T> void intersectMultiCells(std::vector<LINT> *activecells, LINT cell0, LINT cell1,
  TTriangles<T> *tris, TTriangles<T> *other, 
  std::vector<std::set<LINT>> *celltris, std::vector<std::set<LINT>> *ocelltris,
  std::vector<TPoint<T>> *centres, std::vector<TPoint<T>> *ocentres,
  // output :
  std::vector<std::pair<LINT,LINT>> *intersected, std::vector<TPoint<T>> *UVintrs,
  std::vector<std::pair<TPoint<T>,TPoint<T>>> *edges,
  T tolerance, T parmtolerance, std::unordered_set<std::pair<LINT,LINT>,pair_hash> *sintersected = nullptr)
{
  for (LINT i = cell0; i <= cell1; i++)
  {
    LINT index = (*activecells)[i];
    intersectCells(tris,other,index,celltris,ocelltris,centres,ocentres,
      intersected,UVintrs,edges,tolerance,parmtolerance,sintersected);
  }
}

/** Intersect tris which are in one cell cellindex. */
template <class T> void intersectCellsNoUV(TTriangles<T> *tris, TTriangles<T> *other, LINT cellindex,
  std::vector<std::set<LINT>> *celltris, std::vector<std::set<LINT>> *ocelltris,
  std::vector<TPoint<T>> *centres, std::vector<TPoint<T>> *ocentres,
  // output :
  std::vector<std::pair<LINT,LINT>> *intersected, std::vector<TPoint<T>> *UVintrs,
  std::vector<std::pair<TPoint<T>,TPoint<T>>> *edges,
  T tolerance, T parmtolerance, std::unordered_set<std::pair<LINT,LINT>,pair_hash> *sintersected = nullptr)
{
  for (LINT index0 : (*celltris)[cellindex])
  {
    TPoint<T> c = (*centres)[index0];

    std::array<TPoint<T>,3> corners = tris->threeCorners(index0);

    for (LINT index1 : (*ocelltris)[cellindex])
    {
      bool done = false;

#ifdef GLOBAL_INTRCHECK
      if (sintersected)
      {
        imutex.lock();

        done = (sintersected && 
          (sintersected->find(std::pair<LINT,LINT>(index0,index1)) != sintersected->end() ||
          sintersected->find(std::pair<LINT,LINT>(index1,index0)) != sintersected->end()));

        imutex.unlock();
      }
#endif

      if (done)
        continue;

      // check by sizes
      TPoint<T> oc = (*ocentres)[index1];
      T d = !(oc - c);
      if (d > oc.W + c.W)
        continue;

      std::array<TPoint<T>,3> ocorners = other->threeCorners(index1);

      std::vector<TPoint<T>> intrs;
      if (intersectTriangleByTriangle(corners,ocorners,intrs,tolerance,parmtolerance))
      {
        if (intrs.size() == 2)
        {
#ifdef GLOBAL_INTRCHECK
          emutex.lock();
#endif

          edges->push_back(std::pair<TPoint<T>,TPoint<T>>(intrs[0],intrs[1]));

#ifdef GLOBAL_INTRCHECK
          emutex.unlock();
#endif

          if (intersected)
            intersected->push_back(std::pair<LINT,LINT>(index0,index1));
        }
      }

#ifdef GLOBAL_INTRCHECK
      if (sintersected)
      {
        imutex.lock();
        sintersected->insert(std::pair<LINT,LINT>(index0,index1));
        imutex.unlock();
      }
#endif

    }
  }
}

/** Intersect range of cells. */
template <class T> void intersectMultiCellsNoUV(std::vector<LINT> *activecells, LINT cell0, LINT cell1,
  TTriangles<T> *tris, TTriangles<T> *other, 
  std::vector<std::set<LINT>> *celltris, std::vector<std::set<LINT>> *ocelltris,
  std::vector<TPoint<T>> *centres, std::vector<TPoint<T>> *ocentres,
  // output :
  std::vector<std::pair<LINT,LINT>> *intersected, std::vector<TPoint<T>> *UVintrs,
  std::vector<std::pair<TPoint<T>,TPoint<T>>> *edges,
  T tolerance, T parmtolerance, std::unordered_set<std::pair<LINT,LINT>,pair_hash> *sintersected = nullptr)
{
  for (LINT i = cell0; i <= cell1; i++)
  {
    LINT index = (*activecells)[i];
    intersectCellsNoUV(tris,other,index,celltris,ocelltris,centres,ocentres,
      intersected,UVintrs,edges,tolerance,parmtolerance,sintersected);
  }
}

template <class T> bool TTriangles<T>::intersect(TTriangles<T> &other, 
  std::vector<std::vector<TPoint<T>>> &lines, T parmtolerance,
  std::vector<std::vector<TPoint<T>>> *boundary0,
  std::vector<std::vector<TPoint<T>>> *boundary1,
  int numthreads)
{
  // no lines
  lines.clear();

  // intersection pairs of points
  std::vector<std::pair<TPoint<T>,TPoint<T>>> edges;

  // parametric values for intersection points, X,Y - U,V for first surface,
  // Z,W - U,V for second surface
  std::vector<TPoint<T>> UVintrs;

  // make parametric boundaries
  bool makeboundaries = boundary0 && boundary1 && !UVcorners.empty() && !other.UVcorners.empty();

  // for collision test
  std::vector<TPoint<T>> centres,ocentres;
  getCentresAndRadii(centres);
  other.getCentresAndRadii(ocentres);

  // max tri size for cell size
  T minedge,maxedge;
  getEdgeMinMax(minedge,maxedge);
  T ominedge,omaxedge;
  other.getEdgeMinMax(ominedge,omaxedge);

  maxedge = std::max(maxedge,omaxedge) * 1.1; //!!!!!!!

  // this tolerance below is only used to remove duplicate intersection points,
  // no other purpose
  T tolerance = maxedge * parmtolerance;

  // use octree background cells for spacial partitioning...
#ifdef USE_SPACEPARTITIONING 

  // min/max of all triangles
  std::pair<TPoint<T>,TPoint<T>> mm0 = minmax();
  std::pair<TPoint<T>,TPoint<T>> mm1 = other.minmax();
  TPoint<T> min = pointMin(mm0.first,mm1.first);
  TPoint<T> max = pointMax(mm0.second,mm1.second);

  // extend
  extendMinMax(min,max,1.01);
  TPoint<T> dmm = max - min;

  // no intersection possible
  if (dmm.X < tolerance || dmm.Y < tolerance || dmm.Z < tolerance)
    return false;

  // cells big enough for spacial partitioning to identify if two triangles
  // may intersect, they can if they have a node in the same cell only
  OBackground<T> cells(min,max,(LINT) (dmm.X / maxedge),(LINT) (dmm.Y / maxedge),(LINT) (dmm.Z / maxedge),8);

  // total number of cells
  LINT numcells = cells.numBackgroundCells();

  // now distribute all tris over cells
  std::vector<std::set<LINT>> celltris(numcells);
  for (LINT i = 0; i < numFaces(); i++)
  {
    //std::array<TPoint<T>,3> corners = threeCorners(i);
   
    //for (int j = 0; j < 3; j++)
    //{
    //  LINT index = cells.findCell(corners[j]);
    //  assert(index >= 0);

    //  celltris[index].insert(i);
    //}

    // important, do not remove
    std::array<TPoint<T>,8> box;
    faceBox(i,box);

    for (int j = 0; j < 8; j++)
    {
      LINT index = cells.findCell(box[j]);
      assert(index >= 0);

      celltris[index].insert(i);
    }
  }

  std::vector<std::set<LINT>> ocelltris(numcells);
  for (LINT i = 0; i < other.numFaces(); i++)
  {
    //std::array<TPoint<T>,3> corners = other.threeCorners(i);
    
    //for (int j = 0; j < 3; j++)
    //{
    //  LINT index = cells.findCell(corners[j]);
    //  assert(index >= 0);

    //  ocelltris[index].insert(i);
    //}

    // important, do not remove
    std::array<TPoint<T>,8> box;
    other.faceBox(i,box);

    for (int j = 0; j < 8; j++)
    {
      LINT index = cells.findCell(box[j]);
      assert(index >= 0);

      ocelltris[index].insert(i);
    }
  }

  // now distinguish cells which are (1) not empty (2) contain triangles from
  // this and other triangles
  std::vector<LINT> activecells;
  for (LINT i = 0; i < numcells; i++)
  {
    if (!celltris[i].empty() && !ocelltris[i].empty())
    {
      activecells.push_back(i);
    }
  }
  int numactive = int(activecells.size());

  if (numactive == 0)
    return false;

  // now we've got two plain lists of triangles for every cell

#ifndef USE_THREADS
  numthreads = 1;
#endif

  if (numactive < numthreads)
    numthreads = numactive;
  LIMIT_MIN(numthreads,1);

#ifdef DEBUG_TRIS
  outputDebugString(
    std::string("num threads = ") + to_string(numthreads) +
    std::string(" num cells = ") + to_string(numcells) +
    std::string(" active cells = ") + to_string(numactive));
#endif

  if (numthreads <= 1)
  {
#ifdef DEBUG_TRIS
  outputDebugString(
    std::string(" active cell size = ") + 
    to_string(int(celltris[activecells[0]].size())) + " " +
    to_string(int(ocelltris[activecells[0]].size())));
#endif

    std::vector<std::vector<std::pair<TPoint<T>,TPoint<T>>>> tedges(numactive);
    std::vector<std::vector<TPoint<T>>> tUVintrs(numactive);
    std::vector<std::vector<std::pair<LINT,LINT>>> intersected(numactive);

    std::unordered_set<std::pair<LINT,LINT>,pair_hash> sintersected;

    if (makeboundaries)
    {
      for (LINT i = 0; i < numactive; i++)
      {
        LINT index = activecells[i];

        intersectCells(this,&other,index,&celltris,&ocelltris,&centres,&ocentres,
#ifdef GLOBAL_INTRCHECK
          nullptr,&UVintrs,&edges,tolerance,parmtolerance,&sintersected);
#else
          &intersected[i],&tUVintrs[i],&tedges[i],tolerance,parmtolerance,nullptr);
#endif
      }
    } else
    {
      for (LINT i = 0; i < numactive; i++)
      {
        LINT index = activecells[i];
        intersectCellsNoUV(this,&other,index,&celltris,&ocelltris,&centres,&ocentres,
#ifdef GLOBAL_INTRCHECK
          nullptr,&UVintrs,&edges,tolerance,parmtolerance,&sintersected);
#else
          &intersected[i],&tUVintrs[i],&tedges[i],tolerance,parmtolerance,nullptr);
#endif
      }
    }

#ifndef GLOBAL_INTRCHECK
    int count = 0;

    std::set<std::pair<LINT,LINT>,EdgeCompare> intrset;

    for (int t = 0; t < numactive; t++)
    {
      for (int k = 0; k < int(tedges[t].size()); k++)
      {
        // avoid inserting same intersected tris
        auto res = intrset.insert(intersected[t][k]);

        if (res.second)
        {
          if (makeboundaries)
          {
            UVintrs.push_back(tUVintrs[t][k * 2]);
            UVintrs.push_back(tUVintrs[t][k * 2 + 1]);
          }

          std::pair<TPoint<T>,TPoint<T>> edge = tedges[t][k];
          edge.first.W = T(count++);
          edge.second.W = T(count++);
          edges.push_back(edge);
        }
      }
    }
#endif

  } else
  {
    std::vector<std::vector<std::pair<TPoint<T>,TPoint<T>>>> tedges(numthreads);
    std::vector<std::vector<TPoint<T>>> tUVintrs(numthreads);
    std::vector<std::vector<std::pair<LINT,LINT>>> intersected(numthreads);

    std::unordered_set<std::pair<LINT,LINT>,pair_hash> sintersected;

    // distribute over threads
    std::vector<std::thread> threads;
    int dt = int(numactive) / numthreads;
    for (int t = 0; t < numthreads; t++)
    {
      int j1 = t * dt;
      int j2 = (t + 1) * dt - 1;

      if (t == numthreads - 1)
        j2 = int(numactive) - 1;

      std::thread th;
      if (makeboundaries)
      {
        th = std::thread(intersectMultiCells<T>,&activecells,j1,j2,this,&other,
          &celltris,&ocelltris,&centres,&ocentres,
#ifdef GLOBAL_INTRCHECK
          nullptr,&UVintrs,&edges,tolerance,parmtolerance,&sintersected);
#else
          &intersected[t],&tUVintrs[t],&tedges[t],tolerance,parmtolerance,nullptr);
#endif
      } else
      {
        th = std::thread(intersectMultiCellsNoUV<T>,&activecells,j1,j2,this,&other,
          &celltris,&ocelltris,&centres,&ocentres,
#ifdef GLOBAL_INTRCHECK          
          nullptr,&UVintrs,&edges,tolerance,parmtolerance,&sintersected);
#else
          &intersected[t],&tUVintrs[t],&tedges[t],tolerance,parmtolerance,nullptr);
#endif
      }

      threads.push_back(std::move(th));
    }

    for (auto &th : threads) {     
      th.join();
    }

#ifndef GLOBAL_INTRCHECK
    int count = 0;

    std::set<std::pair<LINT,LINT>,EdgeCompare> intrset;

    for (int t = 0; t < numthreads; t++)
    {
      for (int k = 0; k < int(tedges[t].size()); k++)
      {
        // avoid inserting same intersected tris
#ifndef GLOBAL_INTRCHECK
        auto res = intrset.insert(intersected[t][k]);

        if (res.second)
#endif
        {
          if (makeboundaries)
          {
            UVintrs.push_back(tUVintrs[t][k * 2]);
            UVintrs.push_back(tUVintrs[t][k * 2 + 1]);
          }

          std::pair<TPoint<T>,TPoint<T>> edge = tedges[t][k];
          edge.first.W = T(count++);
          edge.second.W = T(count++);
          edges.push_back(edge);
        }
      }
    }
#endif
  }

#else
  
  if (makeboundaries)
  {
    // take every face, intersect with 3 other edges
    for (LINT i = 0; i < numFaces(); i++)
    {
      TPoint<T> c = centres[i];

      std::array<TPoint<T>,3> corners = threeCorners(i);

      // UVcorners must be filled in TBaseSurface::createTriangles()
      std::array<TPoint<T>,3> uvcorners = UVcorners[i];

      for (LINT j = 0; j < other.numFaces(); j++)
      {
        TPoint<T> oc = ocentres[j];
        T d = !(oc - c);
        if (d > oc.W + c.W)
          continue;

        std::array<TPoint<T>,3> ocorners = other.threeCorners(j);
        // UVcorners must be filled in TBaseSurface::createTriangles()
        std::array<TPoint<T>,3> ouvcorners = other.UVcorners[j];

        std::vector<TPoint<T>> intrs;
        if (intersectTriangleByTriangle(corners,ocorners,intrs,tolerance,parmtolerance))
        {
          if (intrs.size() == 2)
          {
            TPoint<T> coord0 = barycentricCoord(corners,intrs[0]);
            TPoint<T> coord1 = barycentricCoord(corners,intrs[1]);
            TPoint<T> ocoord0 = barycentricCoord(ocorners,intrs[0]);
            TPoint<T> ocoord1 = barycentricCoord(ocorners,intrs[1]);

            // these will contain U,V in X,Y for every intersection point
            TPoint<T> uv0,uv1,ouv0,ouv1;
            for (int k = 0; k < 3; k++)
            {
              uv0 += uvcorners[k] * coord0.XYZW[k];
              uv1 += uvcorners[k] * coord1.XYZW[k];
              ouv0 += ouvcorners[k] * ocoord0.XYZW[k];
              ouv1 += ouvcorners[k] * ocoord1.XYZW[k];
            }

            // keep U,V for every intersection on both surfaces
            UVintrs.push_back(TPoint<T>(uv0.X,uv0.Y,ouv0.X,ouv0.Y));
            intrs[0].W = T(UVintrs.size() - 1);
            UVintrs.push_back(TPoint<T>(uv1.X,uv1.Y,ouv1.X,ouv1.Y));
            intrs[1].W = T(UVintrs.size() - 1);

            edges.push_back(std::pair<TPoint<T>,TPoint<T>>(intrs[0],intrs[1]));
          }
        }
      }
    }
  } else
  {
    // take every face, intersect with 3 other edges
    for (LINT i = 0; i < numFaces(); i++)
    {
      TPoint<T> c = centres[i];

      std::array<TPoint<T>,3> corners = threeCorners(i);

      for (LINT j = 0; j < other.numFaces(); j++)
      {
        TPoint<T> oc = ocentres[j];
        T d = !(oc - c);
        if (d > oc.W + c.W)
          continue;

        std::array<TPoint<T>,3> ocorners = other.threeCorners(j);

        std::vector<TPoint<T>> intrs;
        if (intersectTriangleByTriangle(corners,ocorners,intrs,tolerance,parmtolerance))
        {
          if (intrs.size() == 2)
          {
            edges.push_back(std::pair<TPoint<T>,TPoint<T>>(intrs[0],intrs[1]));
          }
        }
      }
    }  
  }

#endif

  // this tolerance must be big
  T bigtolerance = maxedge * 0.1; //!!!!!!

//outputDebugString(std::string("edges ") + to_string(edges.size()) +
//std::string(" bigtolerance ") + to_string(bigtolerance,18)); 

  bool ok = makeUpCurves(edges,bigtolerance,lines,true); 

#if 0
  bool print = false;

  if (print)
  {
    other.saveFaceSTL(19388,"19388.STL","TDCAD",false);
    saveFaceSTL(824,"824.STL","TDCAD",false);
    saveFaceSTL(826,"826.STL","TDCAD",false);
    saveFaceSTL(827,"827.STL","TDCAD",false);
    saveFaceSTL(828,"828.STL","TDCAD",false);
    saveFaceSTL(829,"829.STL","TDCAD",false);
    saveFaceSTL(834,"834.STL","TDCAD",false);

    for (int i = 0; i < int(activecells.size()); i++)
    {
      int cellindex = int(activecells[i]);

      IPosition ipos = cells.backCellIndexToPosition(cellindex);
      std::array<TPoint<T>,8> cellcorners = cells.cell8RealCoordinates(ipos);

      saveBoxSTL(cellcorners,"cell" + to_string(cellindex) + ".STL","TDCAD",false);
    }
  }
#endif

  if (makeboundaries)
  {
    for (int i = 0; i < int(lines.size()); i++)
    {
      boundary0->push_back(std::vector<TPoint<T>>());
      boundary1->push_back(std::vector<TPoint<T>>());
      for (int j = 0; j < int(lines[i].size()); j++)
      {
        int pointno = ROUND(lines[i][j].W);
        boundary0->back().push_back(TPoint<T>(UVintrs[pointno].X,UVintrs[pointno].Y));
        boundary1->back().push_back(TPoint<T>(UVintrs[pointno].Z,UVintrs[pointno].W));
      }
    }
  }

  return ok;
}

template <class T> bool TTriangles<T>::intersect(TTriangles<T> &other, 
  std::vector<std::array<TPoint<T>,8>> &boxes,
  std::vector<std::array<TPoint<T>,8>> &oboxes,
  std::vector<std::vector<TPoint<T>>> &lines, T parmtolerance,
  std::vector<std::vector<TPoint<T>>> *boundary0,
  std::vector<std::vector<TPoint<T>>> *boundary1,
  int numthreads)
{
  // no lines
  lines.clear();

  // intersection pairs of points
  std::vector<std::pair<TPoint<T>,TPoint<T>>> edges;

  // parametric values for intersection points, X,Y - U,V for first surface,
  // Z,W - U,V for second surface
  std::vector<TPoint<T>> UVintrs;

  // make parametric boundaries
  bool makeboundaries = boundary0 && boundary1 && !UVcorners.empty() && !other.UVcorners.empty();

  // for collision test
  std::vector<TPoint<T>> centres,ocentres;
  getCentresAndRadii(centres);
  other.getCentresAndRadii(ocentres);

  // max tri size for cell size
  T minedge,maxedge;
  getEdgeMinMax(minedge,maxedge);
  T ominedge,omaxedge;
  other.getEdgeMinMax(ominedge,omaxedge);

  maxedge = std::max(maxedge,omaxedge) * 1.1; //!!!!!!!

  // this tolerance below is only used to remove duplicate intersection points,
  // no other purpose
  T tolerance = maxedge * parmtolerance;

  // use octree background cells for spacial partitioning...
#ifdef USE_SPACEPARTITIONING 

  // min/max of all triangles
  std::pair<TPoint<T>,TPoint<T>> mm0 = minmax();
  std::pair<TPoint<T>,TPoint<T>> mm1 = other.minmax();
  TPoint<T> min = pointMin(mm0.first,mm1.first);
  TPoint<T> max = pointMax(mm0.second,mm1.second);

  // extend
  extendMinMax(min,max,1.01);
  TPoint<T> dmm = max - min;

  // no intersection possible
  if (dmm.X < tolerance || dmm.Y < tolerance || dmm.Z < tolerance)
    return false;

  // cells big enough for spacial partitioning to identify if two triangles
  // may intersect, they can if they have a node in the same cell only
  OBackground<T> cells(min,max,(LINT) (dmm.X / maxedge),(LINT) (dmm.Y / maxedge),(LINT) (dmm.Z / maxedge),8);

  // total number of cells
  LINT numcells = cells.numBackgroundCells();

  // now distribute all tris over cells
  std::vector<std::set<LINT>> celltris(numcells);

  for (LINT i = 0; i < numFaces(); i++)
  {
    // important, do not remove
    std::array<TPoint<T>,8> box = boxes[i];

    for (int j = 0; j < 8; j++)
    {
      LINT index = cells.findCell(box[j]);
      assert(index >= 0);

      celltris[index].insert(i);
    }
  }

  std::vector<std::set<LINT>> ocelltris(numcells);
  for (LINT i = 0; i < other.numFaces(); i++)
  {
    // important, do not remove
    std::array<TPoint<T>,8> box = oboxes[i];

    for (int j = 0; j < 8; j++)
    {
      LINT index = cells.findCell(box[j]);
      assert(index >= 0);

      ocelltris[index].insert(i);
    }
  }

  // now distinguish cells which are (1) not empty (2) contain triangles from
  // this and other triangles
  std::vector<LINT> activecells;
  for (LINT i = 0; i < numcells; i++)
  {
    if (!celltris[i].empty() && !ocelltris[i].empty())
    {
      activecells.push_back(i);
    }
  }
  int numactive = int(activecells.size());

  if (numactive == 0)
    return false;

  // now we've got two plain lists of triangles for every cell

#ifndef USE_THREADS
  numthreads = 1;
#endif

  if (numactive < numthreads)
    numthreads = numactive;
  LIMIT_MIN(numthreads,1);

#ifdef DEBUG_TRIS
  outputDebugString(
    std::string("num threads = ") + to_string(numthreads) +
    std::string(" num cells = ") + to_string(numcells) +
    std::string(" active cells = ") + to_string(numactive));
#endif

  if (numthreads <= 1)
  {
#ifdef DEBUG_TRIS
  outputDebugString(
    std::string(" active cell size = ") + 
    to_string(int(celltris[activecells[0]].size())) + " " +
    to_string(int(ocelltris[activecells[0]].size())));
#endif

    std::vector<std::vector<std::pair<TPoint<T>,TPoint<T>>>> tedges(numactive);
    std::vector<std::vector<TPoint<T>>> tUVintrs(numactive);
    std::vector<std::vector<std::pair<LINT,LINT>>> intersected(numactive);

    std::unordered_set<std::pair<LINT,LINT>,pair_hash> sintersected;

    if (makeboundaries)
    {
      for (LINT i = 0; i < numactive; i++)
      {
        LINT index = activecells[i];

        intersectCells(this,&other,index,&celltris,&ocelltris,&centres,&ocentres,
#ifdef GLOBAL_INTRCHECK
          nullptr,&UVintrs,&edges,tolerance,parmtolerance,&sintersected);
#else
          &intersected[i],&tUVintrs[i],&tedges[i],tolerance,parmtolerance,nullptr);
#endif
      }
    } else
    {
      for (LINT i = 0; i < numactive; i++)
      {
        LINT index = activecells[i];
        intersectCellsNoUV(this,&other,index,&celltris,&ocelltris,&centres,&ocentres,
#ifdef GLOBAL_INTRCHECK
          nullptr,&UVintrs,&edges,tolerance,parmtolerance,&sintersected);
#else
          &intersected[i],&tUVintrs[i],&tedges[i],tolerance,parmtolerance,nullptr);
#endif
      }
    }

#ifndef GLOBAL_INTRCHECK
    int count = 0;

    std::set<std::pair<LINT,LINT>,EdgeCompare> intrset;

    for (int t = 0; t < numactive; t++)
    {
      for (int k = 0; k < int(tedges[t].size()); k++)
      {
        // avoid inserting same intersected tris
        auto res = intrset.insert(intersected[t][k]);

        if (res.second)
        {
          if (makeboundaries)
          {
            UVintrs.push_back(tUVintrs[t][k * 2]);
            UVintrs.push_back(tUVintrs[t][k * 2 + 1]);
          }

          std::pair<TPoint<T>,TPoint<T>> edge = tedges[t][k];
          edge.first.W = T(count++);
          edge.second.W = T(count++);
          edges.push_back(edge);
        }
      }
    }
#endif

  } else
  {
    std::vector<std::vector<std::pair<TPoint<T>,TPoint<T>>>> tedges(numthreads);
    std::vector<std::vector<TPoint<T>>> tUVintrs(numthreads);
    std::vector<std::vector<std::pair<LINT,LINT>>> intersected(numthreads);

    std::unordered_set<std::pair<LINT,LINT>,pair_hash> sintersected;

    // distribute over threads
    std::vector<std::thread> threads;
    int dt = int(numactive) / numthreads;
    for (int t = 0; t < numthreads; t++)
    {
      int j1 = t * dt;
      int j2 = (t + 1) * dt - 1;

      if (t == numthreads - 1)
        j2 = int(numactive) - 1;

      std::thread th;
      if (makeboundaries)
      {
        th = std::thread(intersectMultiCells<T>,&activecells,j1,j2,this,&other,
          &celltris,&ocelltris,&centres,&ocentres,
#ifdef GLOBAL_INTRCHECK
          nullptr,&UVintrs,&edges,tolerance,parmtolerance,&sintersected);
#else
          &intersected[t],&tUVintrs[t],&tedges[t],tolerance,parmtolerance,nullptr);
#endif
      } else
      {
        th = std::thread(intersectMultiCellsNoUV<T>,&activecells,j1,j2,this,&other,
          &celltris,&ocelltris,&centres,&ocentres,
#ifdef GLOBAL_INTRCHECK          
          nullptr,&UVintrs,&edges,tolerance,parmtolerance,&sintersected);
#else
          &intersected[t],&tUVintrs[t],&tedges[t],tolerance,parmtolerance,nullptr);
#endif
      }

      threads.push_back(std::move(th));
    }

    for (auto &th : threads) {     
      th.join();
    }

#ifndef GLOBAL_INTRCHECK
    int count = 0;

    std::set<std::pair<LINT,LINT>,EdgeCompare> intrset;

    for (int t = 0; t < numthreads; t++)
    {
      for (int k = 0; k < int(tedges[t].size()); k++)
      {
        // avoid inserting same intersected tris
#ifndef GLOBAL_INTRCHECK
        auto res = intrset.insert(intersected[t][k]);

        if (res.second)
#endif
        {
          if (makeboundaries)
          {
            UVintrs.push_back(tUVintrs[t][k * 2]);
            UVintrs.push_back(tUVintrs[t][k * 2 + 1]);
          }

          std::pair<TPoint<T>,TPoint<T>> edge = tedges[t][k];
          edge.first.W = T(count++);
          edge.second.W = T(count++);
          edges.push_back(edge);
        }
      }
    }
#endif
  }

#else
  
  if (makeboundaries)
  {
    // take every face, intersect with 3 other edges
    for (LINT i = 0; i < numFaces(); i++)
    {
      TPoint<T> c = centres[i];

      std::array<TPoint<T>,3> corners = threeCorners(i);

      // UVcorners must be filled in TBaseSurface::createTriangles()
      std::array<TPoint<T>,3> uvcorners = UVcorners[i];

      for (LINT j = 0; j < other.numFaces(); j++)
      {
        TPoint<T> oc = ocentres[j];
        T d = !(oc - c);
        if (d > oc.W + c.W)
          continue;

        std::array<TPoint<T>,3> ocorners = other.threeCorners(j);
        // UVcorners must be filled in TBaseSurface::createTriangles()
        std::array<TPoint<T>,3> ouvcorners = other.UVcorners[j];

        std::vector<TPoint<T>> intrs;
        if (intersectTriangleByTriangle(corners,ocorners,intrs,tolerance,parmtolerance))
        {
          if (intrs.size() == 2)
          {
            TPoint<T> coord0 = barycentricCoord(corners,intrs[0]);
            TPoint<T> coord1 = barycentricCoord(corners,intrs[1]);
            TPoint<T> ocoord0 = barycentricCoord(ocorners,intrs[0]);
            TPoint<T> ocoord1 = barycentricCoord(ocorners,intrs[1]);

            // these will contain U,V in X,Y for every intersection point
            TPoint<T> uv0,uv1,ouv0,ouv1;
            for (int k = 0; k < 3; k++)
            {
              uv0 += uvcorners[k] * coord0.XYZW[k];
              uv1 += uvcorners[k] * coord1.XYZW[k];
              ouv0 += ouvcorners[k] * ocoord0.XYZW[k];
              ouv1 += ouvcorners[k] * ocoord1.XYZW[k];
            }

            // keep U,V for every intersection on both surfaces
            UVintrs.push_back(TPoint<T>(uv0.X,uv0.Y,ouv0.X,ouv0.Y));
            intrs[0].W = T(UVintrs.size() - 1);
            UVintrs.push_back(TPoint<T>(uv1.X,uv1.Y,ouv1.X,ouv1.Y));
            intrs[1].W = T(UVintrs.size() - 1);

            edges.push_back(std::pair<TPoint<T>,TPoint<T>>(intrs[0],intrs[1]));
          }
        }
      }
    }
  } else
  {
    // take every face, intersect with 3 other edges
    for (LINT i = 0; i < numFaces(); i++)
    {
      TPoint<T> c = centres[i];

      std::array<TPoint<T>,3> corners = threeCorners(i);

      for (LINT j = 0; j < other.numFaces(); j++)
      {
        TPoint<T> oc = ocentres[j];
        T d = !(oc - c);
        if (d > oc.W + c.W)
          continue;

        std::array<TPoint<T>,3> ocorners = other.threeCorners(j);

        std::vector<TPoint<T>> intrs;
        if (intersectTriangleByTriangle(corners,ocorners,intrs,tolerance,parmtolerance))
        {
          if (intrs.size() == 2)
          {
            edges.push_back(std::pair<TPoint<T>,TPoint<T>>(intrs[0],intrs[1]));
          }
        }
      }
    }  
  }

#endif

  // this tolerance must be big
  T bigtolerance = maxedge * 0.1; //!!!!!!

//outputDebugString(std::string("edges ") + to_string(edges.size()) +
//std::string(" bigtolerance ") + to_string(bigtolerance,18)); 

  bool ok = makeUpCurves(edges,bigtolerance,lines,true); 

#if 0
  bool print = false;

  if (print)
  {
    other.saveFaceSTL(19388,"19388.STL","TDCAD",false);
    saveFaceSTL(824,"824.STL","TDCAD",false);
    saveFaceSTL(826,"826.STL","TDCAD",false);
    saveFaceSTL(827,"827.STL","TDCAD",false);
    saveFaceSTL(828,"828.STL","TDCAD",false);
    saveFaceSTL(829,"829.STL","TDCAD",false);
    saveFaceSTL(834,"834.STL","TDCAD",false);

    for (int i = 0; i < int(activecells.size()); i++)
    {
      int cellindex = int(activecells[i]);

      IPosition ipos = cells.backCellIndexToPosition(cellindex);
      std::array<TPoint<T>,8> cellcorners = cells.cell8RealCoordinates(ipos);

      saveBoxSTL(cellcorners,"cell" + to_string(cellindex) + ".STL","TDCAD",false);
    }
  }
#endif

  if (makeboundaries)
  {
    for (int i = 0; i < int(lines.size()); i++)
    {
      boundary0->push_back(std::vector<TPoint<T>>());
      boundary1->push_back(std::vector<TPoint<T>>());
      for (int j = 0; j < int(lines[i].size()); j++)
      {
        int pointno = ROUND(lines[i][j].W);
        boundary0->back().push_back(TPoint<T>(UVintrs[pointno].X,UVintrs[pointno].Y));
        boundary1->back().push_back(TPoint<T>(UVintrs[pointno].Z,UVintrs[pointno].W));
      }
    }
  }

  return ok;
}

template <class T> TPoint<T> TTriangles<T>::getCentre() const
{
  TPoint<T> centre(0,0,0);

  if (numFaces())
  {
    for (int i = 0; i < numFaces(); i++)
    {
      centre += faceCentre(i);
    }

    centre /= T(numFaces());
  }

  return centre;
}

template <class T> T TTriangles<T>::faceArea(LINT faceNo) const
{
  std::array<TPoint<T>,3> c = threeCorners(faceNo);

  TPoint<T> v01 = c[1] - c[0];
  TPoint<T> v12 = c[2] - c[1];
                           
  T area = !(v01 ^ v12) * T(0.5);
  return area;
}

template <class T> T TTriangles<T>::area() const
{
  T value = 0.0;

  for (int i = 0; i < numFaces(); i++)
  {
    value += faceArea(i);
  }

  return value;
}

template <class T> TPoint<T> TTriangles<T>::centre() const
{
  TPoint<T> value = 0.0;
  T totalarea = 0.0;

  for (int i = 0; i < numFaces(); i++)
  {
    T area = faceArea(i);
    totalarea += area;

    value += faceCentre(i) * area;
  }

  if (totalarea > TOLERANCE(T))
  {
    value = value / totalarea;
  }

  return value;
}

template <class T> bool TTriangles<T>::getBoundaryNodes(std::vector<std::vector<LINT>> &ilines, T tolerance,
  LINT maxlinelength)
{
  // return value
  bool ok = false;

  // lines are cleared
  ilines.clear();

  // make unique node numbering
  if (!duplicatesExcluded())
  {
    buildConnectivityArray(tolerance); 
  }

  // triangles around edges
  std::map<std::pair<LINT,LINT>,std::vector<LINT>,EdgeCompare> edgeTris;
  getEdgeTris(edgeTris);

  // triangles around every corner
  std::vector<std::vector<LINT>> cornerTris;
  getCornerTris(cornerTris);

  // busy edges
  std::set<std::pair<LINT,LINT>,EdgeCompare> busy;

  // find starting hanging edge
  std::pair<LINT,LINT> startedge;
  if (!findFreeEdge(edgeTris,busy,startedge))
    return false;

  // make line of node numbers
  std::vector<LINT> iline;
  iline.push_back(startedge.first);
  iline.push_back(startedge.second);

  busy.insert(startedge);

  do {
    // take tris around node iline.back()
    LINT node0,node1;
    if (getBoundaryNodesAround(iline.back(),cornerTris,edgeTris,node0,node1))
    {
      if (node0 == iline[iline.size() - 2])
      {
        busy.insert(std::pair<LINT,LINT>(iline.back(),node1));
        if (iline.size() >= maxlinelength) 
        {
          ilines.push_back(iline);
          break;
        }
        iline.push_back(node1);
      } else
      {
        busy.insert(std::pair<LINT,LINT>(iline.back(),node0));
        if (iline.size() >= maxlinelength)
        {
          ilines.push_back(iline);
          break;
        }
        iline.push_back(node0); 
      }
    } else
    {
      // restart
      ilines.push_back(iline);

      if (findFreeEdge(edgeTris,busy,startedge))
      {
        iline.clear();

        iline.push_back(startedge.first);
        iline.push_back(startedge.second);
        busy.insert(startedge);
      } else
      {
        // ok not set
        break;
      }
    }

    if (iline.back() == iline.front())
    {
      ilines.push_back(iline);
      ok = true;
      break;
    }
  } while (true);

  return ok;
}

template <class T> bool TTriangles<T>::getBoundaryNodesAround(LINT node, std::vector<std::vector<LINT>> &cornerTris, 
  std::map<std::pair<LINT,LINT>,std::vector<LINT>,EdgeCompare> &edgeTris, 
  LINT &node0, LINT &node1)
{
  // tris around the node
  std::vector<LINT> tris = cornerTris[node];

  node0 = node1 = -1;
  for (int i = 0; i < int(tris.size()); i++)
  {
    std::pair<LINT,LINT> edges[3];
    getFaceEdges(tris[i],edges);

    bool done = false;

    for (int j = 0; j < 3; j++)
    {
      // hanging edge
      if (edgeTris[edges[j]].size() == 1 && (edges[j].first == node || edges[j].second == node))
      {
        if (node0 < 0)
        {
          node0 = (edges[j].first == node) ? edges[j].second : edges[j].first;
        } else
        {
          node1 = (edges[j].first == node) ? edges[j].second : edges[j].first;
          if (node1 != node0)
            return true;
        }
      }
    }
  }

  return (node0 >= 0) && (node1 >= 0) && (node1 != node0);
}

template <class T> void TTriangles<T>::indexLinesToNodes(std::vector<std::vector<LINT>> &ilines,
  std::vector<std::vector<TPoint<T>>> &lines)
{
  lines.clear();

  for (int i = 0; i < int(ilines.size()); i++)
  {
    lines.push_back(std::vector<TPoint<T>>());
    for (int j = 0; j < int(ilines[i].size()); j++)
    {
      lines.back().push_back(coords[ilines[i][j]]);
    }
  }
}

template <class T> void TTriangles<T>::deleteFace(LINT faceNo)
{
  corners.erase(corners.begin() + faceNo * 3,corners.begin() + (faceNo + 1) * 3);
}

template <class T> bool TTriangles<T>::checkSliver(bool fix)
{
  // test for correct node numbering
  bool ok = true;

  std::vector<int> badtris;
  for (int i = 0; i < corners.size(); i += 3)
  {
    if (corners[i] == corners[i + 1] || corners[i] == corners[i + 2] || corners[i + 1] == corners[i + 2])
    {
      badtris.push_back(i / 3);
      ok = false;
    }
  }

  if (fix && !badtris.empty())
  {
    for (int i = int(badtris.size()) - 1; i >= 0; i--)
    {
      deleteFace(badtris[i]);
    }
  }

  return ok;
}

/** Connect two-point pieces into a closed line. */
template <class T> bool TTriangles<T>::edgesIntoCurves(std::vector<std::pair<TPoint<T>,TPoint<T>>> &edges,
  std::vector<std::vector<TPoint<T>>> &lines, T tolerance)
{
  if (edges.empty())
    return false;

  bool ok = makeUpCurves(edges,tolerance,lines,true);

  return ok;
}

template <class T> bool TTriangles<T>::getBoundaryPoints(std::vector<TPoint<T>> &boundary, T tolerance)
{
  // lines are cleared
  boundary.clear();

  // make unique node numbering
  if (!duplicatesExcluded())
  {
    buildConnectivityArray(tolerance); 
  }

  // triangles around edges
  std::map<std::pair<LINT,LINT>,std::vector<LINT>,EdgeCompare> edgeTris;
  getEdgeTris(edgeTris);

  // free edges
  for (auto &e : edgeTris)
  {
    if (e.second.size() == 1)
    {
      TPoint<T> v0 = coords[e.first.first];
      TPoint<T> v1 = coords[e.first.second];

      boundary.push_back(v0);
      boundary.push_back(v1);
    }
  }

  return true;
}

template <class T> bool TTriangles<T>::getEdgeMinMax(T &min, T &max)
{
  if (numFaces() < 1)
    return false;

  min = std::numeric_limits<T>::max();
  max = -std::numeric_limits<T>::max();
  for (int i = 0; i < numFaces(); i++)
  {
    std::array<TPoint<T>,3> c = threeCorners(i);

    T e0 = !(c[1] - c[0]);
    T e1 = !(c[2] - c[1]);
    T e2 = !(c[0] - c[2]);

    min = std::min<T>(min,std::min<T>(e0,std::min<T>(e1,e2)));
    max = std::max<T>(max,std::max<T>(e0,std::max<T>(e1,e2)));
  }

  return true;
}

template <class T> bool TTriangles<T>::getEdgeMinMax(TPoint<T> &dmin, TPoint<T> &dmax)
{
  if (numFaces() < 1)
    return false;

  dmin.X = dmin.Y = dmin.Z = std::numeric_limits<T>::max();
  dmax.X = dmax.Y = dmax.Z = -std::numeric_limits<T>::max();
  for (int i = 0; i < numFaces(); i++)
  {
    std::array<TPoint<T>,3> c = threeCorners(i);

    TPoint<T> c01min = pointMin<T>(c[0],c[1]);
    TPoint<T> c01max = pointMax<T>(c[0],c[1]);
    TPoint<T> mmin = pointMin<T>(c01min,c[2]);
    TPoint<T> mmax = pointMax<T>(c01max,c[2]);
    TPoint<T> d = mmax - mmin;

    dmin = pointMin<T>(dmin,d);
    dmax = pointMax<T>(dmax,d);
  }

  return true;
}

template <class T> bool TTriangles<T>::getBoundary(std::vector<std::vector<TPoint<T>>> &boundary, T tolerance)
{
  // return value
  bool ok = false;

  // lines are cleared
  boundary.clear();

  // make unique node numbering
  if (!duplicatesExcluded())
  {
    buildConnectivityArray(tolerance); 
  }

  // triangles around edges
  std::map<std::pair<LINT,LINT>,std::vector<LINT>,EdgeCompare> edgeTris;
  getEdgeTris(edgeTris);

  // free edges
  std::vector<std::pair<TPoint<T>,TPoint<T>>> edges;
  for (auto &e : edgeTris)
  {
    std::vector<TPoint<T>> edge;
    TPoint<T> v0 = coords[e.first.first];
    TPoint<T> v1 = coords[e.first.second];

    if (e.second.size() == 1)
    {
      edges.push_back(std::pair<TPoint<T>,TPoint<T>>(v0,v1));
    }
  }

  ok = edgesIntoCurves(edges,boundary,tolerance);

  return ok;
}

template <class T> inline void TTriangles<T>::getFaceCorners(LINT faceNo, LINT &i0, LINT &i1, LINT &i2)
{
  i0 = corners[faceNo * 3 + 0];
  i1 = corners[faceNo * 3 + 1];
  i2 = corners[faceNo * 3 + 2];
}

template <class T> inline void TTriangles<T>::getFaceEdges(LINT faceNo, std::pair<LINT,LINT> edges[3])
{
  LINT i0,i1,i2;
  getFaceCorners(faceNo,i0,i1,i2);
  edges[0].first = i0;
  edges[0].second = i1;
  edges[1].first = i1;
  edges[1].second = i2;
  edges[2].first = i2;
  edges[2].second = i0;
}

template <class T> void TTriangles<T>::getFaceEdges(LINT faceNo, std::pair<TPoint<T>,TPoint<T>> pedges[3])
{
  std::pair<LINT,LINT> edges[3];
  getFaceEdges(faceNo,edges);

  for (int i = 0; i < 3; i++)
  {
    pedges[i].first = coords[edges[i].first];
    pedges[i].second = coords[edges[i].second];
  }
}

template <class T> int TTriangles<T>::cutFaceByPlane(TPlane<T> &plane, LINT faceNo, 
  std::vector<TPoint<T>> &intrs, T tolerance)
{
  intrs.clear();

  std::pair<TPoint<T>,TPoint<T>> pedges[3];
  getFaceEdges(faceNo,pedges);

  for (int i = 0; i < 3; i++)
  {
    TPoint<T> p0 = pedges[i].first;
    TPoint<T> p1 = pedges[i].second;

    TPoint<T> I;
    T U = 0.5;
    if (plane.segmentIntersect(p0,p1,&I,&U,tolerance)) 
    {
      intrs.push_back(I);
    }
  }

  removeDuplicates(intrs,true,tolerance); // true is correct here

  // important : set correct orientation of intersection curve : part to the right is
  // thrown off
  if (intrs.size() == 2)
  {
    TPoint<T> normal = faceNormal(faceNo);
    TPoint<T> d = intrs[1] - intrs[0];
    TPoint<T> cross = d ^ normal;

    if (!(cross > plane.normal))
    {
      std::reverse(intrs.begin(),intrs.end());
    }
  }

  return int(intrs.size());
}

template <class T> void TTriangles<T>::getCentresAndRadii(std::vector<TPoint<T>> &centres, T Rcoef)
{
  for (LINT i = 0; i < numFaces(); i++)
  {
    T R = 0.0;
    TPoint<T> centre = faceCentre(i,&R);
    centre.W = R * Rcoef;

    centres.push_back(centre);
  }
}

/** Delete triangles. */
template <class T> void deleteTriangles(std::vector<TTriangles<T> *> &triangles)
{
  for (int i = 0; i < int(triangles.size()); i++)
  {
    DELETE_CLASS(triangles[i]);
  }

  triangles.clear();
}

}
