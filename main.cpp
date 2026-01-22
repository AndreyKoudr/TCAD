
/*******************************************************************************

  Templated 3D CAD

  All the code is based on STL (standard templates, not stereolithography);
  it is templated what means that it is easy to use, like std::vector, you do 
  not need dlls, all is in front of you, no first-, second- or third-party libs,
  easy for debugging and well commented.

*******************************************************************************/

// you can only use high-level operations
#include "toperations.h"
#include "tblocks.h"

#include "ttriangles.h"
#include "tbasesurface.h"
#include "tbezierpatch.h" 
#include "tbeziersurface.h"
#include "tpointsurface.h" 
#include "tsplinesurface.h" 
#include "tbeziervolume.h" 
#include "tsplinevolume.h" 
#include "tdata.h" 
#include "FFD.h" 
#include "tbrep.h" 

// this stuff is for export and debugging
#include "strings.h"
#include "export.h"
#include "import.h"

#include <iostream>

using namespace tcad;

//==============================================================================
// SET WHAT TO DEBUG (below)
//==============================================================================
#define DEBUG_COMMON 
#define DEBUG_SUBMARINE
#define DEBUG_BOOLEANS

// save debugging files here, end this with slash "/"
#define DEBUG_DIR std::string("C:/AndrewK/MyProjects2/temp/") //!!!!!!!

// Let us choose double as our basic real type
#define T double

//===== Random generator for running examples ==================================

#include <time.h>

// random generator started? 
bool RandGenStarted = false;

// start random generator
void randomize()
{
  if (!RandGenStarted)
  {
    srand((unsigned) time(nullptr));
    RandGenStarted = true;
  }
}

// get random value in range 0..range
int random(int range)
{
  // start random generator if not already started
  if (!RandGenStarted) randomize();

  // return random value
  if (range == RAND_MAX) return rand();
  else return rand() * range / RAND_MAX;
}

// get double random value in range 0..1
double random()
{
  // start random generator if not already started
  if (!RandGenStarted) randomize();

  // return random value
  return ((double) rand()) / ((double) RAND_MAX);
}

/** Make random neighbour duplicates in points. */
void makeRandomNeighbourDuplicates(std::vector<TPoint<T>> &points, std::vector<TPoint<T>> &spoiltpoints, 
  T noisesize)
{
  spoiltpoints = points;

  // insert random duplicates
  int nduplicates = int(points.size());
  for (int i = 0; i < nduplicates; i++)
  {
    // random node, we duplicate it
    int index = random(int(spoiltpoints.size()) - 1);
    TPoint<T> p = spoiltpoints[index];

    // spoil p
    p.X += (random() - 0.5) * noisesize;
    p.Y += (random() - 0.5) * noisesize;
    p.Z += (random() - 0.5) * noisesize;

    spoiltpoints.insert(spoiltpoints.begin() + index,p);
  }
}

/** Make random noise in point coordinates. */
void makeRandomNoise(std::vector<TPoint<T>> &points, std::vector<TPoint<T>> &spoiltpoints, 
  T noisesize)
{
  spoiltpoints = points;

  // insert random duplicates
  for (int i = 0; i < int(spoiltpoints.size()); i++)
  {
    TPoint<T> p = spoiltpoints[i];

    // spoil p
    p.X += (random() - 0.5) * noisesize;
    p.Y += (random() - 0.5) * noisesize;
    p.Z += (random() - 0.5) * noisesize;

    spoiltpoints[i] = p;
  }
}

//===== High-resolution timer ==================================================

bool highrestimer = false;
LARGE_INTEGER timerfrequency = { 1000 };

/** Start mesurement, check if high-res timer is supported. */
void TestTimer()
{
  // which counter to use?
  if (QueryPerformanceFrequency(&timerfrequency))
  {
    highrestimer = true;
  } else
  {
    timerfrequency.QuadPart = 1000;
    highrestimer = false;
  }
}

// get time in seconds
double GetTime()
{
  if (highrestimer)
  {
    LARGE_INTEGER count;
    QueryPerformanceCounter(&count);
    return (double) count.QuadPart / (double) timerfrequency.QuadPart;
  } else
  { 
    return (double) (GetTickCount()) / 1000.0;
  }
}

/** Randomly swap points. */
void makeRandomSwap(int numswaps, std::vector<TPoint<T>> &points, std::vector<TPoint<T>> &spoiltpoints)
{
  spoiltpoints = points;

  // randomly swap
  for (int i = 0; i < numswaps; i++)
  {
    int index0 = random(int(spoiltpoints.size()) - 1);
    int index1 = random(int(spoiltpoints.size()) - 1);

    TPoint<T> temp = spoiltpoints[index0];
    spoiltpoints[index0] = spoiltpoints[index1];
    spoiltpoints[index1] = temp;
  }
}

/** STL file test.
  (1) load name.stl
  (2) is it manifold and solid?
  (3) output boundary (is there a boundary?)
  (4) cut by plane and output line
*/
bool checkTopoCutAndBoundary(const std::string &name, TPlane<T> &plane,
  bool mustbemanifold, bool mustbesolid, T sharpangledeg)
{
  std::string partname;
  bool binary;

  // load file
  TTriangles<T> tris;
  if (!loadTrianglesStl(tris,name + ".stl",partname,binary,0.0))
    return false;

  // get model size to set tolerance
  std::pair<TPoint<T>,TPoint<T>> minmax = tris.minmax();
  T size = !(minmax.second - minmax.first);
  T tolerance = size * PARM_TOLERANCE;

  // is manifold?
  std::vector<std::pair<LINT,LINT>> badedges;
  bool manifold = tris.manifold(tolerance,badedges);
  ASSERT(manifold == mustbemanifold);

  // is solid?
  bool solid = tris.solid(tolerance,badedges);
  ASSERT(solid == mustbesolid);

  // get boundary
  std::vector<std::vector<TPoint<T>>> boundary;
  bool bok = tris.getBoundary(boundary,tolerance);

  // output boundary, the should be none if solid
  if (solid)
  {
    ASSERT(!bok);
  } else
  {
    ASSERT(bok);

    // redivide point curves by sharp corners to make it look as in the 
    // original STL
    std::vector<std::vector<TPoint<T>>> bpoints;
    redividePoints(boundary,bpoints,tolerance,sharpangledeg);

    saveLinesIges(bpoints,DEBUG_DIR + name + " boundary.iges");
  }

  // cut it by plane
  std::vector<std::vector<TPoint<T>>> cutpoints;
  tris.intersectByPlane(plane,cutpoints,tolerance); 

  // redivide point curves by sharp corners to make it look good
  std::vector<std::vector<TPoint<T>>> cpoints;
  redividePoints(cutpoints,cpoints,tolerance);

  saveLinesIges(cpoints,DEBUG_DIR + name + " cut curve.iges");

  return true;
}

/** Rewrite an STL file as binary. */
bool rewriteSTLAsBinary(const std::string &filename)
{
  TTriangles<T> tris;

  std::string partname;
  bool binary;

  if (loadTrianglesStl(tris,filename,partname,binary,0.0))
  {
    return saveTrianglesStl(tris,filename,"TCAD",true);
  } else
  {
    return false;
  }
}

//===== TCAD examples ==========================================================

int main(int argc, char* argv[])
{
  // start console output (for errorMessage("..."))
  startConsole();

  // start time measurement
  TestTimer();
  double starttime = GetTime();

#ifdef DEBUG_COMMON

  /*****************************************************************************

    Part 1 : CURVES. What can we do with them?

  *****************************************************************************/

  cout << "    Part 1 : CURVES" << endl;

  // Let us take points of NACA0012 airfoil, they are in XY but of course
  // all your curves are 3-dimensional

  const std::vector<std::pair<T,T> > NACAxy = {
    {0.000000, 0.000000},
    {0.000685, 0.004611},
    {0.002739, 0.009114},
    {0.006156, 0.013503},
    {0.010926, 0.017770},
    {0.017037, 0.021904},
    {0.024472, 0.025893},
    {0.033210, 0.029726},
    {0.043227, 0.033389},
    {0.054497, 0.036867},
    {0.066987, 0.040145},
    {0.080665, 0.043211},
    {0.095492, 0.046049},
    {0.111427, 0.048648},
    {0.128428, 0.050996},
    {0.146447, 0.053083},
    {0.165435, 0.054902},
    {0.185340, 0.056447},
    {0.206107, 0.057714},
    {0.227680, 0.058702},
    {0.250000, 0.059412},
    {0.273005, 0.059848},
    {0.296632, 0.060015},
    {0.320816, 0.059921},
    {0.345492, 0.059575},
    {0.370590, 0.058989},
    {0.396044, 0.058175},
    {0.421783, 0.057148},
    {0.447736, 0.055923},
    {0.473832, 0.054515},
    {0.500000, 0.052940},
    {0.526168, 0.051216},
    {0.552264, 0.049358},
    {0.578217, 0.047383},
    {0.603956, 0.045307},
    {0.629410, 0.043147},
    {0.654508, 0.040917},
    {0.679184, 0.038634},
    {0.703368, 0.036311},
    {0.726995, 0.033962},
    {0.750000, 0.031603},
    {0.772320, 0.029246},
    {0.793893, 0.026905},
    {0.814660, 0.024593},
    {0.834565, 0.022323},
    {0.853553, 0.020107},
    {0.871572, 0.017959},
    {0.888573, 0.015891},
    {0.904508, 0.013914},
    {0.919335, 0.012042},
    {0.933013, 0.010286},
    {0.945503, 0.008658},
    {0.956773, 0.007168},
    {0.966790, 0.005826},
    {0.975528, 0.004642},
    {0.982963, 0.003626},
    {0.989074, 0.002783},
    {0.993844, 0.002120},
    {0.997261, 0.001644},
    {0.999315, 0.001356},
    {1.000000, 0.001260}
  };

  // TPoint<T> is a template for 3D point with all math operators like +,-,dot 
  // and cross products... defined see (tpoint.h). Let us convert them into TPoint<double>
  #define Point TPoint<T>

  std::vector<Point> points;
  for (auto p : NACAxy)
  {
    // Z of points is zero
    points.push_back(Point(p.first,p.second));
  }

  saveLinesIges(points,DEBUG_DIR + "NACApoints.iges");

  // make point curve of original points
  TPointCurve<T> curve(points);

  //TPoint<T> pos1 = curve.derivative(0.53,0);
  //curve.reverse();
  //TPoint<T> pos2 = curve.derivative(0.47,0);
  //T d = !(pos2 - pos1);

  saveCurveIges(curve,DEBUG_DIR + "NACA.iges");

  /*****************************************************************************
    1.1 Curves : how to exclude neighbour duplicate nodes
  *****************************************************************************/

  cout << "1.1 Curves : how to exclude neighbour duplicate nodes" << endl;

  // make a copy
  std::vector<Point> points1 = points;

  // calculate curve length and curvetolerance to make node duplicates
  T len = calculateLength(points);
  T curvetolerance = len * 0.00001;

  // spoil points by inserting node duplicates
  makeRandomNeighbourDuplicates(points,points1,curvetolerance);

  // do not sort coordinates (false), exclude only neighbours
  bool ok = removeDuplicates(points1,false,curvetolerance * 10.0);

  ASSERT(ok);

  TPointCurve<T> curve1(points1);

  // save to compare two curves
  saveTwoCurvesIges(curve,curve1,DEBUG_DIR + "Compare_removed_duplicates.iges");

  /*****************************************************************************
    1.2 Curves : compare original points and low-power LSQ segment, do not expect
      a good correspondence - the lsq is a simple poly of power 3
  *****************************************************************************/

  cout << "1.2 Curves : compare original points and low-power LSQ segment, do not expect a good correspondence - the lsq is a simple poly of power 3" << endl;

  TLSQSegment<T> lsqsegment(points,3);
  // save to compare two curves
  saveTwoCurvesIges(curve,lsqsegment,DEBUG_DIR + "Compare_points_and_lsq_segment.iges");

  /*****************************************************************************
    1.3 Curves : compare points and Bezier segment, do not expect
      a good correspondence - Bezier segment is qubic
  *****************************************************************************/

  cout << "1.3 Curves : compare points and Bezier segment, do not expect a good correspondence - Bezier segment is qubic" << endl;

  TBezierSegment<T> beziersegment(points,END_FIXED,END_FIXED);
  // save to compare two curves
  saveTwoCurvesIges(curve,beziersegment,DEBUG_DIR + "Compare_points_and_bezier_segment.iges");

  /*****************************************************************************
    1.4 Curves : compare points and segment of orthogonal polynomials of
      degree 4. Starting END_ROUNDED end is set to handle the round LE
  *****************************************************************************/

  cout << "1.4 Curves : compare points and segment of orthogonal polynomials of degree 4. Starting END_ROUNDED end is set to handle the round LE" << endl;

  TOrthoSegment<T> orthosegment(points,END_ROUNDED,END_FREE,8,GAUSSINT_20);
  // save to compare two curves
  saveTwoCurvesIges(curve,orthosegment,DEBUG_DIR + "Compare_points_and_orthopoly_segment.iges");

  /*****************************************************************************
    1.5 Curves : compare points and Bezier curve (collection of Bezier segments)
  *****************************************************************************/

  cout << "1.5 Curves : compare points and Bezier curve (collection of Bezier segments)" << endl;

  TBezierCurve<T> beziercurve(points,10,END_FIXED,END_FIXED);
  // save to compare two curves
  saveTwoCurvesIges(curve,beziercurve,DEBUG_DIR + "Compare_points_and_bezier_curve.iges");

  /*****************************************************************************
    1.6 Curves : compare points and approximated b-spline
  *****************************************************************************/

  cout << "1.6 Curves : compare points and approximated b-spline" << endl;

  TSplineCurve<T> splinecurve(points,int(points.size() - 1),SPLINE_DEGREE,END_CLAMPED,END_CLAMPED); 
  // save to compare two curves
  saveTwoCurvesIges(curve,splinecurve,DEBUG_DIR + "Compare_points_and_spline_curve.iges");

  /*****************************************************************************
    1.7 Curves : compare points and approximated b-spline with less number 
    of control points (10 here)
  *****************************************************************************/

  cout << "1.7 Curves : compare points and approximated b-spline with less number of control points (10 here)" << endl;

  TSplineCurve<T> splinecurve10(points,10,SPLINE_DEGREE,END_CLAMPED,END_CLAMPED); 
  // save to compare two curves
  saveTwoCurvesIges(curve,splinecurve10,DEBUG_DIR + "Compare_points_and_spline_curve_K10.iges");

  /*****************************************************************************
    1.8 Curves : how to smooth a curve (curve2 below) with various approximants
  *****************************************************************************/

  cout << "1.8 Curves : how to smooth a curve (curve2 below) with variuos approximants" << endl;

  // make random noise, points2 is a spoilt set
  std::vector<TPoint<T>> points2;
  makeRandomNoise(points,points2,len * 0.005);

  // do not spoil points near start/end : apply "hat" function to diminish
  // the changes near ends with maximum at the middle = 0.5
  // and set poly power as 0.5
  std::vector<T> coefs;
  makeHat(int(points2.size()),coefs,0.5,0.5);

  points2 = points + (points2 - points) * coefs;

  //===== smooth by orthogonal polynomials =====================================
  // these will be smoothed points
  std::vector<TPoint<T>> points3 = points2;
  smoothPointsByOrtho(points3,END_ROUNDED,END_FREE);

  TPointCurve<T> curve2(points2);
  TPointCurve<T> curve3(points3);

  // compare curves : spoilt and smoothed
  saveTwoCurvesIges(curve2,curve3,DEBUG_DIR + "Compare_points_orthosmooth.iges");

  //===== smooth by single Bezier segment ======================================
  // these will be smoothed points
  std::vector<TPoint<T>> points4 = points2;
  smoothPointsByBezier(points4,END_FIXED,END_FIXED);

  TPointCurve<T> curve4(points4);

  // compare curves : spoilt and smoothed
  saveTwoCurvesIges(curve2,curve4,DEBUG_DIR + "Compare_points_beziersegmentsmooth.iges");

  //===== smooth by Bezier curve ===============================================
  // these will be smoothed points
  std::vector<TPoint<T>> points5 = points2;
  smoothPointsByBezierCurve(points5,4,END_CLAMPED,END_CLAMPED);

  TPointCurve<T> curve5(points5);

  // compare curves : spoilt and smoothed
  saveTwoCurvesIges(curve2,curve5,DEBUG_DIR + "Compare_points_beziercurvesmooth.iges");

  //===== smooth by spline curve ===============================================
  // these will be smoothed points
  std::vector<TPoint<T>> points6 = points2;
  smoothPointsBySplineCurve(points6,10,END_CLAMPED,END_FREE);

  TPointCurve<T> curve6(points6);

  // compare curves : spoilt and smoothed
  saveTwoCurvesIges(curve2,curve6,DEBUG_DIR + "Compare_points_splinecurvesmooth.iges");

  /*****************************************************************************
    1.9 Curves : how to diminish effect of smoothing near curve ends (much 
    needed when handling surface control points)
  *****************************************************************************/

  cout << "1.9 Curves : how to diminish effect of smoothing near curve ends (much needed when handling surface control points)" << endl;

  // do not spoil points near start/end : apply "hat" function to diminish
  // the changes near ends with maximum at the middle = 0.5
  // and set poly power as 0.5
  std::vector<T> coefs6;
  makeHat(int(points6.size()),coefs6,0.5,0.5);

  // diminish effect between points2 (spoilt) and the smooth points6 
  std::vector<TPoint<T>> points6mod = points6;

  points6mod = points2 + (points6 - points2) * coefs6;

  TPointCurve<T> curve6mod(points6mod);

  // compare curves : spoilt and smoothed
  saveTwoCurvesIges(curve2,curve6mod,DEBUG_DIR + "Compare_points_splinecurvesmooth_mitigated.iges");

  /*****************************************************************************
    1.10 Curves : how to find parameter U for a point on curve
  *****************************************************************************/

  cout << "1.10 Curves : how to find parameter U for a point on curve" << endl;

  // take a random point
  int index = random(int(points.size()));
  // we want a parameter for this point
  TPoint<T> p = points[index];

  // take spline curve of 10 intervals
  T U = splinecurve.findUforPoint(p);
  // check this value by comparison with calculated position on the curve
  TPoint<T> pactual = splinecurve.derivative(U,0);
  T diff = !(pactual - p);

  ASSERT(diff < 0.005);

  // take Bezier curve
  U = beziercurve.findUforPoint(p);
  // check this value by comparison with calculated position on the curve
  pactual = beziercurve.derivative(U,0);
  diff = !(pactual - p);

  ASSERT(diff < 0.005);

  // take ortho segment
  U = orthosegment.findUforPoint(p);
  // check this value by comparison with calculated position on the curve
  pactual = orthosegment.derivative(U,0);
  diff = !(pactual - p);

 //!!!!!!! ASSERT(diff < 0.005);

  /*****************************************************************************
    1.11 Curves : how to intersect by plane
  *****************************************************************************/

  cout << "1.11 Curves : how to intersect by plane" << endl;

  // horizontal plane at Y = 0.03
  TPlane<T> plane(TPoint<T>(0.0,0.03,0.0),TPoint<T>(1.0,0.03,0.0),TPoint<T>(1.0,0.03,1.0),ok);

  // two intersections are expected
  std::vector<T> Upoints;
  int numintrs = splinecurve.intersectByPlane(plane,Upoints,curvetolerance);

  ASSERT(numintrs == 2);

  // find intersection points
  std::vector<TPoint<T>> intrs;
  for (int i = 0; i < numintrs; i++)
  {
    intrs.push_back(splinecurve.derivative(Upoints[i],0));
  }

  // are they on plane?
  for (int i = 0; i < int(intrs.size()); i++)
  {
    T dist = std::abs(plane.distance(intrs[i]));
    ASSERT(dist < curvetolerance);
  }

  /*****************************************************************************
    1.12 Curves : how to intersect two curves
  *****************************************************************************/

  cout << "1.12 Curves : how to intersect two curves" << endl;

  // make another curve : ellipse to intersect NACA points
  std::vector<TPoint<T>> epoints;
  makeEllipseXY(100,TPoint<T>(0.5,0.053),0.5,0.053,epoints);

  TSplineCurve<T> ecurve(epoints,SPLINE_DEGREE);

  // compare curves : spoilt and smoothed
  saveTwoCurvesIges(curve,ecurve,DEBUG_DIR + "Intersecting points with ellipse.iges");

  // these are pairs of parameters in X,Y for the two curves for
  // every intersection point
  std::vector<TPoint<T>> UV;

  // two intersections are expected
  numintrs = curve.findIntersections(ecurve,UV,curvetolerance); 

  ASSERT(numintrs == 2);

  // check coincidence
  for (int i = 0; i < int(UV.size()); i++)
  {
    TPoint<T> p0 = curve.derivative(UV[i].X,0);
    TPoint<T> p1 = ecurve.derivative(UV[i].Y,0);
    T dist = !(p1 - p0);
    ASSERT(dist < curvetolerance);
  }

  /*****************************************************************************
    1.13 Curves : cut a piece of curve from the previous example
  *****************************************************************************/

  cout << "1.13 Curves : cut a piece of curve from the previous example" << endl;

  // curve was cut from UV[0].X to UV[1].X
  std::vector<TPoint<T>> cutpoints;
  curve.cutPiece(51,UV[0].X,UV[1].X,cutpoints);

  // we order points here
  TSplineCurve<T> cutcurve(cutpoints,40,SPLINE_DEGREE);

  // compare curves : unordered and ordered
  saveCurveIges(cutcurve,DEBUG_DIR + "Piece of curve.iges");

  /*****************************************************************************
    1.14 Curves : order unordered points
  *****************************************************************************/

  cout << "1.14 Curves : order unordered points" << endl;

  std::vector<TPoint<T>> unorderedpoints;
  makeRandomSwap(int(points.size()) / 2,points,unorderedpoints);

  // unordered curve
  TPointCurve<T> unorderedcurve(unorderedpoints);

  // we order points here
  TPointCurve<T> orderedcurve(unorderedpoints);
  orderedcurve.order(curvetolerance);

  // compare curves : unordered and ordered
  saveTwoCurvesIges(unorderedcurve,orderedcurve,DEBUG_DIR + "Unordered and ordered curves.iges");

  /*****************************************************************************

    Part 2 : SURFACES

  *****************************************************************************/

  cout << "    Part 2 : SURFACES" << endl;

  /*****************************************************************************
    2.1 Surfaces : triangles : save in STL
  *****************************************************************************/

  cout << "2.1 Surfaces : triangles : save in STL" << endl;

  TTriangles<T> NACAtris;
  NACAtris.makeNACA0012(50,10,2.0);

  // tolerance for surfaces
  T NACAsize = NACAtris.maxSize();
  //T NACAtolerance = curvetolerance;
  T NACAtolerance = NACAsize * PARM_TOLERANCE;

  saveTrianglesStl(NACAtris,DEBUG_DIR + "NACA.stl");

  /*****************************************************************************
    2.2 Surfaces : triangles : manifold? solid?
  *****************************************************************************/

  cout << "2.2 Surfaces : triangles : manifold? solid?" << endl;

  // is manifold?
  std::vector<std::pair<LINT,LINT>> badedges;
  bool manifold = NACAtris.manifold(NACAtolerance,badedges);
  ASSERT(manifold);

  // is solid?
  bool solid = NACAtris.solid(NACAtolerance,badedges);
  ASSERT(solid);

  /*****************************************************************************
    2.3 Surfaces : triangles : cut by plane
  *****************************************************************************/

  cout << "2.3 Surfaces : triangles : cut by plane" << endl;

  // cut triangles by a plane (made by three points)
  TPlane<T> cutplane(TPoint<T>(-0.5,0.0,0.9),TPoint<T>(0.5,0.0,-0.7),TPoint<T>(0.5,1.0,-0.7),ok); 
  std::vector<std::vector<TPoint<T>>> NACAcuts;
  NACAtris.intersectByPlane(cutplane,NACAcuts,NACAtolerance); 

  /*****************************************************************************
    2.4 Surfaces : triangles : how to find a sharpest point in intersection curve
  *****************************************************************************/

  cout << "2.4 Surfaces : triangles : how to find a sharpest point in intersection curve" << endl;

  if (NACAcuts.size() == 1)
  {
    std::vector<TPoint<T>> NACAcutpoints = NACAcuts[0];
  
    // the curve (NACAcutpoints) starts from an arbitrary point after intersection, 
    // let's make it start from a point with highest curvature like TE in airfoils
    TPointCurve<T> temp(NACAcutpoints);
    int sharpindex = temp.findPointOfMaxCurvature();
    temp.shiftClosed(sharpindex,NACAtolerance);
    NACAcutpoints = temp.controlPoints();

    // make, for example, a spline curve on these points
    TSplineCurve<T> NACAcutcurve(NACAcutpoints,200,SPLINE_DEGREE,END_CLAMPED,END_CLAMPED);

    saveCurveIges(NACAcutcurve,DEBUG_DIR + "NACA cut curve.iges");
  } else
  {
    ASSERT(false);
  }

  /*****************************************************************************
    2.5 Surfaces : triangles : how to get triangulation boundary
  *****************************************************************************/

  cout << "2.5 Surfaces : triangles : how to get triangulation boundary" << endl;

  // NACA case is solid, no boundary points are expected
  std::vector<std::vector<TPoint<T>>> NACAtrisboundary;
  bool bok = NACAtris.getBoundary(NACAtrisboundary,NACAtolerance);

  ASSERT(!bok);

  // now take more complicated boundary
  std::string partname;
  bool binary;

  TTriangles<T> fuselage;
  loadTrianglesStl(fuselage,"fuselage.stl",partname,binary,0.0);

  std::pair<TPoint<T>,TPoint<T>> fminmax = fuselage.minmax();
  T fsize = !(fminmax.second - fminmax.first);
  T ftolerance = fsize * PARM_TOLERANCE;

  std::vector<std::vector<TPoint<T>>> fboundary;
  bool fok = fuselage.getBoundary(fboundary,ftolerance);

  ASSERT(fok);

  // redivide point curves by sharp corners to make it look good
  std::vector<std::vector<TPoint<T>>> fnewpoints;
  redividePoints(fboundary,fnewpoints,ftolerance,10.0);

  saveLinesIges(fnewpoints,DEBUG_DIR + "fuselage boundary.iges");

  /*****************************************************************************
    2.6 Surfaces : triangles : more complicated cases, load STL, check if 
    manifold and solid, generate boundary and cut by plane
  *****************************************************************************/

  cout << "2.6 Surfaces : triangles : more complicated cases, load STL, check if manifold and solid, generate boundary and cut by plane" << endl;

  // cut it by plane, define by normal and one point on plane
  TPlane<T> scutplane(TPoint<T>(1.0,0.0,0.0),TPoint<T>(5.8,0.0,0.0)); 

  bool sok1 = checkTopoCutAndBoundary("shuttle",scutplane,true,true,10.0);
  ASSERT(sok1);

  // already binary
  //rewriteSTLAsBinary("wing.stl");

  // cut it by plane, define by normal and one point on plane
  TPlane<T> wcutplane(TPoint<T>(0.0,0.0,1.0),TPoint<T>(0.0,0.0,0.0)); 

  bool wok1 = checkTopoCutAndBoundary("wing",wcutplane,true,true,45.0);
  //!!!!!!! ASSERT(wok1);

  /*****************************************************************************
    2.7 Surfaces : triangles : find intersection of one set of triangles with 
    another
  *****************************************************************************/

  cout << "2.7 Surfaces : triangles : find intersection of one set of triangles with another" << endl;

  // load another fuselage
  TTriangles<T> fuselage1;
  loadTrianglesStl(fuselage1,"fuselage.stl",partname,binary,0.0);

  // rotate the second one around X
  TTransform<T> t;
  t.Rotate(TPoint<T>(1.0,0.0,0.0),90.0 * CPI);
  fuselage1.makeTransform(t);

  // save tranformed
  saveTrianglesStl(fuselage1,DEBUG_DIR + "fuselage rotated.stl");

  std::vector<std::vector<TPoint<T>>> fcutpoints;
  bool iok = fuselage.intersect(fuselage1,true,fcutpoints,PARM_TOLERANCE); 

  ASSERT(iok);

  // redivide point curves by sharp corners to make it look good
  std::vector<std::vector<TPoint<T>>> fcpoints;
  // only 1 degree for sharp edges is here to divide the curve into
  // straight-line segments
  redividePoints(fcutpoints,fcpoints,NACAtolerance,1.0);

  saveLinesIges(fcpoints,DEBUG_DIR + "fuselage-fuselage intersection curve.iges");

  /*****************************************************************************
    2.8 Surfaces : Bezier patch
  *****************************************************************************/

  cout << "2.8 Surfaces : Bezier patch" << endl;

  TBezierSegment<T> SU0(TPoint<T>(0,0,0),TPoint<T>(1,0.0,0.5),TPoint<T>(2,0,0),TPoint<T>(3,0,0));
  TBezierSegment<T> S1V(TPoint<T>(3,0,0),TPoint<T>(3,1,0),TPoint<T>(3,2,0.0),TPoint<T>(3,3,0));
  TBezierSegment<T> SU1(TPoint<T>(0,3,0),TPoint<T>(1,3,0),TPoint<T>(2,3,0),TPoint<T>(3,3,0));
  TBezierSegment<T> S0V(TPoint<T>(0,0,0),TPoint<T>(0,1,0),TPoint<T>(0,2,0),TPoint<T>(0,3,0));

  TBezierPatch<T> bpatch(&SU0,&S1V,&SU1,&S0V);

  TTriangles<T> btris;
  bpatch.createTriangles(btris,20,20); 

  // save tranformed
  saveTrianglesStl(btris,DEBUG_DIR + "bezier patch.stl");

  /*****************************************************************************
    2.9 Surfaces : Bezier surface, a composite of Bezier patches
  *****************************************************************************/

  cout << "2.9 Surfaces : Bezier surface, a composite of Bezier patches" << endl;

  // we make a twisted airfoil surface here
  std::vector<std::vector<TPoint<T>>> NACAsurfpoints;
  makeNACASurface(NACAsurfpoints,51);

  // create a Bezier surface of Bezier patches from these points
  TBezierSurface<T> bsurface(NACAsurfpoints,20,10);

  // these triangles are to display them in STL
  TTriangles<T> bstris;
  bsurface.createTriangles(bstris,51,26,0.5,1.0); 

  // save
  saveTrianglesStl(bstris,DEBUG_DIR + "bezier surface.stl");

  /*****************************************************************************
    2.10 Surfaces : point surface, a regilar net of points with linear 
      interpolation between
  *****************************************************************************/

  cout << "2.10 Surfaces : point surface, a regular net of points with linear interpolation between" << endl;

  // create a point surface
  TPointSurface<T> psurface(NACAsurfpoints);

  // these triangles are to display them in STL
  TTriangles<T> pstris;
  psurface.createTriangles(pstris,51,26,0.5,1.0); 

  // save
  saveTrianglesStl(pstris,DEBUG_DIR + "point surface.stl");

  /*****************************************************************************
    2.11 Surfaces : B-spline surface, interpolated and approximated, with and 
    without clamping
  *****************************************************************************/

  cout << "2.11 Surfaces : B-spline surface, interpolated and approximated, with and without clamping" << endl;

  //TSplineCurve<T> splinecurve0(NACAsurfpoints.back(),int(NACAsurfpoints.back().size()) - 1,
  //  SPLINE_DEGREE,END_CLAMPED,END_CLAMPED); 
  //// save to compare two curves
  //saveCurveIges(splinecurve0,DEBUG_DIR + "curve.iges");

  // interpolated...

  // create a point surface
  TSplineSurface<T> issurface(NACAsurfpoints,SPLINE_DEGREE,SPLINE_DEGREE,
    END_CLAMPED,END_CLAMPED,END_FREE,END_FREE);

  saveSurfaceIges(&issurface,DEBUG_DIR + "spline surface interpolated.iges");

  // these triangles are to display them in STL
  TTriangles<T> isstris;
  issurface.createTriangles(isstris,51,26,0.5,1.0); 

  // save
  saveTrianglesStl(isstris,DEBUG_DIR + "spline surface interpolated.stl");

  // approximated...

  // create a point surface
  TSplineSurface<T> apsurface(NACAsurfpoints,10,SPLINE_DEGREE,30,SPLINE_DEGREE,
    END_CLAMPED,END_FREE,END_FREE,END_FREE);

  saveSurfaceIges(&apsurface,DEBUG_DIR + "spline surface approximated.iges");

  // these triangles are to display them in STL
  TTriangles<T> apstris;
  apsurface.createTriangles(apstris,51,26,0.5,1.0); 

  // save
  saveTrianglesStl(apstris,DEBUG_DIR + "spline surface approximated.stl");

  // test all derivatives
  T Ustest = 0.5;
  T Vstest = 0.5;

  // note : first derivatives by parameters U,V are normalised for comparison
  // as they may be different in length due to typical non-uniform parameterisation
  // in approximated splines
  TPoint<T> bpos = bsurface.derivative(Ustest,Vstest,PARAMETER_ANY,0);
  TPoint<T> bderU = +bsurface.derivative(Ustest,Vstest,PARAMETER_U,1);
  TPoint<T> bderV = +bsurface.derivative(Ustest,Vstest,PARAMETER_V,1);

  TPoint<T> ppos = psurface.derivative(Ustest,Vstest,PARAMETER_ANY,0);
  TPoint<T> pderU = +psurface.derivative(Ustest,Vstest,PARAMETER_U,1);
  TPoint<T> pderV = +psurface.derivative(Ustest,Vstest,PARAMETER_V,1);

  TPoint<T> ispos = issurface.derivative(Ustest,Vstest,PARAMETER_ANY,0);
  TPoint<T> isderU = +issurface.derivative(Ustest,Vstest,PARAMETER_U,1);
  TPoint<T> isderV = +issurface.derivative(Ustest,Vstest,PARAMETER_V,1);

  TPoint<T> appos = apsurface.derivative(Ustest,Vstest,PARAMETER_ANY,0);
  TPoint<T> apderU = +apsurface.derivative(Ustest,Vstest,PARAMETER_U,1);
  TPoint<T> apderV = +apsurface.derivative(Ustest,Vstest,PARAMETER_V,1);

  /*****************************************************************************
    2.12 Surfaces : intersection of two, get intersection curve and two 
    parametric curves for trimming
  *****************************************************************************/

  cout << "2.12 Surfaces : intersection of two, get intersection curve and two parametric curves for trimming" << endl;

  // make a cylinder
  std::vector<std::vector<TPoint<T>>> cylpoints;
  makeCylinder(11,-2.0,+2.0,21,0.2,0.2,cylpoints,0.0,180.0);

  // create cylindrical surface
  TSplineSurface<T> cylsurface(cylpoints,SPLINE_DEGREE,SPLINE_DEGREE,
    END_CLAMPED,END_CLAMPED,END_FREE,END_FREE);

  TTransform<T> cylt;
  cylt.Rotate(TPoint<T>(0.0,1.0,0.0),90.0 * CPI);
  cylt.Rotate(TPoint<T>(1.0,0.0,0.0),-90.0 * CPI);

  cylsurface.makeTransform(&cylt);

  saveSurfaceIges(&cylsurface,DEBUG_DIR + "cylindrical surface.iges");

  std::vector<std::vector<TPoint<T>>> wcintersections; 
  std::vector<std::vector<TPoint<T>>> boundary0,boundary1;

  // intersect is here
  bool wcok = apsurface.intersect(cylsurface,true,wcintersections,boundary0,boundary1,
    PARM_TOLERANCE,
    MANY_POINTS2D,MANY_POINTS2D,0.5,1.0,1.0,1.0, // round leading edge at U = 0
    MANY_POINTS2D,MANY_POINTS2D,1.0,1.0,1.0,1.0);

  ASSERT(wcok);

  // intersected tris triangles to display them in STL
  TTriangles<T> apsurfacetris;
  apsurface.createTriangles(apsurfacetris,MANY_POINTS2D,MANY_POINTS2D,0.5); // round leading edge at U = 0
  saveTrianglesStl(apsurfacetris,DEBUG_DIR + "wing surface tris intersected.stl");

  // intersected tris triangles to display them in STL
  TTriangles<T> cylsurfacetris;
  cylsurface.createTriangles(cylsurfacetris,100,100); 
  saveTrianglesStl(cylsurfacetris,DEBUG_DIR + "cylinder surface tris intersected.stl");

  std::vector<std::vector<TPoint<T>>> cylboundarypoints;
  cylsurface.boundaryIntoPoints(boundary1,cylboundarypoints);

  std::vector<std::vector<TPoint<T>>> wcboundarypoints;
  apsurface.boundaryIntoPoints(boundary0,wcboundarypoints);

  saveLinesIges(wcintersections,DEBUG_DIR + "wing-cylinder intersection curve.iges");
  saveLinesIges(cylboundarypoints,DEBUG_DIR + "wing-cylinder intersection curve from cylinder.iges");
  saveLinesIges(wcboundarypoints,DEBUG_DIR + "wing-cylinder intersection curve from wing.iges");

  /*****************************************************************************
    2.13 Surfaces : intersect surface by plane, get intersection line and
    parametric curve for trimming
  *****************************************************************************/

  cout << "2.13 Surfaces : intersect surface by plane, get intersection line and parametric curve for trimming" << endl;

  // make plane by normal and point
  TPlane<T> plplane(TPoint<T>(0.7071,0.0,0.7071),TPoint<T>(0.0,0.0,-0.8));

  std::vector<std::vector<TPoint<T>>> plintersections; 
  std::vector<std::vector<TPoint<T>>> plboundary;

  // intersect is here
  bool plok = apsurface.intersectByPlane(plplane,plintersections,plboundary,
    NACAtolerance,PARM_TOLERANCE,
    MANY_POINTS2D,MANY_POINTS2D,0.5,1.0,1.0,1.0); // round leading edge at U = 0

  ASSERT(plok);

  std::vector<std::vector<TPoint<T>>> plboundarypoints;
  apsurface.boundaryIntoPoints(plboundary,plboundarypoints);

  saveLinesIges(plintersections,DEBUG_DIR + "wing-plane intersection curve.iges");
  saveLinesIges(plboundarypoints,DEBUG_DIR + "wing-plane intersection curve from surface.iges");

  /*****************************************************************************
    2.14 Surfaces : create B-spline surface from any other surface
  *****************************************************************************/

  cout << "2.14 Surfaces : create B-spline surface from any other surface" << endl;

  // create cylindrical Bezier surface
  TBezierSurface<T> bcsurface(cylpoints,20,20);

  // create spline surface from it
  TSplineSurface<T> bcspline(bcsurface,10,SPLINE_DEGREE,10,SPLINE_DEGREE,
    END_CLAMPED,END_CLAMPED,END_FREE,END_FREE,0.5,0.5,1.0,1.0);

  // these triangles are to display them in STL
  TTriangles<T> bcstris;
  bcsurface.createTriangles(bcstris,21,21,0.5,0.5); 

  // save
  saveTrianglesStl(bcstris,DEBUG_DIR + "bezier cylinder.stl");
  saveSurfaceIges(&bcspline,DEBUG_DIR + "cylindrical spline from Bezier.iges");

  /*****************************************************************************
    2.15 Surfaces : make a trimmed B-spline surface from intersections by plane. 
    Take the cut curve, close boundary and save this trimmed surface
  *****************************************************************************/

  cout << "2.15 Surfaces : make a trimmed B-spline surface from intersections by plane. Take the cut curve, close boundary and save this trimmed surface" << endl;

  // apsurface (B-spline approximated) is cut again by plane. We need to close 
  // this boundary and save this trimmed surface.

  saveSurfaceIges(&apsurface,DEBUG_DIR + "spline surface for cut by plane and trimming.iges");

  // make plane by normal and point, the cut part is gonna be a triangle
  TPlane<T> tr1plane(TPoint<T>(0.7071,0.0,0.7071),TPoint<T>(0.0,0.0,-0.8));

  std::vector<std::vector<TPoint<T>>> tr1intersections; 
  std::vector<std::vector<TPoint<T>>> tr1boundary;

  // intersect is here
  bool tr1ok = apsurface.intersectByPlane(tr1plane,tr1intersections,tr1boundary,
    NACAtolerance,PARM_TOLERANCE,
    MANY_POINTS2D,MANY_POINTS2D,0.5,1.0,1.0,1.0); // round leading edge at U = 0

  std::vector<std::vector<std::vector<TPoint<T>>>> tr1closedboundary;
  bool tr1ok1 = apsurface.closeBoundaryLoop(tr1boundary,tr1closedboundary);

  ASSERT(tr1ok1);

  saveTrimmedSurfaceIges(&apsurface,tr1closedboundary,DEBUG_DIR + "spline surface cut by plane 1 and trimmed.iges");

  // similar more complicated case for loop construction
  {
    // make plane by normal and point, the cut part is gonna be a triangle
    TPlane<T> tr1plane(TPoint<T>(0.7071,0.0,-0.7071),TPoint<T>(0.0,0.0,-0.2));

    std::vector<std::vector<TPoint<T>>> tr1intersections; 
    std::vector<std::vector<TPoint<T>>> tr1boundary;

    // intersect is here
    bool tr1ok = apsurface.intersectByPlane(tr1plane,tr1intersections,tr1boundary,
      NACAtolerance,PARM_TOLERANCE,
      MANY_POINTS2D,MANY_POINTS2D,0.5,1.0,1.0,1.0); // round leading edge at U = 0

    std::vector<std::vector<std::vector<TPoint<T>>>> tr1closedboundary;
    bool tr1ok1 = apsurface.closeBoundaryLoop(tr1boundary,tr1closedboundary);

    ASSERT(tr1ok1);

    saveTrimmedSurfaceIges(&apsurface,tr1closedboundary,DEBUG_DIR + "spline surface cut by plane 1a and trimmed.iges");
  }

  // now we flip the plane normal and try again 

  // make plane by normal and point
  TPlane<T> tr2plane(TPoint<T>(-0.7071,0.0,-0.7071),TPoint<T>(0.0,0.0,-0.8));

  std::vector<std::vector<TPoint<T>>> tr2intersections; 
  std::vector<std::vector<TPoint<T>>> tr2boundary;

  // intersect is here
  bool tr2ok = apsurface.intersectByPlane(tr2plane,tr2intersections,tr2boundary,
    NACAtolerance,PARM_TOLERANCE,
    MANY_POINTS2D,MANY_POINTS2D,0.5,1.0,1.0,1.0); // round leading edge at U = 0

  std::vector<std::vector<std::vector<TPoint<T>>>> tr2closedboundary;
  bool tr2ok1 = apsurface.closeBoundaryLoop(tr2boundary,tr2closedboundary);

  ASSERT(tr2ok1);

  saveTrimmedSurfaceIges(&apsurface,tr2closedboundary,DEBUG_DIR + "spline surface cut by plane 2 and trimmed.iges");

  /*****************************************************************************
    2.16 Surfaces : make a trimmed B-spline surface from surface-surface 
    intersection. Take the cut curve, close boundary and save this trimmed 
    surface
  *****************************************************************************/

  cout << "2.16 Surfaces : make a trimmed B-spline surface from surface-surface intersection. Take the cut curve, close boundary and save this trimmed surface" << endl;

  saveSurfaceIges(&apsurface,DEBUG_DIR + "this spline surface to be cut by another cylindrical surface.iges");
  saveSurfaceIges(&cylsurface,DEBUG_DIR + "cylindrical surface to be used in cutting.iges");

  std::vector<std::vector<TPoint<T>>> tr3intersections; 
  std::vector<std::vector<TPoint<T>>> tr3boundary0,tr3boundary1;

  // intersect is here
  bool tr3ok = apsurface.intersect(cylsurface,true,tr3intersections,tr3boundary0,tr3boundary1,
    PARM_TOLERANCE,
    MANY_POINTS2D,MANY_POINTS2D,0.5,1.0,1.0,1.0, // round leading edge at U = 0
    MANY_POINTS2D,MANY_POINTS2D,1.0,1.0,1.0,1.0);

  ASSERT(tr3ok);

  // close boundary on first surface
  std::vector<std::vector<std::vector<TPoint<T>>>> tr3closedboundary0;
  bool tr3ok1 = apsurface.closeBoundaryLoop(tr3boundary0,tr3closedboundary0);

  ASSERT(tr3ok1);

  saveTrimmedSurfaceIges(&apsurface,tr3closedboundary0,DEBUG_DIR + "spline surface cut by another surface and trimmed.iges");

  /*****************************************************************************
    2.17 Surfaces : make a hole (two loops) with trimmed B-spline surface
  *****************************************************************************/

  cout << "2.17 Surfaces : make a hole (two loops) with trimmed B-spline surface" << endl;

  // make a lower surface
  std::vector<std::vector<TPoint<T>>> NACAsurfpointslower;
  makeNACASurface(NACAsurfpointslower,51,TPoint<T>(1.0,-1.0,1.0));

  // to make correct normal
  reverseRows(NACAsurfpointslower);

  // create a point surface : free and clamped : LE/TE swapped!
  TSplineSurface<T> apsurfacelower(NACAsurfpointslower,10,SPLINE_DEGREE,30,SPLINE_DEGREE,
    END_FREE,END_CLAMPED,END_FREE,END_FREE);

  saveSurfaceIges(&apsurfacelower,DEBUG_DIR + "spline surface approximated (lower).iges");

  // now we've got apsurface and apsurfacelower for an airfoil

  // cylindrical surface cylsurface, combine them all into trimmed surfaces

  std::vector<std::vector<TPoint<T>>> tr4intersections; 
  std::vector<std::vector<TPoint<T>>> tr4boundary0,tr4boundary1;

  bool tr4ok = apsurface.intersect(cylsurface,true,tr4intersections,tr4boundary0,tr4boundary1,
    PARM_TOLERANCE,
    MANY_POINTS2D,MANY_POINTS2D,0.5,1.0,1.0,1.0, // round leading edge at U = 0
    MANY_POINTS2D,MANY_POINTS2D,1.0,1.0,1.0,1.0);

  ASSERT(tr4ok);

  // close boundary on first surface
  std::vector<std::vector<std::vector<TPoint<T>>>> tr4closedboundary0;
  bool tr4ok1 = apsurface.closeBoundaryLoop(tr4boundary0,tr4closedboundary0);

  ASSERT(tr4ok1);

  saveTrimmedSurfaceIges(&apsurface,tr4closedboundary0,DEBUG_DIR + "upper airfoil surface cut by another surface and trimmed.iges");

  std::vector<std::vector<TPoint<T>>> tr5intersections; 
  std::vector<std::vector<TPoint<T>>> tr5boundary0,tr5boundary1;

  bool tr5ok = apsurfacelower.intersect(cylsurface,true,tr5intersections,tr5boundary0,tr5boundary1,
    PARM_TOLERANCE,
    MANY_POINTS2D,MANY_POINTS2D,0.5,1.0,1.0,1.0, // round leading edge at U = 0
    MANY_POINTS2D,MANY_POINTS2D,1.0,1.0,1.0,1.0);

  ASSERT(tr5ok);

  // close boundary on first surface
  std::vector<std::vector<std::vector<TPoint<T>>>> tr5closedboundary0;
  bool tr5ok1 = apsurfacelower.closeBoundaryLoop(tr5boundary0,tr5closedboundary0);

  ASSERT(tr5ok1);

  saveTrimmedSurfaceIges(&apsurfacelower,tr5closedboundary0,DEBUG_DIR + "lower airfoil surface cut by another surface and trimmed.iges");

  // cut cylinder, its trimming curve contains two parts from upper and lower surfaces
  std::vector<std::vector<TPoint<T>>> tr6intersections; 
  std::vector<std::vector<TPoint<T>>> tr6boundary0,tr6boundary1;

  // intersect with upper surface...
  bool tr6ok = cylsurface.intersect(apsurface,true,tr6intersections,tr6boundary0,tr6boundary1,
    PARM_TOLERANCE,
    MANY_POINTS2D,MANY_POINTS2D,0.5,1.0,1.0,1.0, // round leading edge at U = 0
    MANY_POINTS2D,MANY_POINTS2D,1.0,1.0,1.0,1.0);

  ASSERT(tr6ok);

  // intersect with lower surface...
  bool tr7ok = cylsurface.intersect(apsurfacelower,true,tr6intersections,tr6boundary0,tr6boundary1,
    PARM_TOLERANCE,
    MANY_POINTS2D,MANY_POINTS2D,0.5,1.0,1.0,1.0, // round leading edge at U = 0
    MANY_POINTS2D,MANY_POINTS2D,1.0,1.0,1.0,1.0);

  ASSERT(tr7ok);

  // close boundary (hole) on first surface
  std::vector<std::vector<std::vector<TPoint<T>>>> tr6closedboundary0;
  bool tr6ok1 = cylsurface.closeBoundaryLoop(tr6boundary0,tr6closedboundary0);

  ASSERT(tr6ok1);

  // add outer boundary
  std::vector<std::vector<TPoint<T>>> tr6outerloop;

  saveTrimmedSurfaceIges(&cylsurface,tr6closedboundary0,DEBUG_DIR + "cylindrical surface cut by two airfoil parts and trimmed.iges");

  /*****************************************************************************
    2.18 Surfaces : find parameters U,V for a point on (or close to) surface 
  *****************************************************************************/

  cout << "2.18 Surfaces : find parameters U,V for a point on (or close to) surface" << endl;

  {
    // prepare (fine) mesh of points with nodal parameters
    // create finer mesh by createPoints() for more accurate results
    std::vector<TPoint<T>> points,UVpoints;
    int k1 = 0;
    int k2 = 0;
    apsurface.createPoints(points,&UVpoints,&k1,&k2);

    for (int i = 0; i < 10; i++)
    {
      TPoint<T> UV(random(),random());

      TPoint<T> pos = apsurface.position(UV.X,UV.Y);
      TPoint<T> UVap = apsurface.findUVforPoint(points,UVpoints,k1,k2,pos);

      T distap = !(UVap - UV);
      ASSERT(distap < 0.015);
    }
  }

  /*****************************************************************************

    Part 3 : VOLUMES

  *****************************************************************************/

  cout << "    Part 3 : VOLUMES" << endl;

  /*****************************************************************************
    3.1 Volumes : create a Bezier volume around our wing, check definition of
    parametric values U,V,W for inner X,Y,Z points
  *****************************************************************************/

  {
    // min/max for wing
    TPoint<T> min,max;
    apsurface.calculateMinMax(&min,&max);

    // extend min/max
    extendMinMax(min,max,1.01);

    // create box as Bezier volume around the wing
    TBezierVolume<T> volume(min,max,5,2,6);

    // create fine mesh of points
    std::vector<TPoint<T>> points,UVWpoints;
    int k1 = 0;
    int k2 = 0;
    int k3 = 0;
    volume.createPoints(points,&UVWpoints,&k1,&k2,&k3,11,11,11); 

    // test UVW definitions
    for (int i = 0; i < 10; i++)
    {
      TPoint<T> UVW(random(),random(),random());

      TPoint<T> pos = volume.position(UVW.X,UVW.Y,UVW.Z);
      TPoint<T> UVWap = volume.findUVWforPoint(points,UVWpoints,k1,k2,k3,pos);

      T distap = !(UVWap - UVW);
      ASSERT(distap < 0.00001);
    }

    // test first derivatives
    for (int i = 0; i < 10; i++)
    {
      TPoint<T> UVW(random(),random(),random());

      TPoint<T> derU = volume.derivative(UVW.X,UVW.Y,UVW.Z,PARAMETER_U,1);
      TPoint<T> derV = volume.derivative(UVW.X,UVW.Y,UVW.Z,PARAMETER_V,1);
      TPoint<T> derW = volume.derivative(UVW.X,UVW.Y,UVW.Z,PARAMETER_W,1);
    }
  }

  /*****************************************************************************
    3.2 Volumes : compare Bezier and B-spline based volume derivatives
  *****************************************************************************/

  {
    // min/max for wing
    TPoint<T> min,max;
    apsurface.calculateMinMax(&min,&max);

    // extend min/max
    extendMinMax(min,max,1.01);

    // create box as Bezier volume around the wing
    TBezierVolume<T> bvolume(min,max,5,2,6);
    TSplineVolume<T> svolume(min,max,SPLINE_DEGREE,SPLINE_DEGREE,SPLINE_DEGREE,8,8,8);

    // compare first derivatives after calling makeSplineControlPoints() with uniformparms;
    // results become identical
    for (int i = 0; i < 10; i++)
    {
      TPoint<T> UVW(random(),random(),random());

      TPoint<T> bderU = bvolume.derivative(UVW.X,UVW.Y,UVW.Z,PARAMETER_U,1);
      TPoint<T> bderV = bvolume.derivative(UVW.X,UVW.Y,UVW.Z,PARAMETER_V,1);
      TPoint<T> bderW = bvolume.derivative(UVW.X,UVW.Y,UVW.Z,PARAMETER_W,1);

      TPoint<T> sderU = svolume.derivative(UVW.X,UVW.Y,UVW.Z,PARAMETER_U,1);
      TPoint<T> sderV = svolume.derivative(UVW.X,UVW.Y,UVW.Z,PARAMETER_V,1);
      TPoint<T> sderW = svolume.derivative(UVW.X,UVW.Y,UVW.Z,PARAMETER_W,1);

      T diff0 = !(bderU - sderU);
      T diff1 = !(bderV - sderV);
      T diff2 = !(bderW - sderW);
      ASSERT(diff0 < 0.0001);
      ASSERT(diff1 < 0.0001);
      ASSERT(diff2 < 0.0001);
    }

    std::vector<TSplineSurface<T> *> faces;
    for (int i = 0; i < 6; i++)
    {
      TSplineSurface<T> *face = svolume.getFace(i);
      faces.push_back(face);
    }

    saveSurfacesIges(faces,DEBUG_DIR + "box spline volume.iges");

    for (int i = 0; i < 6; i++)
    {
      DELETE_CLASS(faces[i]);
    }
  }

  /*****************************************************************************

    Part 4 : FFD

  *****************************************************************************/

  cout << "    Part 4 : FFD" << endl;

  /*****************************************************************************
    4.1 FFD : distort wing shape by displacement of two points
  *****************************************************************************/

  cout << "4.1 FFD : distort wing shape by displacement of two points" << endl;

  {
    // make a upper/lower airfoil surfaces
    std::vector<std::vector<TPoint<T>>> upper,lower;
    makeNACASurface(upper,51,TPoint<T>(1.0,+1.0,1.0),1.0,1.0,0.02,0.0,0.0); 
    makeNACASurface(lower,51,TPoint<T>(1.0,-1.0,1.0),1.0,1.0,0.02,0.0,0.0); 

    // to make correct normal
    reverseRows(upper);

    // make two spline surfaces, LE/TE swapped for one surface - see FREE and CLAMPED
    TSplineSurface<T> supper(upper,30,SPLINE_DEGREE,30,SPLINE_DEGREE,
      END_FREE,END_CLAMPED,END_FREE,END_FREE);
    TSplineSurface<T> slower(lower,30,SPLINE_DEGREE,30,SPLINE_DEGREE,
      END_CLAMPED,END_FREE,END_FREE,END_FREE);

    // bwing and wind are identical
    std::vector<tcad::TSplineSurface<T> *> wing;
    std::vector<tcad::TBaseSurface<T> *> bwing;
    wing.push_back(&supper);
    wing.push_back(&slower);
    bwing.push_back(&supper);
    bwing.push_back(&slower);

    saveSurfacesIges(wing,DEBUG_DIR + "wing before FFD.iges");

    // apply FFD : we want to distort the wing by moving these points to new positions:
    std::vector<TPoint<T>> oldpositions = {TPoint<T>(-0.5,0.0,1.0),TPoint<T>(0.5,0.0,1.0)};
    std::vector<TPoint<T>> newpositions = {TPoint<T>(-0.4,0.0,1.0),TPoint<T>(0.6,0.2,0.9)};
    
    FFD<T> ffd(bwing,oldpositions,newpositions);

    saveSurfacesIges(wing,DEBUG_DIR + "wing after FFD.iges");

    // make spline volume from distored ffd box
    TSplineVolume<T> swing(ffd,10,10,10);

    std::vector<TSplineSurface<T> *> faces;
    for (int i = 0; i < 6; i++)
    {
      TSplineSurface<T> *face = swing.getFace(i);
      faces.push_back(face);
    }

    saveSurfacesIges(faces,DEBUG_DIR + "ffd box spline volume.iges");

    for (int i = 0; i < 6; i++)
    {
      DELETE_CLASS(faces[i]);
    }

  }

  /*****************************************************************************

    Part 5 : blocks, various useful geometries

  *****************************************************************************/

  cout << "    Part 5 : BLOCKS" << endl;

  /*****************************************************************************
    5.1 Blocks : airfoild : preparing data to make an airfoil
  *****************************************************************************/

  cout << "5.1 Blocks : airfoil : preparing data to make an airfoil" << endl;

  {
    // step 1 : camber surface, cylindrical

    // first, we need to generate camber surface with thickness, they are all
    // represented as regular nets of points
    std::vector<std::vector<TPoint<T>>> camberpoints;

    // suppose we shall make 11 (chord) x 6 (span) points for a mesh
    int N1 = 11;
    int N2 = 6;

    // cylinder is starting geometry
    makeCylinder(N2,0.01,0.09,N1,0.18,0.18,camberpoints,10.0,190.0);

    // reverse rows to make LE to the left
    reverseRows(camberpoints);

    // to make normal point to upper surface
    reverseColumns(camberpoints);

    // create camber surface
    TSplineSurface<T> cambersurface(camberpoints,SPLINE_DEGREE,SPLINE_DEGREE,
      END_CLAMPED,END_CLAMPED,END_FREE,END_FREE);

    saveSurfaceIges(&cambersurface,DEBUG_DIR + "camber surface step 1.iges");

    // step 2 : camber surface, distort by FFD
    // bwing and wind are identical
    std::vector<tcad::TSplineSurface<T> *> wing;
    std::vector<tcad::TBaseSurface<T> *> bwing;
    wing.push_back(&cambersurface);
    bwing.push_back(&cambersurface);

    // apply FFD : we want to distort the wing by moving these points to new positions:
    std::vector<TPoint<T>> oldpositions = {
      TPoint<T>(-0.177,-0.031,0.01)
    };

    std::vector<TPoint<T>> newpositions = {
      TPoint<T>(-0.23,-0.056,-0.023)
    };
    
    FFD<T> ffd(bwing,oldpositions,newpositions);

    saveSurfacesIges(wing,DEBUG_DIR + "camber surface step 2.iges");

    // these triangles are to save them in STL, this mesh has TE to the right
    // 
    TTriangles<T> tris;
    cambersurface.createTriangles(tris,N1,N2,1.0,1.0); 

    // save
    saveTrianglesStl(tris,DEBUG_DIR + "camber surface step 2.stl");

    // now the camber surface is defined by N1 x N2 mesh, let's make a thickness array
    // remembering that U goes from left to right, U = 0 being LE, V increases with
    // increasing Z
    T chord = cambersurface.Usize() * 2.0;
    std::vector<std::vector<T>> thickness = {
      {0.0, 0.046 * chord, 0.057 * chord, 0.060 * chord, 0.058 * chord, 0.053 * chord, 0.045 * chord, 0.036 * chord, 0.025 * chord, 0.014 * chord, 0.0 * chord},
      {0.0, 0.046 * chord, 0.057 * chord, 0.060 * chord, 0.058 * chord, 0.053 * chord, 0.045 * chord, 0.036 * chord, 0.025 * chord, 0.014 * chord, 0.0 * chord},
      {0.0, 0.046 * chord, 0.057 * chord, 0.060 * chord, 0.058 * chord, 0.053 * chord, 0.045 * chord, 0.036 * chord, 0.025 * chord, 0.014 * chord, 0.0 * chord},
      {0.0, 0.046 * chord, 0.057 * chord, 0.060 * chord, 0.058 * chord, 0.053 * chord, 0.045 * chord, 0.036 * chord, 0.025 * chord, 0.014 * chord, 0.0 * chord},
      {0.0, 0.046 * chord, 0.057 * chord, 0.060 * chord, 0.058 * chord, 0.053 * chord, 0.045 * chord, 0.036 * chord, 0.025 * chord, 0.014 * chord, 0.0 * chord},
      {0.0, 0.046 * chord, 0.057 * chord, 0.060 * chord, 0.058 * chord, 0.053 * chord, 0.045 * chord, 0.036 * chord, 0.025 * chord, 0.014 * chord, 0.0 * chord}
    };

  /*****************************************************************************
    5.2 Blocks : airfoil : how to make an airfoil from camber surface and 
    thickness
  *****************************************************************************/

  cout << "5.2 Blocks : airfoil : how to make an airfoil from camber surface and thickness" << endl;

    // make camber points from FFD-distorted wing
    std::vector<TPoint<T>> temp;
    std::vector<std::vector<TPoint<T>>> newcamberpoints;
    int k1 = 0;
    int k2 = 0;
    cambersurface.createPoints(temp,nullptr,&k1,&k2,N1,N2);
    ASSERT(k1 + 1 == N1);
    ASSERT(k2 + 1 == N2);
    points1Dto2D(temp,k1,k2,newcamberpoints);

    // we have all the data : camber + thickness N1 x N2 arrays, make spline surfaces
    std::vector<TSplineSurface<T> *> surfaces;
    makeAirfoil(newcamberpoints,thickness,surfaces); 

    saveSurfacesIges(surfaces,DEBUG_DIR + "airfoil surfaces.iges");

    deleteSurfaces(surfaces);
  }

  /*****************************************************************************
    5.3 Blocks : box : solid output
  *****************************************************************************/

  cout << "5.3 Blocks : box : solid output" << endl;

  {
    TPoint<T> min(0.0,0.0,0.0);
    TPoint<T> max(1.0,1.0,1.0);
    TSplineVolume<T> svolume(min,max,SPLINE_DEGREE,SPLINE_DEGREE,SPLINE_DEGREE,8,8,8);

    std::vector<TSplineSurface<T> *> surfaces;
    for (int i = 0; i < 6; i++)
    {
      TSplineSurface<T> *face = svolume.getFace(i);
      surfaces.push_back(face);
    }

    // name for debugging
    nameSurfaces<T>(surfaces,"box");

    // make boundaries
    std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> boundariesUV;
    closeOuterBoundary<T>(surfaces,boundariesUV);

    saveSurfacesIges(surfaces,DEBUG_DIR + "solid box surfaces.iges");

    saveTrimmedSurfacesIges(surfaces,boundariesUV,DEBUG_DIR + "solid box surfaces trimmed.iges");

    bool ok = saveSolidIges(surfaces,boundariesUV,DEBUG_DIR + "solid box.iges",NACAtolerance);

    ASSERT(ok);

    deleteSurfaces(surfaces);
  }

#endif

#ifdef DEBUG_SUBMARINE

#if 1 

  /*****************************************************************************
    5.4 Blocks : blade
  *****************************************************************************/

  cout << "5.4 Blocks : blade" << endl;

  {
    T subL = 74.0;
    T tolerance = subL * PARM_TOLERANCE;

    std::vector<TPoint<T>> upperlower; 
    std::pair<T,T> res = makeAirfoilPointsXY<T>(E178<T>,true,false,50,upperlower);

    std::vector<std::vector<TPoint<T>>> camberpoints; 
    std::vector<std::vector<T>> thickness;
    makeBladeCamberThickness<T>(upperlower,KiloBlade<T>,50,20,camberpoints,thickness,tolerance); 

    std::vector<TSplineSurface<T> *> surfaces;
    makeAirfoil<T>(camberpoints,thickness,surfaces,SPLINE_DEGREE,SPLINE_DEGREE,END_FREE,END_FREE,END_FREE,END_FREE);

    // rotate around Z
    TTransform<T> t;
    t.Rotate(TPoint<T>(0.0,0.0,1.0),-90.0 * CPI); 
    makeTransform<T>(surfaces,&t);

    saveSurfacesIges(surfaces,DEBUG_DIR + "Kilo propeller surfaces.iges");

    deleteSurfaces(surfaces);

  /*****************************************************************************
    5.5 Blocks : surfaces of revolution : propeller hub
  *****************************************************************************/

    cout << "5.5 Blocks : surfaces of revolution : propeller hub" << endl;

    makeSurfacesOfRevolution<T>(KiloPropHub<T>,AxisZ,AxisX,7,9,8,8,surfaces,0.0,360.0,SPLINE_DEGREE,SPLINE_DEGREE,
      END_CLAMPED,END_CLAMPED,END_FREE,END_FREE); 
    makeSurfacesOfRevolution<T>(KiloPropHubEnd<T>,AxisZ,AxisX,7,9,8,8,surfaces,0.0,360.0,SPLINE_DEGREE,SPLINE_DEGREE,
      END_CLAMPED,END_CLAMPED,END_FREE,END_FREE); 

    t.LoadIdentity();
    t.Translate(TPoint<T>(0.0,0.0,0.3)); 
    t.Rotate(TPoint<T>(0.0,1.0,0.0),+90.0 * CPI); 
    makeTransform<T>(surfaces,&t);

    saveSurfacesIges(surfaces,DEBUG_DIR + "Kilo propeller hub surfaces.iges");

  /*****************************************************************************
    5.6 Blocks : surfaces of revolution : propeller hub solid
  *****************************************************************************/

    cout << "5.6 Blocks : surfaces of revolution : propeller hub solid" << endl;

    // name for debugging
    nameSurfaces<T>(surfaces,"hub");

    // make boundaries
    std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> boundariesUV;
    closeOuterBoundary<T>(surfaces,boundariesUV);

    bool ok = saveSolidIges(surfaces,boundariesUV,DEBUG_DIR + "Kilo propeller hub solid.iges",
      tolerance);

    ASSERT(ok);

    deleteSurfaces(surfaces);
  }

  /*****************************************************************************
    5.7 Blocks : propeller, trimmed and solid
  *****************************************************************************/

  cout << "5.7 Blocks : propeller, trimmed and solid" << endl;

#ifdef _DEBUG
  errorMessage("This thing (makeTrimming()) may be slow in Debug, try Release instead.");
#endif

  // keep them for submarine, propeller surfaces
  std::vector<TSplineSurface<T> *> propsurfaces;
  // propeller trimming curves
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> propboundariesUV;

  {
    T propR = 1.3;
    T tolerance = propR * PARM_TOLERANCE;

    // blade

    std::vector<TPoint<T>> upperlower; 
    std::pair<T,T> res = makeAirfoilPointsXY<T>(E178<T>,true,false,50,upperlower);

    std::vector<std::vector<TPoint<T>>> camberpoints; 
    std::vector<std::vector<T>> thickness;
    makeBladeCamberThickness<T>(upperlower,KiloBlade<T>,50,20,camberpoints,thickness,tolerance); 

    std::vector<TSplineSurface<T> *> surfaces;
    makeAirfoil<T>(camberpoints,thickness,surfaces,SPLINE_DEGREE,SPLINE_DEGREE,END_FREE,END_FREE,END_FREE,END_FREE);

    // rotate around Z
    TTransform<T> t;
    t.Rotate(TPoint<T>(0.0,0.0,1.0),-90.0 * CPI); 
    makeTransform<T>(surfaces,&t);

    int numblades = 7;
    T da = 360.0 / T(numblades);

    nameSurfaces(surfaces,"blade0");

    for (int i = 1; i < numblades; i++)
    {
      TSplineSurface<T> *blade0 = new TSplineSurface<T>(*surfaces[0]);
      TSplineSurface<T> *blade1 = new TSplineSurface<T>(*surfaces[1]);

      T a = T(i) * da;

      TTransform<T> t;
      t.Rotate(TPoint<T>(1.0,0.0,0.0),a * CPI); 
      blade0->makeTransform(&t);
      blade1->makeTransform(&t);

      blade0->name = "blade" + to_string(i) + "0";
      blade1->name = "blade" + to_string(i) + "1";

      surfaces.push_back(blade0);
      surfaces.push_back(blade1);
    }

    // move propeller to hull position
    t.LoadIdentity();
    t.Translate(TPoint<T>(-0.3 - 0.1 - 29.7,0.0,0.0)); 
    makeTransform<T>(surfaces,&t);

    saveSurfacesIges(surfaces,DEBUG_DIR + "Kilo propeller surfaces 1.iges");

    // move propeller to hull position
    t.LoadIdentity();
    t.Translate(TPoint<T>(+0.3 + 0.1 + 29.7,0.0,0.0)); 
    makeTransform<T>(surfaces,&t);

    std::vector<TSplineSurface<T> *> hsurfaces;

    // + hub
#if 0 //!!!
    makeSurfacesOfRevolution<T>(KiloPropHub<T>,AxisZ,AxisX,2,17,16,8,hsurfaces,0.0,360.0,SPLINE_DEGREE,SPLINE_DEGREE,
      END_CLAMPED,END_CLAMPED,END_FREE,END_FREE); 
    makeSurfacesOfRevolution<T>(KiloPropHubEnd<T>,AxisZ,AxisX,2,17,16,8,hsurfaces,0.0,360.0,SPLINE_DEGREE,SPLINE_DEGREE,
      END_CLAMPED,END_CLAMPED,END_FREE,END_FREE); 
#else
    makeSurfacesOfRevolution<T>(KiloPropHub<T>,AxisZ,AxisX,7,9,16,16,hsurfaces,0.0,360.0,SPLINE_DEGREE,SPLINE_DEGREE,
      END_CLAMPED,END_CLAMPED,END_FREE,END_FREE); 
    makeSurfacesOfRevolution<T>(KiloPropHubEnd<T>,AxisZ,AxisX,7,9,16,16,hsurfaces,0.0,360.0,SPLINE_DEGREE,SPLINE_DEGREE,
      END_CLAMPED,END_CLAMPED,END_FREE,END_FREE); 
    //!!! makeSurfacesOfRevolution<T>(KiloPropHub<T>,7,9,8,8,hsurfaces,0.0,360.0,SPLINE_DEGREE,SPLINE_DEGREE,
    //  END_CLAMPED,END_CLAMPED,END_FREE,END_FREE); 
    //makeSurfacesOfRevolution<T>(KiloPropHubEnd<T>,7,9,8,8,hsurfaces,0.0,360.0,SPLINE_DEGREE,SPLINE_DEGREE,
    //  END_CLAMPED,END_CLAMPED,END_FREE,END_FREE); 
#endif

    t.LoadIdentity();
    t.Translate(TPoint<T>(0.0,0.0,0.3)); 
    t.Rotate(TPoint<T>(0.0,1.0,0.0),+90.0 * CPI); 
    makeTransform<T>(hsurfaces,&t);

    // move propeller to hull position
    t.LoadIdentity();
    t.Translate(TPoint<T>(-0.3 - 0.1 - 29.7,0.0,0.0)); 
    makeTransform<T>(hsurfaces,&t);

    nameSurfaces(hsurfaces,"hub");

    saveSurfacesIges(hsurfaces,DEBUG_DIR + "Kilo propeller surfaces 2.iges");

    // move propeller to hull position
    t.LoadIdentity();
    t.Translate(TPoint<T>(+0.3 + 0.1 + 29.7,0.0,0.0)); 
    makeTransform<T>(hsurfaces,&t);

    // make boundaries
    std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> boundariesUV;
    std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> boundariesUV1;

    // mutual intersections between faces, in process, estimate big tolerance as max
    // difference between boundary curves
    makeTrimming(surfaces,hsurfaces,boundariesUV,boundariesUV1,UNITE, 
      propsurfaces,propboundariesUV,
      tolerance,PARM_TOLERANCE,true);

    // move propeller to hull position
    t.LoadIdentity();
    t.Translate(TPoint<T>(-0.3 - 0.1 - 29.7,0.0,0.0)); 
    makeTransform<T>(propsurfaces,&t);

    saveTrimmedSurfacesIges(propsurfaces,propboundariesUV,DEBUG_DIR + "Kilo propeller surfaces trimmed.iges");

    std::vector<std::vector<TPoint<T>>> badedges;
    bool ok = saveSolidIges(propsurfaces,propboundariesUV,DEBUG_DIR + "Kilo propeller solid.iges",
      tolerance,PARM_TOLERANCE,SPLINE_DEGREE,18,&badedges);

    ASSERT(ok);

    if (!ok)
    {
      saveLinesIges<T>(badedges,DEBUG_DIR + "badedges.iges");
    }
  }

  /*****************************************************************************
    5.8 Blocks : submarine hull, axisymmetric solid
  *****************************************************************************/

  cout << "5.8 Blocks : submarine hull, axisymmetric solid" << endl;

  // hull surfaces
  std::vector<TSplineSurface<T> *> hullsurfaces;
  // hull trimming curves
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> hullboundariesUV;

  {
    T hullL = 74.0;
    T tolerance = hullL * PARM_TOLERANCE;

    // make axisymmetric hull
    smoothPointsByBezier(KiloHull<T>[1],END_FIXED,END_CLAMPED,50);
    smoothPointsByBezier(KiloHull<T>[3],END_CLAMPED,END_CLAMPED,50);
    makeSurfacesOfRevolution<T>(KiloHull<T>,AxisZ,AxisX,8,9,16,32,hullsurfaces,0.0,360.0,SPLINE_DEGREE,SPLINE_DEGREE,
      END_CLAMPED,END_CLAMPED,END_FREE,END_CLAMPED); 

    TTransform<T> t;
    t.Rotate(TPoint<T>(0.0,1.0,0.0),+90.0 * CPI); 
    makeTransform<T>(hullsurfaces,&t);

    saveSurfacesIges(hullsurfaces,DEBUG_DIR + "Kilo sub hull surfaces.iges");

    // name for debugging
    nameSurfaces<T>(hullsurfaces,"hull");

    // make boundaries
    closeOuterBoundary<T>(hullsurfaces,hullboundariesUV);

    saveTrimmedSurfacesIges(hullsurfaces,hullboundariesUV,DEBUG_DIR + "Kilo sub hull trimmed.iges");

    std::vector<std::vector<TPoint<T>>> badedges;
    bool ok = saveSolidIges(hullsurfaces,hullboundariesUV,DEBUG_DIR + "Kilo sub hull solid.iges",
      tolerance,PARM_TOLERANCE,SPLINE_DEGREE,18,&badedges);

    ASSERT(ok);

    if (!ok)
    {
      saveLinesIges<T>(badedges,DEBUG_DIR + "badedges.iges");
    }
  }

  /*****************************************************************************
    5.9 Blocks : attach propeller to submarine hull by shaft
  *****************************************************************************/

  //!!! The process is not very fast, making it faster !!!

  cout << "5.9 Blocks : attach propeller to submarine hull by shaft" << endl;

  // shaft : axisymmetric faces, not solid, to connect propeller to hull
  // and make single solid
  std::vector<TSplineSurface<T> *> shaftsurfaces;
  // shaft trimming curves
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> shaftboundariesUV;

  // this is the whole sub
  std::vector<TSplineSurface<T> *> subsurfaces;
  // hull trimming curves
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> subboundariesUV;

  {
    T hullL = 74.0;
    T tolerance = hullL * PARM_TOLERANCE;

    std::vector<TPoint<T>> shaftcontour;
    shaftcontour.push_back(TPoint<T>(0.2,0.0,-29.9));
    shaftcontour.push_back(TPoint<T>(0.2,0.0,-29.6));

    // make axisymmetric shaft
    makeSurfacesOfRevolution<T>(shaftcontour,AxisZ,AxisX,8,9,8,8,shaftsurfaces,0.0,360.0,SPLINE_DEGREE,SPLINE_DEGREE,
      END_FREE,END_FREE,END_FREE,END_FREE); 

    TTransform<T> t;
    t.Rotate(TPoint<T>(0.0,1.0,0.0),+90.0 * CPI); 
    makeTransform<T>(shaftsurfaces,&t);

    saveSurfacesIges(shaftsurfaces,DEBUG_DIR + "Kilo sub shaft.iges");

    // name for debugging
    nameSurfaces<T>(shaftsurfaces,"shaft");

    // make boundaries
    closeOuterBoundary<T>(shaftsurfaces,shaftboundariesUV);

    // whole sub
    subsurfaces = propsurfaces;
    subboundariesUV = propboundariesUV;

    subsurfaces.insert(subsurfaces.end(),hullsurfaces.begin(),hullsurfaces.end());
    subboundariesUV.insert(subboundariesUV.end(),hullboundariesUV.begin(),hullboundariesUV.end());

    saveTrimmedSurfacesIges(subsurfaces,subboundariesUV,DEBUG_DIR + "Kilo sub surfaces trimmed before connection with shaft.iges");

    // make a solid from shaft and propeller+hull

    // make solid, do not clear old boundaries
    makeTrimming(shaftsurfaces,subsurfaces,shaftboundariesUV,subboundariesUV,UNITE,
      subsurfaces,subboundariesUV,
      tolerance,PARM_TOLERANCE,false);

    saveTrimmedSurfacesIges(subsurfaces,subboundariesUV,DEBUG_DIR + "Kilo sub surfaces trimmed.iges");

    std::vector<std::vector<TPoint<T>>> badedges;
    bool ok = saveSolidIges(subsurfaces,subboundariesUV,DEBUG_DIR + "Kilo sub solid with propeller.iges",
      tolerance,PARM_TOLERANCE,SPLINE_DEGREE,18,&badedges);

    ASSERT(ok);

    if (!ok)
    {
      saveLinesIges<T>(badedges,DEBUG_DIR + "badedges.iges");
    }
  }

  /*****************************************************************************
    5.10 Blocks : submarine fin
  *****************************************************************************/

  cout << "5.10 Blocks : sub fin" << endl;

  // sub fin
  std::vector<TSplineSurface<T> *> subfinsurfaces;
  // sub fin trimming curves
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> subfinboundariesUV;

  {
    // base (lowest) fin contour
    std::vector<TPoint<T>> basepoints = KiloFinBase<T>;
    smoothPointsByOrtho(basepoints,END_FIXED,END_ROUNDED,4,GAUSSINT_8,true,50,0.5,0.5);

    TPointCurve<T> basecurve(basepoints);
    saveCurveIges(basecurve,DEBUG_DIR + "Sub fin base curve.iges");

    std::vector<TPoint<T>> upperpoints = basepoints;

    // upper fin contour
    T XCOEF = 42.0 / 44.0;
    TTransform<T> t;
    t.Resize(TPoint<T>(XCOEF,0.9,1.0));
    // the stem wall is vertical at height Z = 17.0 * KHCOEF
    t.Translate(TPoint<T>(0.272,0.0,17.0 * KHCOEF));

    makeTransform(upperpoints,&t);

    TPointCurve<T> uppercurve(upperpoints);
    saveCurveIges(uppercurve,DEBUG_DIR + "Sub fin upper curve.iges");

    // top contour ("rounding top")
    T radius = KHCOEF;
    std::vector<TPoint<T>> toppoints = upperpoints;

    TPoint<T> min,max;
    calculateMinMax(toppoints,&min,&max);
    TPoint<T> d = max - min;
    TPoint<T> centre = (min + max) * 0.5;
    T xcoef = (d.X - radius * 2.0) / d.X;
    T ycoef = (d.Y - radius) / d.Y;

    t.LoadIdentity();
    t.Translate(-centre);
    t.Resize(TPoint<T>(xcoef,ycoef,1.0));
    t.Translate(centre + TPoint<T>(0.0,-radius * 0.5,radius));

    makeTransform(toppoints,&t);

    TPointCurve<T> topcurve(toppoints);
    saveCurveIges(topcurve,DEBUG_DIR + "Sub fin top curve.iges");

    // make surfaces now
    std::vector<std::vector<TPoint<T>>> points;
    points.push_back(basepoints);
    points.push_back(upperpoints);
    points.push_back(toppoints);

    // make two symmetric (around XZ) walls and ceiling
    makeTwoWallsAndFlatBottomTop<T>(points,subfinsurfaces,false,true,40,40,8,SPLINE_DEGREE,
      END_CLAMPED,END_CLAMPED,END_FREE,END_FREE);

    // move walls and ceiling to its position
    t.LoadIdentity();
    t.Translate(TPoint<T>(50.0 * KHCOEF,0.0,13.0 * KHCOEF));
    makeTransform(subfinsurfaces,&t);

    saveSurfacesIges(subfinsurfaces,DEBUG_DIR + "Sub fin surfaces.iges");
  }

  /*****************************************************************************
    5.11 Blocks : submarine hull + propeller + fin, all single solid
  *****************************************************************************/

  cout << "5.11 Blocks : submarine hull + propeller + fin, all single solid" << endl;

  {
    T hullL = 74.0;
    T tolerance = hullL * PARM_TOLERANCE;

    // name for debugging
    nameSurfaces<T>(subfinsurfaces,"fin");

    // make boundaries
    closeOuterBoundary<T>(subfinsurfaces,subfinboundariesUV);

    // make a solid from fin and propeller + hull

    // make solid, do not clear old boundaries
    makeTrimming(subfinsurfaces,subsurfaces,subfinboundariesUV,subboundariesUV,UNITE,
      subsurfaces,subboundariesUV,
      tolerance,PARM_TOLERANCE,false);

    saveTrimmedSurfacesIges(subsurfaces,subboundariesUV,DEBUG_DIR + "Kilo sub+fin+propeller surfaces trimmed.iges");

    std::vector<std::vector<TPoint<T>>> badedges;
    bool ok = saveSolidIges(subsurfaces,subboundariesUV,DEBUG_DIR + "Kilo sub+fin+propeller solid.iges",
      tolerance,PARM_TOLERANCE,SPLINE_DEGREE,18,&badedges);

    ASSERT(ok);

    if (!ok)
    {
      saveLinesIges<T>(badedges,DEBUG_DIR + "badedges.iges");
    }
  }

  /*****************************************************************************
    5.12 Blocks : submarine hump
  *****************************************************************************/

  cout << "5.12 Blocks : sub hump" << endl;

  // sub hump
  std::vector<TSplineSurface<T> *> subhumpsurfaces;
  // sub hump trimming curves
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> subhumpboundariesUV;

  {
    // base (upper) hump contour
    std::vector<TPoint<T>> upperpoints = KiloHumpBase<T>;
    smoothPointsByOrtho(upperpoints,END_ROUNDED,END_ROUNDED,8,GAUSSINT_8,true,50,0.5,0.5);

    // move it up to position of upper hump line
    TTransform<T> t;
    t.Resize(TPoint<T>(1.0,0.75,1.0));
    t.Translate(TPoint<T>(0.0,0.0,8.0 * KHCOEF));
    makeTransform(upperpoints,&t);

    TPointCurve<T> uppercurve(upperpoints);
    saveCurveIges(uppercurve,DEBUG_DIR + "Sub hump upper curve.iges");

    // this is lower contour
    std::vector<TPoint<T>> bottompoints = upperpoints;
    t.LoadIdentity();
    // widen hum to the bottom
 //   t.Resize(TPoint<T>(1.0,1.5,1.0));
    t.Translate(TPoint<T>(0.0,0.0,-6.0 * KHCOEF));
    makeTransform(bottompoints,&t);

    TPointCurve<T> bottomcurve(bottompoints);
    saveCurveIges(bottomcurve,DEBUG_DIR + "Sub hump bottom curve.iges");

    // top contour ("rounding top")
    T radius = KHCOEF;
    std::vector<TPoint<T>> toppoints = upperpoints;

    TPoint<T> min,max;
    calculateMinMax(toppoints,&min,&max);
    TPoint<T> d = max - min;
    TPoint<T> centre = (min + max) * 0.5;
    T xcoef = (d.X - radius * 2.0) / d.X;
    T ycoef = (d.Y - radius) / d.Y;

    t.LoadIdentity();
    t.Translate(-centre);
    t.Resize(TPoint<T>(xcoef,ycoef,1.0));
    t.Translate(centre + TPoint<T>(0.0,-radius * 0.5,radius));

    makeTransform(toppoints,&t);

    TPointCurve<T> topcurve(toppoints);
    saveCurveIges(topcurve,DEBUG_DIR + "Sub hump top curve.iges");

    // make surfaces now
    std::vector<std::vector<TPoint<T>>> points;
    points.push_back(bottompoints);
    points.push_back(upperpoints);
    points.push_back(toppoints);

    // make two symmetric (around XZ) walls and ceiling
    makeTwoWallsAndFlatBottomTop<T>(points,subhumpsurfaces,false,true,40,40,8,SPLINE_DEGREE, 
      END_CLAMPED,END_CLAMPED,END_FREE,END_FREE);

    saveSurfacesIges(subhumpsurfaces,DEBUG_DIR + "Sub hump surfaces.iges");
  }

  /*****************************************************************************
    5.13 Blocks : submarine hull + propeller + fin + hump, all single solid
  *****************************************************************************/

  cout << "5.13 Blocks : submarine hull + propeller + fin + hump, all single solid" << endl;

  {
    T hullL = 74.0;
    T tolerance = hullL * PARM_TOLERANCE;

    // name for debugging
    nameSurfaces<T>(subhumpsurfaces,"hump");

    // make boundaries
    closeOuterBoundary<T>(subhumpsurfaces,subhumpboundariesUV);

    // make a solid from hump and propeller + hull

    // make solid, do not clear old boundaries
    makeTrimming(subhumpsurfaces,subsurfaces,subhumpboundariesUV,subboundariesUV,UNITE,
      subsurfaces,subboundariesUV,
      tolerance,PARM_TOLERANCE,false);

    saveTrimmedSurfacesIges(subsurfaces,subboundariesUV,DEBUG_DIR + "Kilo sub+fin+hump+propeller surfaces trimmed.iges");

    std::vector<std::vector<TPoint<T>>> badedges;
    bool ok = saveSolidIges(subsurfaces,subboundariesUV,DEBUG_DIR + "Kilo sub+fin+hump+propeller solid.iges",
      tolerance,PARM_TOLERANCE,SPLINE_DEGREE,18,&badedges);

    ASSERT(ok);

    if (!ok)
    {
      saveLinesIges<T>(badedges,DEBUG_DIR + "badedges.iges");
    }
  }

#endif

  /*****************************************************************************
    5.14 Blocks : filleted wing with rudder cut out
  *****************************************************************************/

  cout << "5.14 Blocks : filleted wing with rudder cut out" << endl;

  // rudder
  TBrep<T> brudder;

  {
    T subL = 74.0;
    T tolerance = subL * PARM_TOLERANCE;

    //===== rudder =====

    // numbers of spline intervals along U and V
    int K1 = 40;
    int K2 = 20;

    // wing 1 : single stern bottom
    std::vector<TPoint<T>> upperlower; 
    std::pair<T,T> res = makeAirfoilPointsXY<T>(NACA0009<T>,true,false,K1,upperlower,true);

    std::vector<std::vector<TPoint<T>>> camberpoints; 
    std::vector<std::vector<T>> thickness;
    makeBladeCamberThickness<T>(upperlower,KiloWing1Contour<T>,K1,K2,camberpoints,thickness,tolerance); 

    // create camber surface
    TSplineSurface<T> cambersurface(camberpoints,K1,SPLINE_DEGREE,K2,SPLINE_DEGREE,
      END_CLAMPED,END_CLAMPED,END_FREE,END_FREE);

    saveSurfaceIges(&cambersurface,DEBUG_DIR + "Kilo wing 1 (rudder) camber surface.iges");

    std::vector<TSplineSurface<T> *> surfaces;
    makeAirfoil<T>(camberpoints,thickness,surfaces,SPLINE_DEGREE,SPLINE_DEGREE,
      END_CLAMPED,END_CLAMPED,END_FREE,END_FREE);

    makeFlatRoundedTop(surfaces[0],surfaces[1],0.1,surfaces);

    // rotate around
    TTransform<T> t;
    t.Rotate(TPoint<T>(1,0,0),180.0 * CPI);
    // move to place
    t.Translate(TPoint<T>(-84.0 * KHCOEF));

    makeTransform(surfaces,&t);

    nameSurfaces<T>(surfaces,"rudder");

    brudder = TBrep<T>(surfaces,tolerance);

    //==== rudder cutter : cylinder + box =====

    // contains two boxes + half cylinder with NO FACES INSIDE

    TBrep<T> bruddercut(tolerance);

    bruddercut.makeBox(TPoint<T>(-28.0,-0.3,-4.4),TPoint<T>(-24.0,0.0,-2.0));
    bruddercut.makeBox(TPoint<T>(-28.0,0.0,-4.4),TPoint<T>(-24.0,0.3,-2.0));

    bruddercut.deleteFace(8);
    bruddercut.deleteFace(7);
    bruddercut.deleteFace(3);
    bruddercut.deleteFace(1);

    // half of cylinder is part of cutter body
    TBrep<T> bcylinder(tolerance);

    // this is hollow shaft to hinge plate
    TBrep<T> bshaft0(tolerance);
    TBrep<T> bshaft1(tolerance);

    bcylinder.makeCylinder(2.4,0.3,0.3,"",4,4,1,1,64,MANY_POINTS2D,-90.0,90.0);

    bshaft0.makeCylinder(0.2,0.05,0.05,"",4,4,0,0,64,MANY_POINTS2D,0.0,360.0);
    bshaft1.makeCylinder(0.2,0.05,0.05,"",4,4,0,0,64,MANY_POINTS2D,0.0,360.0);

    // move half cylinder to stern
    t.LoadIdentity();
    t.Translate(TPoint<T>(-24.0,0.0,-3.2));
    bcylinder.makeTransform(&t);

    t.LoadIdentity();
    t.Translate(TPoint<T>(-24.0,0.0,-2.005));
    bshaft0.makeTransform(&t);

    t.LoadIdentity();
    t.Translate(TPoint<T>(-24.0,0.0,-4.416));
    bshaft1.makeTransform(&t);

    // two boxes + half cylinder
    bruddercut.addFaces(bcylinder);

    bruddercut.closeOuterBoundary();
    bruddercut.saveSurfacesIges(DEBUG_DIR + "Rudder cutter surfaces trimmed.iges");

    bshaft0.closeOuterBoundary();
    bshaft0.saveSurfacesIges(DEBUG_DIR + "Rudder shaft 0 surfaces trimmed.iges");
    bshaft1.closeOuterBoundary();
    bshaft1.saveSurfacesIges(DEBUG_DIR + "Rudder shaft 1 surfaces trimmed.iges");

    //===== Intersect rudder with cutter body to get rudder plane shape =====

    TBrep<T> bplane = brudder ^ bruddercut;

    // make it a bit smaller to fit the cut
    t.LoadIdentity();
    t.Translate(-TPoint<T>(-24.0,0.0,-3.2));
    t.Resize(TPoint<T>(0.95,0.95,0.95));
    t.Rotate(TPoint<T>(0.0,0.0,1.0),10.0 * CPI);
    t.Translate(TPoint<T>(-24.0,0.0,-3.2));
    bplane.makeTransform(&t);

    bplane.saveSurfacesIges(DEBUG_DIR + "Rudder plane surfaces trimmed.iges");

    //===== Subtract cutter body from rudder to make space for rudder plane =====

    brudder = brudder - bruddercut;

    brudder.clearNames();
    brudder.saveSurfacesIges(DEBUG_DIR + "Rudder surfaces trimmed.iges");

    bplane.nameSurfaces("plane");
    bshaft0.nameSurfaces("shaft0");
    bplane = bplane + bshaft0;
    bplane = bplane + bshaft1;

    bplane.saveSurfacesIges(DEBUG_DIR + "Rudder plane surfaces trimmed.iges");

    brudder.nameSurfaces("rudder");
    bplane.nameSurfaces("plane");
    brudder = brudder + bplane;

    brudder.saveSurfacesIges(DEBUG_DIR + "Rudder surfaces trimmed.iges");
  }

#if 1

  /*****************************************************************************
    5.15 Blocks : submarine with rudder
  *****************************************************************************/

  cout << "5.15 Blocks : sub with rudder" << endl;

  TBrep<T> bsub = TBrep<T>(subsurfaces,subboundariesUV);

  {
    // plus rudder
    bsub = bsub + brudder;

    // save surfaces
    bsub.saveSurfacesIges(DEBUG_DIR + "Kilo sub+fin+hump+propeller+rudder surfaces trimmed.iges");

    // save solid
    std::vector<std::vector<TPoint<T>>> badedges;
    bool ok = bsub.saveSolidIges(DEBUG_DIR + "Kilo sub+fin+hump+propeller+rudder solid.iges",
      bsub.tolerance,PARM_TOLERANCE,SPLINE_DEGREE,18,&badedges);

    ASSERT(ok);

    if (!ok)
    {
      saveLinesIges<T>(badedges,DEBUG_DIR + "badedges.iges");
    }
  }
#endif

#endif

#ifdef DEBUG_BOOLEANS

  /*****************************************************************************

    Part 6 : booleans

  *****************************************************************************/

  cout << "    Part 6 : BOOLEANS" << endl;

  /*****************************************************************************
    6.0 B-reps : booleans, box and sphere
  *****************************************************************************/

  cout << "6.0 B-reps : booleans, box and sphere" << endl;

  {
    T tolerance = 1.0 * PARM_TOLERANCE;

    // two box positions, second one is hard due to touching face borders
    for (int t = 0; t < 2; t++)
    {
      if (t == 0)
      {
        cout << "6.1 B-reps : booleans, box and sphere, union, subtraction and intersection" << endl;
      } else
      {
        cout << "6.2 B-reps : booleans, box and sphere, hard case (cutting along face edges)" << endl;
      }

      // loop on number of sphere faces
      for (int f = 0; f < 3; f++)
      {
        int numfaces = 2 << f;

        cout << "  " + to_string(numfaces) + " x " + to_string(numfaces * 2) + " sphere faces" << endl;

        //===== Box =====

        TBrep<T> box(tolerance);

        if (t == 0)
        {
          box.makeBox(TPoint<T>(-0.01,0.01,0.01),TPoint<T>(1.01,1.01,1.01));
        } else
        {
          box.makeBox(TPoint<T>(0.0,0.0,0.0),TPoint<T>(1.0,1.0,1.0));
        }

        bool ok0 = box.saveSurfacesIges(DEBUG_DIR + "Brep box surfaces.iges");

        ASSERT(ok0);

        std::vector<std::vector<TPoint<T>>> badedges;
        bool ok1 = box.saveSolidIges(DEBUG_DIR + "Brep box solid.iges",
          tolerance,PARM_TOLERANCE,SPLINE_DEGREE,18,&badedges);

        ASSERT(ok1);

        if (!ok1)
        {
          saveLinesIges<T>(badedges,DEBUG_DIR + "badedges.iges");
        }

        //===== Sphere =====

        TBrep<T> sphere(tolerance);
        sphere.makeSphere(0.5,"",numfaces);

        bool ok2 = sphere.saveSurfacesIges(DEBUG_DIR + "Brep sphere surfaces.iges");

        ASSERT(ok2);

        bool ok3 = sphere.saveSolidIges(DEBUG_DIR + "Brep sphere solid.iges",
          tolerance,PARM_TOLERANCE,SPLINE_DEGREE,18,&badedges);

        ASSERT(ok3);

        if (!ok3)
        {
          saveLinesIges<T>(badedges,DEBUG_DIR + "badedges.iges");
        }

        //===== Union =====

        cout << "box + sphere" << endl;

        TBrep<T> boxplussphere = box + sphere; 

        bool ok4 = boxplussphere.saveSurfacesIges(DEBUG_DIR + "Box+sphere surfaces.iges");

        ASSERT(ok4);

        bool ok5 = boxplussphere.saveSolidIges(DEBUG_DIR + "Box+sphere solid.iges",
          tolerance,PARM_TOLERANCE,SPLINE_DEGREE,18,&badedges); 

        if (!ok5)
        {
          saveLinesIges<T>(badedges,DEBUG_DIR + "badedges.iges");
        }

        ASSERT(ok5);

        //===== Subtraction =====

        {
          cout << "box - sphere" << endl;

          TBrep<T> boxminussphere = box - sphere; 

          bool ok6 = boxminussphere.saveSurfacesIges(DEBUG_DIR + "Box-sphere surfaces.iges");

          ASSERT(ok6);

          bool ok7 = boxminussphere.saveSolidIges(DEBUG_DIR + "Box-sphere solid.iges",
            tolerance,PARM_TOLERANCE,SPLINE_DEGREE,18,&badedges);

          if (!ok7)
          {
            saveLinesIges<T>(badedges,DEBUG_DIR + "badedges.iges");
          }

          ASSERT(ok7);
        }

        {
          cout << "sphere - box" << endl;

          TBrep<T> sphereminusbox = sphere - box;

          bool ok6 = sphereminusbox.saveSurfacesIges(DEBUG_DIR + "Sphere-box surfaces.iges");

          ASSERT(ok6);

          bool ok7 = sphereminusbox.saveSolidIges(DEBUG_DIR + "Sphere-box solid.iges",
            tolerance,PARM_TOLERANCE,SPLINE_DEGREE,18,&badedges);

          if (!ok7)
          {
            saveLinesIges<T>(badedges,DEBUG_DIR + "badedges.iges");
          }

          ASSERT(ok7);
        }

        //===== Intersection =====

        cout << "box ^ sphere" << endl;

        TBrep<T> boxintrsphere = box ^ sphere;

        bool ok8 = boxintrsphere.saveSurfacesIges(DEBUG_DIR + "Box^sphere surfaces.iges");

        ASSERT(ok8);

        bool ok9 = boxintrsphere.saveSolidIges(DEBUG_DIR + "Box^sphere solid.iges",
          tolerance,PARM_TOLERANCE,SPLINE_DEGREE,18,&badedges);

        if (!ok9)
        {
          saveLinesIges<T>(badedges,DEBUG_DIR + "badedges.iges");
        }

        ASSERT(ok9);
      } // # sphere faces
    } // box positions
  }

#endif

  double endtime = GetTime();

  cout << "Time elapsed " << (endtime - starttime) << " sec" << endl;

  return 0;
}