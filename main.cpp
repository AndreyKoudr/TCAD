
/*******************************************************************************

  Templated 3D CAD

  All the code is based on STL (standard templates, not stereolithography);
  it is templated what means that it is easy to use, like std::vector, you do 
  not need dlls, all is in front of you, easy for debugging and well commented.

  tbasics.h         - basics
    MATHEMATICS :
      tmisc         - miscellaneous
      tplane        - plane
      tmatrix       - matrix
      ttransform    - coordinate transformations
      tsystems      - systems of linear equations
      tvalues       - list of values std::vector<T>
      tpoints       - list of 3D points std::vector<std::vector<T>>
      tlsqfitting   - LSQ fitting with regular poly of low power
      tjacobipoly   - fitting with orthogonal polynomials
      tgaussint     - Gauss/Hammer numerical integration
    3D OBJECTS :
      tpoint            - 3-dimensional (0 parametric dimensions) point
      tbasecurve        - 3-dimensional (1 parametric dimension) non-composite curve (abstract class)
        tlinearsegment      - 2-point linear segment
        tlsqsegment         - lsq segment of low poly power
        torthosegment       - segment based on Jacobi ortho poly
        tbeziersegment      - Bezier segment
        tpointcurve         - piecewise linear curve of points
        tsplinecurve        - B-spline curve

  Auxilliary for debugging and export
    strings.h/cpp       - strings
    export.h/cpp        - output curve/surface/solid to CAD file

*******************************************************************************/

// you can only use high-level operations
#include "toperations.h"

#include "ttriangles.h"
#include "tbasesurface.h" //!!!!!!!
#include "tbezierpatch.h" //!!!!!!!

// this stuff is for export and debugging
#include "strings.h"
#include "export.h"
#include "import.h"

using namespace tcad;

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

/** Randomly swap points. */
void makeRandomSwap(int numswaps, std::vector<TPoint<T>> &points, std::vector<TPoint<T>> &spoiltpoints)
{
  spoiltpoints = points;

  // randomly swap
  for (int i = 0; i < numswaps; i++)
  {
    int index0 = random(int(spoiltpoints.size()));
    int index1 = random(int(spoiltpoints.size()));

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
  assert(manifold == mustbemanifold);

  // is solid?
  bool solid = tris.solid(tolerance,badedges);
  assert(solid == mustbesolid);

  // get boundary
  std::vector<std::vector<TPoint<T>>> boundary;
  bool bok = tris.getBoundary(boundary,tolerance);

  // output boundary, the should be none if solid
  if (solid)
  {
    assert (!bok);
  } else
  {
    assert(bok);

    // redivide point curves by sharp corners to make it look as in the 
    // original STL
    std::vector<std::vector<TPoint<T>>> bpoints;
    redividePoints(boundary,bpoints,tolerance,sharpangledeg);

    saveLinesIges(bpoints,DEBUG_DIR + name + " boundary.iges");
  }

  // cut it by plane
  std::vector<std::vector<TPoint<T>>> cutpoints;
  tris.cutByPlane(plane,cutpoints,tolerance); 

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

  /*****************************************************************************

    Part 1 : CURVES. What can we do with them?

  *****************************************************************************/

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

  saveCurveIges(curve,DEBUG_DIR + "NACA.iges");

  /*****************************************************************************
    1.1 Curves : how to exclude neighbour duplicate nodes
  *****************************************************************************/

  // make a copy
  std::vector<Point> points1 = points;

  // calculate curve length and tolerance to make node duplicates
  T len = calculateLength(points);
  T tolerance = len * 0.00001;

  // spoil points by inserting node duplicates
  makeRandomNeighbourDuplicates(points,points1,tolerance);

  // do not sort coordinates (false), exclude only neighbours
  bool ok = removeDuplicates(points1,false,tolerance * 10.0);

  TPointCurve<T> curve1(points1);

  // save to compare two curves
  saveTwoCurvesIges(curve,curve1,DEBUG_DIR + "Compare_removed_duplicates.iges");

  /*****************************************************************************
    1.2 Curves : compare original points and low-power LSQ segment, do not expect
      a good correspondence - the lsq is a simple poly of power 3
  *****************************************************************************/

  TLSQSegment<T> lsqsegment(points,3);
  // save to compare two curves
  saveTwoCurvesIges(curve,lsqsegment,DEBUG_DIR + "Compare_points_and_lsq_segment.iges");

  /*****************************************************************************
    1.3 Curves : compare points and Bezier segment, do not expect
      a good correspondence - Bezier segment is qubic
  *****************************************************************************/

  TBezierSegment<T> beziersegment(points,END_FIXED,END_FIXED);
  // save to compare two curves
  saveTwoCurvesIges(curve,beziersegment,DEBUG_DIR + "Compare_points_and_bezier_segment.iges");

  /*****************************************************************************
    1.4 Curves : compare points and segment of orthogonal polynomials of
      degree 4. Starting END_ROUNDED end is set to handle the round LE.
  *****************************************************************************/

  TOrthoSegment<T> orthosegment(points,END_ROUNDED,END_FREE,8,GAUSSINT_20);
  // save to compare two curves
  saveTwoCurvesIges(curve,orthosegment,DEBUG_DIR + "Compare_points_and_orthopoly_segment.iges");

  /*****************************************************************************
    1.5 Curves : compare points and Bezier curve (collection of Bezier segments).
  *****************************************************************************/

  TBezierCurve<T> beziercurve(points,10,END_FIXED,END_FIXED);
  // save to compare two curves
  saveTwoCurvesIges(curve,beziercurve,DEBUG_DIR + "Compare_points_and_bezier_curve.iges");

  /*****************************************************************************
    1.6 Curves : compare points and approximated b-spline.
  *****************************************************************************/

  TSplineCurve<T> splinecurve(points,int(points.size() - 1),SPLINE_DEGREE,END_CLAMPED,END_CLAMPED); 
  // save to compare two curves
  saveTwoCurvesIges(curve,splinecurve,DEBUG_DIR + "Compare_points_and_spline_curve.iges");

  /*****************************************************************************
    1.7 Curves : compare points and approximated b-spline with less number 
    of control points (10 here).
  *****************************************************************************/

  TSplineCurve<T> splinecurve10(points,10,SPLINE_DEGREE,END_CLAMPED,END_CLAMPED); 
  // save to compare two curves
  saveTwoCurvesIges(curve,splinecurve10,DEBUG_DIR + "Compare_points_and_spline_curve_K10.iges");

  /*****************************************************************************
    1.8 Curves : how to smooth a curve (curve2 below) with variuos approximants.
  *****************************************************************************/

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
  smoothPointsBySplineCurve(points6,10,END_CLAMPED,END_CLAMPED);

  TPointCurve<T> curve6(points6);

  // compare curves : spoilt and smoothed
  saveTwoCurvesIges(curve2,curve6,DEBUG_DIR + "Compare_points_splinecurvesmooth.iges");

  /*****************************************************************************
    1.9 Curves : how to diminish effect of smoothing near curve ends (much 
    needed when handling surface control points)
  *****************************************************************************/

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
    1.10 Curves : how to find parameter U for a point on curve.
  *****************************************************************************/

  // take a random point
  int index = random(int(points.size()));
  // we want a parameter for this point
  TPoint<T> p = points[index];

  // take spline curve of 10 intervals
  T U = splinecurve.findUforPoint(p);
  // check this value by comparison with calculated position on the curve
  TPoint<T> pactual = splinecurve.derivative(U,0);
  T diff = !(pactual - p);

  assert(diff < 0.005);

  // take Bezier curve
  U = beziercurve.findUforPoint(p);
  // check this value by comparison with calculated position on the curve
  pactual = beziercurve.derivative(U,0);
  diff = !(pactual - p);

  assert(diff < 0.005);

  // take ortho segment
  U = orthosegment.findUforPoint(p);
  // check this value by comparison with calculated position on the curve
  pactual = orthosegment.derivative(U,0);
  diff = !(pactual - p);

  assert(diff < 0.005);

  /*****************************************************************************
    1.11 Curves : how to intersect by plane.
  *****************************************************************************/

  // horizontal plane at Y = 0.03
  TPlane<T> plane(TPoint<T>(0.0,0.03,0.0),TPoint<T>(1.0,0.03,0.0),TPoint<T>(1.0,0.03,1.0),ok);

  // two intersections are expected
  std::vector<T> Upoints;
  int numintrs = splinecurve.intersectByPlane(plane,Upoints,tolerance);

  assert(numintrs == 2);

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
    assert(dist < tolerance);
  }

  /*****************************************************************************
    1.12 Curves : how to intersect two curves.
  *****************************************************************************/

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
  numintrs = curve.findIntersections(ecurve,UV,tolerance); 

  assert(numintrs == 2);

  // check coincidence
  for (int i = 0; i < int(UV.size()); i++)
  {
    TPoint<T> p0 = curve.derivative(UV[i].X,0);
    TPoint<T> p1 = ecurve.derivative(UV[i].Y,0);
    T dist = !(p1 - p0);
    assert(dist < tolerance);
  }

  /*****************************************************************************
    1.13 Curves : cut a piece of curve from the previous example
  *****************************************************************************/

  // curve was cut from UV[0].X to UV[1].X
  std::vector<TPoint<T>> cutpoints;
  curve.cutPiece(51,UV[0].X,UV[1].X,cutpoints);

  // we order points here
  TSplineCurve<T> cutcurve(cutpoints,40,SPLINE_DEGREE);

  // compare curves : unordered and ordered
  saveCurveIges(cutcurve,DEBUG_DIR + "Piece of curve.iges");

  /*****************************************************************************
    1.14 Curves : order unordered points.
  *****************************************************************************/

  std::vector<TPoint<T>> unorderedpoints;
  makeRandomSwap(int(points.size()) / 2,points,unorderedpoints);

  // unordered curve
  TPointCurve<T> unorderedcurve(unorderedpoints);

  // we order points here
  TPointCurve<T> orderedcurve(unorderedpoints);
  orderedcurve.order(tolerance);

  // compare curves : unordered and ordered
  saveTwoCurvesIges(unorderedcurve,orderedcurve,DEBUG_DIR + "Unordered and ordered curves.iges");

  /*****************************************************************************

    Part 2 : SURFACES

  *****************************************************************************/

  /*****************************************************************************
    2.1 Surfaces : triangles : save in STL
  *****************************************************************************/

  TTriangles<T> NACAtris;
  NACAtris.makeNACA0012(50,10,2.0);

  saveTrianglesStl(NACAtris,DEBUG_DIR + "NACA.stl");

  /*****************************************************************************
    2.2 Surfaces : triangles : manifold? solid?
  *****************************************************************************/

  // is manifold?
  std::vector<std::pair<LINT,LINT>> badedges;
  bool manifold = NACAtris.manifold(tolerance,badedges);
  assert(manifold);

  // is solid?
  bool solid = NACAtris.solid(tolerance,badedges);
  assert(solid);

  /*****************************************************************************
    2.3 Surfaces : triangles : cut by plane
  *****************************************************************************/

  // cut triangles by a plane (made by three points)
  TPlane<T> cutplane(TPoint<T>(-0.5,0.0,0.9),TPoint<T>(0.5,0.0,-0.7),TPoint<T>(0.5,1.0,-0.7),ok); 
  std::vector<std::vector<TPoint<T>>> NACAcuts;
  NACAtris.cutByPlane(cutplane,NACAcuts,tolerance); 

  /*****************************************************************************
    2.4 Surfaces : triangles : how to find a sharpest point in intersection curve
  *****************************************************************************/

  if (NACAcuts.size() == 1)
  {
    std::vector<TPoint<T>> NACAcutpoints = NACAcuts[0];
  
    // the curve (NACAcutpoints) starts from an arbitrary point after intersection, 
    // let's make it start from a point with highest curvature like TE in airfoils
    TPointCurve<T> temp(NACAcutpoints);
    int sharpindex = temp.findPointOfMaxCurvature();
    temp.shiftClosed(sharpindex,tolerance);
    NACAcutpoints = temp.controlPoints();

    // make, for example, a spline curve on these points
    TSplineCurve<T> NACAcutcurve(NACAcutpoints,200,SPLINE_DEGREE,END_CLAMPED,END_CLAMPED);

    saveCurveIges(NACAcutcurve,DEBUG_DIR + "NACA cut curve.iges");
  } else
  {
    assert(false);
  }

  /*****************************************************************************
    2.5 Surfaces : triangles : how to get triangulation boundary
  *****************************************************************************/

  // NACA case is solid, no boundary points are expected
  std::vector<std::vector<TPoint<T>>> NACAtrisboundary;
  bool bok = NACAtris.getBoundary(NACAtrisboundary,tolerance);

  assert(!bok);

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

  // redivide point curves by sharp corners to make it look good
  std::vector<std::vector<TPoint<T>>> fnewpoints;
  redividePoints(fboundary,fnewpoints,ftolerance,10.0);

  saveLinesIges(fnewpoints,DEBUG_DIR + "fuselage boundary.iges");

  /*****************************************************************************
    2.6 Surfaces : triangles : more complicated cases, load STL, check if 
    manifold and solid, generate boundary and cut by plane
  *****************************************************************************/

  // cut it by plane, define by normal and one point on plane
  TPlane<T> scutplane(TPoint<T>(1.0,0.0,0.0),TPoint<T>(5.8,0.0,0.0)); 

  bool sok1 = checkTopoCutAndBoundary("shuttle",scutplane,true,true,10.0);
  assert(sok1);


  rewriteSTLAsBinary("wing.stl");

  // cut it by plane, define by normal and one point on plane
  TPlane<T> wcutplane(TPoint<T>(0.0,0.0,1.0),TPoint<T>(0.0,0.0,0.0)); 

  bool wok1 = checkTopoCutAndBoundary("wing",wcutplane,true,true,45.0);
  assert(wok1);

  /*****************************************************************************
    2.7 Surfaces : triangles : find intersection of one set of triangles with 
    another
  *****************************************************************************/

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
  bool iok = fuselage.intersect(fuselage1,fcutpoints,ftolerance); 

  assert(iok);

  // redivide point curves by sharp corners to make it look good
  std::vector<std::vector<TPoint<T>>> fcpoints;
  // only 1 degree for sharp edges is here to divide the curve into
  // straight-line segments
  redividePoints(fcutpoints,fcpoints,tolerance,1.0);

  saveLinesIges(fcpoints,DEBUG_DIR + "fuselage-fuselage intersection curve.iges");

  /*****************************************************************************
    2.8 Surfaces : Bezier patch
  *****************************************************************************/

  TBezierSegment<T> SU0(TPoint<T>(0,0,0),TPoint<T>(1,0.0,0.5),TPoint<T>(2,0,0),TPoint<T>(3,0,0));
  TBezierSegment<T> S1V(TPoint<T>(3,0,0),TPoint<T>(3,1,0),TPoint<T>(3,2,0),TPoint<T>(3,3,0));
  TBezierSegment<T> SU1(TPoint<T>(0,3,0),TPoint<T>(1,3,0),TPoint<T>(2,3,0),TPoint<T>(3,3,0));
  TBezierSegment<T> S0V(TPoint<T>(0,0,0),TPoint<T>(0,1,0),TPoint<T>(0,2,0),TPoint<T>(0,3,0));

  TBezierPatch<T> bpatch(&SU0,&S1V,&SU1,&S0V);

  TTriangles<T> btris;
  bpatch.createTriangles(btris,100,100); //!!!!!!!

  // save tranformed
  saveTrianglesStl(btris,DEBUG_DIR + "bezier patch.stl");

  return 0;
}