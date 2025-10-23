
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

// this stuff is for export and debugging
#include "strings.h"
#include "export.h"

using namespace tcad;

// save debugging files here, end this with slash "/"
#define DEBUG_DIR std::string("")

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


//===== TCAD examples ==========================================================

int main(int argc, char* argv[])
{

  /*****************************************************************************

    Part 1 : Curves. What can we do with them?

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
    1.1 Curves : compare points and low-power LSQ segment, do not expect
      a good correspondence - the lsq is a simple poly of power 3
  *****************************************************************************/

  TLSQSegment<T> lsqsegment(points,3);
  // save to compare two curves
  saveTwoCurvesIges(curve,lsqsegment,DEBUG_DIR + "Compare_points_and_lsq_segment.iges");

  /*****************************************************************************
    1.2 Curves : compare points and Bezier segment, do not expect
      a good correspondence - Bezier segment is qubic
  *****************************************************************************/

  TBezierSegment<T> beziersegment(points,END_FIXED,END_FIXED);
  // save to compare two curves
  saveTwoCurvesIges(curve,beziersegment,DEBUG_DIR + "Compare_points_and_bezier_segment.iges");

  /*****************************************************************************
    1.3 Curves : compare points and segment of orthogonal polynomials of
      degree 4. "beta" 0.5 is set to handle the round LE.
  *****************************************************************************/

  TOrthoSegment<T> orthosegment(points,END_ROUNDED,END_UNCORRECTED,8,GAUSSINT_20);
  // save to compare two curves
  saveTwoCurvesIges(curve,orthosegment,DEBUG_DIR + "Compare_points_and_orthopoly_segment.iges");

  /*****************************************************************************
    1.4 Curves : compare points and Bezier curve.
  *****************************************************************************/

  TBezierCurve<T> beziercurve(points,10,END_FIXED,END_FIXED);
  // save to compare two curves
  saveTwoCurvesIges(curve,beziercurve,DEBUG_DIR + "Compare_points_and_bezier_curve.iges");

  /*****************************************************************************
    1.5 Curves : compare points and b-spline.
  *****************************************************************************/

  TSplineCurve<T> splinecurve(points,int(points.size() - 1),3,END_CLAMPED,END_CLAMPED); 
  // save to compare two curves
  saveTwoCurvesIges(curve,splinecurve,DEBUG_DIR + "Compare_points_and_spline_curve.iges");

  /*****************************************************************************
    1.6 Curves : compare points and b-spline of degree 2. You can also select 
    less b-spline control points (10 here).
  *****************************************************************************/

  TSplineCurve<T> splinecurve10(points,10,3,END_CLAMPED,END_CLAMPED); 
  // save to compare two curves
  saveTwoCurvesIges(curve,splinecurve10,DEBUG_DIR + "Compare_points_and_spline_curve_K10.iges");

  /*****************************************************************************
    1.7 Curves : how to smooth a curve (curve2 below).
  *****************************************************************************/

  // make random noise, points2 is a spoilt set
  std::vector<TPoint<T>> points2;
  makeRandomNoise(points,points2,len * 0.005);

  // do not spoil points near start/end : apply "hat" function to diminish
  // the changes near ends with maximum (1.0) at the middle = 0.5
  // and set poly power as 0.5
  std::vector<T> coefs;
  makeHat(int(points2.size()),coefs,0.5,0.5);

  points2 = points + (points2 - points) * coefs;

  //===== smooth by orthogonal polynomials =====================================
  // these will be smoothed points
  std::vector<TPoint<T>> points3 = points2;
  smoothPointsByOrtho(points3,END_ROUNDED,END_UNCORRECTED);

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
    1.8 Curves : how to find parameter U for a point on curve.
  *****************************************************************************/

  // take a random point
  int index = random(int(points.size()));
  // we want a parameter for this point
  TPoint<T> p = points[index];

  // take spline curve of 10 intervals
  T U = splinecurve10.findUforPoint(p);
  // check this value by comparison with calculated position on the curve
  TPoint<T> pactual = splinecurve10.derivative(U,0);
  T diff = !(pactual - p);

  assert(diff < 0.001);


  return 0;
}