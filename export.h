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

/**
  Export into CAD files.
*/

#pragma once

#include "tpoint.h"
#include "tpoints.h"
#include "tsplinecurve.h"
#include "tsplinesurface.h"
#include "ttriangles.h"
#include "strings.h"

#include <vector>
#include <string>

// 13/10
#define CRLF "\r\n"

//===== Common =================================================================

/** Save lines into file. */
bool writeLines(const std::vector<std::string> &lines, const std::string &filename);

/** Write a line in header. */
std::string writeHeaderLine(int line, double val);

//===== IGES - specific ========================================================

// first part of IGES header
extern const char *IgesHeader0[18];
extern const char *IgesHeader1[3];
extern const char *dirline1260;
extern const char *dirline1260l;
extern const char *dirline1261;
extern const char *dirline1280;
extern const char *dirline1281;

/** Make Iges directory line. */
std::string makeIgesDirectoryLine0(const char *dirline, int parmdata, int parmcount, int *dirlineno);

/** Make Iges second directory line. */
std::string makeIgesDirectoryLine1(const char *dirline, int parmdata, int parmcount, int *dirlineno, const char *name);

/** Add string to fall inside positions 1..64 */
void addIgesString(std::vector<std::string> &lines, std::string substr, int dirline, int *count, std::string &igesstr);

/** Finalize lines. */
void finalize(std::vector<std::string> &lines, int dirline, int *count, std::string &igesstr);

/** Add lines with spline curve parameters. */
template <class T> void addIgesCurveLines(std::vector<std::string> &lines, tcad::TSplineCurve<T> *curve, 
  int dirline, int *count, int numdigits = 18)
{
  std::string igesstr = "";

  // save all info
  addIgesString(lines,std::to_string(126),dirline,count,igesstr);
  addIgesString(lines,std::to_string(curve->K1),dirline,count,igesstr);
  addIgesString(lines,std::to_string(curve->M1),dirline,count,igesstr);
  addIgesString(lines,std::to_string(0),dirline,count,igesstr); //this is difference with surface - one additional parm
  addIgesString(lines,std::to_string(0),dirline,count,igesstr);
  addIgesString(lines,std::to_string(1),dirline,count,igesstr); // polynomial
  addIgesString(lines,std::to_string(0),dirline,count,igesstr);

  int i;
  for (i = 0; i < curve->Uknots.size(); i++)
  {
    addIgesString(lines,to_string(double(curve->Uknots[i]),numdigits),dirline,count,igesstr);
  }
  for (i = 0; i < curve->K1 + 1; i++)
  {
    addIgesString(lines,to_string(1.0, 7),dirline,count,igesstr);
  }
  for (i = 0; i < curve->K1 + 1; i++)
  {
    addIgesString(lines,to_string(double(curve->controlPoints()[i].X),numdigits),dirline,count,igesstr);
    addIgesString(lines,to_string(double(curve->controlPoints()[i].Y),numdigits),dirline,count,igesstr);
    addIgesString(lines,to_string(double(curve->controlPoints()[i].Z),numdigits),dirline,count,igesstr);
  }
  addIgesString(lines,to_string(0.0,numdigits),dirline,count,igesstr);
  addIgesString(lines,to_string(1.0,numdigits),dirline,count,igesstr);

  finalize(lines,dirline,count,igesstr);
}

/** Add lines for a spline surface. */
template <class T> void addIgesSurfaceLines(std::vector<std::string> &lines, tcad::TSplineSurface<T> *surface, int dirline, int *count,
  int numdigits = 18)
{
  std::string igesstr = "";

  // save all info
  addIgesString(lines,std::to_string(128),dirline,count,igesstr);
  addIgesString(lines,std::to_string(surface->K1),dirline,count,igesstr);
  addIgesString(lines,std::to_string(surface->K2),dirline,count,igesstr);
  addIgesString(lines,std::to_string(surface->M1),dirline,count,igesstr);
  addIgesString(lines,std::to_string(surface->M2),dirline,count,igesstr);
  addIgesString(lines,std::to_string((int) 0),dirline,count,igesstr);
  addIgesString(lines,std::to_string((int) 0),dirline,count,igesstr);
  addIgesString(lines,std::to_string((int) 1),dirline,count,igesstr); // polynomial
  addIgesString(lines,std::to_string((int) 0),dirline,count,igesstr);
  addIgesString(lines,std::to_string((int) 0),dirline,count,igesstr);

  // knots
  for (int i = 0; i < (surface->K1 + surface->M1 + 2); i++)
  {
    addIgesString(lines,to_string(surface->Uknots[i],numdigits),dirline,count,igesstr);
  }
  for (int i = 0; i < (surface->K2 + surface->M2 + 2); i++)
  {
    addIgesString(lines,to_string(surface->Vknots[i],numdigits),dirline,count,igesstr);
  }

  // weights
  for (int i = 0; i < int(surface->controlPoints().size()); i++)
  {
    addIgesString(lines,to_string(1.0,numdigits),dirline,count,igesstr);
  }

  for (int i = 0; i < int(surface->controlPoints().size()); i++)
  {
    addIgesString(lines,to_string(surface->controlPoints()[i].X,numdigits),dirline,count,igesstr);
    addIgesString(lines,to_string(surface->controlPoints()[i].Y,numdigits),dirline,count,igesstr);
    addIgesString(lines,to_string(surface->controlPoints()[i].Z,numdigits),dirline,count,igesstr);
  }
  addIgesString(lines,to_string(0.0,numdigits),dirline,count,igesstr);
  addIgesString(lines,to_string(1.0,numdigits),dirline,count,igesstr);
  addIgesString(lines,to_string(0.0,numdigits),dirline,count,igesstr);
  addIgesString(lines,to_string(1.0,numdigits),dirline,count,igesstr);

  finalize(lines,dirline,count,igesstr);
}

/** Make the whole collection of lines to save them. */
template <class T> bool makeCurveLinesIges(std::vector<std::vector<tcad::TPoint<T>>> &curves, 
  std::vector<std::string> &lines, int splinedegree = 2, int minpoints = 6)
{
  lines.clear();

  if (curves.empty())
    return false;

  // calculate size
  double msize = 0.0;
  for (auto curve : curves)
  {
    double len = tcad::calculateLength(curve);
    msize = std::max<double>(msize,len);
  }

  // make header
  std::string s;
  for (int i = 0; i < 18; i++)
  {
    s = IgesHeader0[i];
    if (i == 7)
    {
      s = writeHeaderLine(7,TOLERANCE(double));
    } else if (i == 8)
    {
      s = writeHeaderLine(8,msize);
    }
    lines.push_back(s);
  }

  int slines = 1;
  int glines = 11;
  int dlines = (3 + int(curves.size())) * 2;

  // save double directory lines for each surface
  int dcount = 7;
  int dir126line = int(lines.size());
  for (int i = 0; i < int(curves.size()); i++)
  {
    // for entity 126
    s = makeIgesDirectoryLine0(dirline1260l, -1, -1, &dcount);
    lines.push_back(s);
    s = makeIgesDirectoryLine1(dirline1261, -1, -1, &dcount, "");
    lines.push_back(s);
  }

  // save 3 lines  
  // save header
  for (int i = 0; i < 3; i++)
  {
    s = IgesHeader1[i];
    lines.push_back(s);
  }

  // save surfaces
  dcount = 7;
  int pcount = 4;

  for (int i = 0; i < int(curves.size()); i++)
  {
    std::vector<tcad::TPoint<T>> &curve = curves[i];

    // redivide to make spline available
    if (curve.size() < minpoints) 
    {
      tcad::TPointCurve<T> pointcurve(curve);
      // there must be at least two non-identical points
      if (!pointcurve.ok())
        continue;

      pointcurve.redivide(minpoints);
      curve = pointcurve.controlPoints();
    }

    // modify directory
    lines[dir126line] = makeIgesDirectoryLine0(dirline1260l,pcount,-1,&dcount);
    dir126line++;

    tcad::TSplineCurve<T> scurve(curve,int(curve.size()) - 1,splinedegree); 

    // add lines
    int before = int(lines.size());
    addIgesCurveLines(lines,&scurve,dcount - 1,&pcount);
    int after = int(lines.size());

    // modify directory
    lines[dir126line] = makeIgesDirectoryLine1(dirline1261, -1,(after - before),&dcount,"");
    dir126line++;
  }

  int plines = pcount - 1;
  int tlines = 1;

  s = std::string("S") + to_string(slines,7) +
    std::string("G") + to_string(glines,7) +
    std::string("D") + to_string(dlines,7) +
    std::string("P") + to_string(plines,7) +
    "                                        " +
    std::string("T") + to_string(tlines,7) + CRLF;
  lines.push_back(s);

  return true;
}

/** Make the whole collection of lines to save them. */
template <class T> bool makeSurfaceLinesIges(std::vector<tcad::TSplineSurface<T> *> &surfaces, 
  std::vector<std::string> &lines)
{
  lines.clear();

  if (surfaces.empty())
    return false;

  // calculate size
  double msize = 0.0;
  for (auto &surface : surfaces)
  {
    double len = surface->maxSize();
    msize = std::max<double>(msize,len);
  }

  // make header
  std::string s;
  for (int i = 0; i < 18; i++)
  {
    s = IgesHeader0[i];
    if (i == 7)
    {
      s = writeHeaderLine(7,TOLERANCE(double));
    } else if (i == 8)
    {
      s = writeHeaderLine(8,msize);
    }
    lines.push_back(s);
  }

  int slines = 1;
  int glines = 11;
  int dlines = (3 + int(surfaces.size())) * 2;

  // save REAL directory lines for each surface
  int dcount = 7;
  int dir128line = int(lines.size());
  for (int i = 0; i < int(surfaces.size()); i++)
  {
    // for entity 128
    s = makeIgesDirectoryLine0(dirline1280, -1, -1, &dcount);
    lines.push_back(s);
    s = makeIgesDirectoryLine1(dirline1281, -1, -1, &dcount, "");
    lines.push_back(s);
  }

  // save 3 lines  
  // save header
  for (int i = 0; i < 3; i++)
  {
    s = IgesHeader1[i];
    lines.push_back(s);
  }

  // save surfaces
  dcount = 7;
  int pcount = 4;

  for (int i = 0; i < int(surfaces.size()); i++)
  {
    tcad::TSplineSurface<T> *surface = surfaces[i];

    // modify directory
    lines[dir128line] = makeIgesDirectoryLine0(dirline1280,pcount,-1,&dcount);
    dir128line++;

    // add lines
    int before = int(lines.size());
    addIgesSurfaceLines(lines,surface,dcount - 1,&pcount);
    int after = int(lines.size());

    // modify directory
    lines[dir128line] = makeIgesDirectoryLine1(dirline1281,-1,(after - before),&dcount,"");
    dir128line++;
  }

  int plines = pcount - 1;
  int tlines = 1;

  s = std::string("S") + to_string(slines,7) +
    std::string("G") + to_string(glines,7) +
    std::string("D") + to_string(dlines,7) +
    std::string("P") + to_string(plines,7) +
    "                                        " +
    std::string("T") + to_string(tlines,7) + CRLF;
  lines.push_back(s);

  return true;
}


//===== Curves =================================================================

/** Save curves as points. */
template <class T> bool saveLinesIges(std::vector<std::vector<tcad::TPoint<T>>> &curves, const std::string &filename)
{
  std::vector<std::string> lines;

  if (makeCurveLinesIges(curves,lines))
  {
    bool ok = writeLines(lines,filename);
    return ok;
  } else
  {
    return false;
  }
}

/** Save curve as points. */
template <class T> bool saveLinesIges(std::vector<tcad::TPoint<T>> &curve, const std::string &filename)
{
  std::vector<std::vector<tcad::TPoint<T>>> curves;
  curves.push_back(curve);

  return saveLinesIges(curves,filename);
}

/** Save any curve (based on TBaseCurve) as points. */
template <class T> bool saveCurveIges(tcad::TBaseCurve<T> &curve, const std::string &filename)
{
  std::vector<tcad::TPoint<T>> points;
  curve.createPoints(points);

  return saveLinesIges(points,filename);
}

/** Save any TWO curve (based on TBaseCurve) as points e,g, to compare. */
template <class T> bool saveTwoCurvesIges(tcad::TBaseCurve<T> &curve0, tcad::TBaseCurve<T> &curve1, 
  const std::string &filename)
{
  std::vector<tcad::TPoint<T>> points0,points1;
  curve0.createPoints(points0);
  curve1.createPoints(points1);

  std::vector<std::vector<tcad::TPoint<T>>> points;
  points.push_back(points0);
  points.push_back(points1);

  return saveLinesIges(points,filename);
}

/** Save any THREE curve (based on TBaseCurve) as points e,g, to compare. */
template <class T> bool saveThreeCurvesIges(tcad::TBaseCurve<T> &curve0, tcad::TBaseCurve<T> &curve1, 
  tcad::TBaseCurve<T> &curve2, const std::string &filename)
{
  std::vector<tcad::TPoint<T>> points0,points1,points2;
  curve0.createPoints(points0);
  curve1.createPoints(points1);
  curve2.createPoints(points2);

  std::vector<std::vector<tcad::TPoint<T>>> points;
  points.push_back(points0);
  points.push_back(points1);
  points.push_back(points2);

  return saveLinesIges(points,filename);
}

//===== Surfaces ===============================================================

/** Save surfaces in IGES. */
template <class T> bool saveSurfacesIges(std::vector<tcad::TSplineSurface<T> *> &surfaces, const std::string &filename)
{
  std::vector<std::string> lines;

  if (makeSurfaceLinesIges(surfaces,lines))
  {
    bool ok = writeLines(lines,filename);
    return ok;
  } else
  {
    return false;
  }
}

/** Save surface in IGES. */
template <class T> bool saveSurfaceIges(tcad::TSplineSurface<T> *surface, const std::string &filename)
{
  std::vector<tcad::TSplineSurface<T> *> surfaces;
  surfaces.push_back(surface);

  return saveSurfacesIges(surfaces,filename);
}

//===== STL ====================================================================

template <class T> bool saveTrianglesStl(tcad::TTriangles<T> &triangles, const std::string &filename, 
  const std::string &partname = "TCAD", const bool binary = true)
{
  return triangles.saveSTL(filename,partname,binary);
}

//===== OBJ ====================================================================

template <class T> bool saveTrianglesObj(tcad::TTriangles<T> &triangles, const std::string &filename, 
  const std::string &partname = "TCAD")
{
  return triangles.saveOBJ(filename,partname);
}



//===== Step-specific ==========================================================

