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
#include "tsolid.h"
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
extern const char *dirline1410;
extern const char *dirline1411;
extern const char *dirline1430;
extern const char *dirline1431;
extern const char *dirline1260parametric;
extern const char *dirline1260;
extern const char *dirline1261;

extern const char *sdirline5020;
extern const char *sdirline5021;
extern const char *sdirline5140;
extern const char *sdirline5141;
extern const char *sdirline1260;
extern const char *sdirline1261;
extern const char *sdirline5040;
extern const char *sdirline5041;
extern const char *sdirline1280;
extern const char *sdirline1281;
extern const char *sdirline5080;
extern const char *sdirline5081;
extern const char *sdirline5100;
extern const char *sdirline5101;
extern const char *sdirline1860;
extern const char *sdirline1860final;
extern const char *sdirline1861;


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

/** Entity 141 : add to lines. */
template <class T> void addIges141(std::vector<std::string> &lines, int surfaceDE, int curveDE1,
  int curveDE2, int dirline, int *count, int numcurves = 1)
{
  std::string igesstr = "";

  // save all info
  addIgesString(lines,to_string(141),dirline,count,igesstr);
  addIgesString(lines,to_string(1),dirline,count,igesstr);

  // preferences : 0 - unspecified, 1 - model, 2 - parameter, 3 - equal pref
  addIgesString(lines,to_string(1),dirline,count,igesstr);

  addIgesString(lines,to_string(surfaceDE),dirline,count,igesstr);

  addIgesString(lines,to_string(numcurves),dirline,count,igesstr);

  for (int i = 0; i < numcurves; i++)
  {
    addIgesString(lines,to_string(curveDE2 + i * 2),dirline,count,igesstr);
    addIgesString(lines,to_string(1),dirline,count,igesstr);
    addIgesString(lines,to_string(1),dirline,count,igesstr);
    addIgesString(lines,to_string(curveDE1 + i * 2),dirline,count,igesstr);
  }

  finalize(lines,dirline,count,igesstr);
}

/** Entity 143 : add to lines. */
template <class T> void addIges143(std::vector<std::string> &lines, int surfaceDE, 
  std::vector<int> curve141DEs, int dirline, int *count)
{
  std::string igesstr = "";

  addIgesString(lines,to_string(143),dirline,count,igesstr);
  addIgesString(lines,to_string(1),dirline,count,igesstr);
  addIgesString(lines,to_string(surfaceDE),dirline,count,igesstr);
  addIgesString(lines,to_string(int(curve141DEs.size())),dirline,count,igesstr);
  for (int i = 0; i < int(curve141DEs.size()); i++)
  {
    addIgesString(lines,to_string(curve141DEs[i]),dirline,count,igesstr);
  }

  finalize(lines,dirline,count,igesstr);
}

/** Entity 502. */
template <class T> void addIges502(std::vector<std::string> &lines, std::vector<TPoint<T>> &vertices, int dirline, 
  int *count, int numdigits)
{
  std::string igesstr = "";

  // save all info
  addIgesString(lines,to_string(502),dirline,count,igesstr);
  addIgesString(lines,to_string((int) vertices.size()),dirline,count,igesstr);

  for (int i = 0; i < vertices.size(); i++)
  {
    addIgesString(lines,to_string(vertices[i].X,numdigits),dirline,count,igesstr);
    addIgesString(lines,to_string(vertices[i].Y,numdigits),dirline,count,igesstr);
    addIgesString(lines,to_string(vertices[i].Z,numdigits),dirline,count,igesstr);
  }

  finalize(lines,dirline,count,igesstr);
}

/** Entity 504. */
template <class T> void addIges504(std::vector<std::string> &lines, 
  std::vector<std::array<LINT,11>> &edges, 
  int verticesDE, std::vector<int> &boundaryDEs, int dirline, int *count)
{
  std::string igesstr = "";

  // save all info
  addIgesString(lines,to_string(504),dirline,count,igesstr);
  addIgesString(lines,to_string((int) edges.size()),dirline,count,igesstr);

  int count1 = 0;
  for (auto &e : edges)
  {
    addIgesString(lines,to_string(boundaryDEs[count1]),dirline,count,igesstr);
    addIgesString(lines,to_string(verticesDE),dirline,count,igesstr);
    addIgesString(lines,to_string(e[0] + 1),dirline,count,igesstr);
    addIgesString(lines,to_string(verticesDE),dirline,count,igesstr);
    addIgesString(lines,to_string(e[1] + 1),dirline,count,igesstr);
    count1++;
  }

  finalize(lines,dirline,count,igesstr);
}

/** Make the whole collection of lines to save them. */
template <class T> bool makeCurveLinesIges(std::vector<std::vector<tcad::TPoint<T>>> &curves, 
  std::vector<std::string> &lines, int splinedegree = SPLINE_DEGREE, int minpoints = 6)
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

/** Save curves as points. */
template <class T> bool saveLinesIges(std::vector<std::vector<tcad::TPoint<T>>> &curves, 
  const std::string &filename, int splinedegree = SPLINE_DEGREE, int minpoints = 6)
{
  std::vector<std::string> lines;

  if (makeCurveLinesIges(curves,lines,splinedegree,minpoints))
  {
    bool ok = writeLines(lines,filename);
    return ok;
  } else
  {
    return false;
  }
}

/** Entity 508. */
template <class T> void addIges508(std::vector<tcad::TSplineSurface<T> *> &surfaces, 
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV,
  std::vector<TPoint<T>> &vertices, std::vector<TPoint<T>> &middlevertices, T tolerance,
  std::vector<std::string> &lines, 
  std::vector<std::array<LINT,11>> &edges, 
  int surface, int loop, int loopsize, int edgeDE, 
  int dirline, int *count)
{
  std::string igesstr = "";

  // save all info
  addIgesString(lines,to_string(508),dirline,count,igesstr);

  // some boundary pieces can be split
  addIgesString(lines,to_string(loopsize),dirline,count,igesstr);

  // find edge numbers for this surface and this loop
  std::vector<std::pair<int,bool>> iedges;

#if 1 //!!!!!!!
  std::vector<std::vector<TPoint<T>>> dloop0,dloop1;
  loopOK(surfaces,boundariesUV,vertices,middlevertices,edges,surface,loop,loopsize,
    iedges,tolerance,dloop0,dloop1);

  #ifdef DEBUG_SOLID
    saveLinesIges<T>(dloop0,std::string("loop0_") + to_string(surface) + "_" + to_string(loop) + ".iges");
    saveLinesIges<T>(dloop1,std::string("loop1_") + to_string(surface) + "_" + to_string(loop) + ".iges");
  #endif
#else
  findLoopEdges<T>(edges,surface,loop,loopsize,iedges);
#endif;

  assert(!iedges.empty());

  // it should be ok even with degenerated edges as see the code :
  // "// remove this degenerated piece of boundary, only ONE!"
  assert(loopsize == iedges.size());

  for (int i = 0; i < iedges.size(); i++)
  {
    int iedge = iedges[i].first;
    bool reversed = iedges[i].second;

    addIgesString(lines,to_string(0),dirline,count,igesstr);
    addIgesString(lines,to_string(edgeDE),dirline,count,igesstr);

    addIgesString(lines,to_string(iedge + 1),dirline,count,igesstr);
    addIgesString(lines,to_string((int) (!reversed)),dirline,count,igesstr);

    addIgesString(lines,to_string((int) 0),dirline,count,igesstr);
  }

  finalize(lines,dirline,count,igesstr);
}

/** Entity 510. */
template <class T> void addIges510(std::vector<std::string> &lines, int surfaceDE, 
  std::vector<int> &loopDEs, int dirline, int *count)
{
  std::string igesstr = "";

  addIgesString(lines,to_string(510),dirline,count,igesstr);
  addIgesString(lines,to_string(surfaceDE),dirline,count,igesstr);
  addIgesString(lines,to_string(int(loopDEs.size())),dirline,count,igesstr);
  addIgesString(lines,to_string(1),dirline,count,igesstr); // the first loop is the outer
//  addIgesString(lines,to_string(0),dirline,count,igesstr); //!!! no outer loop is identified
  for (int i = 0; i < int(loopDEs.size()); i++)
  {
    addIgesString(lines,to_string(loopDEs[i]),dirline,count,igesstr);
  }

  finalize(lines,dirline,count,igesstr);
}

/** Entity 514. */
template <class T> void addIges514(std::vector<std::string> &lines, std::vector<int> &faceDEs, int dirline, int *count)
{
  std::string igesstr = "";

  addIgesString(lines,to_string(514),dirline,count,igesstr);
  addIgesString(lines,to_string((int) faceDEs.size()),dirline,count,igesstr);

  for (int i = 0; i < faceDEs.size(); i++)
  {
    addIgesString(lines,to_string(faceDEs[i]),dirline,count,igesstr);
    addIgesString(lines,to_string((int) 1),dirline,count,igesstr);
  }

  finalize(lines,dirline,count,igesstr);
}

/** Entity 186. */
template <class T> void addIges186(std::vector<std::string> &lines, int shellDE, int dirline, int *count)
{
  std::string igesstr = "";

  addIgesString(lines,to_string(186),dirline,count,igesstr);
  addIgesString(lines,to_string(shellDE),dirline,count,igesstr);
  addIgesString(lines,to_string(1),dirline,count,igesstr);
  addIgesString(lines,to_string(0),dirline,count,igesstr);

  finalize(lines,dirline,count,igesstr);
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

  // save T directory lines for each surface
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

/** Make iges lines for trimmed surfaces. */
template <class T> bool makeTrimmedSurfaceLinesIges(std::vector<tcad::TSplineSurface<T> *> &surfaces, 
//surface     loop        piece       points
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV,
  std::vector<std::string> &lines, int splinedegree = SPLINE_DEGREE)
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
  int dlines = 3;

  // save T directory lines for each surface
  int dcount = 7;
  int dir128line = int(lines.size());
  for (int i = 0; i < int(surfaces.size()); i++)
  {
    char name[32] = { 0 };
    
    s = makeIgesDirectoryLine0(dirline1280,-1,-1,&dcount);
    lines.push_back(s);
    s = makeIgesDirectoryLine1(dirline1281,-1,-1,&dcount,name);
    lines.push_back(s);
    dlines++;
    
    if (!boundariesUV[i].empty())
    {
      // loops
      for (int il = 0; il < int(boundariesUV[i].size()); il++)
      {
        for (int j = 0; j < int(boundariesUV[i][il].size()); j++)
        {
          // curve B (126)
          s = makeIgesDirectoryLine0(dirline1260,-1,-1,&dcount);
          lines.push_back(s);
          s = makeIgesDirectoryLine1(dirline1261,-1,-1,&dcount,name);
          lines.push_back(s);
          dlines++;

          // curve C (126)
          s = makeIgesDirectoryLine0(dirline1260,-1,-1,&dcount);
          lines.push_back(s);
          s = makeIgesDirectoryLine1(dirline1261,-1,-1,&dcount,name);
          lines.push_back(s);
          dlines++;
        }

        // bounding curve (141)
        s = makeIgesDirectoryLine0(dirline1410,-1,-1,&dcount);
        lines.push_back(s);
        s = makeIgesDirectoryLine1(dirline1411,-1,-1,&dcount,name);
        lines.push_back(s);
        dlines++;
      }

      // trimmed surface (143)
      s = makeIgesDirectoryLine0(dirline1430,-1,-1,&dcount);
      lines.push_back(s);
      s = makeIgesDirectoryLine1(dirline1431,-1,-1,&dcount,name);
      lines.push_back(s);
      dlines++;
    }
  }

  // total number of directory lines
  dlines = dlines * 2;

  // save 3 lines  
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

    // patch name
    char name[32] = { 0 };
    
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

    if (!boundariesUV[i].empty())
    {
      // loops
      std::vector<int> list141;
      int surfacedir = 0;
      for (int il = 0; il < int(boundariesUV[i].size()); il++)
      {
        // parametric boundary is input, we need XYZ boundary as well
        std::vector<std::vector<tcad::TPoint<T>>> boundary;
        surface->boundaryIntoPoints(boundariesUV[i][il],boundary);

        // spline boundary curve is combined of boundary parts which 
        // connect to each other with C1 discontinuity
        for (int j = 0; j < int(boundariesUV[i][il].size()); j++)
        {
          tcad::TSplineCurve<T> B(boundariesUV[i][il][j],int(boundary[j].size()) - 1,splinedegree,tcad::END_FREE,tcad::END_FREE); //!!!

          lines[dir128line] = makeIgesDirectoryLine0(dirline1260parametric,pcount,-1,&dcount);
          dir128line++;

          // add lines
          int before = int(lines.size());
          addIgesCurveLines(lines,&B,dcount - 1,&pcount);
          int after = int(lines.size());

          // modify directory
          lines[dir128line] = makeIgesDirectoryLine1(dirline1261,-1,(after - before),&dcount,name);
          dir128line++;
        }

        // curve C in XYZ
        for (int j = 0; j < int(boundary.size()); j++)
        {
          tcad::TSplineCurve<T> C(boundary[j],int(boundary[j].size()) - 1,splinedegree,tcad::END_FREE,tcad::END_FREE); //!!!

          lines[dir128line] = makeIgesDirectoryLine0(dirline1260,pcount,-1,&dcount);
          dir128line++;

          // add lines
          before = int(lines.size());
          addIgesCurveLines(lines,&C,dcount - 1,&pcount);
          after = int(lines.size());

          // modify directory
          lines[dir128line] = makeIgesDirectoryLine1(dirline1261,-1,(after - before),&dcount,name);
          dir128line++;
        }

        // save 141 
        lines[dir128line] = makeIgesDirectoryLine0(dirline1410,pcount,-1,&dcount);
        dir128line++;

        int numcurves = int(boundary.size());
        if (il == 0)
          surfacedir = dcount - (numcurves * 2 + 2) * 2 + 1;

        // add lines B and C
        before = int(lines.size());
        addIges141<T>(lines,surfacedir,dcount - (numcurves * 2 + 1) * 2 + 1,
          dcount - (numcurves + 1) * 2 + 1,
          dcount - 1,&pcount,numcurves);
        after = int(lines.size());

        // modify directory
        lines[dir128line] = makeIgesDirectoryLine1(dirline1411,-1,(after - before),&dcount,name);
        dir128line++;

        list141.push_back(dcount - 1 * 2);
      } // loops

      // save 143 
      lines[dir128line] = makeIgesDirectoryLine0(dirline1430,pcount,-1,&dcount);
      dir128line++;

      // add lines
      before = int(lines.size());
      addIges143<T>(lines,surfacedir,list141,dcount - 1,&pcount);
      after = int(lines.size());

      // modify directory
      lines[dir128line] = makeIgesDirectoryLine1(dirline1431,-1,(after - before),&dcount,name);
      dir128line++;
    } // loops
  } // loop on surfaces

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

/** Make Iges lines to save a solid. */
template <class T> bool makeSolidLinesIges(std::vector<tcad::TSplineSurface<T> *> &surfaces, 
//surface     loop        piece       points
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV,
  std::vector<std::string> &lines, T tolerance, int splinedegree = SPLINE_DEGREE, int numdigits = 18, 
  std::vector<std::vector<TPoint<T>>> *pbadedges = nullptr, int attempts = 40)
{
  lines.clear();

  if (surfaces.empty())
    return false;

  // unique vertices
  std::vector<TPoint<T>> vertices,middlevertices;

  // edges mapping to faces and face boundaries
  // straight list of edges, "compressed" edgemap
  std::vector<std::array<LINT,11>> edges;

  // Create non-manifold solid model
  if (!createSolidEdges(surfaces,boundariesUV,vertices,middlevertices,edges,
    tolerance,pbadedges,attempts))
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
  for (int i = 0; i < 12; i++)
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
  int dlines = 0;
  int dcount = 1;

  //===== Directory ============================================================

  // save T directory lines for each surface
  int sdirline = int(lines.size());

  // for entity 502
  s = makeIgesDirectoryLine0(sdirline5020,-1,-1,&dcount);
  lines.push_back(s);
  s = makeIgesDirectoryLine1(sdirline5021,-1,-1,&dcount,"");
  lines.push_back(s);
  dlines++;

  // 126 curves
  for (int i = 0; i < edges.size(); i++)
  {
    // curve C (126)
    s = makeIgesDirectoryLine0(sdirline1260,-1,-1,&dcount);
    lines.push_back(s);
    s = makeIgesDirectoryLine1(sdirline1261,-1,-1,&dcount,"");
    lines.push_back(s);
    dlines++;
  }

  // for entity 504
  s = makeIgesDirectoryLine0(sdirline5040,-1,-1,&dcount);
  lines.push_back(s);
  s = makeIgesDirectoryLine1(sdirline5041,-1,-1,&dcount,"");
  lines.push_back(s);
  dlines++;

  // 128 curves each with 508 and 510
  for (int i = 0; i < int(surfaces.size()); i++)
  {
    // patch name
    char name[32] = { 0 };

    // for entity 128
    s = makeIgesDirectoryLine0(sdirline1280,-1,-1,&dcount);
    lines.push_back(s);
    s = makeIgesDirectoryLine1(sdirline1281,-1,-1,&dcount,name);
    lines.push_back(s);
    dlines++;

    // 508 loops
    for (int j = 0; j < int(boundariesUV[i].size()); j++) 
    {
      s = makeIgesDirectoryLine0(sdirline5080,-1,-1,&dcount);
      lines.push_back(s);
      s = makeIgesDirectoryLine1(sdirline5081,-1,-1,&dcount,name);
      lines.push_back(s);
      dlines++;
    }

    // 510
    s = makeIgesDirectoryLine0(sdirline5100,-1,-1,&dcount);
    lines.push_back(s);
    s = makeIgesDirectoryLine1(sdirline5101,-1,-1,&dcount,name);
    lines.push_back(s);
    dlines++;
  }

  // for entity 514
  s = makeIgesDirectoryLine0(sdirline5140,-1,-1,&dcount);
  lines.push_back(s);
  s = makeIgesDirectoryLine1(sdirline5141,-1,-1,&dcount,"");
  lines.push_back(s);
  dlines++;

  // for entity 186
  s = makeIgesDirectoryLine0(sdirline1860final,-1,-1,&dcount);
  lines.push_back(s);
  s = makeIgesDirectoryLine1(sdirline1861,-1,-1,&dcount,"");
  lines.push_back(s);
  dlines++;

  //===== Directory over, implementation =======================================

  // total number of directory lines
  dlines = dlines * 2;

  dcount = 1;
  int pcount = 1;

  int verticesDE = dcount;

  lines[sdirline] = makeIgesDirectoryLine0(sdirline5020, pcount, -1, &dcount);
  sdirline++;

  size_t before = lines.size();
  addIges502<T>(lines,vertices,dcount - 1,&pcount,numdigits);
  size_t after = lines.size();

  // modify directory
  lines[sdirline] = makeIgesDirectoryLine1(sdirline5021,-1,int(after - before),&dcount,"");
  sdirline++;

  // 126
  // directory pointers for all boundary pieces one by one
  std::vector<int> boundaryDEs;
  std::vector<std::vector<TPoint<T>>> boundarysegments;
  std::vector<tcad::TSplineCurve<T>> boundarycurves;
  std::vector<int> parmboundaryDEs;

  for (auto &e : edges)
  {
    // prepare XYZ boundary piece
    std::vector<TPoint<T>> points;
    getBoundaryPartXYZ<T>(surfaces,boundariesUV,int(e[3]),int(e[4]),int(e[5]),points);

    // make spline curve
    tcad::TSplineCurve<T> C(points,int(points.size()) - 1,splinedegree,tcad::END_CLAMPED,tcad::END_CLAMPED); //!!!

    lines[sdirline] = makeIgesDirectoryLine0(dirline1260,pcount,-1,&dcount);
    sdirline++;

    // add lines
    before = int(lines.size());
    addIgesCurveLines(lines,&C,dcount - 1,&pcount);

    boundaryDEs.push_back(dcount - 1);
    boundarysegments.push_back(points);
    boundarycurves.push_back(C);

    after = int(lines.size());

    // modify directory
    lines[sdirline] = makeIgesDirectoryLine1(sdirline1261,-1,int(after - before),&dcount,"");
    sdirline++;
  }

  lines[sdirline] = makeIgesDirectoryLine0(sdirline5040,pcount,-1,&dcount);
  sdirline++;

  before = lines.size();

  addIges504<T>(lines,edges,verticesDE,boundaryDEs,dcount - 1,&pcount);

  int edgeDE = dcount - 1;
  after = lines.size();

  // modify directory
  lines[sdirline] = makeIgesDirectoryLine1(sdirline5041,-1,int(after - before),&dcount,"");
  sdirline++;

  // all 128 curves each with 508 and 510
  std::vector<int> faceDEs;

  for (int i = 0; i < int(surfaces.size()); i++)
  {
    // patch name
    char name[32] = { 0 };

    TSplineSurface<T> *surface = surfaces[i];

    int surfaceDE = -1;

    // modify directory
    lines[sdirline] = makeIgesDirectoryLine0(sdirline1280,pcount,-1,&dcount);
    sdirline++;

    // add lines
    int before = int(lines.size());
    addIgesSurfaceLines(lines,surface,dcount - 1,&pcount);
    surfaceDE = dcount - 1;
    int after = int(lines.size());

    // modify directory
    lines[sdirline] = makeIgesDirectoryLine1(sdirline1281,-1,(after - before),&dcount,name);
    sdirline++;

    // 508
    std::vector<int> loopDEs;
    for (int j = 0; j < int(boundariesUV[i].size()); j++) 
    {
      lines[sdirline] = makeIgesDirectoryLine0(sdirline5080,pcount,-1,&dcount);
      sdirline++;

      before = int(lines.size());

      addIges508<T>(surfaces,boundariesUV,vertices,middlevertices,tolerance,
        lines,edges,i,j,int(boundariesUV[i][j].size()),edgeDE,dcount - 1,&pcount);

      loopDEs.push_back(dcount - 1);
      after = int(lines.size());

      // modify directory
      lines[sdirline] = makeIgesDirectoryLine1(sdirline5081,-1,(after - before),&dcount,name);
      sdirline++;
    }

    // 510
    lines[sdirline] = makeIgesDirectoryLine0(sdirline5100,pcount,-1,&dcount);
    sdirline++;

    before = int(lines.size());
    addIges510<T>(lines,surfaceDE,loopDEs,dcount - 1,&pcount);
    faceDEs.push_back(dcount - 1);
    after = int(lines.size());

    // modify directory
    lines[sdirline] = makeIgesDirectoryLine1(sdirline5101,-1,(after - before),&dcount,name);
    sdirline++;
  }

  lines[sdirline] = makeIgesDirectoryLine0(sdirline5140,pcount,-1,&dcount);
  sdirline++;

  before = lines.size();
  addIges514<T>(lines,faceDEs,dcount - 1,&pcount);
  int shellDE = dcount - 1;
  after = lines.size();

  // modify directory
  lines[sdirline] = makeIgesDirectoryLine1(sdirline5141,-1,int(after - before),&dcount,"Shell");
  sdirline++;

  lines[sdirline] = makeIgesDirectoryLine0(sdirline1860final,pcount,-1,&dcount);
  sdirline++;

  before = lines.size();
  addIges186<T>(lines,shellDE,dcount - 1,&pcount);
  int solidDE = dcount - 1;
  after = lines.size();

  // modify directory
  lines[sdirline] = makeIgesDirectoryLine1(sdirline1861,-1,int(after - before),&dcount,"Object");

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

///** Save curves as points. DECLARED ABOVE. */
//template <class T> bool saveLinesIges(std::vector<std::vector<tcad::TPoint<T>>> &curves, 
//  const std::string &filename, int splinedegree = SPLINE_DEGREE, int minpoints = 6)
//{
//  std::vector<std::string> lines;
//
//  if (makeCurveLinesIges(curves,lines,splinedegree,minpoints))
//  {
//    bool ok = writeLines(lines,filename);
//    return ok;
//  } else
//  {
//    return false;
//  }
//}

/** Save curve as points. */
template <class T> bool saveLinesIges(std::vector<tcad::TPoint<T>> &curve, 
  const std::string &filename, int splinedegree = SPLINE_DEGREE, int minpoints = 6)
{
  std::vector<std::vector<tcad::TPoint<T>>> curves;
  curves.push_back(curve);

  return saveLinesIges(curves,filename,splinedegree,minpoints);
}

/** Save any curve (based on TBaseCurve) as points. */
template <class T> bool saveCurveIges(tcad::TBaseCurve<T> &curve, const std::string &filename,
  int splinedegree = SPLINE_DEGREE, int minpoints = 6)
{
  std::vector<tcad::TPoint<T>> points;
  curve.createPoints(points);

  return saveLinesIges(points,filename,splinedegree,minpoints);
}

/** Save any TWO curve (based on TBaseCurve) as points e,g, to compare. */
template <class T> bool saveTwoCurvesIges(tcad::TBaseCurve<T> &curve0, tcad::TBaseCurve<T> &curve1, 
  const std::string &filename, int splinedegree = SPLINE_DEGREE, int minpoints = 6)
{
  std::vector<tcad::TPoint<T>> points0,points1;
  curve0.createPoints(points0);
  curve1.createPoints(points1);

  std::vector<std::vector<tcad::TPoint<T>>> points;
  points.push_back(points0);
  points.push_back(points1);

  return saveLinesIges(points,filename,splinedegree,minpoints);
}

/** Save any THREE curve (based on TBaseCurve) as points e,g, to compare. */
template <class T> bool saveThreeCurvesIges(tcad::TBaseCurve<T> &curve0, tcad::TBaseCurve<T> &curve1, 
  tcad::TBaseCurve<T> &curve2, const std::string &filename, int splinedegree = SPLINE_DEGREE, int minpoints = 6)
{
  std::vector<tcad::TPoint<T>> points0,points1,points2;
  curve0.createPoints(points0);
  curve1.createPoints(points1);
  curve2.createPoints(points2);

  std::vector<std::vector<tcad::TPoint<T>>> points;
  points.push_back(points0);
  points.push_back(points1);
  points.push_back(points2);

  return saveLinesIges(points,filename,splinedegree,minpoints);
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

/** Save trimmed surfaces in IGES. boundaries must be closed. */
template <class T> bool saveTrimmedSurfacesIges(std::vector<tcad::TSplineSurface<T> *> &surfaces, 
//surface     loop        piece       points
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV, const std::string &filename,
  int splinedegree = SPLINE_DEGREE)
{
  std::vector<std::string> lines;

  if (makeTrimmedSurfaceLinesIges(surfaces,boundariesUV,lines,splinedegree))
  {
    bool ok = writeLines(lines,filename);
    return ok;
  } else
  {
    return false;
  }
}

/** Save trimmed surface in IGES. boundary must be closed. */
template <class T> bool saveTrimmedSurfaceIges(tcad::TSplineSurface<T> *surface, 
//loop        piece       points
  std::vector<std::vector<std::vector<tcad::TPoint<T>>>> &boundaryUV, const std::string &filename,
  int splinedegree = SPLINE_DEGREE)
{
  std::vector<tcad::TSplineSurface<T> *> surfaces;
  surfaces.push_back(surface);

  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> boundariesUV;
  boundariesUV.push_back(boundaryUV);

  return saveTrimmedSurfacesIges(surfaces,boundariesUV,filename,splinedegree);
}

/** Save trimmed surfaces as solid in IGES. All surfaces must have closed boundaries. */
template <class T> bool saveSolidIges(std::vector<tcad::TSplineSurface<T> *> &surfaces, 
//surface     loop        piece       points
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &boundariesUV, 
  const std::string &filename,
  T tolerance, int splinedegree = SPLINE_DEGREE, int numdigits = 18, 
  std::vector<std::vector<TPoint<T>>> *pbadedges = nullptr, int attempts = 40)
{
  std::vector<std::string> lines;

  if (makeSolidLinesIges(surfaces,boundariesUV,lines,tolerance,splinedegree,
    numdigits,pbadedges,attempts))
  {
    bool ok = writeLines(lines,filename);
    return ok;
  } else
  {
    return false;
  }
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

