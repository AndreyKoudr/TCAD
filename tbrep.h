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

  tbrep.h

  Common class for breps (surfaces + trimming curves).

*******************************************************************************/

#pragma once

#include "tpoints.h"
#include "tmisc.h"
#include "tsplinesurface.h"

namespace tcad {

/** Common class for breps (surfaces + trimming curves). For boolean operations 
  and ease of use. */
template <class T> struct TBrep {
public:

  // untrimmed surfaces
  std::vector<TSplineSurface<T> *> surfaces;

  // trimming curves for surfaces
  //surface   loop        piece      UV points in X,Y  
  std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> boundariesUV;

  /** Constructor. */
  TBrep()
  {
  }

  /** Destructor. */
  ~TBrep()
  {
    deleteSurfaces(surfaces);
  }

  /** Copy constructor. */
  TBrep(const TBrep &other)  
  {
    deleteSurfaces(surfaces);

    for (int i = 0; i < int(other.surfaces.size()); i++)
    {
      TSplineSurface<T> *surface = new TSplineSurface<T>(*(other.surfaces[i]));
      surfaces.push_back(surface);
    }

    boundariesUV = other.boundariesUV;
  }

  /** Assignment operator. */
  TBrep &operator = (const TBrep &other)  
  {
    deleteSurfaces(surfaces);

    for (int i = 0; i < int(other.surfaces.size()); i++)
    {
      TSplineSurface<T> *surface = new TSplineSurface<T>(*(other.surfaces[i]));
      surfaces.push_back(surface);
    }

    boundariesUV = other.boundariesUV;

    return *this;
  }

  /** Constructor, takes ownership of psurfaces pointers, do not delete them. */
  TBrep(std::vector<TSplineSurface<T> *> &psurfaces, 
    std::vector<std::vector<std::vector<std::vector<tcad::TPoint<T>>>>> &pboundariesUV)
  {
    surfaces = psurfaces;
    boundariesUV = pboundariesUV;
  }


};

}