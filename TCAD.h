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

  Files :
    BASICS          (basic types, macros etc)
    INFRASTRUCTURE  (vectors, matrices, integration, linear systems, planes, points etc)
    CURVES
    SURFACES
    VOLUMES
    OCTREE          (for spacial partitioning only)
    OPERATIONS      (FFD etc)
    BLOCKS          (breps are inside)
    AUXILIARY       (strings, export, import)

*******************************************************************************/

#pragma once

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
