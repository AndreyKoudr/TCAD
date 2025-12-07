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

#include "ttriangles.h"

namespace tcad {

#ifdef USE_MUTEX
  std::mutex imutex;
  std::mutex emutex;
#endif

bool edgeComp(const std::pair<LINT,LINT> &a, const std::pair<LINT,LINT> &b)
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

bool edgesEqual(const std::pair<LINT,LINT> &a, const std::pair<LINT,LINT> &b, bool *reversed)
{
  bool ok0 = (a.first == b.first && a.second == b.second);
  bool ok1 = (a.first == b.second && a.second == b.first);

  if (reversed)
    *reversed = ok1;

  return ok0 || ok1;
}

bool edgeContainsTri(std::pair<LINT,LINT> edge, LINT tri, std::map<std::pair<LINT,LINT>,std::vector<LINT>,EdgeCompare> &edgeTris,
  LINT &othertri)
{
  othertri = -1;
  bool found = false;
  for (int i = 0; i < edgeTris[edge].size(); i++)
  {
    if (edgeTris[edge][i] == tri)
    {
      found = true;
    } else
    {
      othertri = edgeTris[edge][i];
    }
  }
  return found;
}

int findIntersectionByTri(LINT tri, std::vector<std::pair<LINT,LINT>> &iedges, std::vector<bool> &ibusy, 
  std::map<std::pair<LINT,LINT>,std::vector<LINT>,EdgeCompare> &edgeTris, LINT &othertri)
{
  for (int i = 0; i < int(iedges.size()); i++)
  {
    if (!ibusy[i])
    {
      if (edgeContainsTri(iedges[i],tri,edgeTris,othertri))
        return i;
    }
  }

  return -1;
}

bool findFreeEdge(std::map<std::pair<LINT,LINT>,std::vector<LINT>,EdgeCompare> &edgeTris,
  std::set<std::pair<LINT,LINT>,EdgeCompare> &busy, std::pair<LINT,LINT> &startedge)
{ 
  bool found = false;
  for (auto &e : edgeTris)
  {
    if (busy.find(e.first) == busy.end() && e.second.size() == 1)
    {
      startedge = e.first;
      found = true;
      break;
    }
  }

  return found;
}

}