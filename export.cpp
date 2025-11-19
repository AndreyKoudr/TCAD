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

#include "export.h"

bool writeLines(const std::vector<std::string> &lines, const std::string &filename)
{
  FILE *fp = nullptr;
  if (fopen_s(&fp,filename.c_str(),"wb") != 0)
    return false;

  int error = 0;
  for (auto &line : lines)
  {
    error = fprintf(fp,"%s",line.c_str());
    if (error < 0) return false;
  }

  fclose(fp);

  return true;
}

// first part of IGES header
const char *IgesHeader0[18] =
{
  "                                                                        S      1\r\n",
  "1H,,1H;,,                                                               G      1\r\n",
  "08H000.iges,                                                            G      2\r\n",
  "10HCAD writer,10HCAD writer,                                            G      3\r\n",
  "32,38,6,308,15,                                                         G      4\r\n",
  ",                                                                       G      5\r\n",
  "1.0D0,6,1HM,1,0.254D0,13H070623.120000,                                 G      6\r\n",
  "0.01D0,                                                                 G      7\r\n",
  "1.0D0,                                                                  G      8\r\n",
  "03HCAD,                                                                 G      9\r\n",
  "10HCAD writer,                                                          G     10\r\n",
  "10,0,13H070623.120000;                                                  G     11\r\n",
  "     314       1       0       0       0       0       0       000000200D      1\r\n",
  "     314       0       1       1       0       0       0   COLOR       0D      2\r\n",
  "     406       2       0       0       1       0       0       000000300D      3\r\n",
  "     406       0      -1       1       3       0       0LEVELDEF       0D      4\r\n",
  "     406       3       0       0       2       0       0       000000300D      5\r\n",
  "     406       0      -1       1       3       0       0LEVELDEF       0D      6\r\n"
};

// second part of IGES header
const char *IgesHeader1[3] =
{
  "314,0.0,0.0,0.0,20HRGB(   0,   0,   0 );                         0000001P      1\r\n",
  "406,2,1,7HDefault;                                               0000003P      2\r\n",
  "406,2,2,12HIGES level 0;                                         0000005P      3\r\n"
};

/*
 0        1         2         3         4         5         6         7         8
 12345678901234567890123456789012345678901234567890123456789012345678901234567890
*/

const char *dirline1280 =
"     128       4       0       0       2       0       0       000000000D      7\r\n";
const char *dirline1281 =
"     128       0       0    1167       0       0       0 TrimSrf       0D      8\r\n";

const char *dirline1260 =
"     126       4       0       0       2       0       0       000000000D      7\r\n";

// for saving lines
const char *dirline1260l =
"     126       4       0       0       2       0       0       000000000D      7\r\n";

const char *dirline1260parametric =
"     126       4       0       0       2       0       0       000010500D      7\r\n";
const char *dirline1261 =
"     126       0       0    1167       0       0       0 TrimSrf       0D      8\r\n";

const char *dirline1420 =
"     142       4       0       0       2       0       0       000000000D      7\r\n"; // 05 is parametric
const char *dirline1421 =
"     142       0      -1    1167       0       0       0 TrimSrf       0D      8\r\n";

const char *dirline1440 =
"     144       4       0       0       2       0       0       000000000D      7\r\n";
const char *dirline1441 =
"     144       0      -1    1167       0       0       0 TrimSrf       0D      8\r\n";

const char *dirline1410 =
"     141       4       0       0       2       0       0       000000000D      7\r\n"; // 05 is parametric
const char *dirline1411 =
"     141       0      -1    1167       0       0       0 TrimSrf       0D      8\r\n";

const char *dirline1430 =
"     143       4       0       0       2       0       0       000000000D      7\r\n";
const char *dirline1431 =
"     143       0      -1    1167       0       0       0 TrimSrf       0D      8\r\n";

const char *dirline1200 =
"     120       4       0       0       2       0       0       000000000D      7\r\n";
const char *dirline1201 =
"     120       0      -1    1167       0       0       0 Cylindr       0D      8\r\n";

const char *dirline1100 =
"     110       4       0       0       2       0       0       000000000D      7\r\n";
const char *dirline1101 =
"     110       0      -1    1167       1       0       0 Lineseg       0D      8\r\n"; // form number 1

const char *dirline5020 =
"     502       1       0       0       0       0       0       000000000D      7\r\n";
const char *dirline5021 =
"     502       0       0       9       1       0       0               0D      8\r\n";

const char *dirline5040 =
"     504       1       0       0       0       0       0       000000000D      7\r\n";
const char *dirline5041 =
"     504       0       0       3       1       0       0               0D      8\r\n";

const char *dirline5080 =
"     508       1       0       0       0       0       0       000000000D      7\r\n";
const char *dirline5081 =
"     508       0       0       1       1       0       0               0D      8\r\n";

const char *dirline5100 =
"     510       1       0       0       0       0       0       000000000D      7\r\n";
const char *dirline5101 =
"     510       0       0       1       1       0       0               0D      8\r\n";

const char *dirline5140 =
"     514       1       0       0       0       0       0       000000000D      7\r\n";
const char *dirline5141 =
"     514       0       0       1       1       0       0               0D      8\r\n";

const char *dirline1860 =
"     186       1       0       0       0       0       0       000000000D      7\r\n";
const char *dirline1861 =
"     186       0       0       1       0       0       0               0D      8\r\n";

const char *dirline3080 =
"     308       1       0       0       0       0       0       000000000D      7\r\n";
const char *dirline3081 =
"     308       0       0       1       0       0       0               0D      8\r\n";



const char *dirline1240 =
"     124       1       0       0       0       0       0       000010001D      7\r\n";
const char *dirline1241 =
"     124       0       0       1       0       0       0               0D      8\r\n";

const char *dirline4080 =
"     408       1       0       0       0       0       0       000000001D      7\r\n";
const char *dirline4081 =
"     408       0       0       1       0       0       0               0D      8\r\n";

// final line if not addsubfigureentities
const char *dirline1860final =
"     186       1       0       0       0       0       0       000000000D      7\r\n";




/*
 0        1         2         3         4         5         6         7         8
 12345678901234567890123456789012345678901234567890123456789012345678901234567890
*/

// "000010001D" - the first "1" from the left means "dependent and not displayed"

// 128 must be shown for ccm+
const char *sdirline1280 =
"     128       4       0       0       0       0       0       000010000D      7\r\n";
const char *sdirline1281 =
"     128       0       0    1167       0       0       0               0D      8\r\n";

const char *sdirline1260 =
"     126       4       0       0       0       0       0       000010000D      7\r\n";
const char *sdirline1260parametric =
"     126       4       0       0       0       0       0       000010500D      7\r\n";
const char *sdirline1261 =
"     126       0       0    1167       0       0       0               0D      8\r\n";

// for saving lines
const char *sdirline1260l =
"     126       4       0       0       2       0       0       000000000D      7\r\n";


const char *sdirline1420 =
"     142       4       0       0       0       0       0       000010000D      7\r\n"; // 05 is parametric
const char *sdirline1421 =
"     142       0      -1    1167       0       0       0               0D      8\r\n";

const char *sdirline1440 =
"     144       4       0       0       0       0       0       000010000D      7\r\n";
const char *sdirline1441 =
"     144       0      -1    1167       0       0       0               0D      8\r\n";

const char *sdirline1410 =
"     141       4       0       0       0       0       0       000010000D      7\r\n"; // 05 is parametric
const char *sdirline1411 =
"     141       0      -1    1167       0       0       0               0D      8\r\n";

const char *sdirline1430 =
"     143       4       0       0       0       0       0       000010000D      7\r\n";
const char *sdirline1431 =
"     143       0      -1    1167       0       0       0               0D      8\r\n";

const char *sdirline1200 =
"     120       4       0       0       0       0       0       000010000D      7\r\n";
const char *sdirline1201 =
"     120       0      -1    1167       0       0       0               0D      8\r\n";

const char *sdirline1100 =
"     110       4       0       0       0       0       0       000010000D      7\r\n";
const char *sdirline1101 =
"     110       0      -1    1167       1       0       0               0D      8\r\n"; // form number 1

const char *sdirline5020 =
"     502       1       0       0       0       0       0       000010000D      7\r\n";
const char *sdirline5021 =
"     502       0       0       9       1       0       0               0D      8\r\n";

// must be 000010001
const char *sdirline5040 =
"     504       1       0       0       0       0       0       000010000D      7\r\n";
const char *sdirline5041 =
"     504       0       0       3       1       0       0               0D      8\r\n";

const char *sdirline5080 =
"     508       1       0       0       0       0       0       000010000D      7\r\n";
const char *sdirline5081 =
"     508       0       0       1       1       0       0               0D      8\r\n";

const char *sdirline5100 =
"     510       1       0       0       0       0       0       000000000D      7\r\n";
const char *sdirline5101 =
"     510       0       0       1       1       0       0               0D      8\r\n";

const char *sdirline5140 =
"     514       1       0       0       0       0       0       000010000D      7\r\n";
const char *sdirline5141 =
"     514       0       0       1       1       0       0               0D      8\r\n";

// "solid" entity if addsubfigureentities (not last)
const char *sdirline1860 =
"     186       1       0       0       0       0       0       000010000D      7\r\n";
// final line if not addsubfigureentities
const char *sdirline1860final =
"     186       1       0       0       0       0       0       000000000D      7\r\n";
const char *sdirline1861 =
"     186       0       0       1       0       0       0               0D      8\r\n";

const char *sdirline3080 =
"     308       1       0       0       0       0       0       000010000D      7\r\n";
const char *sdirline3081 =
"     308       0       0       1       0       0       0               0D      8\r\n";

const char *sdirline1240 =
"     124       1       0       0       0       0       0       000010000D      7\r\n";
const char *sdirline1241 =
"     124       0       0       1       0       0       0               0D      8\r\n";

const char *sdirline4080 =
"     408       1       0       0       0       0       0       000000000D      7\r\n";
const char *sdirline4081 =
"     408       0       0       1       0       0       0               0D      8\r\n";




std::string writeHeaderLine(int line, double val)
{
  std::string s = IgesHeader0[line];

  std::string str = trimReal(to_string(val,20)) + ",";
  memmove(&s[0],str.c_str(),str.size());
  return s;
}

std::string makeIgesDirectoryLine0(const char *dirline, int parmdata, int parmcount, int *dirlineno)
{
  std::string s = dirline;

  std::string substr = padFromLeft(to_string(parmdata),8);
  memmove(&(s[8]),substr.c_str(),substr.size());

  substr = padFromLeft(to_string(*dirlineno),7);
  memmove(&(s[73]),substr.c_str(),substr.size());

  (*dirlineno)++;

  return s;
}

std::string makeIgesDirectoryLine1(const char *dirline, int parmdata, int parmcount, int *dirlineno, const char *name)
{
  std::string s = dirline;

  std::string substr = padFromLeft(to_string(parmcount),8);
  memmove(&(s[24]),substr.c_str(),substr.size());

  substr = padFromLeft(to_string(*dirlineno),7);
  memmove(&(s[73]),substr.c_str(),substr.size());

  if (strlen(name) > 0)
  {
    std::string substr = name;
    substr = trim(substr);
    substr = padFromLeft(substr,8);
    assert(substr.size() == 8);
    memmove(&(s[56]),substr.c_str(),substr.size());
  }

  (*dirlineno)++;

  return s;
}

void addIgesString(std::vector<std::string> &lines, std::string substr, int dirline, int *count, std::string &igesstr)
{
  substr = trimLeft(substr);
  substr = trimReal(substr);

  // plus comma
  int substrlen = int(substr.size()) + 1;
  int currentlen = int(igesstr.size());
  if ((currentlen + substrlen) <= 64)
  {
    igesstr = igesstr + substr + ",";
  } else
  {
    // complete string and add to the list
    igesstr = padFromRight(igesstr,65,' ');
    igesstr = igesstr + to_string(dirline,7) + "P" + padFromLeft(to_string(*count), 7, ' ') + CRLF;

    // with ending zero
    lines.push_back(igesstr);

    igesstr = substr + ",";

    (*count)++;
  }
}

void finalize(std::vector<std::string> &lines, int dirline, int *count, std::string &igesstr)
{
  int currentlen = int(igesstr.size());
  if (currentlen)
  {
    // get rid of final comma
    igesstr.erase(igesstr.end() - 1);
    igesstr = igesstr + ";";

    // complete string and add to the list
    igesstr = padFromRight(igesstr,65,' ');
    igesstr = igesstr + to_string(dirline,7) + "P" + padFromLeft(to_string(*count),7,' ') + CRLF;

    // with ending zero
    lines.push_back(igesstr);

    (*count)++;
    igesstr = "";
  }
}


