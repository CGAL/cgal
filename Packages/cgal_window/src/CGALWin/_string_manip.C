// Copyright (c) 1999  Martin-Luther-University Halle-Wittenberg (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Matthias Baesken, Algorithmic Solutions


#include <CGAL/LEDA/string_manip.h>
#include <string>
#include <cctype>


namespace CGAL {

using std::string;

#if defined(__BORLANDC__)
using std::isspace;
#endif

std::list<string> break_into_words(string s)
{ 
  std::list<string> L;

  int len = s.length();
  int pos1 = 0;
  while (pos1 < len && isspace(s[pos1])) pos1++;

  while (pos1 < len)
  { int pos2 = pos1;
    char sep = ' ';
    if (s[pos2] == '"') sep = '"';
    if (s[pos2] == '\'') sep = '\'';
    pos2++;
    while (pos2 < len && s[pos2] != sep) pos2++;
    if (pos2 == len || sep == ' ') pos2--;
    L.push_back(s.substr(pos1,pos2-pos1+1));
    pos1 = pos2+1;
    while (pos1 < len && isspace(s[pos1]) && pos1 < len) pos1++;
   }

  return L;
}



std::list<string> break_into_lines(string s)
{ 
  std::list<string> L;

  int len = s.length();
  int pos1 = 0;

  while (pos1 < len)
  { int pos2 = pos1;
    while (pos2 < len && s[pos2] != '\n') pos2++;
    L.push_back(s.substr(pos1,pos2-pos1));
    pos1 = pos2+1;
   }

  if (s[len-1] == '\n') L.push_back("");

  return L;
}


}
  

