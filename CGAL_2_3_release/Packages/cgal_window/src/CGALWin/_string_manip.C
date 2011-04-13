// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-75 $
// release_date  : $CGAL_Date: 2001/06/21 $
//
// file          : src/CGALWin/_string_manip.C
// package       : cgal_window (1.0)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.0.2
// revision_date : 25 June 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================


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
  

