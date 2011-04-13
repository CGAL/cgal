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
// file          : include/CGAL/LEDA/string_manip.h
// package       : cgal_window (1.0.3)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.0.3
// revision_date : 25 June 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================


#ifndef CGAL_WINDOW_STRING_MANIP_H
#define CGAL_WINDOW_STRING_MANIP_H

#if defined(CGAL_USE_CGAL_HEADERS)
#include <CGAL/basic.h>
#endif

#include <string>
#include <list>

namespace CGAL {


extern  std::list<std::string>  break_into_words(std::string s);
/*{\Mfunc   returns the list of words (separated by white space) of |s|. }*/

extern  std::list<std::string>  break_into_lines(std::string s);
/*{\Mfunc   returns the list of lines (separated by newline) of |s|. }*/

}

#endif

