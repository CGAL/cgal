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
// release       : 
// release_date  : 2001, May 23
//
// file          : include/CGAL/LEDA/window_point.h
// package       : cgal_window (0.9.7)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 0.9.7
// revision_date : 23 May 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================

#ifndef CGAL_WINDOW_POINT_H
#define CGAL_WINDOW_POINT_H

#if defined(CGAL_USE_CGAL_HEADERS)
#include <CGAL/basic.h>
#endif

#include <CGAL/LEDA/basic.h>

namespace CGAL {

class __exportC window_point {

  double xw, yw; // x and y coordinates ...

public:

  window_point() : xw(0), yw(0) { }

  window_point(double xc, double yc) : xw(xc), yw(yc) { }
  
  window_point(const window_point& p) : xw(p.xw), yw(p.yw) { }

  window_point& operator=(const window_point& p) { xw=p.xw; yw=p.yw; return *this; }

/*{\Moperations }*/

  // retrieving x and y coordinates ...
  
  double x() const { return xw; }
  
  double y() const { return yw; }
  
  double xcoord() const { return xw; }
  
  double ycoord() const { return yw; }  

  bool operator==(const window_point& q) const { return ((x()==q.x()) && (y()==q.y()) ); }
  bool operator!=(const window_point& q) const { return !operator==(q);}

};



}

#endif
