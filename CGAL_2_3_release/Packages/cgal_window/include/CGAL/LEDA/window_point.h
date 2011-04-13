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
// file          : include/CGAL/LEDA/window_point.h
// package       : cgal_window (1.0.3)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.0.3
// revision_date : 25 June 2001
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

/*{\Manpage {window_point} {} {Points}}*/

namespace CGAL {

class __exportC window_point {

/*{\Mdefinition
The data type $window\_point$ is a simple point data type used in the interface of
the window class.
}*/

  double xw, yw; // x and y coordinates ...

public:

/*{\Mcreation wpt}*/

  window_point() : xw(0), yw(0) { }
/*{\Mcreate  creates a window point with coordinates $(0,0)$. }*/

  window_point(double xc, double yc) : xw(xc), yw(yc) { }
/*{\Mcreate  creates a window point with coordinates $(xc,yc)$. }*/
  
  window_point(const window_point& p) : xw(p.xw), yw(p.yw) { }

  window_point& operator=(const window_point& p) { xw=p.xw; yw=p.yw; return *this; }

/*{\Moperations }*/

  // retrieving x and y coordinates ...
  
  double x() const { return xw; }
/*{\Mop  returns the $x$-coordinate of |\Mvar|. }*/
  
  double y() const { return yw; }
/*{\Mop  returns the $y$-coordinate of |\Mvar|. }*/
  
  double xcoord() const { return xw; }
/*{\Mop  returns the $x$-coordinate of |\Mvar|. }*/
  
  double ycoord() const { return yw; }  
/*{\Mop  returns the $y$-coordinate of |\Mvar|. }*/  

  bool operator==(const window_point& q) const { return ((x()==q.x()) && (y()==q.y()) ); }
  bool operator!=(const window_point& q) const { return !operator==(q);}

};



}

#endif
