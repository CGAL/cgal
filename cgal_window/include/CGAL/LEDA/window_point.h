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
