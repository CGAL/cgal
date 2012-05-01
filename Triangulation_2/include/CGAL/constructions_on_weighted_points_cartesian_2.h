// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Julia Flototto, Mariette Yvinec

#ifndef CGAL_CONSTRUCTIONS_ON_WEIGHTED_POINTS_CARTESIAN_2_H
#define CGAL_CONSTRUCTIONS_ON_WEIGHTED_POINTS_CARTESIAN_2_H

//#include <CGAL/Cartesian/redefine_names_2.h>

namespace CGAL { 


template < class RT >
void
weighted_circumcenter_translateC2(const RT &dqx, const RT &dqy, const RT &dqw,
				  const RT &drx, const RT &dry, const RT &drw,
				  RT &dcx,       RT &dcy)
{
  // Given 3 points P, Q, R, this function takes as input:
  // qx-px, qy-py,qw-pw,  rx-px, ry-py, rw-pw.  And returns cx-px, cy-py,
  // where (cx, cy) are the coordinates of the circumcenter C.

  // What we do is intersect the radical axis
  RT r2 = CGAL_NTS square(drx) + CGAL_NTS square(dry) - drw;
  RT q2 = CGAL_NTS square(dqx) + CGAL_NTS square(dqy) - dqw;

  RT den = RT(2) * determinant(dqx, dqy, drx, dry);

  // The 3 points aren't collinear.
  // Hopefully, this is already checked at the upper level.
  CGAL_assertion ( den != RT(0) );

  // One possible optimization here is to precompute 1/den, to avoid one
  // division.  However, we loose precision, and it's maybe not worth it (?).
  dcx =   determinant (dry, dqy, r2, q2) / den;
  dcy = - determinant (drx, dqx, r2, q2) / den;
}

//template < class RT >
template < class RT, class We>
void
weighted_circumcenterC2( const RT &px, const RT &py, const We &pw,
			 const RT &qx, const RT &qy, const We &qw,
			 const RT &rx, const RT &ry, const We &rw,
			 RT &x, RT &y )
{
  RT dqw = RT(qw-pw);
  RT drw = RT(rw-pw);
 
  weighted_circumcenter_translateC2(qx-px, qy-py, dqw,rx-px, ry-py,drw,x, y);
  x += px;
  y += py;
}


template < class RT , class We>
void
radical_axisC2(const RT &px, const RT &py, const We &pw,
	       const RT &qx, const RT &qy, const We &qw,
	       RT &a, RT &b, RT& c )
{
  a =  RT(2)*(px - qx);
  b =  RT(2)*(py - qy);
  c = - CGAL_NTS square(px) - CGAL_NTS square(py) 
      + CGAL_NTS square(qx) + CGAL_NTS square(qy) 
       +RT(pw) - RT(qw);
}




} //namespace CGAL
#endif //CGAL_CONSTRUCTIONS_ON_WEIGHTED_POINTS_CARTESIAN_2_H
