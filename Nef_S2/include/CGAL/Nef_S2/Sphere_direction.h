// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Seel  <seel@mpi-sb.mpg.de>

#ifndef CGAL_SPHERE_DIRECTION_H
#define CGAL_SPHERE_DIRECTION_H

#include <CGAL/license/Nef_S2.h>


#include <CGAL/basic.h>
#include <CGAL/Kernel/global_functions_3.h>

namespace CGAL {

/*{\Manpage{Sphere_direction}{R}{Directions on the unit sphere}{c}}*/

template <class R_> class Sphere_direction : public R_::Plane_3 {

/*{\Mdefinition An object |\Mvar| of type |\Mname| is a direction
on the surface of the unit sphere.  Such directions can be used to
describe walks that are part of great circles.}*/

public:

/*{\Mtypes 5}*/
typedef R_ R;
/*{\Mtypemember representation class.}*/
typedef typename R::RT RT;
/*{\Mtypemember ring type.}*/

typedef typename R_::Point_3 Point_3;
typedef typename R_::Plane_3 Plane_3;
typedef typename R_::Plane_3 Base;
typedef Sphere_direction<R_> Self;

/*{\Mcreation 5}*/
Sphere_direction() : Base() {}
/*{\Mcreate creates some direction.}*/

Sphere_direction(const Sphere_circle<R>& c) 
/*{\Mcreate creates the direction corresponding to the circle |c|.}*/
  : Base(c) {} 

Sphere_direction(const Sphere_point<R>& p, const Sphere_point<R>&q) 
  : Base(Point_3(0,0,0),p,q) 
/*{\Mcreate creates a direction that describes the orientation of
the great circle through $p$ and $q$ (oriented such that the segment
$pq$ is the shorter one of the two possible ones. \precond $p$ and $q$
are not opposite on $S_2$.}*/ 
{ CGAL_assertion(p!=q.opposite()); 
  Point_3 p1(0,0,0), p4 = CGAL::ORIGIN + ((Base*) this)->orthogonal_vector();
  if ( CGAL::orientation(p1,p,q,p4) != CGAL::POSITIVE )
    *this = Sphere_direction(opposite());
}

Sphere_direction(const typename R::Plane_3& h) 
/*{\Xcreate creates the direction corresponding to the plane |h|.
\precond |h| contains the origin.}*/
 : Base(h) { CGAL_assertion(h.d() == 0); } 

/*{\Moperations 4 2}*/

Sphere_direction<R> opposite() const 
/*{\Mop returns the opposite of |\Mvar|.}*/
{ return Base::opposite(); }

Plane_3 plane() const { return Base(*this); }
/*{\Xop returns the plane supporting |\Mvar|.}*/

}; // Sphere_direction<R>


/* We have:
   1) all directions fixed at p 
   2) d1==d3 possible
   return true iff d1,d2,d3 are stricly ccw ordered around p
   Note: Sphere_directions are Plane_3
         we therefore compare the normal vectors of the planes
         that underly the directions d1,d2,d3 in the plane 
         through 0 and orthogonal to the vector p-0
 */

template <typename R>
bool strictly_ordered_ccw_at(const Sphere_point<R>& p,
  const Sphere_direction<R>& d1,
  const Sphere_direction<R>& d2,
  const Sphere_direction<R>& d3)
{ CGAL_assertion(d1.has_on(p) && d2.has_on(p) && d3.has_on(p));
  typename R::Point_3 p0(0,0,0); 
  typename R::Point_3 p1(CGAL::ORIGIN + d1.orthogonal_vector()); 
  typename R::Point_3 p2(CGAL::ORIGIN + d2.orthogonal_vector()); 
  typename R::Point_3 p3(CGAL::ORIGIN + d3.orthogonal_vector()); 

  if ( d1 == d3 ) return false;
  if ( CGAL::orientation(p0,p,p1,p3) == CGAL::POSITIVE ) {
    return CGAL::orientation(p0,p,p1,p2) == CGAL::POSITIVE &&
           CGAL::orientation(p0,p,p3,p2) == CGAL::NEGATIVE;
  } else {
    return CGAL::orientation(p0,p,p1,p2) == CGAL::POSITIVE ||
           CGAL::orientation(p0,p,p3,p2) == CGAL::NEGATIVE;
  }
}


} //namespace CGAL
#endif //CGAL_SPHERE_DIRECTION_H
