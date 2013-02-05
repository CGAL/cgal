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

#ifndef CGAL_SPHERE_TRIANGLE_H
#define CGAL_SPHERE_TRIANGLE_H

#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>
#include <CGAL/array.h>

namespace CGAL {

template <class R_> class Sphere_triangle_rep 
{ typedef Sphere_point<R_>  Point;
  typedef Sphere_circle<R_> Circle;
  typedef Sphere_triangle_rep<R_> Rep;

  cpp11::array<Point,3>  points_; 
  cpp11::array<Circle,3> circles_;

  friend class Sphere_triangle<R_>;
 
  Sphere_triangle_rep(const Point& p1, const Point& p2, const Point& p3,
        const Circle& c1, const Circle& c2, const Circle& c3) :
    points_(CGAL::make_array(p1,p2,p3)), circles_(CGAL::make_array(c1,c2,c3)) {}
public:
  Sphere_triangle_rep() {}
};


/*{\Manpage{Sphere_triangle}{R}{Triangles on the unit sphere}{t}}*/
template <class R_> class Sphere_triangle : 
  public Handle_for< Sphere_triangle_rep<R_> > {
/*{\Mdefinition An object |\Mvar| of type |\Mname| is a triangle
on the surface of the unit sphere.}*/

public:
/*{\Mtypes 5}*/
typedef R_ R;
/*{\Mtypemember representation class.}*/
typedef typename R::RT RT;
/*{\Mtypemember ring type.}*/

typedef Sphere_triangle_rep<R_> Rep;
typedef Handle_for<Rep>         Base;

/*{\Mcreation 5}*/
Sphere_triangle() : Base() {}
/*{\Mcreate creates some triangle.}*/

Sphere_triangle(
  const Sphere_point<R>& p0, const Sphere_point<R>& p1, 
  const Sphere_point<R>& p2,
  const Sphere_circle<R>& c0, const Sphere_circle<R>& c1, 
  const Sphere_circle<R>& c2) : Base(Rep(p0,p1,p2,c0,c1,c2)) 
/*{\Mcreate creates a triangle spanned by the three points
|p0|, |p1|, |p2|, where the triangle is left of the three circles
|c0|, |c1|, |c2|. \precond $c_i$ contains $p_i$ and $p_{i+1}$ mod 3.}*/
{ CGAL_assertion( c0.has_on(p0) && c0.has_on(p1) );
  CGAL_assertion( c1.has_on(p1) && c1.has_on(p2) );
  CGAL_assertion( c2.has_on(p2) && c0.has_on(p0) );
}

Sphere_triangle(const Sphere_triangle<R>& t) : Base(t) {} 

/*{\Moperations 4 2}*/

const Sphere_point<R>& point(unsigned i) const 
/*{\Mop returns the ith point of |\Mvar|.}*/
{ return this->ptr()->points_[i%3]; }

const Sphere_circle<R>& circle(unsigned i) const 
/*{\Mop returns the ith circle of |\Mvar|.}*/
{ return this->ptr()->circles_[i%3]; }

Sphere_triangle<R> opposite() const 
/*{\Mop returns the opposite of |\Mvar|.}*/
{ return Sphere_triangle<R>(point(0), point(1), point(2),
    circle(0).opposite(), circle(1).opposite(), circle(2).opposite()); }


}; // Sphere_triangle<R>


template <typename R>
std::ostream& operator<<(std::ostream& os, 
                         const CGAL::Sphere_triangle<R>& t)
{ for (int i=0; i<3; ++i) os << t.point(i);
  for (int i=0; i<3; ++i) os << t.circle(i); return os; }

template <typename R>
std::istream& operator>>(std::istream& is, 
                         CGAL::Sphere_triangle<R>& t)
{ CGAL::Sphere_point<R> p1,p2,p3;
  CGAL::Sphere_circle<R> c1,c2,c3;
  if ( !(is >> p1 >> p2 >> p3 >> c1 >> c2 >> c3) ) return is; 
  t = CGAL::Sphere_triangle<R>(p1,p2,p3,c1,c2,c3);
  return is; }

} //namespace CGAL
#endif //CGAL_SPHERE_TRIANGLE_H
