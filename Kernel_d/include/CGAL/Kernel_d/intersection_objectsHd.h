// Copyright (c) 2002  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : ?
#ifndef CGAL_INTERSECTION_OBJECTSHD_H
#define CGAL_INTERSECTION_OBJECTSHD_H

#include <CGAL/basic.h>

namespace CGAL {

/*{\Manpage{Line_line_intersectionHd}{R}{intersecting two lines}}*/

template <class R> 
class Line_line_intersectionHd {

typedef typename R::RT RT;
typedef typename R::FT FT;
typedef typename R::LA LA;
typedef typename R::Point_d Point_d;

public:
enum Intersection_result { NO_INTERSECTION, POINT, LINE };

Intersection_result operator()(
  const Point_d& s1, const Point_d& t1,
  const Point_d& s2, const Point_d& t2,
  Point_d& p, FT& l1, FT& l2) 
/*{\Mfunop returns |NO_INTERSECTION| if the lines which are represented by |s1t1|
and |s2t2| don't intersect, returns |POINT| if they intersect in a
unique point, and returns LINE if they are identical. In the |POINT|
case the point of intersection is assigned to |p|.  Then |p = s1 + l1
* t1-s1| and |p = s2 + l2 * t2-s2|. \precond none of the point pairs
is degenerate.}*/
{ 
  int d = s1.dimension(),i; 
  CGAL_assertion_msg(d==s2.dimension(),
    "intersection: dimensions disagree!"); 
  typename LA::Matrix M(d,2),S; 
  typename LA::Vector b(d), lambda(2), c; 
  RT D; 
  
  RT s1w = s1.homogeneous(d); 
  RT t1w = t1.homogeneous(d); 
  RT s2w = s2.homogeneous(d); 
  RT t2w = t2.homogeneous(d); 
  RT g1w = s1w*t1w; 
  RT g2w = s2w*t2w; 
  RT t12w = t1w*t2w; 

  /* init $d \times 2$ - matrix |M| and $d$ - vector |b| */
  for (i = 0; i < d; i++) { 
    M(i,0) = g2w * (t1.homogeneous(i) * s1w - s1.homogeneous(i) * t1w); 
    M(i,1) = g1w * (s2.homogeneous(i) * t2w - t2.homogeneous(i) * s2w); 
    b[i]   = t12w * (s2.homogeneous(i) * s1w - s1.homogeneous(i) * s2w); 
  }

  if (LA::linear_solver(M,b,lambda,D,S,c)) { 
    if (S.column_dimension()>0) return LINE;
    l1 = R::make_FT(lambda[0],D);
    l2 = R::make_FT(lambda[1],D);
    p = s1 + l1 * (t1 - s1);
    return POINT; 
  }
  return NO_INTERSECTION; 
}

};

/*{\Manpage {Line_hyperplane_intersectionHd}{R} 
{intersecting a line and a hyperplane}}*/

template <class R> 
class Line_hyperplane_intersectionHd {
typedef typename R::RT RT;
typedef typename R::FT FT;
typedef typename R::LA LA;
typedef typename R::Point_d      Point_d;
typedef typename R::Hyperplane_d Hyperplane_d;

public:
enum Intersection_result { NO_INTERSECTION, POINT, LINE };

Intersection_result operator()(const Point_d& s, const Point_d& t,
  const Hyperplane_d& h, Point_d& p, FT& lambda) 
/*{\Mfunop returns |NO_INTERSECTION| if the line represented by |s1t1| and the
hyperplane |h| don't intersect, returns |POINT| if they intersect in a
unique point, and returns LINE if the line is part of the
hyperplane. In the |POINT| case the point of intersection is assigned
to |p|.  Then |p = s1 + lambda * t1-s1|.  \precond the point pair is
not degenerate.}*/
{ 
  CGAL_assertion_msg((h.dimension()==s.dimension() && 
                      h.dimension()==t.dimension()), 
  "Line_hyperplane_intersection_d: dimensions do not agree.");

  int d = h.dimension(),i;
  RT S(0),T(0);
  for (i=0; i<=d; ++i) {
    S += h[i]*s.homogeneous(i);
    T += h[i]*t.homogeneous(i);
  }
  bool s_contained = CGAL_NTS is_zero(S), 
       t_contained = CGAL_NTS is_zero(T);
  if (s_contained && t_contained) { p = s; return LINE; }
  if (s_contained) { p = s; return POINT; }
  if (t_contained) { p = t; return POINT; }
  // now the simple cases are done 

  RT D = S * t.homogeneous(d) - T * s.homogeneous(d);
  if (CGAL_NTS is_zero(D)) return NO_INTERSECTION;

  typename LA::Vector homog(d + 1);
  for (i = 0; i < d; ++i)
    homog[i] = S * t.homogeneous(i) - T * s.homogeneous(i); 
  homog[d] = D;
  p = Point_d(d,homog.begin(),homog.end());
  lambda = R::make_FT(S * t.homogeneous(d), D);
  return POINT;
}

};


} //namespace CGAL

#include <CGAL/Kernel_d/intersection_objects_d.h>

#endif //CGAL_INTERSECTION_OBJECTSHD_H
