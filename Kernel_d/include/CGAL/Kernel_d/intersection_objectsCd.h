// Copyright (c) 2002
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : ?

#ifndef CGAL_INTERSECTION_OBJECTSCD_H
#define CGAL_INTERSECTION_OBJECTSCD_H

#include <CGAL/basic.h>
#include <CGAL/Kernel_d/debug.h>

namespace CGAL {

/*{\Manpage{Line_line_intersectionCd}{R}{intersecting two lines}}*/

template <class R>
class Line_line_intersectionCd {

typedef typename R::FT FT;
typedef typename R::LA LA;
typedef typename R::Point_d Point_d;
typedef typename R::Line_d Line_d;

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
  FT D;

  /* init $d \times 2$ - matrix |M| and $d$ - vector |b| */
  for (i = 0; i < d; i++) {
    M(i,0) = t1.cartesian(i) - s1.cartesian(i);
    M(i,1) = s2.cartesian(i) - t2.cartesian(i);
    b[i]   = s2.cartesian(i) - s1.cartesian(i);
  }

  if (LA::linear_solver(M,b,lambda,D,S,c)) {
    if ( S.column_dimension()>0 ) return LINE;
    l1 = lambda[0]; l2 = lambda[1];
    p = s1 + l1 * (t1 - s1);
#ifdef CGAL_CHECK_EXACTNESS
    Line_d L1(s1,t1), L2(s2,t2);
    CGAL_assertion(L1.has_on(p)&&L2.has_on(p));
#endif
    return POINT;
  }
  return NO_INTERSECTION;
}
};

/*{\Manpage {Line_hyperplane_intersectionCd}{R}
{intersecting a line and a hyperplane}}*/

template <class R>
class Line_hyperplane_intersectionCd {

typedef typename R::FT FT;
typedef typename R::LA LA;
typedef typename R::Point_d Point_d;
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
  FT S = h.value_at(s), T = h.value_at(t);

  bool s_contained = CGAL_NTS is_zero(S),
       t_contained = CGAL_NTS is_zero(T);
  if (s_contained && t_contained) { p = s; return LINE; }
  if (s_contained) { p = s; return POINT; }
  if (t_contained) { p = t; return POINT; }
  // now the simple cases are done

  FT D = S - T;
  if ( CGAL_NTS is_zero(D) ) return NO_INTERSECTION;

  typename LA::Vector v(d);
  for (i = 0; i < d; ++i)
    v[i] = (S * t.cartesian(i) - T * s.cartesian(i))/D;
  p = Point_d(d,v.begin(),v.end()); lambda = S/D;

#ifdef CGAL_CHECK_EXACTNESS
  Line_d l(s,t);
  CGAL_assertion(h.has_on(p)&&l.has_on(p));
#endif

  return POINT;
}

};

} //namespace CGAL

#include <CGAL/Kernel_d/intersection_objects_d.h>

#endif //CGAL_INTERSECTION_OBJECTSCD_H
