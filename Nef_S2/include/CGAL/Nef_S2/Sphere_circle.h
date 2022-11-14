// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel  <seel@mpi-sb.mpg.de>

#ifndef CGAL_SPHERE_CIRCLE_H
#define CGAL_SPHERE_CIRCLE_H

#include <CGAL/license/Nef_S2.h>


#include <CGAL/basic.h>
#include <CGAL/Nef_S2/Sphere_point.h>

namespace CGAL {

template <class R> class Sphere_segment;

/*{\Manpage{Sphere_circle}{R}{Great circles on the unit sphere}{c}}*/

template <class R_> class Sphere_circle : public R_::Plane_3 {

/*{\Mdefinition An object |\Mvar| of type |\Mname| is an oriented
great circle on the surface of a unit sphere.  Such circles correspond
to the intersection of an oriented plane (that contains the origin)
and the surface of $S_2$. The orientation of the great circle is that
of a counterclockwise walk along the circle as seen from the positive
halfspace of the oriented plane.}*/

public:

/*{\Mtypes 5}*/
typedef R_ R;
/*{\Mtypemember representation class.}*/
typedef typename R::RT RT;
/*{\Mtypemember ring type.}*/
typedef std::pair< Sphere_segment<R>,Sphere_segment<R> >
  Sphere_segment_pair;
/*{\Mtypemember sphere segment pair.}*/

typedef typename R_::Plane_3 Plane_3;
typedef typename R_::Line_3 Line_3;
typedef typename R_::Point_3 Point_3;
typedef typename R_::Vector_3 Vector_3;
typedef Sphere_circle<R_> Self;
typedef typename R_::Plane_3 Base;

/*{\Mcreation 5}*/
Sphere_circle() : Base() {}
/*{\Mcreate creates some great circle.}*/

Sphere_circle(const Sphere_point<R>& p, const Sphere_point<R>&q)
  : Base(Point_3(0,0,0),p,q)
/*{\Mcreate creates a great circle through $p$ and $q$.  If $p$ and
$q$ are not antipodal on $S_2$, then this circle is unique and oriented
such that a walk along |\Mvar| meets $p$ just before the shorter segment
between $p$ and $q$. If $p$ and $q$ are antipodal of each other then we
create any great circle that contains $p$ and $q$.}*/
{
  if ( p != q.antipode() ) {
    Point_3 po = CGAL::ORIGIN + Base::orthogonal_vector();
    if (R_().orientation_3_object()(CGAL::ORIGIN, Point_3(p),
                                    Point_3(q), po) != CGAL::POSITIVE )
      *this = opposite();
  } else {
    /* previous method was: *this = Self(Plane_3((0,0,0),q-p));
       but p, q don't belong to the plane ((0,0,0), q-p) */
    if(!Line_3(p,q).has_on(Point_3(1,0,0)))
      *this = Self(Plane_3(p,q,Point_3(1,0,0)));
    else
      *this = Self(Plane_3(p,q,Point_3(0,1,0)));
    /* take one point that doesn't belong to the line (p, q-p) */
  }
}

Sphere_circle(const Plane_3& h) : Base(h)
/*{\Mcreate creates the circle of $S_2$ corresponding to the plane
|h|. If |h| does not contain the origin, then |\Mvar| becomes the
circle parallel to |h| containing the origin.}*/
{
  if(!is_zero(h.d()))
    *this = Plane_3(h.a(),h.b(),h.c(),RT(0));
}

Sphere_circle(const Origin& o, const Vector_3& v) : Base(o,v) {}

Sphere_circle(const Origin&, const Plane_3& h) : Base(h.a(),h.b(),h.c(),RT(0))
{/* Even if |h| does not contain the origin, the circle will contain the origin
and be parallel to |h| */}

Sphere_circle(const RT& x, const RT& y, const RT& z) : Base(x,y,z,0) {}
/*{\Mcreate creates the circle orthogonal to the vector $(x,y,z)$.}*/

Sphere_circle(Sphere_circle<R> c, const Sphere_point<R>& p)
/*{\Mcreate creates a great circle orthogonal to $c$ that contains $p$.
\precond $p$ is not part of $c$.}*/
{ CGAL_assertion(!c.has_on(p));
  if ( c.has_on_negative_side(p) ) c=c.opposite();
  if ( p == c.orthogonal_pole() )
    *this = Sphere_circle<R>(Base(Point_3(0,0,0),p,CGAL::ORIGIN+c.base1()));
  else
    *this = Sphere_circle<R>(Base(Point_3(0,0,0),p,c.orthogonal_pole()));
}

/*{\Moperations 4 2}*/
Sphere_circle<R> opposite() const
/*{\Mop returns the opposite of |\Mvar|.}*/
{ return Sphere_circle<R>(Base::opposite(),Assume_d_equal_0{}); }

bool has_on(const Sphere_point<R>& p) const
/*{\Mop returns true iff |\Mvar| contains |p|.}*/
{ return Base::has_on(p); }

Plane_3 plane() const { return Base(*this); }
/*{\Mop returns the plane supporting |\Mvar|.}*/

Plane_3 plane_through(const Point_3& p) const
/*{\Mop returns the plane parallel to |\Mvar| that
contains point |p|.}*/
{ return Plane_3(p,Base::orthogonal_vector()); }

Sphere_point<R> orthogonal_pole() const
/*{\Mop returns the point that is the pole of the
hemisphere left of |\Mvar|.}*/
{ return CGAL::ORIGIN+Base::orthogonal_vector(); }

Sphere_segment_pair split_at(const Sphere_point<R>& p) const;
/*{\Mop returns the pair of circle segments that is the result
of splitting |\Mvar| at |p| and |p.antipode()|.}*/

Sphere_segment_pair split_at_xy_plane() const;
/*{\Mop returns the pair of circle segments that is the result
of splitting |\Mvar| at the $x$-$y$-coordinate plane if |\Mvar|
is not part of it. Otherwise |\Mvar| is split at the
$x$-$z$-coordinate plane.}*/

private:
struct Assume_d_equal_0{};
Sphere_circle(const Plane_3& h, Assume_d_equal_0) : Base(h){}

}; // Sphere_circle<R>

/*{\Mtext\headerline{Global functions}}*/

template <class R>
bool equal_as_sets(const CGAL::Sphere_circle<R>& c1,
                   const CGAL::Sphere_circle<R>& c2)
/*{\Mfunc returns true iff |c1| and |c2| are equal as unoriented
circles.}*/
{ return c1==c2 || c1==c2.opposite(); }

template <class R>
bool equal_not_opposite(const CGAL::Sphere_circle<R>& c1,
                        const CGAL::Sphere_circle<R>& c2) {
  // function should be called to decide whether two circles
  // are equal or opposites. returns true iff |c1| and |c2| are equal
  if(c1.a() != 0) return CGAL::sign(c1.a()) == CGAL::sign(c2.a());
  if(c1.b() != 0) return CGAL::sign(c1.b()) == CGAL::sign(c2.b());
  return CGAL::sign(c1.c()) == CGAL::sign(c2.c());
}

template <typename R>
Sphere_point<R> intersection(const Sphere_circle<R>& c1,
                             const Sphere_circle<R>& c2)
/*{\Mfunc returns one of the two intersection points of
|c1| and |c2|. \precond |c1 != c2| as sets.}*/
{
  CGAL_assertion(!equal_as_sets(c1,c2));
  CGAL_NEF_TRACEN("circle_intersection "<<c1<<" "<<c2);
  return R().construct_cross_product_vector_3_object()(
         c1.orthogonal_vector(), c2.orthogonal_vector());
}


} //namespace CGAL
#endif //CGAL_SPHERE_CIRCLE_H
