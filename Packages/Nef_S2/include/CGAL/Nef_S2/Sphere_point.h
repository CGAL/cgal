#ifndef CGAL_SPHERE_POINT_H
#define CGAL_SPHERE_POINT_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

/*{\Mtext \headerline{Restricted Spherical Geometry}

We introduce geometric objects that are part of the
spherical surface $S_2$ and operations on them. We define types
|Sphere_point|, |Sphere_circle|, |Sphere_segment|, and
|Sphere_direction|.  |Sphere_point|s are points on $S_2$,
|Sphere_circle|s are oriented great circles of $S_2$,
|Sphere_segment|s are oriented parts of |Sphere_circles| bounded by a
pair of |Sphere_point|s, and |Sphere_direction|s are directions that
are part of great circles (a direction is usually defined to be a
vector without length, that floats around in its underlying space and
can be used to specify a movement at any point of the underlying
space; in our case we use directions only at points that are part of
the great circle that underlies also the direction.)

Note that we have to consider special geometric properties of the
objects. For example two points that are part of a great circle define
two |Sphere_segment|s, and two arbitrary |Sphere_segment|s can
intersect in two points. 

If we restrict our geometric objects to a so-called perfect hemisphere
of $S_2$\footnote{A perfect hemisphere of $S_2$ is an open half-sphere
plus an open half-circle in the boundary of the open half-sphere plus one
endpoint of the half-circle.} then the restricted objects behave like
in classical geometry, e.g., two points define exactly one segment,
two segments intersect in at most one interior point
(non-degenerately), or three non-cocircular sphere points can be
qualified as being positively or negatively oriented.}*/

/*{\Moptions print_title=yes }*/ 
/*{\Manpage{Sphere_point}{R}{Points on the unit sphere}{p}}*/ 
 
template <class R_> 
class Sphere_point : public R_::Point_3 {

/*{\Mdefinition An object |\Mvar| of type |\Mname| is a point on the
surface of a unit sphere. Such points correspond to the nontrivial
directions in space and similarly to the equivalence classes of all
nontrivial vectors under normalization.}*/ 

public:
/*{\Mtypes 5}*/
typedef R_ R;
/*{\Mtypemember representation class.}*/
typedef typename R_::FT FT;
/*{\Xtypemember ring number type.}*/
typedef typename R_::RT RT;
/*{\Mtypemember field number type.}*/
typedef Sphere_point<R>    Self;
typedef typename R::Point_3     Base;
typedef typename R::Vector_3    Vector_3;
typedef typename R::Direction_3 Direction_3;

/*{\Mcreation 5}*/
Sphere_point() : Base() {}
/*{\Mcreate creates some sphere point.}*/

Sphere_point(int x, int y, int z) :
Base(x,y,z,1) { CGAL_assertion(x!=0 || y!=0 || z!=0); }

Sphere_point(const RT& x, const RT& y, const RT& z) :
/*{\Mcreate creates a sphere point corresponding to the point of
intersection of the ray starting at the origin in direction $(x,y,z)$
and the surface of $S_2$.}*/
  Base(x,y,z,1) { CGAL_assertion(x!=0 || y!=0 || z!=0); }

Sphere_point(const Base& p) : Base(p) {} 
Sphere_point(const Vector_3& v) : Base(CGAL::ORIGIN+v) {} 
Sphere_point(const Direction_3& d) : 
  Base(CGAL::ORIGIN+d.vector()) {} 

/*{\Moperations 4 2}*/

/*{\Mtext Access to the coordinates is provided by the following
operations. Note that the vector $(x,y,z)$ is not normalized.}*/

RT x() const { return hx(); }
/*{\Mop the $x$-coordinate.}*/
RT y() const { return hy(); }
/*{\Mop the $y$-coordinate.}*/
RT z() const { return hz(); }
/*{\Mop the $z$-coordinate.}*/

bool operator==(const Sphere_point<R>& q) const 
/*{\Mbinop Equality.}*/
{ return Direction_3(Base(*this)-ORIGIN)==
         Direction_3(q-ORIGIN); }

bool operator!=(const Sphere_point<R>& q) const
/*{\Mbinop Inequality.}*/
{ return !operator==(q); }

Sphere_point<R> antipode() const 
/*{\Mop returns the antipode of |\Mvar|.}*/
{ return ORIGIN + -(Base(*this)-ORIGIN); }

/*
Sphere_point<R> simplified() const {
  int d = (hx().degree() > hy().degree() ? hx().degree() : hy().degree());
  d = (d > hz().degree() ? d : hz().degree());
  CGAL_assertion(d < 3);
  return Sphere_point(x()(d),y()(d),z()(d));  
}
*/

Sphere_point<R> normalized() const {

  RT g = (hx()==0) ? ((hy()==0) ? ((hz()==0) ? 1: hz()): hy()): hx();
  
  if(hy() != 0) g = gcd(g,hy());
  if(hz() != 0) g = gcd(g,hz());
	  
  g = g.abs();

  RT x = hx()/g;
  RT y = hy()/g;
  RT z = hz()/g;

  return Sphere_point(Base(x,y,z,1));
}

}; // Sphere_point<R>

template <typename R>
CGAL::Point_3<R> operator+(
  const CGAL::Point_3<R>& p, const Sphere_point<R>& q)
{ return p + (q-CGAL::ORIGIN); }

CGAL_END_NAMESPACE
#endif //CGAL_SPHERE_POINT_H
