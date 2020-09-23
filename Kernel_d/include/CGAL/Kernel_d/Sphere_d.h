// Copyright (c) 2000,2001
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
// Author(s)     : Michael Seel

#ifndef CGAL_SPHERE_D_H
#define CGAL_SPHERE_D_H

#include <CGAL/basic.h>
#include <vector>
#include <CGAL/Dimension.h>
#include <CGAL/use.h>
#include <CGAL/Handle.h>
#include <CGAL/Kernel_d/Tuple_d.h>
#include <CGAL/Kernel_d/Aff_transformation_d.h>
#include <CGAL/Kernel_d/Vector_d.h>

namespace CGAL {

template <class R> class Sphere_d;
template <class R> bool equal_as_sets(const Sphere_d<R>&, const Sphere_d<R>&);

template <class R>
class  Sphere_d_rep  {

  typedef typename R::Point_d Point_d;

  friend class Sphere_d<R>;
  friend bool equal_as_sets <>
    (const Sphere_d<R>&, const Sphere_d<R>&);

  std::vector< Point_d > P; // d+1 defining points, index range 0-d
  Orientation orient;       // orientation(P)
  Point_d* cp;              // pointer to center (lazy calculated)

public:
  Sphere_d_rep() : cp(0) {}
  Sphere_d_rep(int d)  : P(d), cp(0) {}

  template <class ForwardIterator>
  Sphere_d_rep(int d, ForwardIterator first, ForwardIterator last) :
     P(first,last), cp(0)
  { TUPLE_DIM_CHECK(P.begin(),P.end(),Sphere_d);
    CGAL_assertion(d+1==int(P.size()));
    CGAL_USE(d);
    typename R::Orientation_d orientation_;
    orient = orientation_(P.begin(),P.end()); }

  ~Sphere_d_rep() { if (cp) delete cp; }

};  // Sphere_d_rep<R>

/*{\Manpage {Sphere_d}{R}{Simple Spheres}{S}}*/

template <class R_>
class Sphere_d : public Handle_for< Sphere_d_rep<R_> > {

/*{\Mdefinition
An instance $S$ of the data type |Sphere_d| is an oriented sphere in
some $d$-dimensional space. A sphere is defined by $d+1$ points with
rational coordinates (class |Point_d<R>|). We use $A$ to denote the
array of the defining points.  A set $A$ of defining points is
\emph{legal} if either the points are affinely independent or if the
points are all equal. Only a legal set of points defines a sphere in
the geometric sense and hence many operations on spheres require the
set of defining points to be legal.  The orientation of $S$ is equal
to the orientation of the defining points, i.e., |orientation(A)|. }*/

public:
  typedef CGAL::Dynamic_dimension_tag Ambient_dimension;
  typedef CGAL::Dynamic_dimension_tag Feature_dimension;

/*{\Mtypes 4}*/

typedef Sphere_d_rep<R_>  Rep;
typedef Handle_for<Rep>   Base;
typedef Sphere_d<R_>      Self;
typedef typename R_::Point_d Point_d;

using Base::ptr;

Sphere_d(const Base& b) : Base(b) {}

typedef R_ R;
/*{\Mtypemember the representation type.}*/

typedef typename R::RT RT;
/*{\Mtypemember the ring type.}*/

typedef typename R::FT FT;
/*{\Mtypemember the field type.}*/

typedef typename R::LA LA;
/*{\Mtypemember the linear algebra layer.}*/

typedef typename std::vector< Point_d >::const_iterator point_iterator;
/*{\Mtypemember a read-only iterator for points defining the sphere.}*/

/*{\Mcreation 4}*/

Sphere_d(int d = 0) : Base( Rep(d+1) )
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname|. |\Mvar|
is initialized to the empty sphere centered at the origin of
$d$-dimensional space. }*/
{
  Point_d p(d);
  for (int i = 0; i <= d; i++) ptr()->P[i] = p;
  ptr()->orient = ZERO;
  ptr()->cp = new Point_d(p);
}

template <class ForwardIterator>
Sphere_d(int d, ForwardIterator first, ForwardIterator last) :
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname|. |\Mvar| is
initialized to the sphere through the points in |A = set [first,last)|.
\precond $A$ consists of $d+1$ $d$-dimensional points.}*/
  Base( Rep(d,first,last) ) {}


/*{\Moperations 4 3}*/

int dimension() const
/*{\Mop returns the dimension of |\Mvar|.}*/
  { return static_cast<int>(ptr()->P.size()) - 1; }

Point_d point(int i) const
/*{\Mop returns the $i$-th defining point. \precond $0 \le i \le |dim|$.}*/
{ CGAL_assertion_msg((0<=i && i<=dimension()),
    "Sphere_d::point(): index out of range.");
  return ptr()->P[i];
}

point_iterator points_begin() const { return ptr()->P.begin(); }
/*{\Mop returns an iterator pointing to the first defining point.}*/

point_iterator points_end() const { return ptr()->P.end(); }
/*{\Mop returns an iterator pointing beyond the last defining point.}*/

bool is_degenerate() const { return (ptr()->orient == CGAL::ZERO); }
/*{\Mop returns true iff the defining points are not full dimenional.}*/

bool is_legal() const
/*{\Mop returns true iff the set of defining points is legal.
A set of defining points is legal iff their orientation is
non-zero or if they are all equal.}*/
{ if (ptr()->orient != ZERO ) return true;
  const std::vector< Point_d >& A = ptr()->P;
  Point_d po = A[0];
  for (int i = 1; i < int(A.size()); ++i)
    if (A[i] != po) return false;
  return true;
}

Point_d center() const
/*{\Mop  returns the center of |\Mvar|. \precond The orientation
of |\Mvar| is non-zero. }*/
{
  if (ptr()->cp == 0) {
    if (ptr()->orient == 0) {
      const std::vector< Point_d >& A = ptr()->P;
      Point_d po = A[0];
      for (int i = 1; i < int(A.size()); ++i)
        if (A[i] != po)
          CGAL_error_msg("Sphere_d::center(): points are illegal.");
      const_cast<Self&>(*this).ptr()->cp = new Point_d(A[0]);
      return *(ptr()->cp);
    }
    typename R::Center_of_sphere_d center_of_sphere_;
    const_cast<Self&>(*this).ptr()->cp =
      new Point_d(center_of_sphere_(points_begin(),points_end()));
  }
  return *(ptr()->cp);
}



FT squared_radius() const
/*{\Mop returns the squared radius of the sphere.}*/
{ if (is_degenerate()) return 0;
  return (point(0)-center()).squared_length();
}

Orientation orientation()  const { return ptr()->orient; }
/*{\Mop returns the orientation of |\Mvar|.}*/

Oriented_side oriented_side(const Point_d& p) const
/*{\Mop returns either the constant |ON_ORIENTED_BOUNDARY|,
|ON_POSITIVE_SIDE|, or |ON_NEGATIVE_SIDE|, iff p lies on the boundary,
properly on the positive side, or properly on the negative side of
circle, resp.}*/
{ typename R::Side_of_oriented_sphere_d side;
  return side(points_begin(),points_end(),p); }

Bounded_side bounded_side(const Point_d& p) const
/*{\Mop returns |ON_BOUNDED_SIDE|, |ON_BOUNDARY|, or
|ON_UNBOUNDED_SIDE| iff p lies properly inside, on
 the boundary, or properly outside of circle, resp.}*/
{ typename R::Side_of_bounded_sphere_d side;
  return side(points_begin(),points_end(),p); }

bool has_on_positive_side (const Point_d& p) const
/*{\Mop returns |\Mvar.oriented_side(p)==ON_POSITIVE_SIDE|.}*/
{ return oriented_side(p) == ON_POSITIVE_SIDE; }

bool has_on_negative_side (const Point_d& p) const
/*{\Mop returns |\Mvar.oriented_side(p)==ON_NEGATIVE_SIDE|.}*/
{ return oriented_side(p) == ON_NEGATIVE_SIDE; }

bool has_on_boundary (const Point_d& p) const
/*{\Mop returns |\Mvar.oriented_side(p)==ON_ORIENTED_BOUNDARY|,
which is the same as |\Mvar.bounded_side(p)==ON_BOUNDARY|.}*/
{ return oriented_side(p) == ON_ORIENTED_BOUNDARY; }

bool has_on_bounded_side (const Point_d& p) const
/*{\Mop returns |\Mvar.bounded_side(p)==ON_BOUNDED_SIDE|.}*/
{ return (int(ptr()->orient) * int(oriented_side(p))) > 0 ; }

bool has_on_unbounded_side (const Point_d& p) const
/*{\Mop returns |\Mvar.bounded_side(p)==ON_UNBOUNDED_SIDE|.}*/
{ return (int(ptr()->orient) * int(oriented_side(p))) < 0; }

Sphere_d<R> opposite() const
/*{\Mop returns the sphere with the same center and squared
  radius as |\Mvar| but with opposite orientation.}*/
{ CGAL_assertion(dimension()>1);
  if (is_degenerate()) return *this;
  std::vector< Point_d > V(points_begin(),points_end());
  std::swap(V[0],V[1]);
  return Sphere_d<R>(dimension(),V.begin(),V.end());
}

Sphere_d<R> transform(const Aff_transformation_d<R>& t) const
/*{\Mopl returns $t(s)$. }*/
{ std::vector< Point_d > B(points_begin(),points_end());
  typename std::vector< Point_d >::iterator it;
  for (it = B.begin(); it != B.end(); ++it)
    *it = it->transform(t);
  return Sphere_d<R>(dimension(),B.begin(),B.end());
}

Sphere_d<R> operator+(const Vector_d<R>& v) const
/*{\Mbinop returns the sphere translated by |v|. }*/
{ std::vector< Point_d > B(points_begin(),points_end());
  typename std::vector< Point_d >::iterator it;
  for (it = B.begin(); it != B.end(); ++it) it->operator+=(v);
  return Sphere_d<R>(dimension(),B.begin(),B.end());
}

bool operator==(const Sphere_d<R>& D) const
{ if (this->identical(D)) return true;
  if (dimension() != D.dimension()) return false;
  return (center()==D.center() &&
          squared_radius() == D.squared_radius() &&
          orientation() == D.orientation());
}

bool operator!=(const Sphere_d<R>& D) const
{ return !operator==(D); }


}; // end of class Sphere_d

/*{\Mtext \headerline{Non-Member Functions} }*/
template <class R>
bool weak_equality(const Sphere_d<R>& s1, const Sphere_d<R>& s2)
/*{\Mfunc Test for equality as unoriented spheres.}*/
{ if (s1.identical(s2)) return true;
  if (s1.dimension() != s2.dimension()) return false;
  return (s1.center()==s2.center() &&
          s1.squared_radius() == s2.squared_radius());
}

/*{\Mimplementation Spheres are implemented by a vector of points as
an item type.  All operations like creation, initialization, tests,
input and output of a sphere $s$ take time
$O(|s.dimension()|)$. |dimension()|, point access take constant time.
The space requirement for spheres is $O(|s.dimension()|)$
neglecting the storage room of the points.}*/

template <class R>
std::ostream& operator<<(std::ostream& os, const CGAL::Sphere_d<R>& s)
{ typedef typename Sphere_d<R>::point_iterator iterator;
  os << s.dimension()+1 << " ";
  for (iterator it=s.points_begin(); it!=s.points_end(); ++it)
    os << *it << " ";
  return os;
}

template <class R> std::istream&
operator>>(std::istream& is, CGAL::Sphere_d<R>& s)
{ int d; is >> d;
  std::vector< Point_d<R> > V(d);
  Point_d<R> p;
  while ( d-- ) {
    if (!(is >> p)) return is;
    V[d] = p;
  }
  s = Sphere_d<R>(static_cast<int>(V.size())-1, V.begin(), V.end() );
  return is;
}


/*
The center is only defined if the set of defining points are
legal. If the defining points are all equal the sphere is trivial. So
assume otherwise. Then the center $c$ is the unique point with equal
distance to all the defining points. A point $c$ has equal distance to
point $p_0$ and $p_i$ if it lies on the perpendicual bisector of $p_d$
and $p_i$, i.e., if it satiesfies the plane equation $(p_i - p_d)^T c
= (p_i - p_d) (p_i + p_d)/2$. Note that $p_i - p_d$ is the normal
vector of the bisector hyperplane and $(p_i + p_d)/2$ is the midpoint
of $p_d$ and $p_i$. The equation above translates into the equation \[
\sum_{0 \le j < d} 2*p_{dd}p_{id}(p_{ij}p_{dd} - p_{dj}p_{id})c_j/c_d
= \sum_{0 \le j < d} (p_{ij}p_{dd} - p_{dj}p_{id})(p_{ij}p_{dd} +
p_{dj}p_{id}) \] for the homogeneous coordinates of the points and the
center. We may tentatively assume that $c_d = 1$, solve the
corresponding linear system, and then define the center.
*/

} //namespace CGAL

#endif // CGAL_SPHERE_D_H
//----------------------- end of file ----------------------------------
