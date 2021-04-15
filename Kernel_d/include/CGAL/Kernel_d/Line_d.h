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

#ifndef CGAL_LINE_D_H
#define CGAL_LINE_D_H

#include <CGAL/Kernel_d/Pair_d.h>
#include <CGAL/Kernel_d/Point_d.h>
#include <CGAL/Kernel_d/Direction_d.h>
#include <CGAL/Kernel_d/Segment_d.h>
#include <CGAL/Kernel_d/Ray_d.h>
#include <CGAL/Kernel_d/Aff_transformation_d.h>
#include <CGAL/Dimension.h>

namespace CGAL {

template <class R>
std::istream& operator>>(std::istream&, Line_d<R>&);
template <class R>
std::ostream& operator<<(std::ostream&, const Line_d<R>&);

/*{\Manpage {Line_d}{R}{Lines in d-space}{l}}*/

template <class p_R>
class Line_d : public Handle_for< Pair_d<p_R> > {
  typedef Pair_d<p_R>      Pair;
  typedef Handle_for<Pair> Base;
  typedef Line_d<p_R>      Self;

  using Base::ptr;

/*{\Mdefinition
An instance of data type |Line_d| is an oriented line in
$d$-dimensional Euclidian space.}*/

public:

  typedef CGAL::Dynamic_dimension_tag Ambient_dimension;
  typedef CGAL::Dimension_tag<1>      Feature_dimension;

/*{\Mtypes 5}*/
typedef p_R R;
/*{\Mtypemember the representation type.}*/
typedef typename p_R::RT RT;
/*{\Mtypemember the ring type.}*/
typedef typename p_R::FT FT;
/*{\Mtypemember the field type.}*/
typedef typename p_R::LA LA;
/*{\Mtypemember the linear algebra layer.}*/

typedef typename Vector_d<R>::Base_vector Base_vector;

friend class Ray_d<R>;
friend class Segment_d<R>;

private:
Line_d(const Base& b) : Base(b) {}
public:
/*{\Mcreation 3}*/

Line_d() : Base( Pair() ) {}
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname| and
initializes it to some line in $d$ - dimensional space }*/

Line_d(const Point_d<R>& p, const Point_d<R>& q)
/*{\Mcreate introduces a line through |p| and |q| and oriented
from |p| to |q|. \precond $p$ and $q$ are distinct and have the same
dimension.}*/
 : Base( Pair(p,q) )
{ CGAL_assertion_msg(!ptr()->is_degenerate(),
    "Line_d::constructor: the two points must be different." );
  CGAL_assertion_msg((p.dimension()==q.dimension()),
    "Line_d::constructor: the two points must have the same dimension." );
}

Line_d(const Point_d<R>& p, const Direction_d<R>& dir)
/*{\Mcreate introduces a line through |p| with direction |dir|.
\precond |p| and |dir| have the same dimension, |dir| is not trivial. }*/
  : Base( Pair(p,p+dir.vector()) )
{ CGAL_assertion_msg((p.dimension()==dir.dimension()),
    "Line_d::constructor: the p and dir must have the same dimension." );
  CGAL_assertion_msg(!dir.is_degenerate(),
    "Line_d::constructor: dir must be non-degenerate." );
}

Line_d(const Segment_d<R>& s)
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname| and
initializes it to the line through |s.source()| and |s.target()|
with direction from |s.source()| to |s.target()|.
\precond $s$ is not degenerate. }*/
  : Base( s )
{ CGAL_assertion_msg((!s.is_degenerate()),
    "Line_d::constructor: segment is trivial.");
}

Line_d(const Ray_d<R>& r) : Base(r) {}
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname| and
initializes it to the line through |r.point(1)| and |r.point(2)|. }*/


/*{\Moperations 3 3}*/

int dimension() const { return (ptr()->_p[0].dimension()); }
/*{\Mop returns the dimension of the underlying space.}*/

Point_d<R> point(int i) const
/*{\Mop returns an arbitrary point on |l|.  It holds that |point(i) ==
point(j)|, iff |i==j|. Furthermore, |l| is directed from |point(i)| to
|point(j)|, for all |i < j|.}*/
{ return (ptr()->_p[i%2]); }

Line_d<R> opposite() const
/*{\Mop returns the line |(point(2),point(1))|. }*/
{ return Line_d<R>(point(1),point(0)); }

Direction_d<R> direction() const
/*{\Mop  returns the direction of |\Mvar|. }*/
{ return ptr()->direction(); }

Line_d<R> transform(const Aff_transformation_d<R> & t) const
/*{\Mop returns $t(l)$. }*/
{ return Line_d<R>(point(0).transform(t),point(1).transform(t)); }

Line_d<R> operator+(const Vector_d<R>& v) const
/*{\Mbinop returns |\Mvar+v|, i.e., |\Mvar| translated by vector $v$.}*/
{ return Line_d<R>(point(0)+v,point(1)+v); }

Point_d<R> projection(const Point_d<R>& p) const
/*{\Mop returns the point of intersection of |\Mvar| with the hyperplane
that is orthogonal to |\Mvar| through |p|. }*/
{ Vector_d<R> v = direction().vector();
  Point_d<R> q = point(0);
  FT lambda = ((p-q) * v) / (v*v);
  Point_d<R> res = q + lambda * v;
  return res;
}

bool has_on(const Point_d<R>& p) const
/*{\Mopl returns true if $p$ lies on |\Mvar| and false otherwise. }*/
{ typename R::Position_on_line_d pos; FT dummy;
  return pos(p,point(0),point(1),dummy); }

bool operator==(const Line_d<R>& l1) const
{ if ( this->identical(l1) ) return true;
  if ( dimension() != l1.dimension() ) return false;
  return has_on(l1.point(0)) &&
         direction() == l1.direction();
}

bool operator!=(const Line_d<R>& l1) const
{ return !operator==(l1); }

friend std::istream& operator>> <>
(std::istream&, Line_d<R>&);
friend std::ostream& operator<< <>
(std::ostream&, const Line_d<R>&);

}; // end of class

/*{\Mtext \headerline{Non-Member Functions} }*/

template <class R>
bool weak_equality(const Line_d<R>& l1, const Line_d<R>& l2)
/*{\Mfunc Test for equality as unoriented lines.}*/
{ if (l1.identical(l2)) return true;
  if (l1.dimension()!=l2.dimension()) return false;
  return (l1.has_on(l2.point(0)) &&
          l1.has_on(l2.point(1)));
}

template <class R>
bool parallel(const Line_d<R>& l1, const Line_d<R>& l2)
/*{\Mfunc returns true if |l1| and |l2| are parallel as unoriented lines
and false otherwise. }*/
{ return (l1.direction() == l2.direction() ||
          l1.direction() == -l2.direction()); }

template <class R>
std::istream& operator>>(std::istream& I, Line_d<R>& l)
{ l.copy_on_write(); l.ptr()->read(I);
  CGAL_assertion_msg(l.point(0)!=l.point(1),
    "Line_d::operator>>: trivial line.");
  CGAL_assertion_msg(l.point(0).dimension()==l.point(1).dimension(),
    "Line_d::operator>>: dimensions disagree.");
  return I;
}

template <class R>
std::ostream& operator<<(std::ostream& O, const Line_d<R>& l)
{ l.ptr()->print(O,"Line_d"); return O; }

/*{\Mimplementation
Lines are implemented by a pair of points as an item type.  All
operations like creation, initialization, tests, direction
calculation, input and output on a line $l$ take time
$O(|l.dimension()|)$. |dimension()|, coordinate and point access, and
identity test take constant time.  The operations for intersection
calculation also take time $O(|l.dimension()|)$. The space requirement
is $O(|l.dimension()|)$.}*/


} //namespace CGAL
#endif // CGAL_LINE_D_H
//----------------------- end of file ----------------------------------
