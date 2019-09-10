// Copyright (c) 2000,2001  
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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Michael Seel

#ifndef CGAL_RAY_D_H
#define CGAL_RAY_D_H

#include <CGAL/Kernel_d/Pair_d.h>
#include <CGAL/Kernel_d/Point_d.h> 
#include <CGAL/Kernel_d/Direction_d.h> 
#include <CGAL/Kernel_d/Segment_d.h> 
#include <CGAL/Kernel_d/Line_d.h>
#include <CGAL/Kernel_d/Aff_transformation_d.h>
#include <CGAL/Dimension.h>

namespace CGAL {

template <class R>
std::istream& operator>>(std::istream&, Ray_d<R>&);
template <class R>
std::ostream& operator<<(std::ostream&, const Ray_d<R>&);


/*{\Manpage {Ray_d}{R}{Rays in d-space}{r}}*/

template <class p_R>
class Ray_d : public Handle_for< Pair_d<p_R> > { 
  typedef Pair_d<p_R>       Pair;
  typedef Handle_for<Pair>  Base;
  typedef Ray_d<p_R>        Self;

  using Base::ptr;

/*{\Mdefinition
An instance of data type |Ray_d| is a ray in $d$-dimensional
Euclidian space. It starts in a point called the source of |\Mvar| and
it goes to infinity.}*/

public: 

  typedef CGAL::Dynamic_dimension_tag Ambient_dimension;
  typedef CGAL::Dimension_tag<1>      Feature_dimension;

/*{\Mtypes 4}*/
typedef p_R R;
/*{\Mtypemember the representation type.}*/
typedef typename p_R::RT RT;
/*{\Mtypemember the ring type.}*/
typedef typename p_R::FT FT;
/*{\Mtypemember the field type.}*/
typedef typename p_R::LA LA;
/*{\Mtypemember the linear algebra layer.}*/

typedef typename Vector_d<R>::Base_vector Base_vector;

friend class Line_d<R>; 
friend class Segment_d<R>; 

private:
Ray_d(const Base& b) : Base(b) {}
public: 
/*{\Mcreation 3}*/

Ray_d() : Base( Pair() ) {}
/*{\Mcreate introduces some ray in $d$-dimensional space }*/
 
Ray_d(const Point_d<R>& p, const Point_d<R>& q)
/*{\Mcreate introduces a ray through |p| and |q| and starting at |p|.
\precond $p$ and $q$ are distinct and have the same dimension. }*/
 : Base( Pair(p,q) )
{ CGAL_assertion_msg(!ptr()->is_degenerate(), 
    "Ray_d::constructor: the two points must be different." );
  CGAL_assertion_msg((p.dimension()==q.dimension()), 
    "Ray_d::constructor: the two points must have the same dimension." );
}

Ray_d(const Point_d<R>& p, const Direction_d<R>& dir)
/*{\Mcreate introduces a ray starting in |p| with direction |dir|.
\precond |p| and |dir| have the same dimension and |dir| is not
trivial.}*/
  : Base( Pair(p,p+dir.vector()) )
{ CGAL_assertion_msg((p.dimension()==dir.dimension()), 
    "Ray_d::constructor: the p and dir must have the same dimension." );
  CGAL_assertion_msg(!dir.is_degenerate(), 
    "Ray_d::constructor: dir must be non-degenerate." );
}

Ray_d(const Segment_d<R>& s) 
/*{\Mcreate introduces a ray through |s.source()| and |s.target()| and 
starting at |s.source()|. \precond $s$ is not trivial. }*/
  : Base( s ) 
{ CGAL_assertion_msg(!s.is_degenerate(), 
    "Ray_d::constructor: segment is trivial.");
}

Ray_d(const Ray_d<R>& r) : Base(r) {}

/*{\Moperations 3 3}*/

int dimension() const { return (ptr()->_p[0].dimension()); }
/*{\Mop returns the dimension of the underlying space.}*/

Point_d<R> source() const { return (ptr()->_p[0]); }
/*{\Mop returns the source point of |\Mvar|. }*/

Point_d<R> point(int i) const 
/*{\Mop returns a point on |\Mvar|. |point(0)| is the source.
|point(i)|, with $i>0$, is different from the source. \precond $i
\geq 0$.}*/ 
{ return (ptr()->_p[i%2]); }

Direction_d<R> direction() const 
/*{\Mop returns the direction of |\Mvar|. }*/
{ return ptr()->direction(); }

inline Line_d<R> supporting_line() const; 
/*{\Mop returns the supporting line of |\Mvar|.}*/

Ray_d<R> opposite() const
/*{\Mop returns the ray with direction opposite to |\Mvar|
and starting in |source|.}*/
{ return Ray_d<R>(source(),-direction()); }

Ray_d<R> transform(const Aff_transformation_d<R>& t) const
/*{\Mop returns $t(l)$. }*/
{ return Ray_d<R>(point(0).transform(t),point(1).transform(t)); }

Ray_d<R> operator+(const Vector_d<R>& v) const
/*{\Mbinop returns |\Mvar+v|, i.e., |\Mvar| translated by vector $v$.}*/ 
{ return Ray_d<R>(point(0)+v, point(1)+v); }

bool has_on(const Point_d<R>& p) const 
/*{\Mop A point is on |r|, iff it is equal to the source of |r|, or if it is
in the interior of |r|.}*/
{ typename R::Position_on_line_d pos; FT l;
  if (pos(p,point(0),point(1),l)) return (FT(0)<=l);
  return false;
}

/*{\Mtext \headerline{Non-Member Functions}}*/

bool operator==(const Ray_d<R>& r1) const
{ if ( this->identical(r1) ) return true;
  if ( dimension() != r1.dimension() ) return false;
  return source() == r1.source() && 
         direction() == r1.direction(); 
}

bool operator!=(const Ray_d<R>& r1)
{ return !operator==(r1); }

friend std::istream& operator>> <> 
(std::istream&, Ray_d<R>&);
friend std::ostream& operator<< <> 
(std::ostream&, const Ray_d<R>&); 

}; // end of class

template <class R>
bool parallel(const Ray_d<R>& r1, const Ray_d<R>& r2)
/*{\Mfunc returns true if the unoriented supporting lines of |r1| and |r2|
are parallel and false otherwise. }*/
{ return (r1.direction() == r2.direction()) || 
         (r1.direction() == -(r2.direction())); 
} 

template <class R>
std::istream& operator>>(std::istream& I, Ray_d<R>& r) 
{ r.copy_on_write(); r.ptr()->read(I); 
  CGAL_assertion_msg(r.point(0)!=r.point(1),
    "Line_d::operator>>: trivial ray.");
  CGAL_assertion_msg(r.point(0).dimension()==r.point(1).dimension(),
  "Ray_d::operator>>: dimensions disagree.");
  return I; 
}

template <class R>
std::ostream& operator<<(std::ostream& O, const Ray_d<R>& r)
{ r.ptr()->print(O,"Ray_d"); return O; }

/*{\Mimplementation 
Rays are implemented by a pair of points as an item type.  All
operations like creation, initialization, tests, direction
calculation, input and output on a ray $r$ take time
$O(|r.dimension()|)$. |dimension()|, coordinate and point access, and
identity test take constant time. The space requirement is
$O(|r.dimension()|)$.}*/


} //namespace CGAL
#endif // CGAL_RAYHD_H
//----------------------- end of file ----------------------------------
