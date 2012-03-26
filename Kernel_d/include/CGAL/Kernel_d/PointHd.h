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
// 
// Author(s)     : Michael Seel

#ifndef CGAL_POINTHDXXX_H
#define CGAL_POINTHDXXX_H 

#include <CGAL/basic.h>
#include <CGAL/Origin.h>
#include <CGAL/Quotient.h>
#include <CGAL/Kernel_d/Tuple_d.h>
#include <CGAL/Kernel_d/VectorHd.h> 
#include <CGAL/Kernel_d/Aff_transformationHd.h>


namespace CGAL {
#define PointHd PointHd2

template <class RT, class LA> class PointHd;
template <class RT, class LA>
std::istream& operator>>(std::istream&, PointHd<RT,LA>&);
template <class RT, class LA>
std::ostream& operator<<(std::ostream&, const PointHd<RT,LA>&);


/*{\Moptions outfile=Point_d.man}*/
/*{\Manpage {Point_d} {R} {Points in d-space} {p}}*/
/*{\Msubst 
Hd<RT,LA>#_d<R>
PointHd#Point_d
Quotient<RT>#FT
}*/

template <class _RT, class _LA > 
class PointHd : public Handle_for< Tuple_d<_RT,_LA> > { 
  typedef Tuple_d<_RT,_LA> Tuple;
  typedef Handle_for<Tuple> Base;
  typedef PointHd<_RT,_LA> Self;

  using Base::ptr;

/*{\Mdefinition 
An instance of data type |\Mname| is a point of Euclidean space in
dimension $d$. A point $p = (p_0,\ldots,p_{ d - 1 })$ in
$d$-dimensional space can be represented by homogeneous coordinates
$(h_0,h_1,\ldots,h_d)$ of number type |RT| such that $p_i = h_i/h_d$,
which is of type |FT|. The homogenizing coordinate $h_d$ is positive.

We call $p_i$, $0 \leq i < d$ the $i$-th Cartesian coordinate and
$h_i$, $0 \le i \le d$, the $i$-th homogeneous coordinate. We call $d$
the dimension of the point.}*/

const typename _LA::Vector& vector_rep() const { return ptr()->v; }
_RT& entry(int i) { return ptr()->v[i]; }
const _RT& entry(int i) const { return ptr()->v[i]; }
void invert_rep() { ptr()->invert(); }
PointHd(const Base& b) : Base(b) {}

public: 
/*{\Mtypes 4}*/

typedef _RT RT;
/*{\Mtypemember the ring type.}*/
typedef Quotient<_RT> FT;
/*{\Mtypemember the field type.}*/
typedef _LA LA;
/*{\Mtypemember the linear algebra layer.}*/
typedef typename Tuple::Cartesian_const_iterator Cartesian_const_iterator;
/*{\Mtypemember a read-only iterator for the cartesian coordinates.}*/
typedef typename Tuple::const_iterator Homogeneous_const_iterator;
/*{\Mtypemember a read-only iterator for the homogeneous coordinates.}*/

friend class VectorHd<RT,LA>;
friend class HyperplaneHd<RT,LA>;

/*{\Mcreation 4}*/

PointHd(int d = 0) 
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname| in 
$d$-dimensional space.}*/
  : Base( Tuple(d+1) ) 
{ if ( d > 0 ) entry(d) = 1; }

PointHd(int d, const Origin&) 
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname| in 
$d$-dimensional space, initialized to the origin.}*/
  : Base( Tuple(d+1) )
{ entry(d) = 1; }

template <class InputIterator>
PointHd(int d, InputIterator first, InputIterator last) 
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname| in
dimension |d|.  If |size [first,last) == d| this creates a point with
Cartesian coordinates |set [first,last)|. If |size [first,last) ==
p+1| the range specifies the homogeneous coordinates $|H = set
[first,last)| = (\pm h_0, \pm h_1, \ldots, \pm h_d)$ where the sign
chosen is the sign of $h_d$.  \precond |d| is nonnegative,
|[first,last)| has |d| or |d+1| elements where the last has to be
non-zero, and the value type of |InputIterator| is |RT|.}*/
  : Base( Tuple(d+1,first,last) )
{ RT D = entry(d);
  if ( D == RT(0) ) entry(d) = 1;
  if ( D < RT(0) ) invert_rep();
}

template <class InputIterator>
PointHd (int d, InputIterator first, InputIterator last, 
         const RT& D) 
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname| in
dimension |d| initialized to the point with homogeneous coordinates as
defined by |H = set [first,last)| and |D|: $(\pm |H[0]|, \pm|H[1]|,
\ldots, \pm|H[d-1]|, \pm|D|)$. The sign chosen is the sign of
$D$. \precond |D| is non-zero, the iterator range defines a $d$-tuple
of |RT|, and the value type of |InputIterator| is |RT|. }*/
  : Base( Tuple(d+1,first,last,D) )
{ CGAL_assertion_msg(D!=RT(0),"PointHd::constructor: D must be nonzero.");
  if (D < RT(0)) invert_rep();
}

PointHd(int x, int y, int w = 1) : Base( Tuple((RT)x,(RT)y,(RT)w) )
{ CGAL_assertion_msg((w != 0),"PointHd::construction: w == 0.");
  if (w < 0) invert_rep();
}

PointHd(const RT& x, const RT& y, const RT& w = 1)
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname| in 
$2$-dimensional space.}*/ 
  : Base( Tuple(x,y,w,MatchHelper()) )
{ CGAL_assertion_msg((w!=0),"PointHd::construction: w == 0.");
  if (w < 0) invert_rep();
}

PointHd(int x, int y, int z, int w) : 
  Base( Tuple((RT)x,(RT)y,(RT)z,(RT)w) )
{ CGAL_assertion_msg((w!=0),"PointHd::construction: w == 0.");
  if (w < 0) invert_rep();
}

PointHd(const RT& x, const RT& y, const RT& z, const RT& w) 
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname| in 
$3$-dimensional space.}*/
  : Base( Tuple(x,y,z,w) )
{ CGAL_assertion_msg((w!=0),"PointHd::construction: w == 0.");
  if (w < 0) invert_rep();
}

PointHd(const PointHd<RT,LA>& p) : Base(p) {}
~PointHd() {}     

/*{\Moperations 4 3}*/

int dimension() const  { return ptr()->size()-1; }
/*{\Mop  returns the dimension of |\Mvar|. }*/

Quotient<RT> cartesian(int i) const
/*{\Mop returns the $i$-th Cartesian coordinate of |\Mvar|. 
   \precond $0 \leq i < d$.}*/
{ CGAL_assertion_msg((0<=i && i<dimension()),"PointHd::cartesian():\
  index out of range.");
  return Quotient<RT>(entry(i), entry(dimension())); 
}

Quotient<RT> operator[](int i) const  { return cartesian(i); }
/*{\Marrop returns the $i$-th Cartesian coordinate of |\Mvar|.
   \precond $0 \leq i < d$.}*/

RT homogeneous(int i) const 
/*{\Mop  returns the $i$-th homogeneous coordinate of |\Mvar|.
   \precond $0 \leq i \leq d$.}*/
{ CGAL_assertion_msg((0<=i && i<=(dimension())), 
    "PointHd::homogeneous():index out of range.");
  return entry(i);
}

Cartesian_const_iterator cartesian_begin() const 
/*{\Mop returns an iterator pointing to the zeroth Cartesian coordinate 
$p_0$ of |\Mvar|. }*/
{ return Cartesian_const_iterator(ptr()->begin(),ptr()->last()); }

Cartesian_const_iterator cartesian_end() const 
/*{\Mop returns an iterator pointing beyond the last Cartesian coordinate 
of |\Mvar|. }*/
{ return Cartesian_const_iterator(ptr()->last(),ptr()->last()); }

Homogeneous_const_iterator homogeneous_begin() const 
/*{\Mop returns an iterator pointing to the zeroth homogeneous coordinate 
$h_0$ of |\Mvar|. }*/
{ return ptr()->begin(); }

Homogeneous_const_iterator homogeneous_end() const 
/*{\Mop returns an iterator pointing beyond the last homogeneous coordinate 
of |\Mvar|. }*/
{ return ptr()->end(); }

PointHd<RT,LA> transform(const Aff_transformationHd<RT,LA>& t) const; 
/*{\Mop returns $t(p)$. }*/

/*{\Mtext \headerline{Arithmetic Operators, Tests and IO}}*/

inline VectorHd<RT,LA> operator-(const Origin& o) const; 
/*{\Mbinop  returns the vector $\vec{0p}$.}*/

VectorHd<RT,LA> operator-(const PointHd<RT,LA>& q) const 
/*{\Mbinop  returns $p - q$. \precond |p.dimension() == q.dimension()|.}*/
{ VectorHd<RT,LA> res(dimension()); 
  res.ptr()->homogeneous_sub(ptr(),q.ptr());
  return res; 
}

PointHd<RT,LA> operator+(const VectorHd<RT,LA>& v) const; 
/*{\Mbinop  returns $p + v$. \precond |p.dimension() == v.dimension()|.}*/

PointHd<RT,LA> operator-(const VectorHd<RT,LA>& v) const; 
/*{\Mbinop  returns $p - v$. \precond |p.dimension() == v.dimension()|.}*/

PointHd<RT,LA>& operator+=(const VectorHd<RT,LA>& v); 
/*{\Mbinop  adds |v| to |p|.\\
\precond |p.dimension() == v.dimension()|. }*/

PointHd<RT,LA>& operator-=(const VectorHd<RT,LA>& v); 
/*{\Mbinop  subtracts |v| from |p|.\\
\precond |p.dimension() == v.dimension()|. }*/

static Comparison_result cmp(
  const PointHd<RT,LA>& p1, const PointHd<RT,LA>& p2)
{ Compare_homogeneously<RT,LA> cmpobj;
  return cmpobj(p1.vector_rep(),p2.vector_rep());
}

bool operator==(const PointHd<RT,LA>& q) const
{ if (this->identical(q)) return true;
  if (dimension()!=q.dimension()) return false;
  return cmp(*this,q) == EQUAL; 
}

bool operator!=(const PointHd<RT,LA>& q) const
{ return !(*this==q); }

bool operator==(const Origin&) const
/*{\Mbinop returns true if |\Mvar| is the origin. }*/
{ for (int i = 0; i < dimension(); i++)
    if (homogeneous(i) != RT(0)) return false;
  return true;
}

friend std::istream& operator>> <>
  (std::istream&, PointHd<RT,LA>&);
friend std::ostream& operator<< <> 
  (std::ostream&, const PointHd<RT,LA>&);
/*{\Mtext \headerline{Downward compatibility}
We provide operations of the lower dimensional interface |x()|, |y()|,
|z()|, |hx()|, |hy()|, |hz()|, |hw()|.}*/

RT hx() const { return homogeneous(0); }
RT hy() const { return homogeneous(1); }
RT hz() const { return homogeneous(2); }
RT hw() const { return homogeneous(dimension()); }
Quotient<RT> x()  const { return Quotient<RT>(hx(),hw()); }
Quotient<RT> y()  const { return Quotient<RT>(hy(),hw()); }
Quotient<RT> z()  const { return Quotient<RT>(hz(),hw()); }


}; // PointHd


/*{\Mimplementation 
Points are implemented by arrays of |RT| items.  All operations like
creation, initialization, tests, point - vector arithmetic, input and
output on a point $p$ take time $O(|p.dimension()|)$. |dimension()|,
coordinate access and conversions take constant time.  The space
requirement for points is $O(|p.dimension()|)$.}*/

#undef PointHd 
} //namespace CGAL
#endif // CGAL_POINTHD_H 
//----------------------- end of file ----------------------------------
