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

#ifndef CGAL_HYPERPLANEHD_H
#define CGAL_HYPERPLANEHD_H

#include <CGAL/basic.h>
#include <CGAL/Quotient.h>
#include <CGAL/Kernel_d/PointHd.h> 
#include <CGAL/Kernel_d/VectorHd.h> 
#include <CGAL/Kernel_d/Aff_transformationHd.h>

namespace CGAL {
#define PointHd PointHd2

template <class RT, class LA>
std::istream& operator>>(std::istream&, HyperplaneHd<RT,LA>&);
template <class RT, class LA>
std::ostream& operator<<(std::ostream&, const HyperplaneHd<RT,LA>&); 

/*{\Manpage{Hyperplane_d}{R}{Hyperplanes in d-space}{h}}*/
/*{\Msubst 
Hd<RT,LA>#_d<R>
HyperplaneHd#Hyperplane_d
Quotient<RT>#FT
}*/

template <class _RT, class _LA>
class HyperplaneHd : public Handle_for< Tuple_d<_RT,_LA> > { 
  typedef Tuple_d<_RT,_LA> Tuple;
  typedef Handle_for<Tuple> Base;
  typedef HyperplaneHd<_RT,_LA> Self;

  using Base::ptr;

/*{\Mdefinition An instance of data type |HyperplaneHd| is an
oriented hyperplane in $d$ - dimensional space. A hyperplane $h$ is
represented by coefficients $(c_0,c_1,\ldots,c_d)$ of type |RT|. At
least one of $c_0$ to $c_{ d - 1 }$ must be non-zero.  The plane
equation is $\sum_{ 0 \le i < d } c_i x_i + c_d = 0$, where $x_0$ to
$x_{d-1}$ are Cartesian point coordinates.  
For a particular $x$ the sign of $\sum_{ 0 \le i < d } c_i x_i +
c_d$ determines the position of a point $x$ with respect to the
hyperplane (on the hyperplane, on the negative side, or on the
positive side).

There are two equality predicates for hyperplanes. The (weak)
equality predicate (|weak_equality|) declares two hyperplanes equal if
they consist of the same set of points, the strong equality predicate
(|operator==|) requires in addition that the negative halfspaces
agree. In other words, two hyperplanes are strongly equal if their
coefficient vectors are positive multiples of each other and they are
(weakly) equal if their coefficient vectors are multiples of each
other.}*/

const typename _LA::Vector& vector_rep() const { return ptr()->v; }
_RT& entry(int i) { return ptr()->v[i]; }
const _RT& entry(int i) const { return ptr()->v[i]; }
void invert_rep() { ptr()->invert(); }

public: 
/*{\Mtypes 4}*/

typedef _RT RT;
/*{\Mtypemember the ring type.}*/
typedef Quotient<_RT> FT;
/*{\Mtypemember the field type.}*/
typedef _LA LA;
/*{\Mtypemember the linear algebra layer.}*/
typedef typename Tuple::const_iterator Coefficient_const_iterator;
/*{\Mtypemember a read-only iterator for the coefficients.}*/

/*{\Mcreation h 4}*/

/*{\Moptions nextwarning=no}*/
HyperplaneHd(int d = 0) : Base( Tuple(d+1) ) {}
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname|
initialized to some hyperplane in $d$ - dimensional space. }*/

template <class InputIterator>
HyperplaneHd(int d, InputIterator first, InputIterator last) 
  : Base( Tuple(d+1,first,last) ) {}
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname|
initialized to the hyperplane with coefficients |set [first,last)|.
\precond |size [first,last) == d+1| and the value type of
InputIterator is |RT|.}*/

template <class InputIterator>
HyperplaneHd(int d, InputIterator first, InputIterator last, const RT& D)
  : Base( Tuple(d+1,first,last,D) ) {}
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname|
initialized to the hyperplane with coefficients |set [first,last)| and
|D|. \precond |size [first,last) == d| and the value type of
InputIterator is |RT|.}*/

/* We want to construct a hyperplane that passes through a set |P =
set [first,last)| of points in $d$-dimensional space and has a
specified point $o$ on a specified side.  We simply have to find a
vector $x$ such that $P^T \cdot x = 0$ for every point in $P$. This
amounts to solving a homogeneous linear system. If the system has only
a trivial solution the task at hand is unsolvable and we report an
error.  So assume that the system has a non-trivial solution. Let
vectors $s_1, \ldots, s_k$ span the solution space. if |side == ZERO|
we may take any $s_j$ as the normal vector of our hyperplane. if
$|side| \neq 0$ and the task at hand is solvable there must be a $j$
such that $o^T \cdot s_j \neq 0$.  We take $s_j$ as the normal vector
of our hyperplane and use |o| to normalize the hyperplane equation. */

template <class ForwardIterator> 
void
construct_from_points(ForwardIterator first, ForwardIterator last, 
		      const PointHd<RT,LA>& o, Oriented_side side)
{ 
  TUPLE_DIM_CHECK(first,last,hyperplane::construction);
  CGAL_assertion_msg((first->dimension()==o.dimension()), 
  "hyperplane::construction: dimensions disagree.");

  int d = first->dimension();   // we are in $d$ - dimensional space
  int m = static_cast<int>(std::distance(first,last)); // |P| has $m$ points
  typename LA::Matrix A(m,d + 1); 

  for (int i = 0; i < m; i++) {  /* define $i$-th equation */
    for (int j = 0; j <= d; j++) 
      A(i,j) = first->homogeneous(j); // $j$ - th coord of $i$-th point
    ++first;
  }
  typename LA::Matrix  spanning_vecs; // columns span solution
  int dim = LA::homogeneous_linear_solver(A,spanning_vecs); 

  if (dim == 0)
    CGAL_error_msg("HyperplaneHd::constructor: \
    set P is full dimensional."); 

  if (side == ON_ORIENTED_BOUNDARY) { 
    ptr()->v = spanning_vecs.column(0); 
    return; 
  }

  RT sum = 0; 
  int j;
  for (j = 0; j < dim; j++) { 
    for (int i = 0; i <= d; i++) 
      sum += spanning_vecs(i,j)*o.homogeneous(i); 
    if (sum != 0) break; 
  }

  if (j == dim)  
    CGAL_error_msg("HyperplaneHd::constructor: \
    cannot use o to determine side.");

  ptr()->v = spanning_vecs.column(j);
  if ( ( CGAL_NTS sign(sum) > 0 && side == ON_NEGATIVE_SIDE ) || 
       ( CGAL_NTS sign(sum) < 0 && side == ON_POSITIVE_SIDE ) )   
    invert_rep();
}


template <class ForwardIterator>
HyperplaneHd(ForwardIterator first, ForwardIterator last, 
             const PointHd<RT,LA>& o, 
             Oriented_side side = ON_ORIENTED_BOUNDARY) 
/*{\Mcreate constructs some hyperplane that passes through the points
in |set [first,last)|. If |side| is |ON_POSITIVE_SIDE| or
|ON_NEGATIVE_SIDE| then |o| is on that side of the constructed
hyperplane.  \precond A hyperplane with the stated properties must
exist.  The value type of |ForwardIterator| is |PointHd<RT,LA>|. }*/
  : Base( Tuple(o.dimension()+1) )
{ construct_from_points(first,last,o,side); }

HyperplaneHd(const PointHd<RT,LA>& p, const DirectionHd<RT,LA>& dir) 
/*{\Mcreate constructs the hyperplane with normal direction |dir|
that passes through $p$. The direction |dir| points into the positive
side.  \precond |dir| is not the trivial direction.}*/
  : Base( Tuple(p.dimension()+1) ) { 
  int d = p.dimension(); 
  CGAL_assertion_msg((dir.dimension() == d), "HyperplaneHd::constructor: \
  parameter dimensions disagree.");
  CGAL_assertion_msg((dir.dimension() == d), "HyperplaneHd::constructor: \
  parameter dimensions disagree.");

  RT sum = 0; 
  for (int i = 0; i < d; i++) { 
    sum += dir.delta(i)*p.homogeneous(i); 
    entry(i) = dir.delta(i)*p.homogeneous(d); 
  }
  entry(d) = -sum; 
}

HyperplaneHd(const RT& a, const RT& b, const RT& c) : 
  Base( Tuple(a,b,c,MatchHelper()) ) {} 
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname| in 
$2$-dimensional space with equation $ax+by+c=0$. }*/

HyperplaneHd(int a, int b, int c) : 
  Base( Tuple(RT(a),RT(b),RT(c),MatchHelper()) ) {} 

HyperplaneHd(const RT& a, const RT& b, const RT& c, const RT& d) :
  Base( Tuple(a,b,c,d) ) {} 
/*{\Mcreate introduces a variable |\Mvar| of type |\Mname| in 
$3$-dimensional space with equation $ax+by+cz+d=0$. }*/

HyperplaneHd(int a, int b, int c, int d) : 
  Base( Tuple(RT(a),RT(b),RT(c),RT(d)) ) {} 

HyperplaneHd(const HyperplaneHd<RT,LA>& h) : Base(h) {}
~HyperplaneHd()  {}    

/*{\Moperations 4 2}*/

int dimension() const { return ptr()->size()-1; }
/*{\Mop returns the dimension of |\Mvar|. }*/

RT operator[](int i) const
/*{\Marrop returns the $i$-th coefficient of |\Mvar|. 
   \precond $0 \leq i \leq d$.}*/
{ CGAL_assertion_msg((0<=i && i<=(dimension())), "HyperplaneHd::op[]:\
  index out of range."); return entry(i); }

RT coefficient(int i) const { return entry(i); }
/*{\Mop returns the $i$-th coefficient of |\Mvar|. 
   \precond $0 \leq i \leq d$.}*/

const typename LA::Vector& coefficient_vector() const
/*{\Xop returns the coefficient vector $(c_0,\ldots,c_d)$ of |\Mvar|. }*/
{ return vector_rep(); }

Coefficient_const_iterator coefficients_begin() const 
/*{\Mop returns an iterator pointing to the first coefficient.}*/
{ return ptr()->begin(); }

Coefficient_const_iterator coefficients_end() const 
/*{\Mop returns an iterator pointing beyond the last coefficient.}*/
{ return ptr()->end(); }

VectorHd<RT,LA> orthogonal_vector() const; 
/*{\Mop returns the orthogonal vector of |\Mvar|. It points from the 
negative halfspace into the positive halfspace and its 
homogeneous coordinates are $(c_0, \ldots, c_{d - 1},1)$. }*/
DirectionHd<RT,LA> orthogonal_direction() const 
/*{\Mop returns the orthogonal direction of |\Mvar|. It points from the 
negative halfspace into the positive halfspace. }*/
{ return orthogonal_vector().direction(); }

RT value_at(const PointHd<RT,LA>& p) const
/*{\Xop returns the value of |\Mvar| at the point |p|, i.e., 
$\sum_{ 0 \le i \le d } h_i p_i$.\\
Warning: this value depends on the particular representation 
of |\Mvar| and |p|. }*/
{ CGAL_assertion_msg((dimension()==p.dimension()),"HyperplaneHd::value_at:\
  dimensions disagree.");
  return vector_rep()*p.vector_rep();
}

Oriented_side  oriented_side(const PointHd<RT,LA>& p) const 
/*{\Mop returns the side of the hyperplane |\Mvar| containing $p$. }*/
/*{\Mtext \setopdims{2cm}{2cm}}*/
{ 
  CGAL_assertion_msg((dimension()==p.dimension()), 
  "HyperplaneHd::oriented_side: dimensions do not agree."); 
  return CGAL_NTS sign(value_at(p));
}

bool has_on(const PointHd<RT,LA>& p) const 
/*{\Mop returns true iff point |p| lies on the hyperplane |\Mvar|. }*/
{ return (oriented_side(p) == ON_ORIENTED_BOUNDARY); }

bool has_on_boundary(const PointHd<RT,LA>& p) const 
/*{\Mop returns true iff point |p| lies on the boundary of 
hyperplane |\Mvar|. }*/
{ return (oriented_side(p) == ON_ORIENTED_BOUNDARY); }

bool has_on_positive_side(const PointHd<RT,LA>& p) const 
/*{\Mop returns true iff point |p| lies on the positive side of 
hyperplane |\Mvar|. }*/
{ return (oriented_side(p) == ON_POSITIVE_SIDE); }

bool has_on_negative_side(const PointHd<RT,LA>& p) const 
/*{\Mop returns true iff point |p| lies on the negative side of 
hyperplane |\Mvar|. }*/
{ return (oriented_side(p) == ON_NEGATIVE_SIDE); }

/*{\Mtext \restoreopdims }*/

HyperplaneHd<RT,LA> transform(const Aff_transformationHd<RT,LA>& t) const
/*{\Mop returns $t(h)$.}*/
{ Aff_transformationHd<RT,LA> t_inv = t.inverse();
  typename LA::Vector res = LA::transpose(t_inv.matrix())*vector_rep();
  if ( t_inv.is_odd() ) res = -res;
  return HyperplaneHd<RT,LA>(dimension(),res.begin(),res.end()); }

/*{\Mtext \headerline{Non-Member Functions}}*/

static Comparison_result weak_cmp(
  const HyperplaneHd<RT,LA>&, const HyperplaneHd<RT,LA>&);

static Comparison_result strong_cmp(
  const HyperplaneHd<RT,LA>&, const HyperplaneHd<RT,LA>&);

bool operator==(const HyperplaneHd<RT,LA>& h2) const
{ if (this->identical(h2)) return true;
  if (dimension()!=h2.dimension()) return false;
  return HyperplaneHd<RT,LA>::strong_cmp(*this,h2) == EQUAL; 
}

bool operator!=(const HyperplaneHd<RT,LA>& h2) const
{ return !operator==(h2); }

friend std::istream& operator>> <> 
  (std::istream&, HyperplaneHd<RT,LA>&);
friend std::ostream& operator<< <> 
  (std::ostream&, const HyperplaneHd<RT,LA>&);

}; // end of class HyperplaneHd

template <class RT, class LA>
bool weak_equality(const HyperplaneHd<RT,LA>& h1, 
                   const HyperplaneHd<RT,LA>& h2)
/*{\Mfunc test for weak equality. }*/
{ if (h1.identical(h2)) return true;
  if (h1.dimension()!=h2.dimension()) return false;
  return HyperplaneHd<RT,LA>::weak_cmp(h1,h2) == EQUAL; 
}

/*{\Mimplementation
Hyperplanes are implemented by arrays of integers as an item type.
All operations like creation, initialization, tests, vector
arithmetic, input and output on a hyperplane $h$ take time
$O(|h.dimension()|)$. coordinate access and |dimension()| take
constant time.  The space requirement is $O(|h.dimension()|)$.  }*/

#undef PointHd
} //namespace CGAL
#endif // CGAL_HYPERPLANEHD_H

//----------------------- end of file ----------------------------------
