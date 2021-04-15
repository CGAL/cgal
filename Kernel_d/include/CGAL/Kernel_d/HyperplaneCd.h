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
// Author(s)     : Michael Seel

#ifndef CGAL_HYPERPLANECD_H
#define CGAL_HYPERPLANECD_H

#include <CGAL/basic.h>

namespace CGAL {
#define PointCd PointCd2

template <class FT, class LA>
std::istream& operator>>(std::istream&, HyperplaneCd<FT,LA>&);
template <class FT, class LA>
std::ostream& operator<<(std::ostream&, const HyperplaneCd<FT,LA>&);

template <class _FT, class _LA>
class HyperplaneCd : public Handle_for< Tuple_d<_FT,_LA> > {
  typedef Tuple_d<_FT,_LA> Tuple;
  typedef Handle_for<Tuple> Base;
  typedef HyperplaneCd<_FT,_LA> Self;

  using Base::ptr;

const typename _LA::Vector& vector_rep() const { return ptr()->v; }
_FT& entry(int i) { return ptr()->v[i]; }
const _FT& entry(int i) const { return ptr()->v[i]; }
void invert_rep() { ptr()->invert(); }

public:
typedef _FT RT;
typedef _FT FT;
typedef _LA LA;
typedef typename Tuple::const_iterator Coefficient_const_iterator;

HyperplaneCd(int d = 0) : Base( Tuple(d+1) ) {}

template <class InputIterator>
HyperplaneCd(int d, InputIterator first, InputIterator last, const FT& D)
  : Base( Tuple(d+1,first,last,D) ) {}

template <class InputIterator>
HyperplaneCd(int d, InputIterator first, InputIterator last)
  : Base( Tuple(d+1,first,last) ) {}

template <class ForwardIterator>
void
construct_from_points(ForwardIterator first, ForwardIterator last,
                      const PointCd<FT,LA>& o, Oriented_side side)
{
  // inline due to template parameter
  TUPLE_DIM_CHECK(first,last,hyperplane::construction);
  CGAL_assertion_msg((first->dimension()==o.dimension()),
  "hyperplane::construction: dimensions disagree.");

  int d = first->dimension(); // we are in $d$ - dimensional space
  int m = static_cast<int>(std::distance(first,last)); // |P| has $m$ points
  typename LA::Matrix A(m,d + 1);

  for (int i = 0; i < m; i++) {  /* define $i$-th equation */
    for (int j = 0; j < d; j++)
      A(i,j) = first->cartesian(j); // $j$ - th coord of $i$-th point
    A(i,d) = 1;
    ++first;
  }
  typename LA::Matrix spanning_vecs; // columns span solution
  int dim = LA::homogeneous_linear_solver(A,spanning_vecs);

  CGAL_assertion_msg(dim != 0,
   "HyperplaneCd::constructor: set P is full dimensional.");

  if (side == ON_ORIENTED_BOUNDARY)
  { ptr()->v = spanning_vecs.column(0); return; }

  FT sum = 0; int j;
  for (j = 0; j < dim; j++) {
    for (int i = 0; i < d; i++)
      sum += spanning_vecs(i,j)*o.cartesian(i);
    sum += spanning_vecs(d,j);
    if (sum != FT(0)) break;
  }

  CGAL_assertion_msg(j != dim,
    "HyperplaneCd::constructor: cannot use o to determine side.");

  ptr()->v = spanning_vecs.column(j);
  if ( ( CGAL_NTS sign(sum) > 0 && side == ON_NEGATIVE_SIDE ) ||
       ( CGAL_NTS sign(sum) < 0 && side == ON_POSITIVE_SIDE ) )
    invert_rep();
}

template <class ForwardIterator>
HyperplaneCd(ForwardIterator first, ForwardIterator last,
             const PointCd<FT,LA>& o,
             Oriented_side side = ON_ORIENTED_BOUNDARY)
  : Base( Tuple(o.dimension()+1) )
{ construct_from_points(first,last,o,side); }

HyperplaneCd(const PointCd<FT,LA>& p, const DirectionCd<FT,LA>& dir)
  : Base( Tuple(p.dimension()+1) )
{
  int d = p.dimension();
  CGAL_assertion_msg((dir.dimension() == d),
    "HyperplaneCd::constructor: parameter dimensions disagree.");

  FT sum = 0;
  for (int i = 0; i < d; i++) {
    sum += dir.delta(i)*p.cartesian(i);
    entry(i) = dir.delta(i);
  }
  entry(d) = -sum;
}

HyperplaneCd(const FT& a, const FT& b, const FT& c) :
  Base( Tuple(a,b,c,MatchHelper()) ) {}

HyperplaneCd(int a, int b, int c) :
  Base( Tuple(FT(a),FT(b),FT(c),MatchHelper()) ) {}

HyperplaneCd(const FT& a, const FT& b, const FT& c, const FT& d) :
  Base( Tuple(a,b,c,d) ) {}

HyperplaneCd(int a, int b, int c, int d) :
  Base( Tuple(FT(a),FT(b),FT(c),FT(d)) ) {}

~HyperplaneCd()  {}

int dimension() const { return ptr()->size()-1; }

FT operator[](int i) const
{ CGAL_assertion_msg((0<=i && i<=(dimension())),
  "HyperplaneCd::op[]: index out of range.");
  return entry(i); }

FT coefficient(int i) const { return entry(i); }

const typename LA::Vector& coefficient_vector() const
{ return vector_rep(); }

Coefficient_const_iterator coefficients_begin() const
{ return ptr()->begin(); }
Coefficient_const_iterator coefficients_end() const
{ return ptr()->end(); }

inline VectorCd<FT,LA> orthogonal_vector() const;

DirectionCd<FT,LA> orthogonal_direction() const
{ return orthogonal_vector().direction(); }

FT value_at(const PointCd<FT,LA>& p) const
{ CGAL_assertion_msg((dimension()==p.dimension()),
    "HyperplaneCd::value_at: dimensions disagree.");
  FT res(0);
  for (int i=0; i<dimension(); ++i)
    res += coefficient(i)*p.cartesian(i);
  res += coefficient(dimension());
  return res;
}

Oriented_side oriented_side(const PointCd<FT,LA>& p) const
{
  CGAL_assertion_msg(dimension()==p.dimension(),
    "HyperplaneCd::oriented_side: dimensions do not agree.");
  return CGAL_NTS sign(value_at(p));
}

bool has_on(const PointCd<FT,LA>& p) const
{ return (oriented_side(p) == ON_ORIENTED_BOUNDARY); }

bool has_on_boundary(const PointCd<FT,LA>& p) const
{ return (oriented_side(p) == ON_ORIENTED_BOUNDARY); }

bool has_on_positive_side(const PointCd<FT,LA>& p) const
{ return (oriented_side(p) == ON_POSITIVE_SIDE); }

bool has_on_negative_side(const PointCd<FT,LA>& p) const
{ return (oriented_side(p) == ON_NEGATIVE_SIDE); }

HyperplaneCd<FT,LA> transform(const Aff_transformationCd<FT,LA>& t) const
{ Aff_transformationCd<FT,LA> t_inv = t.inverse();
  typename LA::Vector res = LA::transpose(t_inv.matrix())*vector_rep();
  if ( t_inv.is_odd() ) res = -res;
  return HyperplaneCd<FT,LA>(dimension(),res.begin(),res.end()); }

static Comparison_result weak_cmp(
  const HyperplaneCd<FT,LA>&, const HyperplaneCd<FT,LA>&);
static Comparison_result strong_cmp(
  const HyperplaneCd<FT,LA>&, const HyperplaneCd<FT,LA>&);

bool operator==(const HyperplaneCd<FT,LA>& h2) const
{ if (this->identical(h2)) return true;
  if (dimension()!=h2.dimension()) return false;
  return HyperplaneCd<FT,LA>::strong_cmp(*this,h2) == EQUAL;
}
bool operator!=(const HyperplaneCd<FT,LA>& h2) const
{ return !operator==(h2); }

friend std::istream& operator>> <>
  (std::istream&, HyperplaneCd<FT,LA>&);
friend std::ostream& operator<< <>
  (std::ostream&, const HyperplaneCd<FT,LA>&);

}; // end of class HyperplaneCd

template <class FT, class LA>
bool weak_equality(const HyperplaneCd<FT,LA>& h1,
                   const HyperplaneCd<FT,LA>& h2)
/*{\Mfunc test for weak equality. }*/
{ if (h1.identical(h2)) return true;
  if (h1.dimension()!=h2.dimension()) return false;
  return HyperplaneCd<FT,LA>::weak_cmp(h1,h2) == EQUAL;
}

#undef PointCd
} //namespace CGAL
#endif // CGAL_HYPERPLANECD_H
//----------------------- end of file ----------------------------------
