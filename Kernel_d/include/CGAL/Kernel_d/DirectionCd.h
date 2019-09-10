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

#ifndef CGAL_DIRECTIONCD_H
#define CGAL_DIRECTIONCD_H

#include <CGAL/basic.h>
#include <CGAL/Kernel_d/Tuple_d.h> 

namespace CGAL {

template <class FT, class LA> class DirectionCd;
template <class FT, class LA>
std::istream& operator>>(std::istream&, DirectionCd<FT,LA>&);
template <class FT, class LA>
std::ostream& operator<<(std::ostream&, const DirectionCd<FT,LA>&);

template <class _FT, class _LA>
class DirectionCd : public Handle_for< Tuple_d<_FT,_LA> > { 
  typedef Tuple_d<_FT,_LA> Tuple;
  typedef Handle_for<Tuple> Base;
  typedef DirectionCd<_FT,_LA> Self;

  using Base::ptr;

const typename _LA::Vector& vector_rep() const { return ptr()->v; }
_FT& entry(int i) { return ptr()->v[i]; }
const _FT& entry(int i) const { return ptr()->v[i]; }

public: 
/*{\Mtypes 4}*/

typedef _FT RT;
typedef _FT FT;
typedef _LA LA;
typedef typename Tuple::const_iterator Delta_const_iterator;

class Base_direction {};

friend class VectorCd<FT,LA>; 

DirectionCd(int d = 0) : Base( Tuple(d) ) {}

DirectionCd(const VectorCd<FT,LA>& v);

template <class InputIterator>
DirectionCd(int d, InputIterator first, InputIterator last) : 
  Base( Tuple(d,first,last) ) {}

DirectionCd(int d, Base_direction, int i) : Base( Tuple(d) )
{ if ( d==0 ) return;
  CGAL_assertion_msg((0<=i&&i<d), "DirectionCd::base: index out of range.");
  entry(i) = 1;
}

DirectionCd(const FT& x, const FT y) : Base( Tuple(x,y) ) {}
DirectionCd(int a, int b) : Base( Tuple(FT(a),FT(b)) ) {}
DirectionCd(const FT& x, const FT& y, const FT& z) : 
  Base( Tuple(x,y,z) ) {}
DirectionCd(int a, int b, int c) : Base( Tuple(FT(a),FT(b),FT(c), MatchHelper()) ) {}

DirectionCd(const DirectionCd<FT,LA>& p) : Base(p)  {}
~DirectionCd() {}     

int dimension() const { return ptr()->size(); }

FT delta(int i) const  
{ CGAL_assertion_msg((0<=i && i<(dimension())), 
    "DirectionCd::delta(): index out of range.");
  return entry(i);
}
FT D() { return FT(1); }

FT operator[](int i) const  
{ return delta(i); }

Delta_const_iterator deltas_begin() const { return ptr()->begin(); }
Delta_const_iterator deltas_end() const { return ptr()->end(); }

VectorCd<FT,LA> vector() const; 

bool is_degenerate() const 
{ return vector_rep().is_zero(); }

DirectionCd<FT,LA> transform(const Aff_transformationCd<FT,LA>& t) const; 

DirectionCd<FT,LA>  opposite() const
{ DirectionCd<FT,LA> result(*this); // creates a copied object!
  result.copy_on_write(); // creates a copied object!
  result.ptr()->invert();
  return result; 
}
DirectionCd<FT,LA> operator- () const
{ return opposite(); }

static Comparison_result cmp(
  const DirectionCd<FT,LA>& h1, const DirectionCd<FT,LA>& h2);

bool operator==(const DirectionCd<FT,LA>& w) const
{ if ( this->identical(w) ) return true;
  if ( dimension()!=w.dimension() ) return false;
  return (DirectionCd<RT,LA>::cmp(*this,w) == EQUAL); 
}
bool operator!=(const DirectionCd<FT,LA>& w) const
{ return !operator==(w); }

FT dx() const { return delta(0); }
FT dy() const { return delta(1); }
FT dz() const { return delta(2); }

friend std::istream& operator>> <> 
  (std::istream& I, DirectionCd<FT,LA>& d);
friend std::ostream& operator<< <> 
  (std::ostream& O, const DirectionCd<FT,LA>& d);

}; // end of class DirectionCd


} //namespace CGAL
#endif // CGAL_DIRECTIONCD_H
//----------------------- end of file ----------------------------------
