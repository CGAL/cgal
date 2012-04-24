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

#ifndef CGAL_VECTORCD_H
#define CGAL_VECTORCD_H 

#include <CGAL/basic.h>
#include <CGAL/Kernel_d/Tuple_d.h> 

namespace CGAL {
#define PointCd PointCd2

template <class FT, class LA> class VectorCd;
template <class FT, class LA>
std::istream& operator>>(std::istream&, VectorCd<FT,LA>&);
template <class FT, class LA>
std::ostream& operator<<(std::ostream&, const VectorCd<FT,LA>&);

template <class _FT, class _LA>
class VectorCd : public Handle_for< Tuple_d<_FT,_LA> > { 
  typedef Tuple_d<_FT,_LA>  Tuple;
  typedef Handle_for<Tuple> Base;
  typedef VectorCd<_FT,_LA> Self;

  using Base::ptr;
  using Base::copy_on_write;

typename _LA::Vector& vector_rep() { return ptr()->v; }
const typename _LA::Vector& vector_rep() const { return ptr()->v; }
_FT& entry(int i) { return ptr()->v[i]; }
const _FT& entry(int i) const { return ptr()->v[i]; }
void invert_rep() { ptr()->invert(); }
VectorCd(const Base& b) : Base(b) {}

public: 

typedef _FT RT;
typedef _FT FT;
typedef _LA LA;
typedef typename Tuple::const_iterator Cartesian_const_iterator;
typedef typename Tuple::Homogeneous_const_iterator Homogeneous_const_iterator;

class Base_vector {};

friend class PointCd<FT,LA>;
friend class DirectionCd<FT,LA>;
friend class HyperplaneCd<FT,LA>; 

VectorCd(int d = 0) : Base( Tuple(d) ) {}
VectorCd(int d, Null_vector) : Base( Tuple(d) ) {}

template <class InputIterator>
VectorCd(int d, InputIterator first, InputIterator last) 
  : Base( Tuple(d,first,last) ) 
{ if ( first == last ) return; 
  // else first specifies common denominator:
  CGAL_assertion_msg(*first!=FT(0),
    "VectorCd::constructor: denominator must be nonzero.");
  for (int i=0; i<d; ++i) entry(i)/=*first;
}

template <class InputIterator>
VectorCd(int d, InputIterator first, InputIterator last, 
         const FT& D) : Base( Tuple(d,first,last) )
{ CGAL_assertion_msg(D!=FT(0), "VectorCd::constructor: D must be nonzero.");
  for (int i=0; i<d; ++i) entry(i)/=D;
}

VectorCd(int d, Base_vector, int i) : Base( Tuple(d) )
{ if ( d == 0 ) return;
  CGAL_assertion_msg((0<=i&&i<d),"VectorCd::base: index out of range.");
  entry(i) = 1;
}

VectorCd(const FT& x, const FT& y, const FT& w = 1) 
  : Base( Tuple(x,y) ) 
{ CGAL_assertion_msg((w!= FT(0)), "VectorCd::construction: w == 0.");
  vector_rep()/=w; }

VectorCd(int x, int y, int w = 1) 
  : Base( Tuple((FT)x,(FT)y) ) 
{ CGAL_assertion_msg((w!=0), "VectorCd::construction: w == 0.");
  vector_rep()/=w; }

VectorCd(const FT& x, const FT& y, const FT& z, const FT& w) 
  : Base( Tuple(x,y,z) ) 
{ CGAL_assertion_msg((w!=FT(0)), "VectorCd::construction: w == 0.");
  vector_rep()/=w; }

VectorCd(int x, int y, int z, int w) :
  Base( Tuple((FT)x,(FT)y,(FT)z, MatchHelper()) )
{ CGAL_assertion_msg((w!=0), "VectorCd::construction: w == 0.");
  vector_rep()/=w; }

VectorCd(const VectorCd<FT,LA>& p) : Base(p)  {}
~VectorCd() {}     

int dimension() const { return ptr()->size(); } 

FT cartesian(int i) const 
{ CGAL_assertion_msg((0<=i && i<(dimension())), 
    "VectorCd::cartesian(): index out of range.");
  return entry(i); 
}
FT operator[](int i) const { return cartesian(i); }
FT homogeneous(int i) const 
{ CGAL_assertion_msg((0<=i && i<=(dimension())), 
    "VectorCd::homogeneous(): index out of range.");
  if (i!=dimension()) return entry(i); else return FT(1);
}

FT squared_length() const
{ return vector_rep()*vector_rep(); }

Cartesian_const_iterator cartesian_begin() const 
{ return ptr()->begin(); }
Cartesian_const_iterator cartesian_end() const 
{ return ptr()->end(); }

Homogeneous_const_iterator homogeneous_begin() const 
{ return Homogeneous_const_iterator(ptr()->begin(),ptr()->end()); }
Homogeneous_const_iterator homogeneous_end() const 
{ return Homogeneous_const_iterator(ptr()->beyondend()); }

inline PointCd<FT,LA> to_point() const;

inline DirectionCd<FT,LA> direction() const; 
/*{\Mop returns the direction of |\Mvar|. }*/

VectorCd<FT,LA> transform(const Aff_transformationCd<FT,LA>& t) const; 

VectorCd<FT,LA> scale(const FT& m) const
{ VectorCd<FT,LA> result(*this);
  result.copy_on_write(); 
  result.vector_rep() *= m;
  return result; 
}

void self_scale(const FT& m)
{ copy_on_write();
  vector_rep() *= m;
}

VectorCd<FT,LA>& operator*=(const FT& n) 
{ self_scale(n); return *this; }
VectorCd<FT,LA>& operator*=(int n) 
{ self_scale(n); return *this; }

VectorCd<FT,LA> operator/(int n) const
{ return scale(FT(1)/FT(n)); }
VectorCd<FT,LA> operator/(const FT& n) const
{ return scale(FT(1)/n); }

VectorCd<FT,LA>& operator/=(const FT& n)
{ self_scale(FT(1)/n); return *this; }
VectorCd<FT,LA>& operator/=(int n) 
{ self_scale(FT(1)/FT(n)); return *this; }

FT operator* (const VectorCd<FT,LA>& w) const
{ return vector_rep()*w.vector_rep(); }

VectorCd<FT,LA> operator+(const VectorCd<FT,LA>& w) const 
{ VectorCd<FT,LA> result(w.dimension()); 
  result.ptr()->cartesian_add(ptr(),w.ptr());
  return result; }

VectorCd<FT,LA> operator-(const VectorCd<FT,LA>& w) const 
{ VectorCd<FT,LA> result(w.dimension());
  result.ptr()->cartesian_sub(ptr(),w.ptr());
  return result; }

VectorCd<FT,LA> operator-() const 
{ VectorCd<FT,LA> result(*this);
  result.copy_on_write(); // creates a copied object!
  result.ptr()->invert();
  return result; 
}

VectorCd<FT,LA>& operator+=(const VectorCd<FT,LA>& w) 
{ copy_on_write(); vector_rep() += w.vector_rep(); 
  return *this; }

VectorCd<FT,LA>& operator-=(const VectorCd<FT,LA>& w) 
{ copy_on_write(); vector_rep() -= w.vector_rep(); 
  return *this; }

static Comparison_result cmp(
  const VectorCd<FT,LA>& x, const VectorCd<FT,LA>& y) 
{ Compare_componentwise<FT,LA> cmpobj;
  return cmpobj(x.vector_rep(),y.vector_rep());
}

bool operator==(const VectorCd<FT,LA>& w) const
{ if ( this->identical(w) ) return true;
  if ( dimension() != w.dimension() ) return false;
  return vector_rep()==w.vector_rep();
}
bool operator!=(const VectorCd<FT,LA>& w) const
{ return !operator==(w); }

bool is_zero() const
{ return vector_rep().is_zero(); }

FT hx() const { return cartesian(0); }
FT hy() const { return cartesian(1); }
FT hz() const { return cartesian(2); }
FT hw() const { return FT(1); }
FT x()  const { return cartesian(0); }
FT y()  const { return cartesian(1); }
FT z()  const { return cartesian(2); }

friend std::istream& operator>> <>
  (std::istream& I, VectorCd<FT,LA>& v);
friend std::ostream& operator<< <>
  (std::ostream& O, const VectorCd<FT,LA>& v);

}; // end of class VectorCd

#undef PointCd
} //namespace CGAL
#endif // CGAL_VECTORCD_H 
//----------------------- end of file ----------------------------------
