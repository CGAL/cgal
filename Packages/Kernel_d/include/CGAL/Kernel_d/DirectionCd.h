// ======================================================================
//
// Copyright (c) 2000,2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Kernel_d/DirectionCd.h
// package       : Kernel_d
// maintainer    : Michael Seel <Michael.Seel@mpi-sb.mpg.de>
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Seel
// coordinator   : MPI Saarbruecken (Susan.Hert@mpi-sb.mpg.de)
//
// ======================================================================
#ifndef CGAL_DIRECTIONCD_H
#define CGAL_DIRECTIONCD_H

#ifndef NOCGALINCL
#include <CGAL/basic.h>
#endif
#include <CGAL/Kernel_d/Tuple_d.h> 

CGAL_BEGIN_NAMESPACE

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

const typename _LA::Vector& vector_rep() const { return ptr->v; }
_FT& entry(int i) const { return ptr->v[i]; }

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

#ifndef CGAL_SIMPLE_INTERFACE

template <class InputIterator>
DirectionCd(int d, InputIterator first, InputIterator last) : 
  Base( Tuple(d,first,last) ) {}

#else
#define FIXDIRCD(I) \
DirectionCd(int d, I first, I last) : Base( Tuple(d,first,last) ) {}
FIXDIRCD(int*)
FIXDIRCD(const int*)
FIXDIRCD(RT*)
FIXDIRCD(const RT*)
#undef FIXDIRCD
#endif

DirectionCd(int d, Base_direction, int i) : Base( Tuple(d) )
{ if ( d==0 ) return;
  CGAL_assertion_msg((0<=i&&i<d), "DirectionCd::base: index out of range.");
  entry(i) = 1;
}

DirectionCd(const FT& x, const FT& y) : Base( Tuple(x,y) ) {}
DirectionCd(int a, int b) : Base( Tuple(FT(a),FT(b)) ) {}
DirectionCd(const FT& x, const FT& y, const FT& z) : 
  Base( Tuple(x,y,z) ) {}
DirectionCd(int a, int b, int c) : Base( Tuple(FT(a),FT(b),FT(c)) ) {}

DirectionCd(const DirectionCd<FT,LA>& p) : Base(p)  {}
~DirectionCd() {}     

int dimension() const { return ptr->size(); }

FT delta(int i) const  
{ CGAL_assertion_msg((0<=i && i<(dimension())), 
    "DirectionCd::delta(): index out of range.");
  return entry(i);
}
FT D() { return FT(1); }

FT operator[](int i) const  
{ return delta(i); }

Delta_const_iterator deltas_begin() const { return ptr->begin(); }
Delta_const_iterator deltas_end() const { return ptr->end(); }

VectorCd<FT,LA> vector() const; 

bool is_degenerate() const 
{ return vector_rep().is_zero(); }

DirectionCd<FT,LA> transform(const Aff_transformationCd<FT,LA>& t) const; 

DirectionCd<FT,LA>  opposite() const
{ DirectionCd<FT,LA> result(*this); // creates a copied object!
  result.copy_on_write(); // creates a copied object!
  result.ptr->invert();
  return result; 
}
DirectionCd<FT,LA> operator- () const
{ return opposite(); }

static Comparison_result cmp(
  const DirectionCd<FT,LA>& h1, const DirectionCd<FT,LA>& h2);

bool operator==(const DirectionCd<FT,LA>& w) const
{ if ( identical(w) ) return true;
  if ( dimension()!=w.dimension() ) return false;
  return (DirectionCd<RT,LA>::cmp(*this,w) == EQUAL); 
}
bool operator!=(const DirectionCd<FT,LA>& w) const
{ return !operator==(w); }

FT dx() const { return delta(0); }
FT dy() const { return delta(1); }
FT dz() const { return delta(2); }

friend std::istream& operator>> CGAL_NULL_TMPL_ARGS 
  (std::istream& I, DirectionCd<FT,LA>& d);
friend std::ostream& operator<< CGAL_NULL_TMPL_ARGS 
  (std::ostream& O, const DirectionCd<FT,LA>& d);

}; // end of class DirectionCd


CGAL_END_NAMESPACE
#endif // CGAL_DIRECTIONCD_H
//----------------------- end of file ----------------------------------

