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
// file          : include/CGAL/Kernel_d/PointCd.h
// package       : Kernel_d
// maintainer    : Michael Seel <Michael.Seel@mpi-sb.mpg.de>
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Seel
// coordinator   : MPI Saarbruecken (Susan.Hert@mpi-sb.mpg.de)
//
// ======================================================================
#ifndef CGAL_POINTCDXXX_H
#define CGAL_POINTCDXXX_H 

#ifndef NOCGALINCL
#include <CGAL/basic.h>
#include <CGAL/Origin.h>
#endif
#include <CGAL/Kernel_d/Tuple_d.h>

CGAL_BEGIN_NAMESPACE
#define PointCd PointCd2

template <class FT, class LA> class PointCd;
template <class FT, class LA>
std::istream& operator>>(std::istream&, PointCd<FT,LA>&);
template <class FT, class LA>
std::ostream& operator<<(std::ostream&, const PointCd<FT,LA>&);

template <class _FT, class _LA > 
class PointCd : public Handle_for< Tuple_d<_FT,_LA> > { 
  typedef Tuple_d<_FT,_LA> Tuple;
  typedef Handle_for<Tuple> Base;
  typedef PointCd<_FT,_LA> Self;

typename _LA::Vector& vector_rep() { return ptr->v; }
const typename _LA::Vector& vector_rep() const { return ptr->v; }
_FT& entry(int i) const { return ptr->v[i]; }
PointCd(const Base& b) : Base(b) {}

public: 
/*{\Mtypes 4}*/

typedef _FT RT;
typedef _FT FT;
typedef _LA LA;
typedef typename Tuple::const_iterator Cartesian_const_iterator;
typedef typename Tuple::Homogeneous_const_iterator Homogeneous_const_iterator;

friend class VectorCd<FT,LA>;
friend class HyperplaneCd<FT,LA>;

/*{\Mcreation 4}*/

PointCd(int d = 0) : Base( Tuple(d) ) {}
PointCd(int d, const Origin&) : Base( Tuple(d) ) {}

#ifndef CGAL_SIMPLE_INTERFACE

template <class InputIterator>
PointCd(int d, InputIterator first, InputIterator last) 
  : Base( Tuple(d,first,last) ) 
{ if ( first == last ) return; 
  // else first specifies common denominator:
  CGAL_assertion_msg(FT(*first)!=FT(0),
    "PointCd::constructor: denominator must be nonzero.");
  for (register int i=0; i<d; ++i) entry(i)/=FT(*first);
}

template <class InputIterator>
PointCd (int d, InputIterator first, InputIterator last, 
  const FT& D) : Base( Tuple(d,first,last) )
{ CGAL_assertion_msg(D!=FT(0),"PointCd::constructor: D must be nonzero.");
  for (register int i=0; i<d; ++i) entry(i)/=D;
}

#else
#define FIXPNTCD(I) \
PointCd(int d, I first, I last) : Base( Tuple(d,first,last) ) \
{ if ( first == last ) return; \
  CGAL_assertion_msg(FT(*first)!=FT(0), \
    "PointCd::constructor: denominator must be nonzero."); \
  for (register int i=0; i<d; ++i) entry(i)/=FT(*first); \
} \
PointCd (int d, I first, I last, const FT& D) : Base( Tuple(d,first,last) ) \
{ CGAL_assertion_msg(D!=FT(0),"PointCd::constructor: D must be nonzero."); \
  for (register int i=0; i<d; ++i) entry(i)/=D; \
}

FIXPNTCD(int*)
FIXPNTCD(const int*)
FIXPNTCD(RT*)
FIXPNTCD(const RT*)
#undef FIXPNTCD
#endif 

PointCd(int x, int y, int w = 1) : Base( Tuple((FT)x,(FT)y) )
{ CGAL_assertion_msg(w!=0,"PointCd::construction: w == 0."); 
  vector_rep()/=w; }

PointCd(const FT& x, const FT& y, const FT& w = 1) 
  : Base( Tuple(x,y) )
{ CGAL_assertion_msg(w!=FT(0),"PointCd::construction: w == 0."); 
  vector_rep()/=w; }

PointCd(int x, int y, int z, int w) : 
  Base( Tuple(FT(x),FT(y),FT(z)) )
{ CGAL_assertion_msg(w!=0,"PointCd::construction: w == 0."); 
  vector_rep()/=w; }

PointCd(const FT& x, const FT& y, const FT& z, const FT& w) 
  : Base( Tuple(x,y,z) )
{ CGAL_assertion_msg(w!=FT(0),"PointCd::construction: w == 0.");
  vector_rep()/=w; }

PointCd(const PointCd<FT,LA>& p) : Base(p) {}
~PointCd() {}     

int dimension() const  { return ptr->size(); }

FT cartesian(int i) const
{ CGAL_assertion_msg((0<=i && i<dimension()),
    "PointCd::cartesian(): index out of range.");
  return entry(i); 
}
FT operator[](int i) const { return cartesian(i); }
FT homogeneous(int i) const 
{ CGAL_assertion_msg((0<=i && i<=(dimension())), 
    "PointCd::homogeneous(): index out of range.");
  if (i!=dimension()) return entry(i); else return FT(1);
}

Cartesian_const_iterator cartesian_begin() const 
{ return ptr->begin(); }
Cartesian_const_iterator cartesian_end() const 
{ return ptr->end(); }

Homogeneous_const_iterator homogeneous_begin() const 
{ return Homogeneous_const_iterator(ptr->begin(),ptr->end()); }
Homogeneous_const_iterator homogeneous_end() const 
{ return Homogeneous_const_iterator(ptr->beyondend()); }

PointCd<FT,LA> transform(const Aff_transformationCd<FT,LA>& t) const;

inline VectorCd<FT,LA> operator-(const Origin& o) const; 

VectorCd<FT,LA> operator-(const PointCd<FT,LA>& q) const 
{ VectorCd<FT,LA> res(dimension()); 
  res.ptr->cartesian_sub(ptr,q.ptr);
  return res; 
}

PointCd<FT,LA> operator+(const VectorCd<FT,LA>& v) const;
PointCd<FT,LA> operator-(const VectorCd<FT,LA>& v) const;
PointCd<FT,LA>& operator+=(const VectorCd<FT,LA>& v); 
PointCd<FT,LA>& operator-=(const VectorCd<FT,LA>& v); 

static Comparison_result cmp(
  const PointCd<FT,LA>& p1, const PointCd<FT,LA>& p2)
{ Compare_componentwise<FT,LA> cmpobj;
  return cmpobj(p1.vector_rep(),p2.vector_rep());
}

bool operator==(const PointCd<FT,LA>& q) const
{ if (identical(q)) return true;
  if (dimension()!=q.dimension()) return false;
  return vector_rep()==q.vector_rep();
}

bool operator!=(const PointCd<FT,LA>& q) const
{ return !(*this==q); }

bool operator==(const Origin&) const
{ for (int i = 0; i < dimension(); i++)
    if (cartesian(i) != FT(0)) return false;
  return true;
}

friend std::istream& operator>> CGAL_NULL_TMPL_ARGS
  (std::istream&, PointCd<FT,LA>&);
friend std::ostream& operator<< CGAL_NULL_TMPL_ARGS
  (std::ostream&, const PointCd<FT,LA>&);

FT hx() const { return cartesian(0); }
FT hy() const { return cartesian(1); }
FT hz() const { return cartesian(2); }
FT hw() const { return FT(1); }
FT x()  const { return cartesian(0); }
FT y()  const { return cartesian(1); }
FT z()  const { return cartesian(2); }

}; // PointCd

#undef PointCd2
CGAL_END_NAMESPACE
#endif // CGAL_POINTCDXXX_H 
//----------------------- end of file ----------------------------------

