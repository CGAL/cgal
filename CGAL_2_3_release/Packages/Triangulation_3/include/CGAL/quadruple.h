// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/quadruple.h
// revision      : $Revision$
// author(s)     : Andreas Fabri <Andreas.Fabri@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_QUADRUPLE_H
#define CGAL_QUADRUPLE_H

#include <CGAL/basic.h>
CGAL_BEGIN_NAMESPACE

template <class T1, class T2, class T3, class T4>
struct quadruple 
{
  T1 first;
  T2 second;
  T3 third;
  T4 fourth;

  quadruple() {}

  quadruple(const T1& a, const T2& b, const T3& c, const T4& d)
    : first(a), second(b), third(c), fourth(d) 
    {}
};

template <class T1, class T2, class T3, class T4>
inline 
quadruple<T1, T2, T3, T4>
make_quadruple(const T1& x, const T2& y, const T3& z, const T4& zz)
{
  return quadruple<T1, T2, T3, T4>(x, y, z, zz);
}

template <class T1, class T2, class T3, class T4>
inline
bool
operator==(const quadruple<T1, T2, T3, T4>& x,
           const quadruple<T1, T2, T3, T4>& y) 
{ 
  return ( (x.first == y.first) && 
	   (x.second == y.second) && 
	   (x.third == y.third) && 
	   (x.fourth == y.fourth) ); 
}

template <class T1, class T2, class T3, class T4>
inline
bool
operator<(const quadruple<T1, T2, T3, T4>& x,
	  const quadruple<T1, T2, T3, T4>& y)
{ 
  return ( x.first < y.first || 
	   ( (x.first == y.first) && (x.second < y.second) ) ||
	   ( (x.first == y.first) && (x.second == y.second) && 
	                             (x.third < y.third) ) ||
	   ( (x.first == y.first) && (x.second == y.second) && 
	                             (x.third == y.third) ) &&
	   (x.fourth < y.fourth) );
}

CGAL_END_NAMESPACE

#endif // CGAL_QUADRUPLE_H
