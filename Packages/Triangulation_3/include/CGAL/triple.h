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
// file          : include/CGAL/triple.h
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_TRIPLE_H
#define CGAL_TRIPLE_H

#include <CGAL/basic.h>
CGAL_BEGIN_NAMESPACE

template <class T1, class T2, class T3>
struct triple 
{
  T1 first;
  T2 second;
  T3 third;

  triple() {}

  triple(const T1& a, const T2& b, const T3& c)
    : first(a), second(b), third(c) 
    {}
};

template <class T1, class T2, class T3>
inline 
triple<T1, T2, T3> make_triple(const T1& x, const T2& y, const T3& z)
{
  return triple<T1, T2, T3>(x, y, z);
}

template <class T1, class T2, class T3>
inline bool operator==(const triple<T1, T2, T3>& x,
		       const triple<T1, T2, T3>& y) 
{ 
  return ( (x.first == y.first) && 
	   (x.second == y.second) && 
	   (x.third == y.third) ); 
}

template <class T1, class T2, class T3>
inline
bool operator<(const triple<T1, T2, T3>& x,
	       const triple<T1, T2, T3>& y)
{ 
  return ( x.first < y.first || 
	   ( (x.first == y.first) && (x.second < y.second) ) ||
	   ( (x.first == y.first) && (x.second == y.second) && 
	                             (x.third < y.third) ) );
}

CGAL_END_NAMESPACE

#endif // CGAL_TRIPLE_H
