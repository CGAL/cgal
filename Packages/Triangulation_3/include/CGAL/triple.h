// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// coordinator   : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIPLE_H
#define CGAL_TRIPLE_H

//#include <CGAL/circulator.h>

// #ifdef __GNUG__
// #define CGAL_NULL_TYPE const void*
// #else // __GNUG__ //
// #define CGAL_NULL_TYPE int
// #endif // __GNUG__ //

template <class T1, class T2, class T3>
struct CGAL_triple {
  T1 first;
  T2 second;
  T3 third;

  CGAL_triple() {}

  CGAL_triple(const T1& a, const T2& b, const T3& c)
    : first(a), second(b), third(c) 
    {}

//   bool operator==( CGAL_NULL_TYPE n ) const
//   {
//     CGAL_triangulation_assertion( n == NULL );
//     return (first == NULL);
//   }

//   bool operator!=( CGAL_NULL_TYPE n ) const
//   {
//     return !(this == n);
//   }

};

template <class T1, class T2, class T3>
inline 
CGAL_triple<T1, T2, T3> CGAL_make_triple(const T1& x, const T2& y, const T3& z)
{
  return CGAL_triple<T1, T2, T3>(x, y, z);
}

template <class T1, class T2, class T3>
inline bool operator==(const CGAL_triple<T1, T2, T3>& x,
		       const CGAL_triple<T1, T2, T3>& y) 
{ 
  return ( (x.first == y.first) && 
	   (x.second == y.second) && 
	   (x.third == y.third) ); 
}

template <class T1, class T2, class T3>
inline
bool operator<(const CGAL_triple<T1, T2, T3>& x,
	       const CGAL_triple<T1, T2, T3>& y)
{ 
  return ( x.first < y.first || 
	   ( (x.first == y.first) && (x.second < y.second) ) ||
	   ( (x.first == y.first) && (x.second == y.second) && 
	                             (x.third < y.third) ) );
}

#endif CGAL_TRIPLE_H
