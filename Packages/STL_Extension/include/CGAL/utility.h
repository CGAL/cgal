// ============================================================================
//
// Copyright (c) 1997, 1998, 1999, 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : utility.h
// chapter       : $CGAL_Chapter: STL Extensions for CGAL $
// package       : $CGAL_Package: STL_Extension $
// source        : stl_extension.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@cs.unc.edu>
//
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
//
// STL like utilities (triple and such)
// ============================================================================

#ifndef CGAL_UTILITY_H
#define CGAL_UTILITY_H 1

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

//+---------------------------------------------------------------------+
//| triple class                                                        |
//+---------------------------------------------------------------------+

template <class T1, class T2, class T3>
struct triple
{
  typedef T1 first_type;
  typedef T2 second_type;
  typedef T3 third_type;

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
//+---------------------------------------------------------------------+
//| quadruple class                                                     |
//+---------------------------------------------------------------------+

template <class T1, class T2, class T3, class T4>
struct quadruple
{
  typedef T1 first_type;
  typedef T2 second_type;
  typedef T3 third_type;
  typedef T4 fourth_type;

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

#endif // CGAL_UTILITY_H //
// EOF //
