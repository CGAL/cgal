// Copyright (c) 2003  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
//                 Sylvain Pion <Sylvain.Pion@sophia.inria.fr>

#ifndef CGAL_UTILITY_H
#define CGAL_UTILITY_H 1

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

//+---------------------------------------------------------------------+
//| Triple class                                                        |
//+---------------------------------------------------------------------+

template <class T1, class T2, class T3>
struct Triple
{
  typedef T1 first_type;
  typedef T2 second_type;
  typedef T3 third_type;

  T1 first;
  T2 second;
  T3 third;

  Triple() {}

  Triple(const T1& a, const T2& b, const T3& c)
  : first(a), second(b), third(c)
  {}

  template <class U, class V, class W>
  Triple(const U& a, const V& b, const W& c)
  : first(a), second(b), third(c)
  {}

  template <class U, class V, class W>
  Triple& operator=(const Triple<U, V, W> &t) {
    first = t.first;
    second = t.second;
    third = t.third;
    return *this;
  }
};

template <class T1, class T2, class T3>
inline
Triple<T1, T2, T3> make_triple(const T1& x, const T2& y, const T3& z)
{
  return Triple<T1, T2, T3>(x, y, z);
}

template <class T1, class T2, class T3>
inline bool operator==(const Triple<T1, T2, T3>& x,
                       const Triple<T1, T2, T3>& y)
{
  return ( (x.first == y.first) &&
           (x.second == y.second) &&
           (x.third == y.third) );
}

template <class T1, class T2, class T3>
inline
bool operator<(const Triple<T1, T2, T3>& x,
               const Triple<T1, T2, T3>& y)
{
  return ( x.first < y.first ||
           ( !(y.first < x.first) &&
             ( x.second < y.second ||
               ( !(y.second < x.second) && x.third < y.third ) ) ) );
}
//+---------------------------------------------------------------------+
//| Quadruple class                                                     |
//+---------------------------------------------------------------------+

template <class T1, class T2, class T3, class T4>
struct Quadruple
{
  typedef T1 first_type;
  typedef T2 second_type;
  typedef T3 third_type;
  typedef T4 fourth_type;

  T1 first;
  T2 second;
  T3 third;
  T4 fourth;

  Quadruple() {}

  Quadruple(const T1& a, const T2& b, const T3& c, const T4& d)
  : first(a), second(b), third(c), fourth(d)
  {}

  template <class U, class V, class W, class X>
  Quadruple(const U& a, const V& b, const W& c, const X& d)
  : first(a), second(b), third(c), fourth(d)
  {}

  template <class U, class V, class W, class X>
  Quadruple& operator=(const Quadruple<U, V, W, X> &q) {
    first = q.first;
    second = q.second;
    third = q.third;
    fourth = q.fourth;
    return *this;
  }
};

template <class T1, class T2, class T3, class T4>
inline
Quadruple<T1, T2, T3, T4>
make_quadruple(const T1& x, const T2& y, const T3& z, const T4& zz)
{
  return Quadruple<T1, T2, T3, T4>(x, y, z, zz);
}

template <class T1, class T2, class T3, class T4>
inline
bool
operator==(const Quadruple<T1, T2, T3, T4>& x,
           const Quadruple<T1, T2, T3, T4>& y)
{
  return ( (x.first == y.first) &&
           (x.second == y.second) &&
           (x.third == y.third) &&
           (x.fourth == y.fourth) );
}

template <class T1, class T2, class T3, class T4>
inline
bool
operator<(const Quadruple<T1, T2, T3, T4>& x,
          const Quadruple<T1, T2, T3, T4>& y)
{
  return ( x.first < y.first ||
           ( !(y.first < x.first) &&
             ( x.second < y.second ||
               ( !(y.second < x.second) &&
                 ( x.third < y.third ||
                   !(y.third < x.third) && x.fourth < y.fourth) ) ) ) );
}


CGAL_END_NAMESPACE

#endif // CGAL_UTILITY_H //
// EOF //
