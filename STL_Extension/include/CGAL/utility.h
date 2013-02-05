// Copyright (c) 2003  
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
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
//                 Sylvain Pion

#ifndef CGAL_UTILITY_H
#define CGAL_UTILITY_H 1

#include <CGAL/basic.h>

// The Triple and Quadruple classes are NOT RECOMMENDED anymore.
// We recommend that you use cpp11::tuple or cpp11::array instead
// for new uses.

namespace CGAL {

namespace internal {

template <int i, typename T>
struct Tuple_get;

template <typename T>
struct Tuple_get<0, T>
{
  typedef typename T::first_type result_type;
  static result_type       & get(T       & t) { return t.first; }
  static result_type const & get(T const & t) { return t.first; }
};

template <typename T>
struct Tuple_get<1, T>
{
  typedef typename T::second_type result_type;
  static result_type       & get(T       & t) { return t.second; }
  static result_type const & get(T const & t) { return t.second; }
};

template <typename T>
struct Tuple_get<2, T>
{
  typedef typename T::third_type result_type;
  static result_type       & get(T       & t) { return t.third; }
  static result_type const & get(T const & t) { return t.third; }
};

template <typename T>
struct Tuple_get<3, T>
{
  typedef typename T::fourth_type result_type;
  static result_type       & get(T       & t) { return t.fourth; }
  static result_type const & get(T const & t) { return t.fourth; }
};

}

//+---------------------------------------------------------------------+
//| Triple class                                                        |
//+---------------------------------------------------------------------+

template <class T1, class T2, class T3>
class Triple
{
  typedef Triple<T1, T2, T3> Self;

public:

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

  template < int i >
  typename internal::Tuple_get<i, Self>::result_type const &
  get() const
  {
    return internal::Tuple_get<i, Self>::get(*this);
  }

  template < int i >
  typename internal::Tuple_get<i, Self>::result_type &
  get()
  {
    return internal::Tuple_get<i, Self>::get(*this);
  }

};

template <class T1, class T2, class T3>
inline
Triple<T1, T2, T3> make_triple(const T1& x, const T2& y, const T3& z)
{
  return Triple<T1, T2, T3>(x, y, z);
}

template <class T1, class T2, class T3>
inline
Triple<T1, T2, T3> make_tuple(const T1& x, const T2& y, const T3& z)
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
inline bool operator!=(const Triple<T1, T2, T3>& x,
                       const Triple<T1, T2, T3>& y)
{
  return !(x == y);
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
class Quadruple
{
  typedef Quadruple<T1, T2, T3, T4>  Self;

public:

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

  template < int i >
  typename internal::Tuple_get<i, Self>::result_type const &
  get() const
  {
    return internal::Tuple_get<i, Self>::get(*this);
  }

  template < int i >
  typename internal::Tuple_get<i, Self>::result_type &
  get()
  {
    return internal::Tuple_get<i, Self>::get(*this);
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
Quadruple<T1, T2, T3, T4>
make_tuple(const T1& x, const T2& y, const T3& z, const T4& zz)
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
operator!=(const Quadruple<T1, T2, T3, T4>& x,
           const Quadruple<T1, T2, T3, T4>& y)
{
  return ! (x == y);
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
                   (!(y.third < x.third) && x.fourth < y.fourth)) ) ) ) );
}


} //namespace CGAL

#endif // CGAL_UTILITY_H
