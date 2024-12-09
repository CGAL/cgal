// Copyright (c) 2003
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
//                 Sylvain Pion

#ifndef CGAL_UTILITY_H
#define CGAL_UTILITY_H 1

#include <CGAL/config.h>
#include <utility>
#include <functional>
#include <type_traits>
#include <boost/functional/hash.hpp>

// The Triple and Quadruple classes are NOT RECOMMENDED anymore.
// We recommend that you use std::tuple or std::array instead
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


  friend std::size_t hash_value(Triple<T1,T2,T3> const& t)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, t.first);
        boost::hash_combine(seed, t.second);
        boost::hash_combine(seed, t.third);

        return seed;
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

struct Default_using_type
{
  template <typename Argument, typename Value>
  struct Get {
      typedef Argument type;
  };

  template <typename Value>
  struct Get<Default_using_type, Value> {
      typedef typename Value::type type;
  };
};

template <class T = Default_using_type,
          class Compare = std::less<>,
          class T1, class T2,
          class A = typename Default_using_type::Get<T,
            typename std::common_type<
              typename std::decay<T1>::type,
              typename std::decay<T2>::type > >::type,
          class P = std::pair<A, A> >
inline P make_sorted_pair(T1&& t1, T2&& t2, Compare comp = Compare())
{
  return comp(t1, t2) ? P(std::forward<T1>(t1), std::forward<T2>(t2))
                      : P(std::forward<T2>(t2), std::forward<T1>(t1));
}

template <class Pair>
auto make_sorted_pair(Pair&& pair) {
  auto&& [a, b] = std::forward<Pair>(pair);
  return make_sorted_pair(std::forward<decltype(a)>(a), std::forward<decltype(b)>(b));
}

template <class F>
class [[nodiscard]] Scope_exit {
  CGAL_NO_UNIQUE_ADDRESS F exit_function;

public:
  template <typename G>
  explicit Scope_exit(G&& g, std::enable_if_t<!std::is_same_v<std::decay_t<G>, Scope_exit>>* = nullptr)
      : exit_function(std::forward<G>(g))
  {
  }

  Scope_exit(const Scope_exit&) = delete;
  Scope_exit& operator=(const Scope_exit&) = delete;
  Scope_exit(Scope_exit&&) = delete;
  Scope_exit& operator=(Scope_exit&&) = delete;

  ~Scope_exit() {
    exit_function();
  }
};

template<typename F> Scope_exit(F) -> Scope_exit<F>;

template <typename F>
Scope_exit<F> make_scope_exit(F&& f) {
  return Scope_exit<F>(std::forward<F>(f));
}

} //namespace CGAL

namespace std {

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4099) // For VC10 it is class hash
#endif

#ifndef CGAL_CFG_NO_STD_HASH
template <class T1, class T2, class T3>
struct hash<CGAL::Triple<T1,T2,T3>>
{
  std::size_t operator()(const CGAL::Triple<T1,T2,T3>& t) const
  {
    return hash_value(t);
  }
};
#endif // CGAL_CFG_NO_STD_HASH

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

} // std namespace

#endif // CGAL_UTILITY_H
