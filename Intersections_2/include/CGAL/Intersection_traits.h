// Copyright (c) 2011 GeometryFactory (France). All rights reserved.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Philipp Möller

#ifndef CGAL_INTERSECTION_TRAITS_H
#define CGAL_INTERSECTION_TRAITS_H

#include <CGAL/Kernel_traits.h>
#include <CGAL/Object.h>
#include <CGAL/assertions.h>
#include <CGAL/Dimension.h>

#include <variant>

#include <type_traits>

#define CGAL_INTERSECTION_TRAITS_2(A, B, R1, R2)                \
  template<typename K>     \
  struct Intersection_traits<K, typename K::A, typename K::B>  { \
    typedef typename std::variant<typename K::R1, typename K::R2 >    \
                     variant_type;                                      \
    typedef typename std::optional< variant_type > result_type;       \
  };

#define CGAL_INTERSECTION_TRAITS_3(A, B, R1, R2, R3)            \
  template<typename K>     \
  struct Intersection_traits<K, typename K::A, typename K::B>  { \
    typedef typename std::variant<typename K::R1, typename K::R2,     \
                                    typename K::R3> variant_type;       \
    typedef typename std::optional< variant_type > result_type;       \
  };

#define CGAL_INTERSECTION_FUNCTION(A, B, DIM)                           \
  template<typename K>                                                  \
  inline                                                                \
  decltype(auto) \
  intersection(const A<K>& a, const B<K>& b) {                          \
    return BOOST_PP_CAT(K().intersect_, BOOST_PP_CAT(DIM, _object()(a, b))); \
  }                                                                     \
  template<typename K>                                                  \
  inline                                                                \
  decltype(auto) \
  intersection(const B<K>& b, const A<K>& a) {                          \
    return BOOST_PP_CAT(K().intersect_, BOOST_PP_CAT(DIM, _object()(b, a))); \
  }

#define CGAL_INTERSECTION_FUNCTION_SELF(A, DIM)                         \
  template<typename K>                                                  \
  inline                                                                \
  decltype(auto) \
  intersection(const A<K> & a, const A<K> & b) {                          \
    return BOOST_PP_CAT(K().intersect_, BOOST_PP_CAT(DIM, _object()(a, b))); \
  }

#define CGAL_DO_INTERSECT_FUNCTION(A, B, DIM)              \
  template<typename K>                                     \
  inline typename K::Boolean                               \
  do_intersect(const A<K>& a, const B<K>& b) {             \
    return BOOST_PP_CAT(K().do_intersect_, BOOST_PP_CAT(DIM, _object()(a, b))); \
  }                                                        \
  template<typename K>                                     \
  inline typename K::Boolean                               \
  do_intersect(const B<K>& b, const A<K>& a) {             \
    return BOOST_PP_CAT(K().do_intersect_, BOOST_PP_CAT(DIM, _object()(b, a))); \
  }

#define CGAL_DO_INTERSECT_FUNCTION_SELF(A, DIM)                         \
  template<typename K>                                                  \
  inline typename K::Boolean                                            \
  do_intersect(const A<K> & a, const A<K> & b) {                          \
    return BOOST_PP_CAT(K().do_intersect_, BOOST_PP_CAT(DIM, _object()(a, b))); \
  }

namespace CGAL {

// only declarationn
template<typename, typename, typename>
struct Intersection_traits {
  // This defaults to Object, if we use VERSION < 2 and do nothing
  // otherwise.
};


namespace Intersections {
namespace internal {

// this function is used to call either make_object or a
// Intersection_traits::result_type constructor to create return
// values. The Object version takes some dummy template arguments
// that are needed for the return of the Intersection_traits. In
// theory a one parameter variant could be returned, but this
// _could_ come with conversion overhead and so we rather go for
// the real type.
// Overloads for empty returns are also provided.
  template<typename F, typename A, typename B, typename T>
  decltype(auto)
  intersection_return(T&& t) { return decltype(std::declval<F>()(std::declval<A>(), std::declval<B>()))(std::forward<T>(t)); }
  template<typename F, typename A, typename B>
  decltype(auto)
  intersection_return() { return decltype(std::declval<F>()(std::declval<A>(), std::declval<B>()))(); }

// Something similar to wrap around boost::get and object_cast to
// prevent ifdefing too much. Another way could be to introduce an
// overload of boost::get for Object.  We only provide the pointer
// casts here. But use references to const as parameters. This makes
// it somewhat nicer.
template<typename T>
inline
const T* intersect_get(const CGAL::Object& o) {
  return CGAL::object_cast<T>(&o);
}

template<typename T, typename ... U>
inline
const T* intersect_get(const std::optional< std::variant<U...> >& v) {
  return std::get_if<T>(&*v);
}

template<typename T, typename ... U>
inline
const T* intersect_get(const std::variant<U...> & v) {
  return std::get_if<T>(&v);
}

template<typename A, typename B>
decltype(auto)
intersection_impl(const A& a, const B& b, CGAL::Dimension_tag<2>) {
  typedef typename CGAL::Kernel_traits<A>::Kernel Kernel;
  return Kernel().intersect_2_object()(a, b);
}

template<typename A, typename B>
decltype(auto)
intersection_impl(const A& a, const B& b, Dimension_tag<3>) {
  typedef typename CGAL::Kernel_traits<A>::Kernel Kernel;
  return Kernel().intersect_3_object()(a, b);
}

template<typename A, typename B>
typename Intersection_traits< typename CGAL::Kernel_traits<A>::Kernel, A, B>::result_type
intersection_impl(const A& a, const B& b, Dynamic_dimension_tag) {
  typedef typename CGAL::Kernel_traits<A>::Kernel Kernel;
  return Kernel().intersect_d_object()(a, b);
}

template<typename A, typename B>
inline auto // K::Boolean
do_intersect_impl(const A& a, const B& b, CGAL::Dimension_tag<2>) {
  typedef typename CGAL::Kernel_traits<A>::Kernel Kernel;
  return Kernel().do_intersect_2_object()(a, b);
}

template<typename A, typename B>
inline auto // K::Boolean
do_intersect_impl(const A& a, const B& b, Dimension_tag<3>) {
  typedef typename CGAL::Kernel_traits<A>::Kernel Kernel;
  return Kernel().do_intersect_3_object()(a, b);
}

template<typename A, typename B>
inline auto // K::Boolean
do_intersect_impl(const A& a, const B& b, Dynamic_dimension_tag) {
  typedef typename CGAL::Kernel_traits<A>::Kernel Kernel;
  return Kernel().do_intersect_d_object()(a, b);
}

} // namespace internal
} // namespace Intersections

// See overloads in the respective header files

// template<typename A, typename B>
// inline
// typename Intersection_traits< typename Kernel_traits<A>::Kernel, A, B>::result_type >::type
// intersection(const A& a, const B& b) {
//   static_assert(std::is_same<typename A::Ambient_dimension, typename B::Ambient_dimension>::value),
//                               "intersection with objects of different dimensions not supported";
//   return internal::intersection_impl(a, b, typename A::Ambient_dimension());
// }

// template<typename A, typename B>
// inline
// auto // K::Boolean
// do_intersect(const A& a, const B& b) {
//   static_assert(std::is_same<typename A::Ambient_dimension, typename B::Ambient_dimension>::value,
//                         "do_intersect with objects of different dimensions not supported");
//   return internal::do_intersect_impl(a, b, typename A::Ambient_dimension());
// }

} // CGAL

#endif /* CGAL_INTERSECTION_TRAITS_H */
