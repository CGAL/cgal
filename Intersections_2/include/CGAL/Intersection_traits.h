// Copyright (c) 2011 GeometryFactory (France). All rights reserved.
// All rights reserved.
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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Philipp Möller

#ifndef CGAL_INTERSECTION_TRAITS_H
#define CGAL_INTERSECTION_TRAITS_H

#include <CGAL/Kernel_traits.h>
#include <CGAL/Object.h>
#include <CGAL/assertions.h>
#include <CGAL/Dimension.h>
#include <CGAL/result_of.h>

#include <boost/type_traits/is_same.hpp>
#include <boost/variant.hpp>

// The macro CGAL_INTERSECTION_VERSION controls which version of the
// intersection is used.
// Currently two values are supported:
// - 1, which means intersections with CGAL::Object
// - 2, which means intersections with Intersection_traits and the 
//      corresponding APIs in other modules
// The default value is 2.

#if !defined(CGAL_INTERSECTION_VERSION)
#define CGAL_INTERSECTION_VERSION 2
#endif

#if CGAL_INTERSECTION_VERSION < 2

#define CGAL_INTERSECTION_TRAITS_2(A, B, R1, R2)
#define CGAL_INTERSECTION_TRAITS_3(A, B, R1, R2, R3)

#else

#define CGAL_INTERSECTION_TRAITS_2(A, B, R1, R2)                \
  template<typename K>     \
  struct Intersection_traits<K, typename K::A, typename K::B>  { \
    typedef typename boost::variant<typename K::R1, typename K::R2 >    \
                     variant_type;                                      \
    typedef typename boost::optional< variant_type > result_type;       \
  };  

#define CGAL_INTERSECTION_TRAITS_3(A, B, R1, R2, R3)            \
  template<typename K>     \
  struct Intersection_traits<K, typename K::A, typename K::B>  { \
    typedef typename boost::variant<typename K::R1, typename K::R2,     \
                                    typename K::R3> variant_type;       \
    typedef typename boost::optional< variant_type > result_type;       \
  };                                                                    

#endif

#define CGAL_INTERSECTION_FUNCTION(A, B, DIM)                           \
  template<typename K>                                                  \
  inline                                                                \
  typename cpp11::result_of<BOOST_PP_CAT(typename K::Intersect_, DIM)(typename K::A, typename K::B)>::type \
  intersection(const A<K>& a, const B<K>& b) {                          \
    return BOOST_PP_CAT(K().intersect_, BOOST_PP_CAT(DIM, _object()(a, b))); \
  }                                                                     \
  template<typename K>                                                  \
  inline                                                                \
  typename cpp11::result_of<BOOST_PP_CAT(typename K::Intersect_, DIM)(typename K::A, typename K::B)>::type \
  intersection(const B<K>& a, const A<K>& b) {                          \
    return BOOST_PP_CAT(K().intersect_, BOOST_PP_CAT(DIM, _object()(a, b))); \
  }

#define CGAL_INTERSECTION_FUNCTION_SELF(A, DIM)                         \
  template<typename K>                                                  \
  inline                                                                \
  typename cpp11::result_of<BOOST_PP_CAT(typename K::Intersect_, DIM)(typename K::A, typename K::A)>::type \
  intersection(const A<K> & a, const A<K> & b) {                          \
    return BOOST_PP_CAT(K().intersect_, BOOST_PP_CAT(DIM, _object()(a, b))); \
  }

#define CGAL_DO_INTERSECT_FUNCTION(A, B, DIM)              \
  template<typename K>                                     \
  inline bool                                              \
  do_intersect(const A<K>& a, const B<K>& b) {             \
    return BOOST_PP_CAT(K().do_intersect_, BOOST_PP_CAT(DIM, _object()(a, b))); \
  }                                                        \
  template<typename K>                                     \
  inline bool                                              \
  do_intersect(const B<K>& a, const A<K>& b) {             \
    return BOOST_PP_CAT(K().do_intersect_, BOOST_PP_CAT(DIM, _object()(a, b))); \
  }

#define CGAL_DO_INTERSECT_FUNCTION_SELF(A, DIM)                         \
  template<typename K>                                                  \
  inline bool                                                           \
  do_intersect(const A<K> & a, const A<K> & b) {                          \
    return BOOST_PP_CAT(K().do_intersect_, BOOST_PP_CAT(DIM, _object()(a, b))); \
  }

namespace CGAL {

// only declarationn
template<typename, typename, typename>
struct Intersection_traits {
  // This defaults to Object, if we use VERSION < 2 and do nothing
  // otherwise.
  #if CGAL_INTERSECTION_VERSION < 2
  typedef CGAL::Object result_type;
  #endif
};


// Alias that gets the Kernel automatically and does some error checking.
// Including corresponding specialization for Bbox, as it has no Kernel.
template<typename A, typename B>
class IT : public Intersection_traits< typename Kernel_traits<A>::Kernel, A, B > {
  typedef typename Kernel_traits<A>::Kernel A_Kernel;
  typedef typename Kernel_traits<B>::Kernel B_Kernel;
  // CGAL_static_assertion_msg( (boost::is_same< A_Kernel, B_Kernel>::value),
  //                            "IT instantiated with objects from two different Kernels");
};

class Bbox_2;
class Bbox_3;

template<typename B>
class IT<Bbox_2, B> : public Intersection_traits< typename Kernel_traits<B>::Kernel, CGAL::Bbox_2, B >
{ };

template<typename B>
class IT<Bbox_3, B> : public Intersection_traits< typename Kernel_traits<B>::Kernel, CGAL::Bbox_3, B >
{ };

namespace internal {

// this function is used to call either make_object or a
// Intersection_traits::result_type constructor to create return
// values. The Object version takes some dummy template arguments
// that are needed for the return of the Intersection_traits. In
// theory a one parameter variant could be returned, but this
// _could_ come with conversion overhead and so we rather go for
// the real type.
// Overloads for empty returns are also provided.
#if CGAL_INTERSECTION_VERSION < 2
  #if defined(CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE)
    template<typename, typename, typename, typename T>
    inline
    CGAL::Object intersection_return(const T& t) { return CGAL::make_object(t); }
  #else
    template<typename, typename, typename, typename T>
    inline
    CGAL::Object intersection_return(T&& t) { return CGAL::make_object(std::forward<T>(t)); }
  #endif // CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE
  template<typename, typename, typename>
  inline
  CGAL::Object intersection_return() { return CGAL::Object(); }
#else
#if defined(CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE)
    template<typename F, typename A, typename B, typename T>
    inline typename cpp11::result_of<F(A, B)>::type
    intersection_return(const T& t) { return typename cpp11::result_of<F(A, B)>::type(t); }
  #else
    template<typename F, typename A, typename B, typename T>
    inline typename cpp11::result_of<F(A, B)>::type
    intersection_return(T&& t) { return typename cpp11::result_of<F(A, B)>::type(std::forward<T>(t)); }
  #endif // CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE
  template<typename F, typename A, typename B>
  inline typename cpp11::result_of<F(A, B)>::type
  intersection_return() { return typename cpp11::result_of<F(A, B)>::type(); }
#endif // CGAL_INTERSECTION_VERSION < 2

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

template<typename T, BOOST_VARIANT_ENUM_PARAMS(typename U)>
inline
const T* intersect_get(const boost::optional< boost::variant<BOOST_VARIANT_ENUM_PARAMS(U)> >& v) {
  return boost::get<T>(&*v);
}

template<typename T, BOOST_VARIANT_ENUM_PARAMS(typename U)>
inline
const T* intersect_get(const boost::variant<BOOST_VARIANT_ENUM_PARAMS(U)> & v) {
  return boost::get<T>(&v);
}

template<typename A, typename B>
typename cpp11::result_of<typename CGAL::Kernel_traits<A>::Kernel::Intersect_2(A, B)>::type
intersection_impl(const A& a, const B& b, CGAL::Dimension_tag<2>) {
  typedef typename CGAL::Kernel_traits<A>::Kernel Kernel;
  return Kernel().intersect_2_object()(a, b);
}

template<typename A, typename B>
typename cpp11::result_of<typename CGAL::Kernel_traits<A>::Kernel::Intersect_3(A, B)>::type
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
inline bool
do_intersect_impl(const A& a, const B& b, CGAL::Dimension_tag<2>) {
  typedef typename CGAL::Kernel_traits<A>::Kernel Kernel;
  return Kernel().do_intersect_2_object()(a, b);
}

template<typename A, typename B>
inline bool
do_intersect_impl(const A& a, const B& b, Dimension_tag<3>) {
  typedef typename CGAL::Kernel_traits<A>::Kernel Kernel;
  return Kernel().do_intersect_3_object()(a, b);
}

template<typename A, typename B>
inline bool
do_intersect_impl(const A& a, const B& b, Dynamic_dimension_tag) {
  typedef typename CGAL::Kernel_traits<A>::Kernel Kernel;
  return Kernel().do_intersect_d_object()(a, b);
}

} // internal

// See overloads in the respective header files

// template<typename A, typename B>
// inline
// typename Intersection_traits< typename Kernel_traits<A>::Kernel, A, B>::result_type >::type
// intersection(const A& a, const B& b) {
//   CGAL_static_assertion_msg( (boost::is_same<typename A::Ambient_dimension, typename B::Ambient_dimension>::value),
//                               "intersection with objects of different dimensions not supported");
//   return internal::intersection_impl(a, b, typename A::Ambient_dimension());
// }

// template<typename A, typename B>
// inline
// bool
// do_intersect(const A& a, const B& b) {
//   CGAL_static_assertion_msg((boost::is_same<typename A::Ambient_dimension, typename B::Ambient_dimension>::value), 
//                         "do_intersect with objects of different dimensions not supported");
//   return internal::do_intersect_impl(a, b, typename A::Ambient_dimension());
// }

} // CGAL

#endif /* CGAL_INTERSECTION_TRAITS_H */
