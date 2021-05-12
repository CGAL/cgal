// Copyright (c) 2013 GeometryFactory (France). All rights reserved.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Philipp Möller and Sebastien Loriot

#ifndef CGAL_CIRCULAR_KERNEL_2_INTERSECTION_TRAITS_H
#define CGAL_CIRCULAR_KERNEL_2_INTERSECTION_TRAITS_H

#include <CGAL/Intersection_traits.h>

#include <boost/variant.hpp>

namespace CGAL {

template <typename CK, typename T1, typename T2>
struct CK2_Intersection_traits
{};

// Intersection_traits for the circular kernel

// The additional CGAL_ADDITIONAL_VARIANT_FOR_ICL ( = int) in the variant
// has the only purpose to work around a bug of the Intel compiler,
// which without it produces the error
// /usr/include/boost/type_traits/has_nothrow_copy.hpp(36): internal error: bad pointer
// template struct has_nothrow_copy_constructor : public integral_constant{};
// See also https://github.com/CGAL/cgal/issues/1581

template<typename CK>
struct CK2_Intersection_traits<CK, typename CK::Circle_2, typename CK::Circle_2>
{
  typedef typename
  boost::variant< typename CK::Circle_2,
                  typename std::pair< typename CK::Circular_arc_point_2,
                                      unsigned int >
                  CGAL_ADDITIONAL_VARIANT_FOR_ICL
                  >
  type;
};

template<typename CK>
struct CK2_Intersection_traits<CK, typename CK::Circular_arc_2, typename CK::Circular_arc_2>
{
  typedef typename
  boost::variant< typename CK::Circular_arc_2,
                  typename std::pair< typename CK::Circular_arc_point_2,
                                      unsigned int >
                  CGAL_ADDITIONAL_VARIANT_FOR_ICL
                  >
  type;
};

template<typename CK>
struct CK2_Intersection_traits<CK, typename CK::Line_arc_2, typename CK::Line_arc_2>
{
  typedef typename
  boost::variant< typename CK::Line_arc_2,
                  typename std::pair< typename CK::Circular_arc_point_2,
                                      unsigned int >
                  CGAL_ADDITIONAL_VARIANT_FOR_ICL
                  >
  type;
};

template<typename CK>
struct CK2_Intersection_traits<CK, typename CK::Line_arc_2, typename CK::Circle_2>
{
  typedef typename
  boost::variant< typename std::pair< typename CK::Circular_arc_point_2,
                                      unsigned int >
                  CGAL_ADDITIONAL_VARIANT_FOR_ICL
                  >
  type;
};

template<typename CK>
struct CK2_Intersection_traits<CK, typename CK::Circle_2, typename CK::Line_arc_2>
  : public CK2_Intersection_traits<CK, typename CK::Line_arc_2, typename CK::Circle_2>
{};

template<typename CK>
struct CK2_Intersection_traits<CK, typename CK::Line_arc_2, typename CK::Circular_arc_2>
{
  typedef typename
  boost::variant< typename std::pair< typename CK::Circular_arc_point_2,
                                      unsigned int >
                  CGAL_ADDITIONAL_VARIANT_FOR_ICL
                  >
  type;
};

template<typename CK>
struct CK2_Intersection_traits<CK, typename CK::Circular_arc_2, typename CK::Line_arc_2>
  : public CK2_Intersection_traits<CK, typename CK::Line_arc_2, typename CK::Circular_arc_2>
{};

template<typename CK>
struct CK2_Intersection_traits<CK, typename CK::Line_arc_2, typename CK::Line_2>
{
  typedef typename
  boost::variant< typename CK::Line_arc_2,
                  typename std::pair< typename CK::Circular_arc_point_2,
                                      unsigned int >
                  CGAL_ADDITIONAL_VARIANT_FOR_ICL
                  >
  type;
};

template<typename CK>
struct CK2_Intersection_traits<CK, typename CK::Line_2, typename CK::Line_arc_2>
  : public CK2_Intersection_traits<CK, typename CK::Line_arc_2, typename CK::Line_2>
{};

template<typename CK>
struct CK2_Intersection_traits<CK, typename CK::Line_2, typename CK::Circular_arc_2>
{
  typedef typename
  boost::variant< typename std::pair< typename CK::Circular_arc_point_2,
                                      unsigned int >
                  CGAL_ADDITIONAL_VARIANT_FOR_ICL
                  >
  type;
};

template<typename CK>
struct CK2_Intersection_traits<CK, typename CK::Circular_arc_2, typename CK::Line_2>
  : public CK2_Intersection_traits<CK, typename CK::Line_2, typename CK::Circular_arc_2>
{};

// Circular_arc_2 Circle_2 simply aliases
template<typename CK>
struct CK2_Intersection_traits<CK, typename CK::Circular_arc_2, typename CK::Circle_2>
  : public CK2_Intersection_traits<CK, typename CK::Circular_arc_2, typename CK::Circular_arc_2>
{};

template<typename CK>
struct CK2_Intersection_traits<CK, typename CK::Circle_2, typename CK::Circular_arc_2>
  : public CK2_Intersection_traits<CK, typename CK::Circular_arc_2, typename CK::Circular_arc_2>
{};

// !!! undocumented !!! //

// Line_2 Circle_2
template<typename CK>
struct CK2_Intersection_traits<CK, typename CK::Line_2, typename CK::Circle_2>
{
  typedef typename
  boost::variant< typename std::pair< typename CK::Circular_arc_point_2,
                                      unsigned int >
                  CGAL_ADDITIONAL_VARIANT_FOR_ICL
                  >
  type;
};

template<typename CK>
struct CK2_Intersection_traits<CK, typename CK::Circle_2, typename CK::Line_2> :
    public CK2_Intersection_traits<CK, typename CK::Line_2, typename CK::Circle_2>
{};

template<typename CK>
struct CK2_Intersection_traits<CK, typename CK::Line_2, typename CK::Line_2>
{
  typedef typename Intersection_traits<CK, typename CK::Line_2, typename CK::Line_2>::result_type type;
};

} //end of namespace CGAL

namespace CGAL{
namespace internal{

// this function is used to call either make_object or a
// CK2_Intersection_traits::result_type constructor to create return
// values. The Object version takes some dummy template arguments
// that are needed for the return of the Intersection_traits. In
// theory a one parameter variant could be returned, but this
// _could_ come with conversion overhead and so we rather go for
// the real type.
// Overloads for empty returns are also provided.
  template<typename RT, typename T>
  inline RT
  ck2_intersection_return(T&& t) { return RT(std::forward<T>(t)); }
  template<typename RT>
  inline RT
  ck2_intersection_return() { return RT(); }

} } //end of namespace CGAL::internal


#endif // CGAL_CIRCULAR_KERNEL_2_INTERSECTION_TRAITS_H
