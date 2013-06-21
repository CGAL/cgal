// Copyright (c) 2013 GeometryFactory (France). All rights reserved.
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
//
// $URL$
// $Id$
//
//
// Author(s)     : Philipp Möller and Sebastien Loriot

#ifndef CGAL_CIRCULAR_KERNEL_2_INTERSECTION_TRAITS_H
#define CGAL_CIRCULAR_KERNEL_2_INTERSECTION_TRAITS_H

//this include is needed to know the value of CGAL_INTERSECTION_VERSION
#include <CGAL/Intersection_traits.h>

#if !(CGAL_INTERSECTION_VERSION < 2)

#include <boost/variant.hpp>

namespace CGAL {

template <typename CK, typename T1, typename T2>
struct CK2_Intersection_traits
{};

// Intersection_traits for the circular kernel
template<typename CK>
struct CK2_Intersection_traits<CK, typename CK::Circle_2, typename CK::Circle_2>
{
  typedef typename
  boost::variant< typename CK::Circle_2,
                  typename std::pair< typename CK::Circular_arc_point_2,
                                      unsigned int > >
  type;
};

template<typename CK>
struct CK2_Intersection_traits<CK, typename CK::Circular_arc_2, typename CK::Circular_arc_2>
{
  typedef typename
  boost::variant< typename CK::Circular_arc_2,
                  typename std::pair< typename CK::Circular_arc_point_2,
                                      unsigned int > >
  type;
};

template<typename CK>
struct CK2_Intersection_traits<CK, typename CK::Line_arc_2, typename CK::Line_arc_2>
{
  typedef typename
  boost::variant< typename CK::Line_arc_2,
                  typename std::pair< typename CK::Circular_arc_point_2,
                                      unsigned int > >
  type;
};

template<typename CK>
struct CK2_Intersection_traits<CK, typename CK::Line_arc_2, typename CK::Circle_2>
{
  typedef typename
  boost::variant< typename std::pair< typename CK::Circular_arc_point_2,
                                      unsigned int > >
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
                                      unsigned int > >
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
                                      unsigned int > >
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
                                      unsigned int > >
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
                                      unsigned int > >
  type;
};

template<typename CK>
struct CK2_Intersection_traits<CK, typename CK::Circle_2, typename CK::Line_2> :
    public CK2_Intersection_traits<CK, typename CK::Line_2, typename CK::Circle_2>
{};

} //end of namespace CGAL

#else

#include <CGAL/Object.h>

template <typename CK, typename T1, typename T2>
struct CK2_Intersection_traits
{ typedef CGAL::Object type; };

#endif

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
#if CGAL_INTERSECTION_VERSION < 2
  #if defined(CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE)
    template<typename, typename T>
    inline
    CGAL::Object ck2_intersection_return(const T& t) { return CGAL::make_object(t); }
  #else
    template<typename, typename T>
    inline
    CGAL::Object ck2_intersection_return(T&& t) { return CGAL::make_object(std::forward<T>(t)); }
  #endif // CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE
  template<typename>
  inline
  CGAL::Object ck2_intersection_return() { return CGAL::Object(); }
#else
  #if defined(CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE)
    template<typename RT, typename T>
    inline RT
    ck2_intersection_return(const T& t) { return RT(t); }
  #else
    template<typename RT, typename T>
    inline RT
    ck2_intersection_return(T&& t) { return RT(std::forward<T>(t)); }
  #endif // CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE
  template<typename RT>
  inline RT
  ck2_intersection_return() { return RT(); }
#endif // CGAL_INTERSECTION_VERSION < 2

} } //end of namespace CGAL::internal


#endif // CGAL_CIRCULAR_KERNEL_2_INTERSECTION_TRAITS_H
