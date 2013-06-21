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
//
// Author(s)     : Philipp Möller and Sebastien Loriot

#ifndef CGAL_CIRCULAR_KERNEL_3_INTERSECTION_TRAITS_H
#define CGAL_CIRCULAR_KERNEL_3_INTERSECTION_TRAITS_H

//this include is needed to know the value of CGAL_INTERSECTION_VERSION
#include <CGAL/Intersection_traits.h>

#if !(CGAL_INTERSECTION_VERSION < 2)

#include <boost/variant.hpp>
#include <utility>

namespace CGAL {

template <typename SK, typename T1, typename T2, typename T3=void*>
struct SK3_Intersection_traits
{};

template <typename SK>
struct SK3_Intersection_traits<SK, typename SK::Sphere_3, typename SK::Line_3>
{ 
  typedef boost::variant< 
            std::pair< typename SK::Circular_arc_point_3, unsigned int > 
          > type; 
};

template <typename SK>
struct SK3_Intersection_traits<SK, typename SK::Line_3, typename SK::Sphere_3>
  : SK3_Intersection_traits<SK, typename SK::Sphere_3, typename SK::Line_3> {};

template <typename SK>
struct SK3_Intersection_traits<SK, typename SK::Circle_3, typename SK::Plane_3>
{
  typedef boost::variant< 
            std::pair< typename SK::Circular_arc_point_3, unsigned int >,
            typename SK::Circle_3 
          > type; };

template <typename SK>
struct SK3_Intersection_traits<SK, typename SK::Plane_3, typename SK::Circle_3>
  : SK3_Intersection_traits<SK, typename SK::Circle_3, typename SK::Plane_3> {};

template <typename SK>
struct SK3_Intersection_traits<SK, typename SK::Circle_3, typename SK::Sphere_3>
{ 
  typedef boost::variant< 
            std::pair< typename SK::Circular_arc_point_3, unsigned int >, 
            typename SK::Circle_3 
          > type;
};

template <typename SK>
struct SK3_Intersection_traits<SK, typename SK::Sphere_3, typename SK::Circle_3>
  : SK3_Intersection_traits<SK, typename SK::Circle_3, typename SK::Sphere_3> {};

template <typename SK>
struct SK3_Intersection_traits<SK, typename SK::Circle_3, typename SK::Circle_3>
{
  typedef boost::variant< 
            std::pair <typename SK::Circular_arc_point_3, unsigned int >,
            typename SK::Circle_3 
          > type; 
};

template <typename SK>
struct SK3_Intersection_traits<SK, typename SK::Circle_3, typename SK::Line_3>
{
  typedef boost::variant<
            std::pair <typename SK::Circular_arc_point_3, unsigned int > 
          > type; 
};

template <typename SK>
struct SK3_Intersection_traits<SK, typename SK::Line_3, typename SK::Circle_3>
  : SK3_Intersection_traits<SK, typename SK::Circle_3, typename SK::Line_3> {};

template <typename SK>
struct SK3_Intersection_traits<SK, typename SK::Circular_arc_3, typename SK::Circular_arc_3>
{ 
  typedef boost::variant< 
            typename SK::Circle_3, 
            std::pair <typename SK::Circular_arc_point_3, unsigned int >,
            typename SK::Circular_arc_3 
          > type;
};

template <typename SK>
struct SK3_Intersection_traits<SK, typename SK::Circular_arc_3, typename SK::Plane_3>
{
  typedef boost::variant<
            std::pair <typename SK::Circular_arc_point_3, unsigned int >,
            typename SK::Circular_arc_3 
          > type;
};

template <typename SK>
struct SK3_Intersection_traits<SK, typename SK::Plane_3, typename SK::Circular_arc_3>
  : SK3_Intersection_traits<SK, typename SK::Circular_arc_3, typename SK::Plane_3> {};

template <typename SK>
struct SK3_Intersection_traits<SK, typename SK::Line_arc_3, typename SK::Line_arc_3>
{ 
  typedef boost::variant< 
            std::pair <typename SK::Circular_arc_point_3, unsigned int >,
            typename SK::Line_arc_3 
          > type; 
};
  
//struct to factorize the following specializations
template <typename SK>
struct SK3_intersect_ternary 
{
  typedef boost::variant< 
            typename SK::Circle_3,
            typename SK::Plane_3,
            typename SK::Sphere_3,
            std::pair< typename SK::Circular_arc_point_3, unsigned >
          > type;
};

template <typename SK>
struct SK3_Intersection_traits<SK, typename SK::Sphere_3, typename SK::Sphere_3, typename SK::Sphere_3> 
  : SK3_intersect_ternary<SK> {};

template <typename SK>
struct SK3_Intersection_traits<SK, typename SK::Sphere_3, typename SK::Sphere_3, typename SK::Plane_3>
  : SK3_intersect_ternary<SK> {};

template <typename SK>
struct SK3_Intersection_traits<SK, typename SK::Plane_3, typename SK::Sphere_3, typename SK::Sphere_3>
  : SK3_intersect_ternary<SK> {};

template <typename SK>
struct SK3_Intersection_traits<SK, typename SK::Plane_3, typename SK::Plane_3, typename SK::Sphere_3>
  : SK3_intersect_ternary<SK> {};

template <typename SK>
struct SK3_Intersection_traits<SK, typename SK::Sphere_3, typename SK::Plane_3, typename SK::Plane_3>
  : SK3_intersect_ternary<SK> {};

} //end of namespace CGAL

#else

#include <CGAL/Object.h>

template <typename CK, typename T1, typename T2, typename T3=void*>
struct SK3_Intersection_traits
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
    CGAL::Object sk3_intersection_return(const T& t) { return CGAL::make_object(t); }
  #else
    template<typename, typename T>
    inline
    CGAL::Object sk3_intersection_return(T&& t) { return CGAL::make_object(std::forward<T>(t)); }
  #endif // CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE
  template<typename>
  inline
  CGAL::Object sk3_intersection_return() { return CGAL::Object(); }
#else
  #if defined(CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE)
    template<typename RT, typename T>
    inline RT
    sk3_intersection_return(const T& t) { return RT(t); }
  #else
    template<typename RT, typename T>
    inline RT
    sk3_intersection_return(T&& t) { return RT(std::forward<T>(t)); }
  #endif // CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE
  template<typename RT>
  inline RT
  sk3_intersection_return() { return RT(); }
#endif // CGAL_INTERSECTION_VERSION < 2

} } //end of namespace CGAL::internal


#endif // CGAL_CIRCULAR_KERNEL_2_INTERSECTION_TRAITS_H
