// Copyright (c) 2010-2013 CNRS and LIRIS' Establishments (France).
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
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_LINEAR_CELL_COMPLEX_COPY_FUNCTORS_H
#define CGAL_LINEAR_CELL_COMPLEX_COPY_FUNCTORS_H 1

#include <CGAL/Dimension.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Cartesian_converter.h>

/* Definition of functors used internally to copy linear cell complex attributes
 * (we need functors as attributes are stored in tuple, thus all the access
 *  must be done at compiling time).
 */
namespace CGAL
{
// ****************************************************************************
namespace internal
{
// ****************************************************************************
// Functor allowing to set the value of a point if point exist, have
// same dimension. For dim>3, if type of points are the same (because no converter).
template<typename Point1, typename Point2>
struct Set_point_d_if_same
{
  static void run(const Point1&, Point2&)
  {}
};

template<typename Point1>
struct Set_point_d_if_same<Point1, Point1>
{
  static void run(const Point1& p1, Point1& p2)
  {
    p2 = p1; // Copy of Point_d having same type
  }
};

template<typename Point1, typename Point2,
         typename T1=typename Ambient_dimension<Point1>::type,
         typename T2=typename Ambient_dimension<Point2>::type>
struct Set_point_if_possible
{
  static void run(const Point1&, Point2&)
  {}
};

template<typename Point1, typename Point2>
struct Set_point_if_possible<Point1, Point2, Dimension_tag<2>, Dimension_tag<2> >
{
  static void run(const Point1& p1, Point2& p2)
  {
    p2 = Cartesian_converter<typename Kernel_traits<Point1>::Kernel,
        typename Kernel_traits<Point2>::Kernel>(p1);
  }
};

template<typename Point1, typename Point2>
struct Set_point_if_possible<Point1, Point2, Dimension_tag<3>, Dimension_tag<3> >
{
  static void run(const Point1& p1, Point2& p2)
  {
    p2 = Cartesian_converter<typename Kernel_traits<Point1>::Kernel,
        typename Kernel_traits<Point2>::Kernel>()(p1);
  }
};

template<typename Point1, typename Point2>
struct Set_point_if_possible<Point1, Point2, Dynamic_dimension_tag,
                             Dynamic_dimension_tag>
{
  static void run(const Point1& p1, Point2& p2)
  {
    if ( p1.dimension()==p2.dimension() )
      Set_point_d_if_same<Point1, Point2>::run(p1, p2);
  }
};

// Set_point_if_exist if Attr1 has a point
template< typename Attr1, typename Attr2,
          typename Point1=typename Attr1::Point,
          typename Point2=typename Attr2::Point >
struct Set_point_if_exist
{
  static void run(const Attr1& a1, Attr2& a2)
  {
    Set_point_if_possible<Point1, Point2>::run(a1.point(), a2.point());
  }
};

template<typename Attr1, typename Attr2, typename Point1 >
struct Set_point_if_exist<Attr1, Attr2, Point1, CGAL::Void>
{
  static void run(const Attr1&, Attr2&)
  {}
};
// ****************************************************************************
// Map1 is the existing map, to convert into map2.
// Case where the two i-attributes are non void.
template< typename Map1, typename Map2, unsigned int i,
          typename Info1, typename Info2,
          typename Point2 >
struct Default_converter_two_non_void_attributes_lcc
{ // Here Info1!=Info2 but Point2!=CGAL::Void (thus Linear_cell_complex)
  static typename Map2::template Attribute_handle<i>::type
  run(Map2& map2, typename Map1::template Attribute_const_handle<i>::type ah)
  {
    typename Map2::template Attribute_handle<i>::type
        res=map2.template create_attribute<i>();
    if ( ah!=NULL )
    {
      // Copy the point of ah if it exists and have same dimension
      Set_point_if_exist<typename Map1::template Attribute_type<i>::type,
          typename Map2::template Attribute_type<i>::type>::run( *ah, *res );
    }
    return res;
  }
};

template< typename Map1, typename Map2, unsigned int i, typename Info, typename Point2 >
struct Default_converter_two_non_void_attributes_lcc<Map1, Map2, i, Info, Info, Point2>
{ // Here Info1==Info2 but Point2!=CGAL::Void (thus Linear_cell_complex)
  static typename Map2::template Attribute_handle<i>::type
  run(Map2& map2, typename Map1::template Attribute_const_handle<i>::type ah)
  {
    typename Map2::template Attribute_handle<i>::type
        res=map2.template create_attribute<i>();
    if ( ah!=NULL )
    {
      res->info()=ah->info();
      // Copy the point of ah if it exists and have same dimension
      Set_point_if_exist<typename Map1::template Attribute_type<i>::type,
          typename Map2::template Attribute_type<i>::type>::run( *ah, *res );
    }
    return res;
  }
};

template< typename Map1, typename Map2, unsigned int i, typename Point2 >
struct Default_converter_two_non_void_attributes_lcc<Map1, Map2, i, void, void, Point2>
{ // Here Info1==Info2==void but Point2!=CGAL::Void (thus Linear_cell_complex)
  static typename Map2::template Attribute_handle<i>::type
  run(Map2& map2, typename Map1::template Attribute_const_handle<i>::type ah)
  {
    typename Map2::template Attribute_handle<i>::type res=
      map2.template create_attribute<i>();
    if ( ah!=NULL )
    {
      // Copy the point of ah if it exists and have same dimension
      Set_point_if_exist<typename Map1::template Attribute_type<i>::type,
          typename Map2::template Attribute_type<i>::type>::
          run( *ah, *res );
    }
    return res;
  }
};

// ****************************************************************************
} // namespace internal
} // namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_COPY_FUNCTORS_H
