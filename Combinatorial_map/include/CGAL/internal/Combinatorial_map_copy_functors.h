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
#ifndef CGAL_COMBINATORIAL_MAP_COPY_FUNCTORS_H
#define CGAL_COMBINATORIAL_MAP_COPY_FUNCTORS_H

#include <CGAL/internal/Combinatorial_map_utility.h>
#include <CGAL/internal/Combinatorial_map_internal_functors.h>
#include <CGAL/Dimension.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Cartesian_converter.h>

/* Definition of functors used internally to copy combinatorial maps attributes
 * (we need functors as attributes are stored in tuple, thus all the access
 *  must be done at compiling time).
 */
namespace CGAL
{
template< typename Map1, typename Map2, unsigned int i,
          typename Attr1=typename Map1::template Attribute_type<i>::type,
          typename Attr2=typename Map2::template Attribute_type<i>::type >
struct Default_converter_cmap_attributes;
// ****************************************************************************
namespace internal
{
// ****************************************************************************
// Map1 is the existing map, to convert into map2.
// Case where the two i-attributes are non void.
template< typename Map1, typename Map2, unsigned int i,
          typename Info1, typename Info2 >
struct Default_converter_two_non_void_attributes_cmap
{
  static typename Map2::template Attribute_handle<i>::type
  run(Map2&, typename Map1::template Attribute_const_handle<i>::type)
  { return NULL; }
};

template< typename Map1, typename Map2, unsigned int i>
struct Default_converter_two_non_void_attributes_cmap<Map1, Map2, i, void, void>
{
  static typename Map2::template Attribute_handle<i>::type
  run(Map2&, typename Map1::template Attribute_const_handle<i>::type)
  { return NULL; }
};

template< typename Map1, typename Map2, unsigned int i, typename Info1 >
struct Default_converter_two_non_void_attributes_cmap<Map1, Map2, i, Info1, void>
{
  static typename Map2::template Attribute_handle<i>::type
  run(Map2&, typename Map1::template Attribute_const_handle<i>::type)
  { return NULL; }
};

template< typename Map1, typename Map2, unsigned int i, typename Info2 >
struct Default_converter_two_non_void_attributes_cmap<Map1, Map2, i, void, Info2>
{
  static typename Map2::template Attribute_handle<i>::type
  run(Map2&, typename Map1::template Attribute_const_handle<i>::type)
  { return NULL; }
};

template< typename Map1, typename Map2, unsigned int i, typename Info >
struct Default_converter_two_non_void_attributes_cmap<Map1, Map2, i, Info, Info>
{
  static typename Map2::template Attribute_handle<i>::type
  run(Map2& map2, typename Map1::template Attribute_const_handle<i>::type ah)
  {
    typename Map2::template Attribute_handle<i>::type
      res = map2.template create_attribute<i>();
    if ( ah!=NULL ) res->info() = ah->info();
    return res;
  }
};

template<typename Map1, typename Map2, typename Converters, unsigned int i,
         bool t=(i>=boost::tuples::length<Converters>::value)>
struct Convert_attribute_functor
{
  static typename Map2::template Attribute_handle<i>::type
  run( const Map1* cmap1, Map2* cmap2, typename Map1::Dart_const_handle dh1,
       typename Map2::Dart_handle dh2, const Converters& converters)
  {
    return
        CGAL::Default_converter_cmap_attributes<Map1, Map2, i>()
        (*cmap1, *cmap2, dh1, dh2);
  }
};

template<typename Map1, typename Map2, typename Converters, unsigned int i>
struct Convert_attribute_functor<Map1,Map2,Converters,i,false>
{
  static typename Map2::template Attribute_handle<i>::type
  run( const Map1* cmap1, Map2* cmap2, typename Map1::Dart_const_handle dh1,
       typename Map2::Dart_handle dh2, const Converters& converters)
  {
    return CGAL::cpp11::get<i>(converters) (*cmap1, *cmap2, dh1, dh2);
  }
};
// ****************************************************************************
// Functor allowing to set the value of a point if point exist, have
// same dimension. For dim>3, if type of points are the same
// (because no converter).
template<typename Point1, typename Point2,
         typename T1=typename Ambient_dimension<Point1>::type,
         typename T2=typename Ambient_dimension<Point2>::type>
struct Set_point_if_possible
{
  static void run(const Point1&, Point2&)
  {}
};

template<typename Point1, typename Point2>
struct Set_point_if_possible<Point1, Point2,
    Dimension_tag<2>, Dimension_tag<2> >
{
  static void run(const Point1& p1, Point2& p2)
  {
    p2 = Cartesian_converter<typename Kernel_traits<Point1>::Kernel,
        typename Kernel_traits<Point2>::Kernel>(p1);
  }
};

template<typename Point1>
struct Set_point_if_possible<Point1, Point1,
    Dimension_tag<2>, Dimension_tag<2> >
{
  static void run(const Point1& p1, Point1& p2)
  {
    p2 = p1;
  }
};

template<typename Point1, typename Point2>
struct Set_point_if_possible<Point1, Point2,
    Dimension_tag<3>, Dimension_tag<3> >
{
  static void run(const Point1& p1, Point2& p2)
  {
    p2 = Cartesian_converter<typename Kernel_traits<Point1>::Kernel,
        typename Kernel_traits<Point2>::Kernel>()(p1);
  }
};

template<typename Point1>
struct Set_point_if_possible<Point1, Point1,
    Dimension_tag<3>, Dimension_tag<3> >
{
  static void run(const Point1& p1, Point1& p2)
  {
    p2 = p1;
  }
};

template<typename Point1>
struct Set_point_if_possible<Point1, Point1,
    Dynamic_dimension_tag, Dynamic_dimension_tag >
{
  static void run(const Point1& p1, Point1& p2)
  {
    p2 = p1;
  }
};
// ****************************************************************************
// Set_point_if_exist if both attribute has a point
template< typename Map1, typename Map2, unsigned int i,
          bool Withpoint1, bool Withpoint2 >
struct Set_point_if_exist
{
  static void run( const Map1* cmap1,
                   Map2* cmap2,
                   typename Map1::Dart_const_handle dh1,
                   typename Map2::Dart_handle dh2 )
  {}
};

template< typename Map1, typename Map2, unsigned int i >
struct Set_point_if_exist<Map1, Map2, i, true, true>
{
  static void run( const Map1* cmap1,
                   Map2* cmap2,
                   typename Map1::Dart_const_handle dh1,
                   typename Map2::Dart_handle dh2 )
  {
    if ( dh1->template attribute<i>()==NULL ) return;

    typename Map2::template Attribute_handle<i>::type
      res = dh2->template attribute<i>();
    if ( res==NULL )
    {
      res = cmap2->template create_attribute<i>();
      cmap2->template set_attribute<i>(dh2, res);
    }
    Set_point_if_possible
        <typename Map1::template Attribute_type<i>::type::Point,
        typename Map2::template Attribute_type<i>::type::Point>::
        run(dh1->template attribute<i>()->point(),
            res->point());
  }
};
// ****************************************************************************
/// Copy enabled attributes from one cmap to other
template<typename Map1, typename Map2, typename Converters>
struct Copy_attributes_functor
{
  template<unsigned int i>
  static void run( const Map1* cmap1,
                   Map2* cmap2,
                   typename Map1::Dart_const_handle dh1,
                   typename Map2::Dart_handle dh2,
                   const Converters& converters)
  {
    if (dh2->template attribute<i>()==NULL)
    {
    typename Map2::template Attribute_handle<i>::type
          res=Convert_attribute_functor<Map1,Map2,Converters,i>::
          run(cmap1, cmap2, dh1, dh2, converters);

    if ( res!=NULL )
      cmap2->template set_attribute<i>(dh2, res);

    Set_point_if_exist<Map1, Map2, i,
        sizeof(has_point<typename Map1::template Attribute_type<i>::type>(NULL))==sizeof(char),
        sizeof(has_point<typename Map2::template Attribute_type<i>::type>(NULL))==sizeof(char)>::
        run(cmap1, cmap2, dh1, dh2);
    }
  }
};
// ****************************************************************************
} // namespace internal
// ****************************************************************************
// Map1 is the existing map, to convert into map2.
// Default converter copy only attributes if they have
// same info types.
template< typename Map1, typename Map2, unsigned int i,
          typename Attr1,typename Attr2 >
struct Default_converter_cmap_attributes
{
  typename Map2::template Attribute_handle<i>::type operator()
  (const Map1& map1, Map2& map2, typename Map1::Dart_const_handle dh1,
   typename Map2::Dart_handle dh2) const
  { return internal::Default_converter_two_non_void_attributes_cmap
      <Map1,Map2,i,typename Attr1::Info,typename Attr2::Info>::
      run(map2, dh1->template attribute<i>()); }
};

template< typename Map1, typename Map2, unsigned int i,
          typename Attr1>
struct Default_converter_cmap_attributes<Map1, Map2, i, Attr1, CGAL::Void>
{
  typename Map2::template Attribute_handle<i>::type operator()
  (const Map1& map1, Map2& map2, typename Map1::Dart_const_handle dh1,
   typename Map2::Dart_handle dh2) const
  { return NULL; }
};

template< typename Map1, typename Map2, unsigned int i,
          typename Attr2>
struct Default_converter_cmap_attributes<Map1, Map2, i, CGAL::Void, Attr2>
{
  typename Map2::template Attribute_handle<i>::type operator()
  (const Map1& map1, Map2& map2, typename Map1::Dart_const_handle dh1,
   typename Map2::Dart_handle dh2) const
  { return NULL; }
};

template< typename Map1, typename Map2, unsigned int i>
struct Default_converter_cmap_attributes<Map1, Map2, i, CGAL::Void, CGAL::Void>
{
  typename Map2::template Attribute_handle<i>::type operator()
  (const Map1& map1, Map2& map2, typename Map1::Dart_const_handle dh1,
   typename Map2::Dart_handle dh2) const
  { return NULL; }
};
// ****************************************************************************
// Map1 is the existing map, to convert into map2.
// Cast converter copy always copy attributes, doing
// a cast. This works only if both types are convertible
// (and this is user responsability to use it only in
//  this case).
template< typename Map1, typename Map2, unsigned int i>
struct Cast_converter_cmap_attributes
{
  typename Map2::template Attribute_handle<i>::type operator()
  (const Map1& map1, Map2& map2, typename Map1::Dart_const_handle dh1,
   typename Map2::Dart_handle dh2) const
  {
    typename Map2::template Attribute_handle<i>::type
      res = map2.template create_attribute<i>();
    if ( dh1->template attribute<i>()!=NULL )
      res->info() = dh1->template attribute<i>()->info();
    return res;
  }
};
// ****************************************************************************
} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_COPY_FUNCTORS_H
