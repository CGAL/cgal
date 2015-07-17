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
#include <CGAL/Combinatorial_map_functors.h>
#include <CGAL/Dimension.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/use.h>
/* Definition of functors used internally to copy combinatorial maps attributes
 * (we need functors as attributes are stored in tuple, thus all the access
 *  must be done at compiling time).
 */
namespace CGAL
{
template< typename Map1, typename Map2, unsigned int i>
struct Default_converter_cmap_attributes;
template< typename Map1, typename Map2>
struct Default_converter_cmap_0attributes_with_point;
// ****************************************************************************
namespace internal
{
// ****************************************************************************
// Map1 is the existing map, to convert into map2.
// Functor called only when both i-attributes have non void info.
// General cases when both info are differents.
template< typename Map1, typename Map2, unsigned int i,
          typename Info1=typename Map1::template
          Attribute_type<i>::type::Info,
          typename Info2=typename Map2::template
          Attribute_type<i>::type::Info >
struct Create_attribute_if_same_info_cmap
{
  static typename Map2::template Attribute_handle<i>::type
  run(Map2&, const Info1& )
  { return Map2::null_handle; }
};

// Special case when both attributes have the same info.
template< typename Map1, typename Map2, unsigned int i, typename Info >
struct Create_attribute_if_same_info_cmap<Map1, Map2, i, Info, Info>
{
  static typename Map2::template Attribute_handle<i>::type
  run(Map2& map2, const Info& info)
  { typename Map2::template Attribute_handle<i>::type
        res = map2.template create_attribute<i>();
    map2.template get_attribute<i>(res).info()=info;
    return res;
  }
};
// ****************************************************************************
// Functor allowing to set the value of a point if point exist, have
// same dimension. For dim>3, if type of points are the same
// (because no converter).
template<typename Point1, typename Point2,
         typename T1=typename Ambient_dimension<Point1>::type,
         typename T2=typename Ambient_dimension<Point2>::type>
struct Set_point_if_possible_cmap
{
  static void run(const Point1&, Point2&)
  {}
};

template<typename Point1, typename Point2>
struct Set_point_if_possible_cmap<Point1, Point2,
    Dimension_tag<2>, Dimension_tag<2> >
{
  static void run(const Point1& p1, Point2& p2)
  {
    p2 = Cartesian_converter<typename Kernel_traits<Point1>::Kernel,
        typename Kernel_traits<Point2>::Kernel>(p1);
  }
};

template<typename Point1>
struct Set_point_if_possible_cmap<Point1, Point1,
    Dimension_tag<2>, Dimension_tag<2> >
{
  static void run(const Point1& p1, Point1& p2)
  {
    p2 = p1;
  }
};

template<typename Point1, typename Point2>
struct Set_point_if_possible_cmap<Point1, Point2,
    Dimension_tag<3>, Dimension_tag<3> >
{
  static void run(const Point1& p1, Point2& p2)
  {
    p2 = Cartesian_converter<typename Kernel_traits<Point1>::Kernel,
        typename Kernel_traits<Point2>::Kernel>()(p1);
  }
};

template<typename Point1>
struct Set_point_if_possible_cmap<Point1, Point1,
    Dimension_tag<3>, Dimension_tag<3> >
{
  static void run(const Point1& p1, Point1& p2)
  {
    p2 = p1;
  }
};

template<typename Point1>
struct Set_point_if_possible_cmap<Point1, Point1,
    Dynamic_dimension_tag, Dynamic_dimension_tag >
{
  static void run(const Point1& p1, Point1& p2)
  {
    p2 = p1;
  }
};
// ****************************************************************************
// Get the ith functor of the converters tuple if i<length of the tuple,
// otherwise get the default converter.
template<typename Map1, typename Map2, unsigned int i, typename Converters,
         bool t=((int)i>=My_length<Converters>::value)>
struct Get_convert_attribute_functor
{
  static typename Map2::template Attribute_handle<i>::type
  run( const Map1* cmap1, Map2* cmap2, typename Map1::Dart_const_handle dh1,
       typename Map2::Dart_handle dh2, const Converters& /*converters*/)
  {
    return
        CGAL::Default_converter_cmap_attributes<Map1, Map2, i>()
        (*cmap1, *cmap2, dh1, dh2);
  }
};

template<typename Map1, typename Map2, unsigned int i, typename Converters>
struct Get_convert_attribute_functor<Map1,Map2,i,Converters,false>
{
  static typename Map2::template Attribute_handle<i>::type
  run( const Map1* cmap1, Map2* cmap2, typename Map1::Dart_const_handle dh1,
       typename Map2::Dart_handle dh2, const Converters& converters)
  {
    return CGAL::cpp11::get<i>(converters) (*cmap1, *cmap2, dh1, dh2);
  }
};
// ****************************************************************************
// Call a given functor if both i-attribute have an non void info
template< typename Map1, typename Map2, unsigned int i,
          typename Converters,
          bool Withinfo1=CGAL::template
          Is_attribute_has_non_void_info
          <typename Map1::template Attribute_type<i>::type>::value,
          bool Withinfo2=CGAL::template
          Is_attribute_has_non_void_info
          <typename Map2::template Attribute_type<i>::type>::value >
struct Call_functor_if_both_attributes_have_info
{
  static typename Map2::template Attribute_handle<i>::type
  run( const Map1*,
       Map2*,
       typename Map1::Dart_const_handle,
       typename Map2::Dart_handle,
       const Converters&)
  { return Map2::null_handle; }
};

template< typename Map1, typename Map2, unsigned int i, typename Converters >
struct Call_functor_if_both_attributes_have_info<Map1, Map2, i,
    Converters, true, true>
{
  static typename Map2::template Attribute_handle<i>::type
  run( const Map1* cmap1,
       Map2* cmap2,
       typename Map1::Dart_const_handle dh1,
       typename Map2::Dart_handle dh2,
       const Converters& converters )
  {
    return Get_convert_attribute_functor<Map1,Map2,i,Converters>::
        run(cmap1, cmap2, dh1, dh2, converters);
  }
};
// ****************************************************************************
// Call a given functor only if both 0-attribute have a point.
// general case i!=0 or one attribute without point.
template< typename Map1, typename Map2, unsigned int i,
          typename Pointconverter,
          bool Withpoint1=CGAL::template Is_attribute_has_point
          <typename Map1::template Attribute_type<i>::type>::value,
          bool Withpoint2=CGAL::template Is_attribute_has_point
          <typename Map2::template Attribute_type<i>::type>::value >
struct Call_functor_if_both_attributes_have_point
{
  static typename Map2::template Attribute_handle<i>::type
  run( const Map1*,
       Map2*,
       typename Map1::Dart_const_handle,
       typename Map2::Dart_handle,
       const Pointconverter&)
  { return Map2::null_handle; }
};
// Specialisation with i==0 and both attributes have points.
template< typename Map1, typename Map2, typename Pointconverter >
struct Call_functor_if_both_attributes_have_point<Map1, Map2, 0,
    Pointconverter, true, true>
{
  static typename Map2::template Attribute_handle<0>::type
  run( const Map1* cmap1,
       Map2* cmap2,
       typename Map1::Dart_const_handle dh1,
       typename Map2::Dart_handle dh2,
       const Pointconverter& pointconverter )
  { return pointconverter(*cmap1, *cmap2, dh1, dh2); }
};
// ****************************************************************************
// Copy attribute when if both i-attributes are non void.
// (note Attr2 could not be Void as copy functor is called only for
// non void attributes)
// General case with both attributes non void.
template<typename Map1, typename Map2, typename Converters,
         typename Pointconverter, unsigned int i,
         typename Attr1=typename Map1::template Attribute_type<i>::type,
         typename Attr2=typename Map2::template Attribute_type<i>::type >
struct Copy_attribute_functor_if_nonvoid
{
  static void run( const Map1* cmap1,
                   Map2* cmap2,
                   typename Map1::Dart_const_handle dh1,
                   typename Map2::Dart_handle dh2,
                   const Converters& converters,
                   const Pointconverter& pointconverter)
  {
    // If dh1 has no i-attribute, nothing to copy.
    if ( cmap1->template attribute<i>(dh1)==Map1::null_handle ) return;

    // If dh2 has already an i-attribute, it was already copied.
    if ( cmap2->template attribute<i>(dh2)!=Map2::null_handle ) return;

    // Otherwise we copy the info if both attribute have non void info.
    typename Map2::template Attribute_handle<i>::type
        res=Call_functor_if_both_attributes_have_info
        <Map1, Map2, i, Converters>::
        run(cmap1, cmap2, dh1, dh2, converters);

    if ( res!=Map2::null_handle )
      cmap2->template set_attribute<i>(dh2, res);

    // And the point if both attributes have points (and only for 0-attributes)
    res=Call_functor_if_both_attributes_have_point
        <Map1, Map2, i, Pointconverter>::
        run(cmap1, cmap2, dh1, dh2, pointconverter);

    if ( res!=Map2::null_handle &&
         cmap2->template attribute<i>(dh2)==Map2::null_handle )
      cmap2->template set_attribute<i>(dh2, res);
  }
};
// Specialisation when attr1 is void, and attr2 is non void i==0. Nothing to
// copy, but if 0-attributes has point and i==0, we need to create
// vertex attributes.
template<typename Map1, typename Map2, typename Converters,
         typename Pointconverter, typename Attr2>
struct Copy_attribute_functor_if_nonvoid<Map1, Map2, Converters,
    Pointconverter, 0, CGAL::Void, Attr2>
{
  static void run( const Map1*,
                   Map2* cmap2,
                   typename Map1::Dart_const_handle,
                   typename Map2::Dart_handle dh2,
                   const Converters&,
                   const Pointconverter&)
  {
    // If dh2 has already an 0-attribute, it was already created.
    if ( cmap2->template attribute<0>(dh2)!=Map2::null_handle ) return;

    // Create the point if 0-attributes has Point.
    if ( CGAL::template Is_attribute_has_point
         <typename Map2::template Attribute_type<0>::type>::value )
      cmap2->template
          set_attribute<0>(dh2, cmap2->template create_attribute<0>());
  }
};
// Specialisation when attr1 is void, and attr2 is non void i!=0.
// Nothing to do.
template<typename Map1, typename Map2, typename Converters, unsigned int i,
         typename Pointconverter, typename Attr2>
struct Copy_attribute_functor_if_nonvoid<Map1, Map2, Converters,
    Pointconverter, i, CGAL::Void, Attr2>
{
  static void run( const Map1*,
                   Map2*,
                   typename Map1::Dart_const_handle,
                   typename Map2::Dart_handle,
                   const Converters&,
                   const Pointconverter&)
  {}
};
// ****************************************************************************
/// Copy enabled attributes from one cmap to other. General case called
/// by copy function in Combinatorial_map on all the non void attributes
/// of Map2. Map1 is the existing map, to convert into map2.
template<typename Map1, typename Map2, typename Converters,
         typename Pointconverter>
struct Copy_attributes_functor
{
  template<unsigned int i>
  static void run( const Map1* cmap1,
                   Map2* cmap2,
                   typename Map1::Dart_const_handle dh1,
                   typename Map2::Dart_handle dh2,
                   const Converters& converters,
                   const Pointconverter& pointconverter)
  { Copy_attribute_functor_if_nonvoid
        <Map1, Map2, Converters, Pointconverter, i>::
        run(cmap1, cmap2, dh1, dh2, converters, pointconverter);
   }
};
// ****************************************************************************
} // namespace internal
// ****************************************************************************
// "Converters" called during the copy of attributes, to copy Info.
// Users can replace them by their own converters.
// Info converter are called only if both i-attributes have non void info,
// if dh1 has an i-attribute and if dh2 does not already has an i-attribute.
// Map1 is the existing map, to convert into map2.
// ****************************************************************************
// Default converter copy only attributes if they have same info types.
template< typename Map1, typename Map2, unsigned int i>
struct Default_converter_cmap_attributes
{
  typename Map2::template Attribute_handle<i>::type operator()
  (const Map1& map1, Map2& map2, typename Map1::Dart_const_handle dh1,
   typename Map2::Dart_handle dh2) const
  {
    CGAL_USE(dh2);
    CGAL_assertion( map1.template attribute<i>(dh1)!=Map1::null_handle );
    CGAL_assertion( map2.template attribute<i>(dh2)==Map2::null_handle );
    return internal::Create_attribute_if_same_info_cmap
      <Map1,Map2,i>::run(map2, map1.template info<i>(dh1));
  }
};
// ****************************************************************************
// Cast converter always copy attributes, doing a cast. This can work only
// if both types are convertible and this is user responsability
// to use it only in this case.
template< typename Map1, typename Map2, unsigned int i>
struct Cast_converter_cmap_attributes
{
  typename Map2::template Attribute_handle<i>::type operator()
  (const Map1& map1, Map2& map2, typename Map1::Dart_const_handle dh1,
   typename Map2::Dart_handle dh2) const
  {
    CGAL_USE(dh2);
    CGAL_assertion( map1.template attribute<i>(dh1)!=Map1::null_handle );
    CGAL_assertion( map2.template attribute<i>(dh2)==Map2::null_handle );
    typename Map2::template Attribute_handle<i>::type
      res = map2.template create_attribute<i>();
    map2.template get_attribute<i>(res).info() =
      (typename Map2::template Attribute_type<i>::type::Info)
        map1.template info<i>(dh1);
    return res;
  }
};
// ****************************************************************************
// "Converters" called during the copy of attributes, to copy Point (for
// attributes having such type defined).
// Users can replace them by their own converters.
// Point converter is called after Info converters; thus it is possible that
// attribute<0> was already created for dh2.
// Point converter is only called if both types of 0-attributes have
// Point type defined, and if dh1 has a 0-attribute.
// Map1 is the existing map, to convert into map2.
// ****************************************************************************
// Default converter for points. Point are copied only if they have same
// types, or in 2D/3D we use Cartesian_converter.
template< typename Map1, typename Map2>
struct Default_converter_cmap_0attributes_with_point
{
  typename Map2::template Attribute_handle<0>::type operator()
  (const Map1& map1, Map2& map2, typename Map1::Dart_const_handle dh1,
   typename Map2::Dart_handle dh2) const
  {
    CGAL_assertion( map1.template attribute<0>(dh1)!=Map1::null_handle );

    typename Map2::template Attribute_handle<0>::type
      res = map2.template attribute<0>(dh2);
    if ( res==Map2::null_handle )
    {
      res = map2.template create_attribute<0>();
    }
    internal::Set_point_if_possible_cmap
        <typename Map1::template Attribute_type<0>::type::Point,
        typename Map2::template Attribute_type<0>::type::Point>::
      run(map1.point(dh1), map2.template get_attribute<0>(res).point());
    return res;
  }
};
// ****************************************************************************
} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_COPY_FUNCTORS_H
