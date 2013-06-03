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

/* Definition of functors used internally to copy combinatorial maps attributes
 * (we need functors as attributes are stored in tuple, thus all the access
 *  must be done at compiling time).
 */
namespace CGAL
{
// ****************************************************************************
namespace internal
{
// ****************************************************************************
// Map1 is the existing map, to convert into map2.
// Case where the two i-attributes are non void.
template< typename Map1, typename Map2, unsigned int i, typename Info1, typename Info2 >
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

// Map1 is the existing map, to convert into map2.
template< typename Map1, typename Map2, unsigned int i,
          typename Attr1=typename Map1::template Attribute_type<i>::type,
          typename Attr2=typename Map2::template Attribute_type<i>::type >
struct Default_converter_cmap_attr
{
  typename Map2::template Attribute_handle<i>::type operator()
  (Map2& map2, typename Map1::Dart_const_handle dh1)
  { return Default_converter_two_non_void_attributes_cmap
      <Map1,Map2,i,typename Attr1::Info,typename Attr2::Info>::
      run(map2, dh1->template attribute<i>()); }
};

template< typename Map1, typename Map2, unsigned int i,
          typename Attr1>
struct Default_converter_cmap_attr<Map1, Map2, i, Attr1, CGAL::Void>
{
  typename Map2::template Attribute_handle<i>::type operator()
  (Map2&, typename Map1::Dart_const_handle)
  { return NULL; }
};

template< typename Map1, typename Map2, unsigned int i,
          typename Attr2>
struct Default_converter_cmap_attr<Map1, Map2, i, CGAL::Void, Attr2>
{
  typename Map2::template Attribute_handle<i>::type operator()
  (Map2&, typename Map1::Dart_const_handle)
  { return NULL; }
};

template< typename Map1, typename Map2, unsigned int i>
struct Default_converter_cmap_attr<Map1, Map2, i, CGAL::Void, CGAL::Void>
{
  typename Map2::template Attribute_handle<i>::type operator()
  (Map2&, typename Map1::Dart_const_handle)
  { return NULL; }
};

template<typename Map1, typename Map2, typename Converters, unsigned int i,
         bool t=(i>=My_length<Converters>::value)>
struct Convert_attribute_functor
{
  static typename Map2::template Attribute_handle<i>::type
  run( Map2* cmap2, typename Map1::Dart_const_handle dh1,
       const Converters& converters)
  {
    return Default_converter_cmap_attr<Map1, Map2, i>() (*cmap2, dh1);
  }
};

template<typename Map1, typename Map2, typename Converters, unsigned int i>
struct Convert_attribute_functor<Map1,Map2,Converters,i,false>
{
  static typename Map2::template Attribute_handle<i>::type
  run( Map2* cmap2, typename Map1::Dart_const_handle dh1,
       const Converters& converters)
  {
    return CGAL::cpp11::get<i>(converters) (*cmap2, dh1);
  }
};


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
    typename Map2::template Attribute_handle<i>::type res=NULL;

 /*   if ( i>=My_length<Converters>::value )
      res = Default_converter_cmap_attr<Map1, Map2, i>()
          (*cmap2, dh1->template attribute<i>());
    else
      res =CGAL::cpp11::get<i>(converters)
          (*cmap2, dh1->template attribute<i>());*/
    res=Convert_attribute_functor<Map1,Map2,Converters,i>::run(cmap2,dh1,converters);

    if ( res!=NULL )
      cmap2->template set_attribute<i>(dh2, res);
  }
};
// ****************************************************************************
} // namespace internal
} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_COPY_FUNCTORS_H
