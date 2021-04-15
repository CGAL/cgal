// Copyright (c) 2010-2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_COMBINATORIAL_MAP_FUNCTORS_H
#define CGAL_COMBINATORIAL_MAP_FUNCTORS_H

#include <CGAL/Dart_const_iterators.h>
#include <CGAL/internal/Combinatorial_map_internal_functors.h>
#include <vector>
#include <boost/mpl/has_xxx.hpp>

/* Definition of functors used to manage attributes (we need functors as
 * attributes are stored in tuple, thus all the access must be done at
 * compiling time). Some of these functors are used with
 * Foreach_enabled_attributes to iterate through all the non void attribs.
 * Functors allowing to group/ungroup attributes are defined in
 * Combinatorial_map_group_functors.h. Some internal functors are defined
 * in internal/Combinatorial_map_internal_functors.h.
 *
 * Reserve_mark_functor<CMap> to reserve one mark, used with
 *   Foreach_enabled_attributes to reserve a mark for each non void attribute.
 *
 * Display_attribute_functor<CMap> to display the address of the i-attribute
 *   of a given dart (can be used with Foreach_enabled_attributes)
 *
 * Set_i_attribute_functor<CMap, i> to set the i-attribute of a given
 *   i-cell.
 *
 * Test_is_valid_attribute_functor<CMap> to test if an attribute is valid
 *    (used with Foreach_enabled_attributes)
 *
 * Is_attribute_has_non_void_info<Attr> to test if the attribute
 *   Attr is non Void and has an non void Info as inner type
 *
 * Is_attribute_has_point<Attr> to test if the attribute
 *   Attr is non Void and has a Point inner type
 *
 */

namespace CGAL
{
/** @file Combinatorial_map_functors.h
 * Definition of functors used for dD Combinatorial map.
 */
// ****************************************************************************
/// Functor used to reserve one mark, used with Foreach_enabled_attributes
/// to reserve a mark for each enabled attribute.
template<typename CMap>
struct Reserve_mark_functor
{
  typedef typename CMap::size_type size_type;

  template <unsigned int i>
  static void run(const CMap& amap, std::vector<size_type>& marks)
  { marks[i] = amap.get_new_mark(); }
};
// ****************************************************************************
/// Functor used to display the address of the i-cell attribute. Can be used
/// with Foreach_enabled_attributes.
template<typename CMap>
struct Display_attribute_functor
{
  template <unsigned int i>
  static void run(const CMap& amap,
                  typename CMap::Dart_const_handle adart)
  {
    if ( amap.template attribute<i>(adart)==nullptr )
      std::cout<<"nullptr";
    else
      amap.template display_attribute<i>(amap.template attribute<i>(adart));
  }
};
// ****************************************************************************
/// Functor used to test if a cell is valid
template<typename CMap>
struct Test_is_valid_attribute_functor
{
  template <unsigned int i>
  static bool run(const CMap& amap,
                  typename CMap::Dart_const_handle adart)
  {
    typedef typename CMap::size_type size_type;

    size_type mark=amap.get_new_mark();
    bool res = true;
    CGAL::internal::Test_is_valid_attribute_functor<CMap>::
        run<i>(amap, adart, mark, &res);

    amap.negate_mark(mark);
    if ( !amap.is_whole_map_marked(mark) )
    {
      for ( CGAL::CMap_dart_const_iterator_basic_of_cell<CMap,i>
            it(*amap, adart, mark); it.cont(); ++it )
        amap.unmark(it, mark);
    }
    CGAL_assertion ( amap.is_whole_map_marked(mark) );
    amap.free_mark(mark);

    return res;
  }
};
// ****************************************************************************
/// Functor used to set the i-attribute of a given i-cell.
template<typename CMap, unsigned int i,
         typename T=typename CMap::template Attribute_type<i>::type>
struct Set_i_attribute_functor
{
  static void run( CMap& amap, typename CMap::Dart_handle dh,
                   typename CMap::template Attribute_handle<i>::type ah )
  {
    amap.template set_attribute<i>(dh, ah);
  }
};
/// Specialization for void attributes.
template<typename CMap, unsigned int i>
struct Set_i_attribute_functor<CMap,i,CGAL::Void>
{
  static void run( CMap&, typename CMap::Dart_handle,
                   typename CMap::template Attribute_handle<i>::type)
  {}
};
// ****************************************************************************
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_point,Point,false)

template<typename Attr, typename Info=typename Attr::Info>
struct Is_nonvoid_attribute_has_non_void_info
{
  static const bool value=true;
};
template<typename Attr>
struct Is_nonvoid_attribute_has_non_void_info<Attr, void>
{
  static const bool value=false;
};

template<typename Attr>
struct Is_attribute_has_non_void_info
{
  static const bool value=Is_nonvoid_attribute_has_non_void_info<Attr>::value;
};
template<>
struct Is_attribute_has_non_void_info<CGAL::Void>
{
  static const bool value=false;
};
// ****************************************************************************
template<typename Attr>
struct Is_attribute_has_point
{ static const bool value=Has_point<Attr>::value; };
// ****************************************************************************
} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_FUNCTORS_H //
// EOF //
