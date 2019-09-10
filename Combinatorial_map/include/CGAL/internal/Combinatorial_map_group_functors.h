// Copyright (c) 2010-2011 CNRS and LIRIS' Establishments (France).
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
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_COMBINATORIAL_MAP_GROUP_FUNCTORS_H
#define CGAL_COMBINATORIAL_MAP_GROUP_FUNCTORS_H

#include <CGAL/Unique_hash_map.h>
#include <CGAL/Combinatorial_map_functors.h>
#include <CGAL/Combinatorial_map_basic_operations.h>

/* Definition of functors used to group/ungroup attributes (we need functors
 * as attributes are stored in tuple, thus all the access must be done at
 * compiling time). Some of these functors are used with
 * Foreach_enabled_attributes to iterate through all the non void attribs.
 * These functors used other functors defined in Combinatorial_map_functors.h
 *
 * Group_attribute_functor_of_dart<CMap> to group the <i>-attributes of two
 *    given darts (except for j-dim). Only the attributes of the two given
 *    darts are possibly modified.
 *
 * Group_attribute_functor_of_dart_run<CMap,i> same than
 *   Group_attribute_functor_of_dart<CMap>::run<i>, with i template argument
 *   given in the struct to enable specialization.
 *
 * Group_attribute_functor<CMap> to group the <i>-attributes of two
 *    given i-cells (except for j-adim). If one i-attribute is NULL, we set the
 *    darts of its i-cell to the second attribute. If both i-attributes are
 *    non NULL, we overide all the i-attribute of the second i-cell to the
 *    first i-attribute.
 *
 * Degroup_attribute_functor_run<CMap> to degroup one i-attributes in two
 *   (except for j-adim).
 *
 * Test_split_attribute_functor<CMap,i> to test if there is some i-attributes
 *   that are split after an operation. Modified darts are given in a
 *   std::deque.
 */
namespace CGAL
{
namespace internal
{
// ************************************************************************
/// Functor used for link_beta to update the i-attributes of
/// dh2 on the attributes of dh1 dart, except if i=j.
///    (j is the dimension of the beta modified between dh1 and dh2,
///     so that after the modification we will have beta_j(dh1)==dh2)
/// Only attributes of dh1 or dh2 can be modified. If one dart as its
/// attribute equal to null, it takes the attributes of the second dart.
/// If both attributes are non null, dh2 takes the attribute of dh1.
template<typename CMap, unsigned int i, unsigned int j=CMap::dimension+1,
         typename T=typename CMap::template Attribute_type<i>::type>
struct Group_nonvoid_attribute_functor_of_dart_run
{
  /// Group the i-attribute of dh1 and dh2.
  static void run(CMap& amap,
                  typename CMap::Dart_handle dh1,
                  typename CMap::Dart_handle dh2)
  {
    CGAL_static_assertion( 1<=i && i<=CMap::dimension );
    CGAL_static_assertion( i!=j && (i!=1 || j!=0) );
    CGAL_static_assertion_msg(CMap::Helper::template
                              Dimension_index<i>::value>=0,
                              "Group_attribute_functor_of_dart_run<i> but "
                              "i-attributes are disabled");
    typename CMap::template Attribute_handle<i>::type
        a1=amap.template attribute<i>(dh1);
    typename CMap::template Attribute_handle<i>::type
        a2=amap.template attribute<i>(dh2);

    // If the two attributes are equal, nothing to do.
    if ( a1==a2 ) return;

    if ( a1==CMap::null_handle ) amap.template set_dart_attribute<i>(dh1, a2);
    else amap.template set_dart_attribute<i>(dh2, a1);
  }
};
// Specialization for i=0 and 2<=j. We modify 0-attribute for beta_j j>=2.
template<typename CMap, unsigned int j, typename T>
struct Group_nonvoid_attribute_functor_of_dart_run<CMap, 0, j, T>
{
  static void run(CMap& amap,
                  typename CMap::Dart_handle dh1,
                  typename CMap::Dart_handle dh2)
  {
    CGAL_static_assertion(j!=0 && j!=1);
    CGAL_static_assertion_msg(CMap::Helper::template
                              Dimension_index<0>::value>=0,
                              "Group_attribute_functor_of_dart_run<0> but "
                              "0-attributes are disabled");
    typename CMap::template Attribute_handle<0>::type
        a1=CMap::null_handle, a2=CMap::null_handle;

    // First extremity
    typename CMap::Dart_handle od = amap.other_extremity(dh2);
    if ( od!=amap.null_handle )
    {
      a1=amap.template attribute<0>(dh1);
      a2=amap.template attribute<0>(od);

      if ( a1==CMap::null_handle && a2!=CMap::null_handle )
        amap.template set_dart_attribute<0>(dh1, a2);
    }

    // Second extremity
    od = amap.other_extremity(dh1);
    if ( od!=amap.null_handle )
    {
      a1=amap.template attribute<0>(od);
      if ( a1!=CMap::null_handle )
        amap.template set_dart_attribute<0>(dh2, a1);
    }
  }
};
// Specialization for i=0 and j=0. We modify 0-attribute for beta_0.
template<typename CMap, typename T>
struct Group_nonvoid_attribute_functor_of_dart_run<CMap, 0, 0, T>
{
  static void run(CMap& amap,
                  typename CMap::Dart_handle dh1,
                  typename CMap::Dart_handle dh2)
  {
    CGAL_static_assertion_msg(CMap::Helper::template
                              Dimension_index<0>::value>=0,
                              "Group_attribute_functor_of_dart_run<0> but "
                              "0-attributes are disabled");
    typename CMap::Dart_handle od = amap.other_extremity(dh2);
    if ( od!=amap.null_handle )
    {
      typename CMap::template Attribute_handle<0>::type
          a1=amap.template attribute<0>(dh1);
      typename CMap::template Attribute_handle<0>::type
          a2=amap.template attribute<0>(od);

      if ( a1==CMap::null_handle && a2!=CMap::null_handle )
        amap.template set_dart_attribute<0>(dh1, a2);
    }
  }
};
// Specialization for i=0 and j=1. We modify 0-attribute for beta_1.
template<typename CMap, typename T>
struct Group_nonvoid_attribute_functor_of_dart_run<CMap, 0, 1, T>
{
  static void run(CMap& amap,
                  typename CMap::Dart_handle dh1,
                  typename CMap::Dart_handle dh2)
  {
    CGAL_static_assertion_msg(CMap::Helper::template
                              Dimension_index<0>::value>=0,
                              "Group_attribute_functor_of_dart_run<0> but "
                              "0-attributes are disabled");
    typename CMap::Dart_handle od = amap.other_extremity(dh1);
    if ( od!=amap.null_handle )
    {
      typename CMap::template Attribute_handle<0>::type
          a1=amap.template attribute<0>(od);

      if ( a1!=CMap::null_handle )
        amap.template set_dart_attribute<0>(dh2, a1);
    }
  }
};
// Specialization for i=j. Do nothing as j is the dimension to not consider.
template<typename CMap, unsigned int i, typename T>
struct Group_nonvoid_attribute_functor_of_dart_run<CMap,i,i,T>
{
  static void run(CMap&,
                  typename CMap::Dart_handle,
                  typename CMap::Dart_handle)
  {}
};
// Specialization for i=1 and j=0. Do nothing as edges attributes are not
// modify when we modify beta_0.
template<typename CMap, typename T>
struct Group_nonvoid_attribute_functor_of_dart_run<CMap,1,0,T>
{
  static void run(CMap&,
                  typename CMap::Dart_handle,
                  typename CMap::Dart_handle)
  {}
};
//------------------------------------------------------------------------------
template<typename CMap, unsigned int i, unsigned int j=CMap::dimension+1,
         typename T=typename CMap::template Attribute_type<i>::type>
struct Group_attribute_functor_of_dart_run
{
  /// Group the i-attribute of dh1 and dh2.
  static void run(CMap& amap,
                  typename CMap::Dart_handle dh1,
                  typename CMap::Dart_handle dh2)
  { Group_nonvoid_attribute_functor_of_dart_run<CMap, i, j, T>::
      run(amap, dh1, dh2); }
};
// Specialization for void attributes.
template<typename CMap, unsigned int i, unsigned int j>
struct Group_attribute_functor_of_dart_run<CMap, i, j, CGAL::Void>
{
  static void run(CMap&,
                  typename CMap::Dart_handle,
                  typename CMap::Dart_handle)
  {}
};
// ************************************************************************
/// Functor used for link_beta to update the attributes of
/// adart2 on the attributes of this dart, except for j-attributes.
///    (j is the dimension of the beta modified between adart1 and adart2,
///     so that after the modification we will have beta_j(adart1)==adart2)
/// We define run<i> to allows to use this functor with
/// Foreach_enabled_attributes.
///   If you know i at compiling time, use directly
///   Group_attribute_functor_of_dart_run.
template<typename CMap, unsigned int j=CMap::dimension+1>
struct Group_attribute_functor_of_dart
{
  template <unsigned int i>
  static void run(CMap& amap,
                  typename CMap::Dart_handle adart1,
                  typename CMap::Dart_handle adart2)
  {
    CGAL::internal::Group_attribute_functor_of_dart_run<CMap,i,j>::
        run(amap,adart1,adart2);
  }
};
// ************************************************************************
// Functor used to group the two i-attributes of the two i-cells, except
// the attribute of j
//    (j is the dimension of the beta modified between adart1 and adart2).
template<typename CMap, unsigned int i, unsigned int j=CMap::dimension+1,
         typename T=typename CMap::template Attribute_type<i>::type>
struct Group_nonvoid_attribute_functor_run
{
  static void run(CMap& amap,
                  typename CMap::Dart_handle adart1,
                  typename CMap::Dart_handle adart2)
  {
    CGAL_static_assertion( 1<=i && i<=CMap::dimension );
    CGAL_static_assertion( i!=j );
    CGAL_static_assertion_msg
        ( CMap::Helper::template Dimension_index<i>::value>=0,
          "Group_attribute_functor_run<i> but i-attributes are disabled" );
    typename CMap::template Attribute_handle<i>::type
        a1=amap.template attribute<i>(adart1);
    typename CMap::template Attribute_handle<i>::type
        a2=amap.template attribute<i>(adart2);

    // If the two attributes are equal, nothing to do.
    if ( a1 == a2 ) return;

    typename CMap::Dart_handle toSet = amap.null_handle;

    // If the attribute associated to adart1 is NULL, set it with
    // the attribute associated to adart2 (necessarily != NULL)
    if (a1 == CMap::null_handle)
    { toSet  = adart1; a1 = a2; }
    else
    {
      toSet = adart2;
      if (a2 != CMap::null_handle)
      {
        CGAL::internal::Call_merge_functor<CMap, i>::run(amap, a1, a2);
      }
    }
    amap.template set_attribute<i>(toSet, a1);
  }
};
// Specialization for i=0 and 2<=j. We update 0-attributes for beta_j j>=2.
// We need to update both extremities of the edge dh1.
template<typename CMap, unsigned int j, typename T>
struct Group_nonvoid_attribute_functor_run<CMap, 0, j, T>
{
  static void run( CMap& amap,
                   typename CMap::Dart_handle dh1,
                   typename CMap::Dart_handle dh2 )
  {
    CGAL_static_assertion_msg
        ( CMap::Helper::template Dimension_index<0>::value>=0,
          "Group_attribute_functor_run<0> but 0-attributes are disabled" );
    CGAL_static_assertion(j!=0 && j!=1);

    typename CMap::template Attribute_handle<0>::type
        a1=CMap::null_handle, a2=CMap::null_handle;
    typename CMap::Dart_handle toSet=amap.null_handle;
    // First extremity
    typename CMap::Dart_handle od=amap.other_extremity(dh2);
    if ( od!=amap.null_handle )
    {
      a1=amap.template attribute<0>(dh1);
      a2=amap.template attribute<0>(od);
      if ( a1!=a2 )
      {
        if ( a1==CMap::null_handle )
        { toSet=dh1; a1=a2; }
        else
        {
          toSet=od;
          if ( a2!=CMap::null_handle )
          {
            CGAL::internal::Call_merge_functor<CMap, 0>::run(amap, a1, a2);
          }
        }
        amap.template set_attribute<0>(toSet, a1);
      }
    }
    // Second extremity
    od = amap.other_extremity(dh1);
    if ( od!=amap.null_handle )
    {
      a1=amap.template attribute<0>(od);
      a2=amap.template attribute<0>(dh2);
      if ( a1!=a2 )
      {
        if ( a1==CMap::null_handle )
        { toSet=od; a1=a2; }
        else
        {
          toSet=dh2;
          if ( a2!=CMap::null_handle )
          {
            CGAL::internal::Call_merge_functor<CMap, 0>::run(amap, a1, a2);
          }
        }
        amap.template set_attribute<0>(toSet, a1);
      }
    }
  }
};
// Specialization for i=0 and j=0. We update 0-attributes for beta_0.
// We need to update the first extremity of the edge dh1.
template<typename CMap, typename T>
struct Group_nonvoid_attribute_functor_run<CMap, 0, 0, T>
{
  static void run( CMap& amap,
                   typename CMap::Dart_handle dh1,
                   typename CMap::Dart_handle dh2 )
  {
    CGAL_static_assertion_msg
        ( CMap::Helper::template Dimension_index<0>::value>=0,
          "Group_attribute_functor_run<0> but 0-attributes are disabled" );
    typename CMap::Dart_handle od=amap.other_extremity(dh2);
    if ( od!=amap.null_handle )
    {
      typename CMap::Dart_handle toSet=amap.null_handle;
      typename CMap::template Attribute_handle<0>::type
          a1=amap.template attribute<0>(dh1);
      typename CMap::template Attribute_handle<0>::type
          a2=amap.template attribute<0>(od);
      if ( a1!=a2 )
      {
        if ( a1==CMap::null_handle )
        { toSet=dh1; a1=a2; }
        else
        {
          toSet=od;
          if ( a2!=CMap::null_handle )
          {
            CGAL::internal::Call_merge_functor<CMap, 0>::run(amap, a1, a2);
          }
        }
        amap.template set_attribute<0>(toSet, a1);
      }
    }
  }
};
// Specialization for i=0 and j=1. We update 0-attributes for beta_1.
// We need to update the second extremity of the edge dh1.
template<typename CMap, typename T>
struct Group_nonvoid_attribute_functor_run<CMap, 0, 1, T>
{
  static void run( CMap& amap,
                   typename CMap::Dart_handle dh1,
                   typename CMap::Dart_handle dh2 )
  {
    CGAL_static_assertion_msg
        ( CMap::Helper::template Dimension_index<0>::value>=0,
          "Group_attribute_functor_run<0> but 0-attributes are disabled" );
    typename CMap::Dart_handle od =amap.other_extremity(dh1);
    if ( od!=amap.null_handle )
    {
      typename CMap::Dart_handle toSet=amap.null_handle;
      typename CMap::template Attribute_handle<0>::type
          a1=amap.template attribute<0>(od);
      typename CMap::template Attribute_handle<0>::type
          a2=amap.template attribute<0>(dh2);
      if ( a1!=a2 )
      {
        if ( a1==CMap::null_handle )
        { toSet=od; a1=a2; }
        else
        {
          toSet=dh2;
          if ( a2!=CMap::null_handle )
          {
            CGAL::internal::Call_merge_functor<CMap, 0>::run(amap, a1, a2);
          }
        }
        amap.template set_attribute<0>(toSet, a1);
      }
    }
  }
};
// Specialization for i=j. Do nothing as j is the dimension to not consider.
template<typename CMap, unsigned int i, typename T>
struct Group_nonvoid_attribute_functor_run<CMap,i,i,T>
{
  static void run(CMap&,
                  typename CMap::Dart_handle,
                  typename CMap::Dart_handle)
  {}
};
// Specialization for i=1 and j=0. Do nothing as edges attributes are not
// modify when we modify beta_0.
template<typename CMap, typename T>
struct Group_nonvoid_attribute_functor_run<CMap,1,0,T>
{
  static void run(CMap&,
                  typename CMap::Dart_handle,
                  typename CMap::Dart_handle)
  {}
};
//------------------------------------------------------------------------------
// Group_attribute_functor
template<typename CMap, unsigned int i, unsigned int j=CMap::dimension+1,
         typename T=typename CMap::template Attribute_type<i>::type>
struct Group_attribute_functor_run
{
  static void run( CMap& amap,
                   typename CMap::Dart_handle d1,
                   typename CMap::Dart_handle d2)
  { Group_nonvoid_attribute_functor_run<CMap, i, j, T>::run(amap, d1, d2); }
};
// Specialization for void attributes.
template<typename CMap, unsigned int i, unsigned int j>
struct Group_attribute_functor_run<CMap, i, j, CGAL::Void>
{
  static void run( CMap&,
                   typename CMap::Dart_handle,
                   typename CMap::Dart_handle )
  {}
};
// ************************************************************************
/// Functor used for sew to update the attributes of
/// adart2 on the attributes of this dart, except for j-attributes.
///    (j is the dimension of the beta modified between adart1 and adart2,
///     so that after the modification we will have beta_j(adart1)==adart2)
/// We define run<i> to allows to use this functor with
/// Foreach_enabled_attributes.
///   If you know i at compiling time, use directly
///   Group_attribute_functor_run.
template<typename CMap, unsigned int j=CMap::dimension+1>
struct Group_attribute_functor
{
  template <unsigned int i>
  static void run(CMap& amap,
                  typename CMap::Dart_handle adart1,
                  typename CMap::Dart_handle adart2)
  { CGAL::internal::Group_attribute_functor_run<CMap,i,j>::
        run(amap,adart1,adart2); }
};
// ************************************************************************
// Functor used to degroup one i-attribute of one i-cell in two, except the
// attribute of j.
template<typename CMap, unsigned int i, unsigned int j=CMap::dimension+1,
         typename T=typename CMap::template Attribute_type<i>::type>
struct Degroup_nonvoid_attribute_functor_run
{
  static void run(CMap& amap,
                  typename CMap::Dart_handle adart1,
                  typename CMap::Dart_handle adart2)
  {
    CGAL_static_assertion( i<=CMap::dimension );
    CGAL_static_assertion( i!=j );
    CGAL_static_assertion_msg
        ( CMap::Helper::template Dimension_index<i>::value>=0,
          "Degroup_attribute_functor_run<i> but i-attributes are disabled" );

    typename CMap::template Attribute_handle<i>::type
        a1=amap.template attribute<i>(adart1);

    // If there is no first attribute, nothing to degroup.
    if ( a1==CMap::null_handle ) return;

    // If the second attribute is non null and already different from a1,
    // nothing to do.
    if ( a1!=amap.template attribute<i>(adart2) &&
         amap.template attribute<i>(adart2)!=CMap::null_handle ) return;

    CGAL_assertion( (!CGAL::belong_to_same_cell<CMap,i>
                     (amap, adart1, adart2)) );

    // As we split, we set the dart of the first attribute to adart1 for which
    // we are sure it belongs to the first i-cell.
    amap.template set_dart_of_attribute<i>(a1, adart1);

    typename CMap::template Attribute_handle<i>::type
      a2 = amap.template copy_attribute<i>(a1);

    amap.template set_attribute<i>(adart2, a2);
    CGAL::internal::Call_split_functor<CMap, i>::run(amap, a1, a2);
  }
};
// Specialization for i==j.
template<typename CMap, unsigned int i, typename T>
struct Degroup_nonvoid_attribute_functor_run<CMap, i, i, T>
{
  static void run(CMap&,
                  typename CMap::Dart_handle,
                  typename CMap::Dart_handle)
  {}
};
//------------------------------------------------------------------------------
template<typename CMap, unsigned int i, unsigned int j=CMap::dimension+1,
         typename T=typename CMap::template Attribute_type<i>::type>
struct Degroup_attribute_functor_run
{
  static void run(CMap& amap,
                  typename CMap::Dart_handle adart1,
                  typename CMap::Dart_handle adart2)
  { Degroup_nonvoid_attribute_functor_run<CMap, i, j, T>::
      run(amap, adart1, adart2); }
};
// Specialization for void attributes.
template<typename CMap, unsigned int i, unsigned int j>
struct Degroup_attribute_functor_run<CMap, i, j, CGAL::Void>
{
  static void run(CMap&,
                  typename CMap::Dart_handle,
                  typename CMap::Dart_handle)
  {}
};
// ************************************************************************
// Function used by Test_split_attribute_functor_run to process one dart.
// Test the split of the i-cell containing the given dart adart.
// When we process a dart, we search in the Unique_hash_map if its
// i-attribute was already found. If yes, it means that we already
// found an i-cell with this attribute, thus this attribute is split.
// We mark (with mark) all the darts of the i-cell containing adart to
// process them exactly once.
template<typename CMap, unsigned int i>
void test_split_attribute_functor_one_dart
( CMap& amap, typename CMap::Dart_handle adart,
  CGAL::Unique_hash_map<typename CMap::template Attribute_handle<i>::type,
                        unsigned int, typename CMap::Hash_function> &
  found_attributes, typename CMap::size_type mark )
{
  CGAL_static_assertion_msg(CMap::Helper::template
                            Dimension_index<i>::value>=0,
                            "Test_split_attribute_functor_one_dart<i> but "
                            "i-attributes are disabled");

  typedef typename CMap::template Attribute_handle<i>::type
      Attribute_handle_i;

  // If the current dart has no attribute, or if it is aldready marked,
  // nothing to do.
  if ( amap.template attribute<i>(adart)==CMap::null_handle ||
       amap.is_marked(adart, mark) )
    return;

  Attribute_handle_i a1 = amap.template attribute<i>(adart);
  if ( found_attributes.is_defined(a1) )
  {  // Here the attribute was already present in the hash_map
    Attribute_handle_i a2 = amap.template
      create_attribute<i>(amap.template get_attribute<i>(a1));

    for ( CGAL::CMap_dart_iterator_basic_of_cell<CMap, i>
          itj(amap, adart, mark); itj.cont(); ++itj )
    {
      amap.template set_dart_attribute<i>(itj, a2);
      amap.mark(itj, mark);
    }
    CGAL::internal::Call_split_functor<CMap, i>::run(amap, a1, a2);
  }
  else
  {
    // Here the attribute was not in the hash_map.
    found_attributes[a1]=1;
    amap.template set_dart_of_attribute<i>(a1, adart);

    for ( CGAL::CMap_dart_iterator_basic_of_cell<CMap, i>
          itj(amap, adart, mark); itj.cont(); ++itj )
    {
      CGAL_assertion( amap.template attribute<i>(itj)==a1 );
      amap.mark(itj, mark);
    }
  }
}
// ************************************************************************
/// Functor used for unsew to test if i-attributes are split after an
/// operation, except for j-attributes.
///   (j is the dimension of the beta modified for darts in modified_darts,
///    if j==0 modified_darts2 are the darts modified for beta_1).
template<typename CMap, unsigned int i, unsigned int j=CMap::dimension+1,
         typename T=typename CMap::template Attribute_type<i>::type>
struct Test_split_nonvoid_attribute_functor_run
{
  // modified_darts is the set of modified darts for beta_j
  static void run( CMap& amap,
                   const std::deque<typename CMap::Dart_handle>
                   &modified_darts,
                   typename CMap::size_type mark_modified_darts=CMap::INVALID_MARK)
  {
    CGAL_static_assertion( 1<=i && i<=CMap::dimension );
    CGAL_assertion( i!=j );
    CGAL_static_assertion_msg(CMap::Helper::template
                              Dimension_index<i>::value>=0,
                              "Test_split_attribute_functor_run<i> but "
                              "i-attributes are disabled");

    typedef typename CMap::template Attribute_handle<i>::type
        Attribute_handle_i;

    CGAL::Unique_hash_map<Attribute_handle_i, unsigned int,
                          typename CMap::Hash_function> found_attributes;

    typename CMap::size_type mark = amap.get_new_mark(); // to mark incident cells.
    typename std::deque<typename CMap::Dart_handle>::const_iterator
        it=modified_darts.begin();
    for ( ; it!=modified_darts.end(); ++it )
    {
      CGAL::internal::test_split_attribute_functor_one_dart<CMap,i>
          (amap, *it, found_attributes, mark);
    }

    // Now we unmark all the marked darts.
    amap.negate_mark(mark);
    for ( it=modified_darts.begin(); it!=modified_darts.end(); ++it )
    {
      if ( mark_modified_darts!=CMap::INVALID_MARK )
        amap.unmark(*it, mark_modified_darts);

      if ( !amap.is_marked(*it, mark) )
        CGAL::mark_cell<CMap, i>(amap, *it, mark);
    }

    CGAL_assertion( amap.is_whole_map_marked(mark) );
    amap.free_mark(mark);
  }
  static void run( CMap& amap,
                   const std::deque<typename CMap::Dart_handle>
                   &modified_darts,
                   const std::deque<typename CMap::Dart_handle>
                   &modified_darts2,
                   typename CMap::size_type mark_modified_darts=CMap::INVALID_MARK)
  {
    CGAL_static_assertion( 1<=i && i<=CMap::dimension );
    CGAL_assertion( i!=j );
    CGAL_static_assertion_msg(CMap::Helper::template
                              Dimension_index<i>::value>=0,
                              "Test_split_attribute_functor_run<i> but "
                              "i-attributes are disabled");

    typedef typename CMap::template Attribute_handle<i>::type
        Attribute_handle_i;

    CGAL::Unique_hash_map<Attribute_handle_i, unsigned int,
                          typename CMap::Hash_function> found_attributes;

    typename CMap::size_type mark = amap.get_new_mark(); // to mark incident cells.
    typename std::deque<typename CMap::Dart_handle>::const_iterator
        it=modified_darts.begin();
    for ( ; it!=modified_darts.end(); ++it )
    {
      CGAL::internal::test_split_attribute_functor_one_dart<CMap,i>
          (amap, *it, found_attributes, mark);
    }
    typename std::deque<typename CMap::Dart_handle>::const_iterator
        it2=modified_darts2.begin();
    for ( ; it2!=modified_darts2.end(); ++it2 )
    {
      CGAL::internal::test_split_attribute_functor_one_dart<CMap,i>
          (amap, *it2, found_attributes, mark);
    }

    // Now we unmark all the marked darts.
    amap.negate_mark(mark);
    for ( it=modified_darts.begin(); it!=modified_darts.end(); ++it )
    {
      if ( mark_modified_darts!=CMap::INVALID_MARK )
        amap.unmark(*it, mark_modified_darts);

      if ( !amap.is_marked(*it, mark) )
        CGAL::mark_cell<CMap, i>(amap, *it, mark);
    }
    for ( it2=modified_darts2.begin(); it2!=modified_darts2.end(); ++it2 )
    {
      if ( mark_modified_darts!=CMap::INVALID_MARK )
        amap.unmark(*it2, mark_modified_darts);

      if ( !amap.is_marked(*it2, mark) )
        CGAL::mark_cell<CMap, i>(amap, *it2, mark);
    }

    CGAL_assertion( amap.is_whole_map_marked(mark) );
    amap.free_mark(mark);
  }
};
// Specialization for i=0 and 2<=j.
template<typename CMap, unsigned int j, typename T>
struct Test_split_nonvoid_attribute_functor_run<CMap, 0, j, T>
{
  static void run( CMap& amap,
                   const std::deque<typename CMap::Dart_handle>
                   &modified_darts,
                   typename CMap::size_type mark_modified_darts=CMap::INVALID_MARK)
  {
    CGAL_assertion( j!=0 && j!=1 );
    CGAL_static_assertion_msg(CMap::Helper::template
                              Dimension_index<0>::value>=0,
                              "Test_split_attribute_functor_run<0> but "
                              "0-attributes are disabled");

    typedef typename CMap::template Attribute_handle<0>::type
        Attribute_handle_0;

    CGAL::Unique_hash_map<Attribute_handle_0, unsigned int,
                          typename CMap::Hash_function> found_attributes;
    typename CMap::Dart_handle od=amap.null_handle;

    typename CMap::size_type mark = amap.get_new_mark(); // to mark incident cells.
    typename std::deque<typename CMap::Dart_handle>::const_iterator
        it=modified_darts.begin();
    for ( ; it!=modified_darts.end(); ++it )
    {
      CGAL::internal::test_split_attribute_functor_one_dart<CMap,0>
          (amap, *it, found_attributes, mark);

      od=amap.other_extremity(*it);
      if ( od!=amap.null_handle )
        CGAL::internal::test_split_attribute_functor_one_dart<CMap,0>
            (amap, od, found_attributes, mark);
    }

    // Now we unmark all the marked darts.
    amap.negate_mark(mark);
    for ( it=modified_darts.begin(); it!=modified_darts.end(); ++it )
    {
      if ( mark_modified_darts!=CMap::INVALID_MARK )
        amap.unmark(*it, mark_modified_darts);

      if ( !amap.is_marked(*it, mark) )
        CGAL::mark_cell<CMap, 0>(amap, *it, mark);

      od=amap.other_extremity(*it);
      if ( od!=amap.null_handle && !amap.is_marked(od, mark) )
        CGAL::mark_cell<CMap, 0>(amap, od, mark);
    }

    CGAL_assertion( amap.is_whole_map_marked(mark) );
    amap.free_mark(mark);
  }
  static void run( CMap& amap,
                   const std::deque<typename CMap::Dart_handle>
                   &modified_darts,
                   const std::deque<typename CMap::Dart_handle>
                   &modified_darts2,
                   typename CMap::size_type mark_modified_darts=CMap::INVALID_MARK)
  {
    CGAL_assertion( j!=0 && j!=1 );
    CGAL_static_assertion_msg(CMap::Helper::template
                              Dimension_index<0>::value>=0,
                              "Test_split_attribute_functor_run<0> but "
                              "0-attributes are disabled");

    typedef typename CMap::template Attribute_handle<0>::type
        Attribute_handle_0;

    CGAL::Unique_hash_map<Attribute_handle_0, unsigned int,
                          typename CMap::Hash_function> found_attributes;
    typename CMap::Dart_handle od=amap.null_handle;

    typename CMap::size_type mark = amap.get_new_mark(); // to mark incident cells.
    typename std::deque<typename CMap::Dart_handle>::const_iterator
        it=modified_darts.begin();
    for ( ; it!=modified_darts.end(); ++it )
    {
      CGAL::internal::test_split_attribute_functor_one_dart<CMap,0>
          (amap, *it, found_attributes, mark);

      od=amap.other_extremity(*it);
      if ( od!=amap.null_handle )
        CGAL::internal::test_split_attribute_functor_one_dart<CMap,0>
            (amap, od, found_attributes, mark);
    }
    typename std::deque<typename CMap::Dart_handle>::const_iterator
        it2=modified_darts2.begin();
    for ( ; it2!=modified_darts2.end(); ++it2 )
    {
      CGAL::internal::test_split_attribute_functor_one_dart<CMap,0>
          (amap, *it2, found_attributes, mark);

      od=amap.other_extremity(*it2);
      if ( od!=amap.null_handle )
        CGAL::internal::test_split_attribute_functor_one_dart<CMap,0>
            (amap, od, found_attributes, mark);
    }

    // Now we unmark all the marked darts.
    amap.negate_mark(mark);
    for ( it=modified_darts.begin(); it!=modified_darts.end(); ++it )
    {
      if ( mark_modified_darts!=CMap::INVALID_MARK )
        amap.unmark(*it, mark_modified_darts);

      if ( !amap.is_marked(*it, mark) )
        CGAL::mark_cell<CMap, 0>(amap, *it, mark);

      od=amap.other_extremity(*it);
      if ( od!=amap.null_handle && !amap.is_marked(od, mark) )
        CGAL::mark_cell<CMap, 0>(amap, od, mark);
    }
    for ( it2=modified_darts2.begin(); it2!=modified_darts2.end(); ++it2 )
    {
      if ( mark_modified_darts!=CMap::INVALID_MARK )
        amap.unmark(*it2, mark_modified_darts);

      if ( !amap.is_marked(*it2, mark) )
        CGAL::mark_cell<CMap, 0>(amap, *it2, mark);

      od=amap.other_extremity(*it2);
      if ( od!=amap.null_handle && !amap.is_marked(od, mark) )
        CGAL::mark_cell<CMap, 0>(amap, od, mark);
    }

    CGAL_assertion( amap.is_whole_map_marked(mark) );
    amap.free_mark(mark);
  }
};
// Specialization for i=0 and j=0.
// For j==0 or j==1, we use only the version with two list of darts,
// modified_darts are darts modified for beta0, and
// modified_darts2 are darts modified for beta1.
template<typename CMap, typename T>
struct Test_split_nonvoid_attribute_functor_run<CMap, 0, 0, T>
{
  static void run( CMap&, const std::deque<typename CMap::Dart_handle>&,
                   typename CMap::size_type =CMap::INVALID_MARK)
  { CGAL_assertion(false); }
  static void run( CMap& amap,
                   const std::deque<typename CMap::Dart_handle>
                   &modified_darts,
                   const std::deque<typename CMap::Dart_handle>
                   &modified_darts2,
                   typename CMap::size_type mark_modified_darts=CMap::INVALID_MARK)
  {
    CGAL_static_assertion_msg(CMap::Helper::template
                              Dimension_index<0>::value>=0,
                              "Test_split_attribute_functor_run<0> but "
                              "0-attributes are disabled");

    typedef typename CMap::template Attribute_handle<0>::type
        Attribute_handle_0;

    CGAL::Unique_hash_map<Attribute_handle_0, unsigned int,
                          typename CMap::Hash_function> found_attributes;
    typename CMap::Dart_handle od=amap.null_handle;

    typename CMap::size_type mark = amap.get_new_mark(); // to mark incident cells.
    typename std::deque<typename CMap::Dart_handle>::const_iterator
        it=modified_darts.begin();
    for ( ; it!=modified_darts.end(); ++it )
    {
      CGAL::internal::test_split_attribute_functor_one_dart<CMap,0>
          (amap, *it, found_attributes, mark);
    }
    typename std::deque<typename CMap::Dart_handle>::const_iterator
        it2=modified_darts2.begin();
    for ( ; it2!=modified_darts2.end(); ++it2 )
    {
      od=amap.other_extremity(*it2);
      if ( od!=amap.null_handle )
        CGAL::internal::test_split_attribute_functor_one_dart<CMap,0>
            (amap, od, found_attributes, mark);
    }

    // Now we unmark all the marked darts.
    amap.negate_mark(mark);
    for ( it=modified_darts.begin(); it!=modified_darts.end(); ++it )
    {
      if ( mark_modified_darts!=CMap::INVALID_MARK )
        amap.unmark(*it, mark_modified_darts);

      if ( !amap.is_marked(*it, mark) )
        CGAL::mark_cell<CMap, 0>(amap, *it, mark);
    }
    for ( it2=modified_darts2.begin(); it2!=modified_darts2.end(); ++it2 )
    {
      if ( mark_modified_darts!=CMap::INVALID_MARK )
        amap.unmark(*it2, mark_modified_darts);

      od=amap.other_extremity(*it2);
      if ( od!=amap.null_handle && !amap.is_marked(od, mark) )
        CGAL::mark_cell<CMap, 0>(amap, od, mark);
    }

    CGAL_assertion( amap.is_whole_map_marked(mark) );
    amap.free_mark(mark);
  }
};
// Specialization for i=0 and j=1.
// Equivalent to i=0 and j=0.
template<typename CMap, typename T>
struct Test_split_nonvoid_attribute_functor_run<CMap, 0, 1, T>
{
  static void run( CMap&, const std::deque<typename CMap::Dart_handle>&,
                   typename CMap::size_type =CMap::INVALID_MARK)
  { CGAL_assertion(false); }
  static void run( CMap& amap, const std::deque<typename CMap::Dart_handle>&
                   modified_darts,
                   const std::deque<typename CMap::Dart_handle>&
                   modified_darts2,
                   typename CMap::size_type mark_modified_darts=CMap::INVALID_MARK)
  { CGAL::internal::Test_split_nonvoid_attribute_functor_run<CMap, 0, 0, T>::
        run(amap, modified_darts, modified_darts2, mark_modified_darts); }
};
// Specialization for i=j.
template<typename CMap, unsigned int i, typename T>
struct Test_split_nonvoid_attribute_functor_run<CMap, i, i, T>
{
  static void run( CMap&, const std::deque<typename CMap::Dart_handle>&,
                   typename CMap::size_type=CMap::INVALID_MARK)
  {}
  static void run( CMap&, const std::deque<typename CMap::Dart_handle>&,
                   const std::deque<typename CMap::Dart_handle>&,
                   typename CMap::size_type=CMap::INVALID_MARK)
  {}
};
// Specialization for i=1 and j=0 (edge attributes are not modified
// when we modify beta_0).
template<typename CMap, typename T>
struct Test_split_nonvoid_attribute_functor_run<CMap, 1, 0, T>
{
  static void run( CMap&, const std::deque<typename CMap::Dart_handle>&,
                   typename CMap::size_type=CMap::INVALID_MARK)
  {}
  static void run( CMap&, const std::deque<typename CMap::Dart_handle>&,
                   const std::deque<typename CMap::Dart_handle>&,
                   typename CMap::size_type=CMap::INVALID_MARK)
  {}
};
//------------------------------------------------------------------------------
template<typename CMap, unsigned int i, unsigned int j=CMap::dimension+1,
         typename T=typename CMap::template Attribute_type<i>::type>
struct Test_split_attribute_functor_run
{
  // modified_darts is the set of modified darts for beta_j
  static void run( CMap& amap,
                   const std::deque<typename CMap::Dart_handle>
                   &modified_darts,
                   typename CMap::size_type mark_modified_darts=CMap::INVALID_MARK)
  { Test_split_nonvoid_attribute_functor_run<CMap, i, j, T>::
      run(amap, modified_darts, mark_modified_darts);
  }
  static void run( CMap& amap,
                   const std::deque<typename CMap::Dart_handle>
                   &modified_darts,
                   const std::deque<typename CMap::Dart_handle>
                   &modified_darts2,
                   typename CMap::size_type mark_modified_darts=CMap::INVALID_MARK)
  { Test_split_nonvoid_attribute_functor_run<CMap, i, j, T>::
      run(amap, modified_darts, modified_darts2, mark_modified_darts);
  }
};
// Specialization for void attributes.
template<typename CMap, unsigned int i, unsigned int j>
struct Test_split_attribute_functor_run<CMap, i, j, CGAL::Void>
{
  static void run( CMap&, const std::deque<typename CMap::Dart_handle>&,
                   typename CMap::size_type=CMap::INVALID_MARK)
  {}
  static void run( CMap&, const std::deque<typename CMap::Dart_handle>&,
                   const std::deque<typename CMap::Dart_handle>&,
                   typename CMap::size_type=CMap::INVALID_MARK)
  {}
};
// ************************************************************************
/// Functor used for unsew to test if i-attributes are split after an
/// operation, except for j-attributes.
/// We define run<i> to allows to use this functor with
/// Foreach_enabled_attributes.
template<typename CMap, unsigned int j=CMap::dimension+1>
struct Test_split_attribute_functor
{
  // Test the split of i-attributes, for all modified darts given in
  // modified_darts, and marked with mark_modified_darts.
  // For each split attribute, create a new i-attribute, associate
  // it with the new i-cell and call onsplit functors.
  template <unsigned int i>
  static void run( CMap& amap,
                   const std::deque<typename CMap::Dart_handle>
                   &modified_darts,
                   typename CMap::size_type mark_modified_darts=CMap::INVALID_MARK)
  {
    CGAL::internal::Test_split_attribute_functor_run<CMap, i, j>::
        run(amap, modified_darts, mark_modified_darts);
  }
  template <unsigned int i>
  static void run( CMap& amap,
                   const std::deque<typename CMap::Dart_handle>
                   &modified_darts,
                   const std::deque<typename CMap::Dart_handle>
                   &modified_darts2,
                   typename CMap::size_type mark_modified_darts=CMap::INVALID_MARK)
  {
    CGAL::internal::Test_split_attribute_functor_run<CMap, i, j>::
        run(amap, modified_darts, modified_darts2, mark_modified_darts);
  }
};
// ************************************************************************
} // namespace internal
} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_GROUP_FUNCTORS_H
