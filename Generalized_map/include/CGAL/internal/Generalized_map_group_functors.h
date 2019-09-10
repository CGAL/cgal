// Copyright (c) 2016 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_GENERALIZED_MAP_GROUP_FUNCTORS_H
#define CGAL_GENERALIZED_MAP_GROUP_FUNCTORS_H 1

#include <CGAL/Unique_hash_map.h>
#include <CGAL/internal/Generalized_map_internal_functors.h>
#include <CGAL/GMap_dart_iterators.h>
#include <CGAL/Combinatorial_map_functors.h>

/* Definition of functors used to group/ungroup attributes (we need functors
 * as attributes are stored in tuple, thus all the access must be done at
 * compiling time). Some of these functors are used with
 * Foreach_enabled_attributes to iterate through all the non void attribs.
 * These functors used other functors defined in Combinatorial_map_functors.h
 *
 * GMap_group_attribute_functor_of_dart<GMap> to group the <i>-attributes of two
 *    given darts (except for j-dim). Only the attributes of the two given
 *    darts are possibly modified.
 *
 * GMap_group_attribute_functor_of_dart_run<GMap,i> same than
 *   GMap_group_attribute_functor_of_dart<GMap>::run<i>, with i template argument
 *   given in the struct to enable specialization.
 *
 * GMap_group_attribute_functor<GMap> to group the <i>-attributes of two
 *    given i-cells (except for j-adim). If one i-attribute is NULL, we set the
 *    darts of its i-cell to the second attribute. If both i-attributes are
 *    non NULL, we overide all the i-attribute of the second i-cell to the
 *    first i-attribute.
 *
 * GMap_degroup_attribute_functor_run<GMap> to degroup one i-attributes in two
 *   (except for j-adim).
 *
 * GMap_test_split_attribute_functor<GMap,i> to test if there is some i-attributes
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
template<typename GMap, unsigned int i, unsigned int j=GMap::dimension+1,
         typename T=typename GMap::template Attribute_type<i>::type>
struct GMap_group_attribute_functor_of_dart_run
{
  /// Group the i-attribute of dh1 and dh2.
  static void run(GMap& amap,
                  typename GMap::Dart_handle dh1,
                  typename GMap::Dart_handle dh2)
  {
    CGAL_static_assertion( i<=GMap::dimension );
    CGAL_static_assertion( i!=j );
    CGAL_static_assertion_msg(GMap::Helper::template
                              Dimension_index<i>::value>=0,
                              "GMap_group_attribute_functor_of_dart_run<i> but "
                              "i-attributes are disabled");
    typename GMap::template Attribute_handle<i>::type
        a1=amap.template attribute<i>(dh1);
    typename GMap::template Attribute_handle<i>::type
        a2=amap.template attribute<i>(dh2);

    // If the two attributes are equal, nothing to do.
    if ( a1==a2 ) return;

    if ( a1==GMap::null_handle ) amap.template set_dart_attribute<i>(dh1, a2);
    else amap.template set_dart_attribute<i>(dh2, a1);
  }
};
// Specialization for void attributes.
template<typename GMap, unsigned int i, unsigned int j>
struct GMap_group_attribute_functor_of_dart_run<GMap, i, j, CGAL::Void>
{
  static void run(GMap&,
                  typename GMap::Dart_handle,
                  typename GMap::Dart_handle)
  {}
};
// Specialization for i=j. Do nothing as j is the dimension to not consider.
template<typename GMap, unsigned int i, typename T>
struct GMap_group_attribute_functor_of_dart_run<GMap,i,i,T>
{
  static void run(GMap&,
                  typename GMap::Dart_handle,
                  typename GMap::Dart_handle)
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
///   GMap_group_attribute_functor_of_dart_run.
template<typename GMap, unsigned int j=GMap::dimension+1>
struct GMap_group_attribute_functor_of_dart
{
  template <unsigned int i>
  static void run(GMap& amap,
                  typename GMap::Dart_handle adart1,
                  typename GMap::Dart_handle adart2)
  {
    CGAL::internal::GMap_group_attribute_functor_of_dart_run<GMap,i,j>::
        run(amap,adart1,adart2);
  }
};
// ************************************************************************
// Functor used to group the two i-attributes of the two i-cells, except
// the attribute of j
//    (j is the dimension of the beta modified between adart1 and adart2).
template<typename GMap, unsigned int i, unsigned int j=GMap::dimension+1,
         typename T=typename GMap::template Attribute_type<i>::type>
struct GMap_group_attribute_functor_run
{
  static void run(GMap& amap,
                  typename GMap::Dart_handle adart1,
                  typename GMap::Dart_handle adart2)
  {
    CGAL_static_assertion( i<=GMap::dimension );
    CGAL_static_assertion( i!=j );
    CGAL_static_assertion_msg
        ( GMap::Helper::template Dimension_index<i>::value>=0,
          "GMap_group_attribute_functor_run<i> but i-attributes are disabled" );
    typename GMap::template Attribute_handle<i>::type
        a1=amap.template attribute<i>(adart1);
    typename GMap::template Attribute_handle<i>::type
        a2=amap.template attribute<i>(adart2);

    // If the two attributes are equal, nothing to do.
    if ( a1 == a2 ) return;

    typename GMap::Dart_handle toSet = amap.null_handle;

    // If the attribute associated to adart1 is NULL, set it with
    // the attribute associated to adart2 (necessarily != NULL)
    if (a1 == GMap::null_handle)
    { toSet  = adart1; a1 = a2; }
    else
    {
      toSet = adart2;
      if (a2 != GMap::null_handle)
      {
        CGAL::internal::Call_merge_functor<GMap, i>::run(amap, a1, a2);
      }
    }
    amap.template set_attribute<i>(toSet, a1);
  }
};
// Specialization for void attributes.
template<typename GMap, unsigned int i, unsigned int j>
struct GMap_group_attribute_functor_run<GMap, i, j, CGAL::Void>
{
  static void run( GMap&,
                   typename GMap::Dart_handle,
                   typename GMap::Dart_handle )
  {}
};
// Specialization for i=j. Do nothing as j is the dimension to not consider.
template<typename GMap, unsigned int i, typename T>
struct GMap_group_attribute_functor_run<GMap,i,i,T>
{
  static void run(GMap&,
                  typename GMap::Dart_handle,
                  typename GMap::Dart_handle)
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
///   GMap_group_attribute_functor_run.
template<typename GMap, unsigned int j=GMap::dimension+1>
struct GMap_group_attribute_functor
{
  template <unsigned int i>
  static void run(GMap& amap,
                  typename GMap::Dart_handle adart1,
                  typename GMap::Dart_handle adart2)
  { CGAL::internal::GMap_group_attribute_functor_run<GMap,i,j>::
        run(amap,adart1,adart2); }
};
// ************************************************************************
// Functor used to degroup one i-attribute of one i-cell in two, except the
// attribute of j.
template<typename GMap, unsigned int i, unsigned int j=GMap::dimension+1,
         typename T=typename GMap::template Attribute_type<i>::type>
struct GMap_degroup_attribute_functor_run
{
  static void run(GMap& amap,
                  typename GMap::Dart_handle adart1,
                  typename GMap::Dart_handle adart2)
  {
    CGAL_static_assertion( i<=GMap::dimension );
    CGAL_static_assertion( i!=j );
    CGAL_static_assertion_msg
        ( GMap::Helper::template Dimension_index<i>::value>=0,
          "GMap_degroup_attribute_functor_run<i> but i-attributes are disabled" );

    typename GMap::template Attribute_handle<i>::type
        a1=amap.template attribute<i>(adart1);

    // If there is no first attribute, nothing to degroup.
    if ( a1==GMap::null_handle ) return;

    // If the second attribute is non null and already different from a1,
    // nothing to do.
    if ( a1!=amap.template attribute<i>(adart2) &&
         amap.template attribute<i>(adart2)!=GMap::null_handle ) return;

    CGAL_assertion( (!CGAL::belong_to_same_cell<GMap,i>
                     (amap, adart1, adart2)) );

    // As we split, we set the dart of the first attribute to adart1 for which
    // we are sure it belongs to the first i-cell.
    amap.template set_dart_of_attribute<i>(a1, adart1);

    typename GMap::template Attribute_handle<i>::type
      a2 = amap.template copy_attribute<i>(a1);

    amap.template set_attribute<i>(adart2, a2);
    CGAL::internal::Call_split_functor<GMap, i>::run(amap, a1, a2);
  }
};
// Specialization for void attributes.
template<typename GMap, unsigned int i, unsigned int j>
struct GMap_degroup_attribute_functor_run<GMap, i, j, CGAL::Void>
{
  static void run(GMap&,
                  typename GMap::Dart_handle,
                  typename GMap::Dart_handle)
  {}
};
// Specialization for i==j.
template<typename GMap, unsigned int i, typename T>
struct GMap_degroup_attribute_functor_run<GMap, i, i, T>
{
  static void run(GMap&,
                  typename GMap::Dart_handle,
                  typename GMap::Dart_handle)
  {}
};
// ************************************************************************
// Function used by GMap_test_split_attribute_functor_run to process one dart.
// Test the split of the i-cell containing the given dart adart.
// When we process a dart, we search in the Unique_hash_map if its
// i-attribute was already found. If yes, it means that we already
// found an i-cell with this attribute, thus this attribute is split.
// We mark (with mark) all the darts of the i-cell containing adart to
// process them exactly once.
template<typename GMap, unsigned int i>
void GMap_test_split_attribute_functor_one_dart
( GMap& amap, typename GMap::Dart_handle adart,
  CGAL::Unique_hash_map<typename GMap::template Attribute_handle<i>::type,
                        unsigned int, typename GMap::Hash_function> &
  found_attributes, typename GMap::size_type mark )
{
  CGAL_static_assertion_msg(GMap::Helper::template
                            Dimension_index<i>::value>=0,
                            "GMap_test_split_attribute_functor_one_dart<i> but "
                            "i-attributes are disabled");

  typedef typename GMap::template Attribute_handle<i>::type
      Attribute_handle_i;

  // If the current dart has no attribute, or if it is aldready marked,
  // nothing to do.
  if ( amap.template attribute<i>(adart)==GMap::null_handle ||
       amap.is_marked(adart, mark) )
    return;

  Attribute_handle_i a1 = amap.template attribute<i>(adart);
  if ( found_attributes.is_defined(a1) )
  {  // Here the attribute was already present in the hash_map
    Attribute_handle_i a2 = amap.template
      create_attribute<i>(amap.template get_attribute<i>(a1));

    for ( CGAL::GMap_dart_iterator_basic_of_cell<GMap, i>
          itj(amap, adart, mark); itj.cont(); ++itj )
    {
      amap.template set_dart_attribute<i>(itj, a2);
      amap.mark(itj, mark);
    }
    CGAL::internal::Call_split_functor<GMap, i>::run(amap, a1, a2);
  }
  else
  {
    // Here the attribute was not in the hash_map.
    found_attributes[a1]=1;
    amap.template set_dart_of_attribute<i>(a1, adart);

    for ( CGAL::GMap_dart_iterator_basic_of_cell<GMap, i>
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
template<typename GMap, unsigned int i, unsigned int j=GMap::dimension+1,
         typename T=typename GMap::template Attribute_type<i>::type>
struct GMap_test_split_attribute_functor_run
{
  // modified_darts is the set of modified darts for beta_j
  static void run( GMap& amap,
                   const std::deque<typename GMap::Dart_handle>
                   &modified_darts,
                   typename GMap::size_type mark_modified_darts=GMap::INVALID_MARK)
  {
    CGAL_static_assertion( i<=GMap::dimension );
    CGAL_assertion( i!=j );
    CGAL_static_assertion_msg(GMap::Helper::template
                              Dimension_index<i>::value>=0,
                              "GMap_test_split_attribute_functor_run<i> but "
                              "i-attributes are disabled");

    typedef typename GMap::template Attribute_handle<i>::type
        Attribute_handle_i;

    CGAL::Unique_hash_map<Attribute_handle_i, unsigned int,
                          typename GMap::Hash_function> found_attributes;

    typename GMap::size_type mark = amap.get_new_mark(); // to mark incident cells.
    typename std::deque<typename GMap::Dart_handle>::const_iterator
        it=modified_darts.begin();
    for ( ; it!=modified_darts.end(); ++it )
    {
      CGAL::internal::GMap_test_split_attribute_functor_one_dart<GMap,i>
          (amap, *it, found_attributes, mark);
    }

    // Now we unmark all the marked darts.
    amap.negate_mark(mark);
    for ( it=modified_darts.begin(); it!=modified_darts.end(); ++it )
    {
      if ( mark_modified_darts!=GMap::INVALID_MARK )
        amap.unmark(*it, mark_modified_darts);

      if ( !amap.is_marked(*it, mark) )
        CGAL::mark_cell<GMap, i>(amap, *it, mark);
    }

    CGAL_assertion( amap.is_whole_map_marked(mark) );
    amap.free_mark(mark);
  }
  static void run( GMap& amap,
                   const std::deque<typename GMap::Dart_handle>
                   &modified_darts,
                   const std::deque<typename GMap::Dart_handle>
                   &modified_darts2,
                   typename GMap::size_type mark_modified_darts=GMap::INVALID_MARK)
  {
    CGAL_static_assertion( i<=GMap::dimension );
    CGAL_assertion( i!=j );
    CGAL_static_assertion_msg(GMap::Helper::template
                              Dimension_index<i>::value>=0,
                              "GMap_test_split_attribute_functor_run<i> but "
                              "i-attributes are disabled");

    typedef typename GMap::template Attribute_handle<i>::type
        Attribute_handle_i;

    CGAL::Unique_hash_map<Attribute_handle_i, unsigned int,
                          typename GMap::Hash_function> found_attributes;

    typename GMap::size_type mark = amap.get_new_mark(); // to mark incident cells.
    typename std::deque<typename GMap::Dart_handle>::const_iterator
        it=modified_darts.begin();
    for ( ; it!=modified_darts.end(); ++it )
    {
      CGAL::internal::GMap_test_split_attribute_functor_one_dart<GMap,i>
          (amap, *it, found_attributes, mark);
    }
    typename std::deque<typename GMap::Dart_handle>::const_iterator
        it2=modified_darts2.begin();
    for ( ; it2!=modified_darts2.end(); ++it2 )
    {
      CGAL::internal::GMap_test_split_attribute_functor_one_dart<GMap,i>
          (amap, *it2, found_attributes, mark);
    }

    // Now we unmark all the marked darts.
    amap.negate_mark(mark);
    for ( it=modified_darts.begin(); it!=modified_darts.end(); ++it )
    {
      if ( mark_modified_darts!=GMap::INVALID_MARK )
        amap.unmark(*it, mark_modified_darts);

      if ( !amap.is_marked(*it, mark) )
        CGAL::mark_cell<GMap, i>(amap, *it, mark);
    }
    for ( it2=modified_darts2.begin(); it2!=modified_darts2.end(); ++it2 )
    {
      if ( mark_modified_darts!=GMap::INVALID_MARK )
        amap.unmark(*it2, mark_modified_darts);

      if ( !amap.is_marked(*it2, mark) )
        CGAL::mark_cell<GMap, i>(amap, *it2, mark);
    }

    CGAL_assertion( amap.is_whole_map_marked(mark) );
    amap.free_mark(mark);
  }
};
// Specialization for void attributes.
template<typename GMap, unsigned int i, unsigned int j>
struct GMap_test_split_attribute_functor_run<GMap, i, j, CGAL::Void>
{
  static void run( GMap&, const std::deque<typename GMap::Dart_handle>&,
                   typename GMap::size_type=GMap::INVALID_MARK)
  {}
  static void run( GMap&, const std::deque<typename GMap::Dart_handle>&,
                   const std::deque<typename GMap::Dart_handle>&,
                   typename GMap::size_type=GMap::INVALID_MARK)
  {}
};
// Specialization for i=j.
template<typename GMap, unsigned int i, typename T>
struct GMap_test_split_attribute_functor_run<GMap, i, i, T>
{
  static void run( GMap&, const std::deque<typename GMap::Dart_handle>&,
                   typename GMap::size_type=GMap::INVALID_MARK)
  {}
  static void run( GMap&, const std::deque<typename GMap::Dart_handle>&,
                   const std::deque<typename GMap::Dart_handle>&,
                   typename GMap::size_type=GMap::INVALID_MARK)
  {}
};
// Specialization for i=1 and j=0 (edge attributes are not modified
// when we modify beta_0).
template<typename GMap, typename T>
struct GMap_test_split_attribute_functor_run<GMap, 1, 0, T>
{
  static void run( GMap&, const std::deque<typename GMap::Dart_handle>&,
                   typename GMap::size_type=GMap::INVALID_MARK)
  {}
  static void run( GMap&, const std::deque<typename GMap::Dart_handle>&,
                   const std::deque<typename GMap::Dart_handle>&,
                   typename GMap::size_type=GMap::INVALID_MARK)
  {}
};
// ************************************************************************
/// Functor used for unsew to test if i-attributes are split after an
/// operation, except for j-attributes.
/// We define run<i> to allows to use this functor with
/// Foreach_enabled_attributes.
template<typename GMap, unsigned int j=GMap::dimension+1>
struct GMap_test_split_attribute_functor
{
  // Test the split of i-attributes, for all modified darts given in
  // modified_darts, and marked with mark_modified_darts.
  // For each split attribute, create a new i-attribute, associate
  // it with the new i-cell and call onsplit functors.
  template <unsigned int i>
  static void run( GMap& amap,
                   const std::deque<typename GMap::Dart_handle>
                   &modified_darts,
                   typename GMap::size_type mark_modified_darts=GMap::INVALID_MARK)
  {
    CGAL::internal::GMap_test_split_attribute_functor_run<GMap, i, j>::
        run(amap, modified_darts, mark_modified_darts);
  }
  template <unsigned int i>
  static void run( GMap& amap,
                   const std::deque<typename GMap::Dart_handle>
                   &modified_darts,
                   const std::deque<typename GMap::Dart_handle>
                   &modified_darts2,
                   typename GMap::size_type mark_modified_darts=GMap::INVALID_MARK)
  {
    CGAL::internal::GMap_test_split_attribute_functor_run<GMap, i, j>::
        run(amap, modified_darts, modified_darts2, mark_modified_darts);
  }
};
// ************************************************************************
} // namespace internal
} // namespace CGAL

#endif // CGAL_GENERALIZED_MAP_GROUP_FUNCTORS_H
