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
#ifndef CGAL_COMBINATORIAL_MAP_INTERNAL_FUNCTORS_H
#define CGAL_COMBINATORIAL_MAP_INTERNAL_FUNCTORS_H

#include <CGAL/Dart_const_iterators.h>
#include <CGAL/Combinatorial_map_basic_operations.h>
#include <CGAL/Dimension.h>
#include <CGAL/Kernel_traits.h>
#include <vector>

/* Definition of functors used internally to manage attributes (we need
 * functors as attributes are stored in tuple, thus all the access must be
 * done at compiling time). Some of these functors are used with
 * Foreach_enabled_attributes to iterate through all the non void attribs.
 * Functors allowing to group/ungroup attributes are defined in
 * internal/Combinatorial_map_group_functors.h. Public functions are defined
 * in Combinatorial_map_functors.h.
 *
 * internal::Swap_attributes_functor<CMap> to swap the i-attributes between
 *    two cmaps having same type
 *
 * internal::Call_split_functor<CMap,i> to call the OnSplit functors on two
 *    given i-attributes.
 *
 * internal::Call_merge_functor<CMap,i> to call the OnMerge functors on two
 *    given i-attributes.
 *
 * internal::Test_is_valid_attribute_functor<CMap> to test if a given i-cell is
 *    valid (all its darts are linked to the same attribute, no other dart is
 *    linked with this attribute).
 *
 * internal::Count_cell_functor<CMap> to count the nuber of i-cells.
 *
 * internal::Count_bytes_one_attribute_functor<CMap> to count the memory
 *    occupied by i-attributes.
 *
 * internal::Decrease_attribute_functor<CMap> to decrease by one the ref
 *    counting of a given i-attribute.
 *
 * internal::Restricted_decrease_attribute_functor<CMap> to decrease by one the
 *    ref counting of a given i-attribute, but without deleting attribute
 *    having nomore dart associated with.
 *
 * internal::Beta_functor<Dart, i...> to call several beta on the given dart.
 *   Indices are given as parameter of the run function.
 *
 * internal::Beta_functor_static<Dart, i...> to call several beta on the given
 *   dart. Indices are given as template arguments.
 *
 * internal::Set_i_attribute_of_dart_functor<CMap, i> to set the i-attribute
 *   of a given dart.
 *
 * internal::Test_is_same_dart_info_functor<Map1, Map2> to test if two
 *   darts have the same info.
 *
 * internal::Test_is_same_attribute_functor<Map1, Map2> to test if two
 *   i-attributes of two darts are isomorphic (ie they have the same info).
 *
 * inernal::Test_is_same_attribute_point_functor<Map1, Map2, i> to test if
 *   the point of two i-attributes are equal.
 *
 * internal::Reverse_orientation_of_map_functor<CMap> to reverse the
 *   orientation of a whole combinatorial map
 *
 * internal::Reverse_orientation_of_connected_component_functor to reverse
 *   the orientation of a connected component in a cmap
 *
 * internal::Init_attribute_functor to initialize all attributes to NULL.
 *
 * internal::Correct_invalid_attributes_functor to correct the i-attribute
 *   associated with a given i-cell
 *
 * internal::Cleanup_useless_attributes to erase all attributes having
 *   no more associated dart
 *
 * internal::Init_id initialize the id of an new element created in a compact
 *   container, when the underlying type defines id (TODO Move in Compact container ?)
 */

namespace CGAL
{
//-----------------------------------------------------------------------------
template<typename Attr>
struct Is_attribute_has_non_void_info;
template<typename Attr>
struct Is_attribute_has_point;
// ****************************************************************************
namespace internal
{
// ****************************************************************************
/// Swap i-attribute between the two maps having same type.
template<typename CMap>
struct Swap_attributes_functor
{
  template<unsigned int i>
  static void run( CMap& cmap1,
                   CMap& cmap2)
  { CGAL::cpp11::get<CMap::Helper::template Dimension_index<i>::value>
        (cmap1.mattribute_containers).swap(
          CGAL::cpp11::get<CMap::Helper::template Dimension_index<i>::value>
          (cmap2.mattribute_containers));
   }
};
// ****************************************************************************
// Functor which call Functor::operator() on the two given cell_attributes
 template<typename CMap, typename Cell_attribute, typename Functor>
struct Apply_cell_functor
{
  static void run(Cell_attribute& acell1, Cell_attribute& acell2)
  {
    Functor() (acell1, acell2);
  }
};
//...except for Null_functor.
template<typename CMap, typename Cell_attribute>
struct Apply_cell_functor<CMap, Cell_attribute, CGAL::Null_functor>
{
  static void run(Cell_attribute&, Cell_attribute&)
  {}
};
// ****************************************************************************
// Functor used to call the On_split functor between the two given darts.
template<typename CMap, unsigned int i,
         typename Enabled=typename CMap::
       #ifndef CGAL_CFG_TEMPLATE_IN_DEFAULT_PARAMETER_BUG
         template
       #endif
         Attribute_type<i>::type>
struct Call_split_functor
{
  typedef typename CMap::template Attribute_type<i>::type Attribute;
  typedef typename Attribute::On_split On_split;

  static void
  run(CMap& amap, typename CMap::template Attribute_handle<i>::type a1,
      typename CMap::template Attribute_handle<i>::type a2)
  {
    // Static version
    CGAL::internal::Apply_cell_functor<CMap, Attribute, On_split>::
      run(amap.template get_attribute<i>(a1),
          amap.template get_attribute<i>(a2));
    // Dynamic version
    if ( CGAL::cpp11::get<CMap::Helper::template Dimension_index<i>::value>
         (amap.m_onsplit_functors) )
      CGAL::cpp11::get<CMap::Helper::template Dimension_index<i>::value>
        (amap.m_onsplit_functors) (amap.template get_attribute<i>(a1),
                                    amap.template get_attribute<i>(a2));
  }
};
// Specialization for disabled attributes.
template<typename CMap,unsigned int i>
struct Call_split_functor<CMap, i, CGAL::Void>
{
  static void run(typename CMap::Dart_handle,
                  typename CMap::Dart_handle)
  {}
};
// ****************************************************************************
// Functor used to call the On_merge functor between the two given darts.
template<typename CMap,unsigned int i,
         typename Enabled=typename CMap::
       #ifndef CGAL_CFG_TEMPLATE_IN_DEFAULT_PARAMETER_BUG
         template
       #endif
         Attribute_type<i>::type>
struct Call_merge_functor
{
  typedef typename CMap::template Attribute_type<i>::type Attribute;
  typedef typename Attribute::On_merge On_merge;

  static void
  run(CMap& amap, typename CMap::template Attribute_handle<i>::type a1,
      typename CMap::template Attribute_handle<i>::type a2)
  {
    // Static version
    CGAL::internal::Apply_cell_functor<CMap, Attribute, On_merge>::
      run(amap.template get_attribute<i>(a1),
          amap.template get_attribute<i>(a2));
    // Dynamic version
    if ( CGAL::cpp11::get<CMap::Helper::template Dimension_index<i>::value>
         (amap.m_onmerge_functors) )
      CGAL::cpp11::get<CMap::Helper::template Dimension_index<i>::value>
        (amap.m_onmerge_functors) (amap.template get_attribute<i>(a1),
                                   amap.template get_attribute<i>(a2));
  }
};
// Specialization for disabled attributes.
template<typename CMap,unsigned int i>
struct Call_merge_functor<CMap, i, CGAL::Void>
{
  static void run(CMap&, typename CMap::Dart_handle,
                  typename CMap::Dart_handle)
  {}
};
// ****************************************************************************
/// Functor used to test if a cell is valid
template<typename CMap>
struct Test_is_valid_attribute_functor
{
  /** Test the validity of a i-cell-attribute.
   * ie all the darts belonging to a i-cell are linked to the same attribute.
   * @param adart a dart.
   * @param amark a mark used to mark darts of the i-cell.
   * @return true iff all the darts of the i-cell link to the same attribute.
   */
  typedef typename CMap::size_type size_type;

  template <unsigned int i>
  static void run(const CMap& amap,
                  typename CMap::Dart_const_handle adart,
                  std::vector<size_type>& marks, bool& ares)
  {
    CGAL_static_assertion_msg(CMap::Helper::template
                              Dimension_index<i>::value>=0,
                              "Test_is_valid_attribute_functor<i> but "
                              " i-attributes are disabled");

    size_type amark = marks[i];
    if ( amap.is_marked(adart, amark) ) return; // dart already test.

    bool valid = true;
    bool found_dart = false;

    typename CMap::template Attribute_const_handle<i>::type
        a=amap.template attribute<i>(adart);

    unsigned int nb = 0;
    for ( typename
            CMap::template Dart_of_cell_basic_const_range<i>::const_iterator
            it=amap.template darts_of_cell_basic<i>(adart, amark).begin();
          it.cont(); ++it )
    {
      if ( amap.template attribute<i>(it) != a )
      {
        std::cout<<"ERROR: an attribute of the "<<i<<"-cell is different. cur:";
        amap.template display_attribute<i>(a);
        std::cout<<" != first:";
        amap.template display_attribute<i>(amap.template attribute<i>(it));
        std::cout<<" for dart ";
        amap.display_dart(it);
        std::cout<<std::endl;
        valid=false;
      }

      if ( a!=amap.null_handle && it==amap.template dart_of_attribute<i>(a) )
        found_dart=true;

      amap.mark(it, amark);
      ++nb;
    }

    if ( a!=amap.null_handle )
    {
      if ( amap.template get_attribute_ref_counting<i>(a)!=nb )
      {
        std::cout<<"ERROR: the number of reference of an "<<i
                <<"-attribute is not correct. Count: "<<nb
               <<" != Store in the attribute: "
              <<amap.template get_attribute_ref_counting<i>(a)
             <<" for dart ";
        amap.display_dart(adart); std::cout<<std::endl;
        valid=false;
      }
      if ( !amap.template is_valid_attribute<i>(a) )
      {
        std::cout<<"ERROR: the dart associated with an "<<i
                <<"-attribute is NULL for dart ";
        amap.display_dart(adart); std::cout<<std::endl;
        valid=false;
      }
      if ( amap.template dart_of_attribute<i>(a)!=amap.null_handle &&
           !found_dart )
      {
        std::cout<<"ERROR: the non NULL dart of an "<<i
                <<"-attribute does not belong to the cell."<<std::endl;
        valid=false;
      }
    }

    if ( !valid ) ares=false;
  }
};
// ****************************************************************************
/// Functor used to correct invalid attributes in an i-cell
template<typename CMap>
struct Correct_invalid_attributes_functor
{
  typedef typename CMap::size_type size_type;

  template <unsigned int i>
  static void run(CMap& amap,
                  typename CMap::Dart_handle adart,
                  std::vector<size_type>& marks)
  {
    CGAL_static_assertion_msg(CMap::Helper::template
                              Dimension_index<i>::value>=0,
                              "Correct_invalid_attributes_functor<i> but "
                              " i-attributes are disabled");

    size_type amark = marks[i];
    typename CMap::template Attribute_handle<i>::type
        a=amap.template attribute<i>(adart);

    // dart already test, or without i-attribute
    if ( amap.is_marked(adart, amark) ) return;
    if ( a==amap.null_handle) { amap.mark(adart, amark); return; }

    // We search if all the darts of the i-cell has the same i-attrib, and we count
    // the number of darts of the i-cell.
    unsigned int nb=0;
    bool found_dart = false;

    for ( CGAL::CMap_dart_iterator_basic_of_cell<CMap,i>
            it(amap, adart, amark); it.cont(); ++it, ++nb )
    {
      if ( a!=amap.template attribute<i>(it) )
      {
        // If two different i-attributes, we could call on_split ?
        amap.template set_dart_attribute<i>(it, a);
      }
      if (it==amap.template dart_of_attribute<i>(a))
      {
        found_dart = true;
      }
      amap.mark(it, amark);
    }

    if (!found_dart)
    {
      // the current i-attrib does not belong to the i-cell
      // so we affect it to the first dart of the i-cell
      amap.template set_dart_of_attribute<i>(a,adart);
    }

    // If the i-cell has less darts than the ref counter of the i-attribute,
    // the i-attribute is shared by different cells => we duplicate it.
    if ( nb!=amap.template get_attribute_ref_counting<i>(a) )
    {
      typename CMap::template Attribute_handle<i>::type
        a2=amap.template copy_attribute<i>(a);
      amap.template set_attribute<i>(adart, a2);
    }
  }
};
// ****************************************************************************
/// Functor used to erase all attributes having no more associated dart
template<typename CMap>
struct Cleanup_useless_attributes
{
  template <unsigned int i>
  static void run(CMap& amap)
  {
    CGAL_static_assertion_msg(CMap::Helper::template
                              Dimension_index<i>::value>=0,
                              "Cleanup_useless_attributes<i> but "
                              " i-attributes are disabled");

    for ( typename CMap::template Attribute_range<i>::type::iterator
            it=amap.template attributes<i>().begin(),
            itend=amap.template attributes<i>().end(); it!=itend; ++it )
    {
      if ( amap.template get_attribute<i>(it).get_nb_refs()==0 )
      {
        amap.template erase_attribute<i>(it);
      }
    }
  }
};
// ****************************************************************************
/// Functor for counting i-cell
template<typename CMap>
struct Count_cell_functor
{
  typedef typename CMap::size_type size_type;

  template <unsigned int i>
  static void run( const CMap& amap,
                   typename CMap::Dart_const_handle adart,
                   std::vector<size_type>& amarks,
                   std::vector<unsigned int>& ares )
  {
    if ( amarks[i]!=CMap::INVALID_MARK && !amap.is_marked(adart, amarks[i]) )
    {
      ++ ares[i];
      CGAL::mark_cell<CMap,i>(amap, adart, amarks[i]);
    }
  }
};
// ****************************************************************************
/// Functor for counting the memory occupation of attributes
/// Be careful not reentrant !!! TODO a  Foreach_enabled_attributes
/// taking an instance of a functor as argument allowing to compute
/// and return values.
template<typename CMap>
struct Count_bytes_one_attribute_functor
{
  template <unsigned int i>
  static void run( const CMap& amap, typename CMap::size_type& res )
  {
    res += amap.template attributes<i>().capacity()*
      sizeof(typename CMap::template Attribute_type<i>::type);
  }
};
template<typename CMap>
struct Count_bytes_all_attributes_functor
{
  static typename CMap::size_type run( const CMap& amap )
  {
    typename CMap::size_type res = 0;
    CMap::Helper::template Foreach_enabled_attributes
      <CGAL::internal::Count_bytes_one_attribute_functor<CMap> >::run(amap, res);
    return res;
  }
};
// ****************************************************************************
/// Decrease the cell attribute reference counting of the given dart.
/// The attribute is removed if there is no more darts linked with it.
template<typename CMap, unsigned int i, typename T=
         typename CMap::template Attribute_type<i>::type>
struct Decrease_attribute_functor_run
{
  static void run(CMap& amap, typename CMap::Dart_handle adart)
  {
    if ( amap.template attribute<i>(adart)!=CMap::null_handle )
    {
      amap.template dec_attribute_ref_counting<i>(amap.template attribute<i>(adart));
      if ( amap.are_attributes_automatically_managed() &&
           amap.template get_attribute_ref_counting<i>
           (amap.template attribute<i>(adart))==0 )
        amap.template erase_attribute<i>(amap.template attribute<i>(adart));
    }
  }
};
/// Specialization for void attributes.
template<typename CMap, unsigned int i>
struct Decrease_attribute_functor_run<CMap, i, CGAL::Void>
{
  static void run(CMap&, typename CMap::Dart_handle)
  {}
};
// ****************************************************************************
/// Functor used to call decrease_attribute_ref_counting<i>
/// on each i-cell attribute enabled
template<typename CMap>
struct Decrease_attribute_functor
{
  template <unsigned int i>
  static void run(CMap& amap, typename CMap::Dart_handle adart)
  { CGAL::internal::Decrease_attribute_functor_run<CMap,i>::run(amap, adart); }
};
// ****************************************************************************
template<typename CMap, unsigned int i, typename T=
         typename CMap::template Attribute_type<i>::type>
struct Restricted_decrease_attribute_functor_run
{
  static void run(CMap& amap, typename CMap::Dart_handle adart)
  {
    if ( amap.template attribute<i>(adart)!=CMap::null_handle )
    {
      amap.template dec_attribute_ref_counting<i>(amap.template attribute<i>(adart));
    }
  }
};
/// Specialization for void attributes.
template<typename CMap, unsigned int i>
struct Restricted_decrease_attribute_functor_run<CMap, i, CGAL::Void>
{
  static void run(CMap&, typename CMap::Dart_handle)
  {}
};
// ****************************************************************************
/// Functor used to call restricted_decrease_attribute_ref_counting<i>
/// on each i-cell attribute enabled
template<typename CMap>
struct Restricted_decrease_attribute_functor
{
  template <unsigned int i>
  static void run(CMap& amap, typename CMap::Dart_handle adart)
  { CGAL::internal::Restricted_decrease_attribute_functor_run<CMap,i>::
        run(amap, adart); }
};
// ****************************************************************************
/// Functor used to initialize all attributes to NULL.
template<typename CMap>
struct Init_attribute_functor
{
  template <int i>
  static void run(CMap& amap, typename CMap::Dart_handle adart)
  { amap.template set_dart_attribute<i>(adart, CMap::null_handle); }
};
// ****************************************************************************
/// Functor used to set the i-attribute of a given dart.
template<typename CMap, unsigned int i, typename T=
         typename CMap::template Attribute_type<i>::type>
struct Set_i_attribute_of_dart_functor
{
  static void run( CMap& amap, typename CMap::Dart_handle dh,
                   typename CMap::template Attribute_handle<i>::type ah )
  {
    amap.template set_dart_attribute<i>(dh, ah);
  }
};
/// Specialization for void attributes.
template<typename CMap, unsigned int i>
struct Set_i_attribute_of_dart_functor<CMap, i, CGAL::Void>
{
  static void run( CMap&, typename CMap::Dart_handle,
                   typename CMap::template Attribute_handle<i>::type)
  {}
};
// ****************************************************************************
// Functor allowing to test if two info are the same or not
template<typename Map1, typename Map2, unsigned int i,
         typename Attr1Info1=typename Map1::template Attribute_type<i>::type::Info,
         typename Info2=typename Map2::template Attribute_type<i>::type::Info>
struct Is_same_info
{
  static bool run(const Map1&, const Map2&,
                  typename Map1::template Attribute_const_handle<i>::type,
                  typename Map2::template Attribute_const_handle<i>::type)
  { return false; }
};

template<typename Map1, typename Map2, unsigned int i, typename Info1>
struct Is_same_info<Map1, Map2, i, Info1, void>
{
  static bool run(const Map1&, const Map2&,
                  typename Map1::template Attribute_const_handle<i>::type,
                  typename Map2::template Attribute_const_handle<i>::type)
  { return false; }
};

template<typename Map1, typename Map2, unsigned int i, typename Info2>
struct Is_same_info<Map1, Map2, i, void, Info2>
{
  static bool run(const Map1&, const Map2&,
                  typename Map1::template Attribute_const_handle<i>::type,
                  typename Map2::template Attribute_const_handle<i>::type)
  { return false; }
};

template<typename Map1, typename Map2, unsigned int i>
struct Is_same_info<Map1, Map2, i, void, void>
{
  static bool run(const Map1&, const Map2&,
                  typename Map1::template Attribute_const_handle<i>::type,
                  typename Map2::template Attribute_const_handle<i>::type)
  { return true; }
};

template<typename Map1, typename Map2, unsigned int i, typename Info>
struct Is_same_info<Map1, Map2, i, Info, Info>
{
  static bool run(const Map1& m1, const Map2& m2,
                  typename Map1::template Attribute_const_handle<i>::type a1,
                  typename Map2::template Attribute_const_handle<i>::type a2)
  { return m1.template info_of_attribute<i>(a1)==
      m2.template info_of_attribute<i>(a2); }
};

// Case of two non void type
template<typename Map1, typename Map2, unsigned int i,
         typename T1=typename Map1::template Attribute_type<i>::type,
         typename T2=typename Map2::template Attribute_type<i>::type>
struct Is_same_attribute_info_functor
{
  static bool run(const Map1& m1, const Map2& m2,
                  typename Map1::Dart_const_handle dh1,
                  typename Map2::Dart_const_handle dh2)
  {
    if (m1.template attribute<i>(dh1)==Map1::null_handle &&
        m2.template attribute<i>(dh2)==Map2::null_handle)
      return true;

    if (m1.template attribute<i>(dh1)==Map1::null_handle ||
        m2.template attribute<i>(dh2)==Map2::null_handle)
      return false;

    return
      Is_same_info<Map1, Map2, i>::run(m1, m2,
                                       m1.template attribute<i>(dh1),
                                       m2.template attribute<i>(dh2));
  }
};

// Case T1==void
template <typename Map1, typename Map2, unsigned int i, typename T2>
struct Is_same_attribute_info_functor<Map1, Map2, i, Void, T2>
{
  static bool run(const Map1&, const Map2& m2,
                  typename Map1::Dart_const_handle,
                  typename Map2::Dart_const_handle dh2)
  { return m2.template attribute<i>(dh2)==Map2::null_handle; }
};

// Case T2==void
template <typename Map1, typename Map2, unsigned int i, typename T1>
struct Is_same_attribute_info_functor<Map1, Map2, i, T1, Void>
{
  static bool run(const Map1& m1, const Map2&,
                  typename Map1::Dart_const_handle dh1,
                  typename Map2::Dart_const_handle)
  { return m1.template attribute<i>(dh1)==Map1::null_handle; }
};

// Case T1==T2==void
template <typename Map1, typename Map2, unsigned int i>
struct Is_same_attribute_info_functor<Map1, Map2, i, Void, Void>
{
  static bool run(const Map1&, const Map2&,
                  typename Map1::Dart_const_handle,
                  typename Map2::Dart_const_handle)
  { return true; }
};
// ****************************************************************************
// Functor allowing to test if two darts have the same info or not.
// Default case, T1!=T2
template<typename Map1, typename Map2,
         typename T1=typename Map1::Dart_info,
         typename T2=typename Map2::Dart_info>
struct Test_is_same_dart_info_functor
{
  static bool run(const Map1&, const Map2&,
                  typename Map1::Dart_const_handle,
                  typename Map2::Dart_const_handle)
  { return false; }
};

// Case T1==T2==CGAL::Void
template<typename Map1, typename Map2>
struct Test_is_same_dart_info_functor<Map1, Map2, CGAL::Void, CGAL::Void>
{
  static bool run(const Map1&, const Map2&,
                  typename Map1::Dart_const_handle,
                  typename Map2::Dart_const_handle)
  { return true; }
};

// Case T1==T2!=CGAL::Void
template<typename Map1, typename Map2, typename T1>
struct Test_is_same_dart_info_functor<Map1, Map2, T1, T1>
{
  static bool run(const Map1& map1, const Map2& map2,
                  typename Map1::Dart_const_handle dh1,
                  typename Map2::Dart_const_handle dh2)
  { return (map1.info(dh1)==map2.info(dh2)); }
};

// ****************************************************************************
// Functor allowing to test if two points are the same or not.
// Here we know both attributes have points.
template< typename Attr1, typename Attr2,
          typename Point1=typename Attr1::Point,
          typename Point2=typename Attr2::Point,
          typename T1=typename Ambient_dimension<Point1>::type >
struct Is_same_point
{
  static bool run(const Attr1&, const Attr2&)
  { return false; }
};

template< typename Attr1, typename Attr2, typename Point>
struct Is_same_point<Attr1, Attr2, Point, Point, Dimension_tag<2> >
{
  static bool run(const Attr1& a1, const Attr2& a2)
  { return typename Kernel_traits<Point>::Kernel::Equal_2()
        (a1.point(),a2.point()); }
};

template< typename Attr1, typename Attr2, typename Point>
struct Is_same_point<Attr1, Attr2, Point, Point, Dimension_tag<3> >
{
  static bool run(const Attr1& a1, const Attr2& a2)
  { return typename Kernel_traits<Point>::Kernel::Equal_3()
        (a1.point(),a2.point()); }
};

template< typename Attr1, typename Attr2, typename Point>
struct Is_same_point<Attr1, Attr2, Point, Point, Dynamic_dimension_tag >
{
  static bool run(const Attr1& a1, const Attr2& a2)
  { return typename Kernel_traits<Point>::Kernel::Equal_d()
        (a1.point(),a2.point()); }
};

// Case of two non void type, with two points
template<typename Map1, typename Map2, int i,
         bool Withpoint1=Is_attribute_has_point
         <typename Map1::template Attribute_type<i>::type>::value,
         bool Withpoint2=Is_attribute_has_point
         <typename Map2::template Attribute_type<i>::type>::value>
struct Test_is_same_attribute_point_functor
{
  static bool run(const Map1& m1, const Map2& m2,
                  typename Map1::Dart_const_handle dh1,
                  typename Map2::Dart_const_handle dh2)
  {
    CGAL_static_assertion( Withpoint1==true && Withpoint2==true );
    if (m1.template attribute<i>(dh1)==Map1::null_handle &&
        m2.template attribute<i>(dh2)==Map2::null_handle)
      return true;

    if (m1.template attribute<i>(dh1)==Map1::null_handle ||
        m2.template attribute<i>(dh2)==Map2::null_handle)
      return false;

    return
      Is_same_point<typename Map1::template Attribute_type<i>::type,
                    typename Map2::template Attribute_type<i>::type>::run
      (m1.template get_attribute<i>(m1.template attribute<i>(dh1)),
       m2.template get_attribute<i>(m2.template attribute<i>(dh2)));
  }
};

// Case of two non void type, first without point
template<typename Map1, typename Map2, int i>
struct Test_is_same_attribute_point_functor<Map1, Map2, i, false, true>
{
  static bool run(const Map1&, const Map2&,
                  typename Map1::Dart_const_handle,
                  typename Map2::Dart_const_handle)
  { return false; }
};

// Case of two non void type, second without point
template<typename Map1, typename Map2, int i>
struct Test_is_same_attribute_point_functor<Map1, Map2, i, true, false>
{
  static bool run(const Map1&, const Map2&,
                  typename Map1::Dart_const_handle,
                  typename Map2::Dart_const_handle)
  { return false; }
};

// Case of two non void type, both without point
template<typename Map1, typename Map2, int i>
struct Test_is_same_attribute_point_functor<Map1, Map2, i, false, false>
{
  static bool run(const Map1&, const Map2&,
                  typename Map1::Dart_const_handle,
                  typename Map2::Dart_const_handle)
  { return true; }
};
// ****************************************************************************
/// Test if the two darts are associated with the same attribute.
template<typename Map1, typename Map2>
struct Test_is_same_attribute_functor
{
  template<unsigned int i>
  static void run(const Map1& m1, const Map2& m2,
                  typename Map1::Dart_const_handle dh1,
                  typename Map2::Dart_const_handle dh2,
                  bool & res)
  {
    if (res)
    {
      res = Is_same_attribute_info_functor
        <Map1, Map2, i>::run(m1, m2, dh1, dh2);
    }
  }
};
// ****************************************************************************
/// Functor to reverse the orientation of a combinatorial map,
/// when 0-attributes are non void.
template <typename CMap, typename Attrib =
          typename CMap::Helper::template Attribute_type<0>::type>
struct Reverse_orientation_of_map_functor
{
  static void run(CMap& amap)
  {
    typename CMap::size_type mymark = amap.get_new_mark();
    CGAL_precondition(amap.is_whole_map_unmarked(mymark));
    CGAL_precondition(amap.is_valid());

    typename CMap::Dart_handle first=NULL, current=NULL, prev=NULL, next=NULL;
    typename CMap::Helper::template Attribute_handle<0>::type
      first_attribute=NULL, next_attribute=NULL;

    for (typename CMap::Dart_range::iterator current_dart=amap.darts().begin(),
           last_dart = amap.darts().end(); current_dart!=last_dart;
         ++current_dart)
    {
      if (amap.is_marked(current_dart, mymark)) continue;

      first=current_dart;
      current=amap.beta(current_dart, 1);
      first_attribute=amap.template attribute<0>(current);
      amap.template inc_attribute_ref_counting<0>(first_attribute);

      do
      {
        amap.mark(current, mymark);
        prev=amap.template beta<0>(current);
        next=amap.template beta<1>(current);
        next_attribute=amap.template attribute<0>(next);
        CGAL::internal::Set_i_attribute_of_dart_functor<CMap, 0>::
          run(amap, current, next_attribute);
        amap.template dart_link_beta<1>(current, prev);
        amap.template dart_link_beta<0>(current, next);
        current = amap.beta(current,0); // the old beta1
      }
      while (current!=first);

      amap.mark(current, mymark);
      prev=amap.template beta<0>(current);
      next=amap.template beta<1>(current);
      CGAL::internal::Set_i_attribute_of_dart_functor<CMap, 0>::
        run(amap, current, first_attribute);
      amap.template dart_link_beta<1>(current, prev);
      amap.template dart_link_beta<0>(current, next);

      amap.template dec_attribute_ref_counting<0>(first_attribute);
    }
    CGAL_postcondition(amap.is_whole_map_marked(mymark));
    CGAL_postcondition(amap.is_valid());
    amap.free_mark(mymark);
  }
};
// ****************************************************************************
// Specialization for void 0-attributes
template <typename CMap>
struct Reverse_orientation_of_map_functor<CMap, CGAL::Void>
{
  static void run(CMap& amap)
  {
    typename CMap::size_type mymark = amap.get_new_mark();
    CGAL_precondition(amap.is_whole_map_unmarked(mymark));
    CGAL_precondition(amap.is_valid());

    for (typename CMap::Dart_range::iterator current_dart=amap.darts().begin(),
         last_dart = amap.darts().end();
         current_dart!=last_dart; ++current_dart)
    {
      if (amap.is_marked(current_dart, mymark)) continue;

      if (amap.is_marked(current_dart, mymark)) continue;
      typename CMap::Dart_handle first_dart_in_cell= current_dart;
      typename CMap::Dart_handle current_dart_in_cell=current_dart;
      do
      {
        amap.mark(current_dart_in_cell, mymark);
        typename CMap::Dart_handle previous_dart_in_cell=
          amap.template beta<0>(current_dart_in_cell);
        typename CMap::Dart_handle next_dart_in_cell=
          amap.template beta<1>(current_dart_in_cell);
        amap.template dart_link_beta<1>(current_dart_in_cell, previous_dart_in_cell);
        amap.template dart_link_beta<0>(current_dart_in_cell, next_dart_in_cell);
        current_dart_in_cell = amap.beta(current_dart_in_cell,0);
      }
      while (current_dart_in_cell!=first_dart_in_cell);
    }
    CGAL_postcondition(amap.is_whole_map_marked(mymark));
    CGAL_postcondition(amap.is_valid());
    amap.free_mark(mymark);
  }
};
// ****************************************************************************
/// Functor to reverse the orientation of a connected component in a given map
template <typename CMap, typename Attrib=
          typename CMap::Helper::template Attribute_type<0>::type>
struct Reverse_orientation_of_connected_component_functor
{
  static void run(CMap& amap, typename CMap::Dart_handle adart,
                  typename CMap::size_type amark=CMap::INVALID_MARK)
  {
    typename CMap::size_type mymark =
      (amark==CMap::INVALID_MARK?amap.get_new_mark():amark);

    typename CMap::Dart_handle first=NULL, current=NULL, prev=NULL, next=NULL;
    typename CMap::Helper::template Attribute_handle<0>::type
      first_attribute=NULL, next_attribute=NULL;

    for (typename CMap::template Dart_of_cell_range<CMap::dimension+1>::iterator
           current_dart=amap.template darts_of_cell<CMap::dimension+1>
           (adart).begin(),
           last_dart=amap.template darts_of_cell<CMap::dimension+1>
           (adart).end();
         current_dart!=last_dart; ++current_dart)
    {
      if (amap.is_marked(current_dart, mymark)) continue;

      first=current_dart;
      current=amap.beta(current_dart, 1);
      first_attribute=amap.template attribute<0>(current);
      amap.template inc_attribute_ref_counting<0>(first_attribute);

      do
      {
        amap.mark(current, mymark);
        prev=amap.template beta<0>(current);
        next=amap.template beta<1>(current);
        next_attribute=amap.template attribute<0>(next);
        CGAL::internal::Set_i_attribute_of_dart_functor<CMap, 0>::
          run(amap, current, next_attribute);
        amap.template dart_link_beta<1>(current, prev);
        amap.template dart_link_beta<0>(current, next);
        current = amap.beta(current,0); // the old beta1
      }
      while (current!=first);

      amap.mark(current, mymark);
      prev=amap.template beta<0>(current);
      next=amap.template beta<1>(current);
      CGAL::internal::Set_i_attribute_of_dart_functor<CMap, 0>::
        run(amap, current, first_attribute);
      amap.template dart_link_beta<1>(current, prev);
      amap.template dart_link_beta<0>(current, next);

      amap.template dec_attribute_ref_counting<0>(first_attribute);
    }

    if (amark==CMap::INVALID_MARK)
    {
      CGAL::unmark_cell<CMap, CMap::dimension+1>(amap, adart, mymark);
      CGAL_postcondition(amap.is_whole_map_unmarked(mymark));
      amap.free_mark(mymark);
    }
  }
};
// ****************************************************************************
// Specialization for void 0-attributes
template <typename CMap>
struct Reverse_orientation_of_connected_component_functor<CMap, CGAL::Void>
{
  static void run(CMap& amap, typename CMap::Dart_handle adart,
                  typename CMap::size_type amark=CMap::INVALID_MARK)
  {
    typename CMap::size_type mymark =
      (amark==CMap::INVALID_MARK?amap.get_new_mark():amark);

    for (typename CMap::template Dart_of_cell_range<CMap::dimension+1>::iterator
           current_dart=amap.template darts_of_cell<CMap::dimension+1>
           (adart).begin(),
           last_dart=amap.template darts_of_cell_basic<CMap::dimension+1>
           (adart).end();
         current_dart!=last_dart; ++current_dart)
    {
      if (amap.is_marked(current_dart, mymark)) continue;

      typename CMap::Dart_handle first_dart_in_cell=current_dart;
      typename CMap::Dart_handle current_dart_in_cell=current_dart;

      do
      {
        amap.mark(current_dart_in_cell, mymark);
        typename CMap::Dart_handle previous_dart_in_cell=
          amap.template beta<0>(current_dart_in_cell);
        typename CMap::Dart_handle next_dart_in_cell=
          amap.template beta<1>(current_dart_in_cell);

        amap.template dart_link_beta<1>(current_dart_in_cell, previous_dart_in_cell);
        amap.template dart_link_beta<0>(current_dart_in_cell, next_dart_in_cell);
        current_dart_in_cell = amap.beta(current_dart_in_cell,0);
      }
      while (current_dart_in_cell != first_dart_in_cell);
    }

    if (amark==CMap::INVALID_MARK)
    {
      CGAL::unmark_cell<CMap::dimension+1>(adart, mymark);
      CGAL_postcondition(amap.is_whole_map_unmarked(mymark));
      amap.free_mark(mymark);
    }
  }
};
// ****************************************************************************
// Beta functor, used to combine several beta.
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
template<typename CMap, typename Dart_handle, typename ... Betas>
struct Beta_functor;

template<typename CMap, typename Dart_handle>
struct Beta_functor<CMap, Dart_handle, int>
{
  static Dart_handle run(CMap& AMap, Dart_handle ADart, int B)
  { return AMap.get_beta(ADart, B); }
};

template<typename CMap, typename Dart_handle>
struct Beta_functor<CMap, Dart_handle, unsigned int>
{
  static Dart_handle run(CMap& AMap, Dart_handle ADart, unsigned int B)
  { return  AMap.get_beta(ADart, B); }
};

template<typename CMap, typename Dart_handle, typename ... Betas>
struct Beta_functor<CMap, Dart_handle, int, Betas...>
{
  static Dart_handle run(CMap& AMap, Dart_handle ADart, int B, Betas... betas)
  { return Beta_functor<CMap, Dart_handle, Betas...>::
      run(AMap, AMap.get_beta(ADart, B), betas...); }
};

template<typename CMap, typename Dart_handle, typename ... Betas>
struct Beta_functor<CMap, Dart_handle, unsigned int, Betas...>
{
  static Dart_handle run(CMap& AMap, Dart_handle ADart, unsigned int B,
                         Betas... betas)
  { return Beta_functor<CMap, Dart_handle, Betas...>::
      run(AMap, AMap.get_beta(ADart, B), betas...); }
};
// ****************************************************************************
template<typename CMap, typename Dart_handle, int ... Betas>
struct Beta_functor_static;

template<typename CMap, typename Dart_handle, int B>
struct Beta_functor_static<CMap, Dart_handle, B>
{
  static Dart_handle run(CMap& AMap, Dart_handle ADart)
  { return AMap.template get_beta<B>(ADart); }
};

template<typename CMap, typename Dart_handle, int B, int ... Betas>
struct Beta_functor_static<CMap, Dart_handle, B, Betas...>
{
  static Dart_handle run(CMap& AMap, Dart_handle ADart)
  { return Beta_functor_static<CMap, Dart_handle, Betas...>::
        run(AMap, AMap.template get_beta<B>(ADart)); }
};
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
// ****************************************************************************
template<typename Container, class WitdId=
         typename Container::value_type::Has_id>
struct Init_id
{
  static void run(Container& c, typename Container::iterator e)
  {
    e->set_id(c.index(e));
  }
};
template<typename Container>
struct Init_id<Container, Tag_false>
{
  static void run(Container&, typename Container::iterator)
  {}
};
// ****************************************************************************
} // namespace internal
} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_INTERNAL_FUNCTORS_H
