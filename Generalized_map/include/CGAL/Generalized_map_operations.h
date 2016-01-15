// Copyright (c) 2014 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_GENERALIZED_MAP_OPERATIONS_H
#define CGAL_GENERALIZED_MAP_OPERATIONS_H 1

#include <CGAL/Generalized_map.h>
#include <CGAL/Generalized_map_insertions.h>
#include <deque>
#include <stack>

namespace CGAL
{
  /** @file Generalized_map_operations.h
   * Some operations to modify a generalized map.
   */

  /** Test if an i-cell can be removed.
   *  An i-cell can be removed if i==CMap::dimension or i==CMap::dimension-1,
   *     or if there are at most two (i+1)-cell incident to it.
   * @param adart a dart of the i-cell.
   * @return true iff the i-cell can be removed.
   */
  template <class CMap, unsigned int i, unsigned int nmi=CMap::dimension-i>
  struct Is_removable_functor
  {
    static bool run(const CMap& amap, typename CMap::Dart_const_handle adart)
    {
      // TODO? Optimisation for dim-2, and to not test all the darts of the cell?
      bool res = true;
      for ( CGAL::GMap_dart_const_iterator_of_cell<CMap,i> it(amap, adart);
            res && it.cont(); ++it )
      {
        if (amap.template alpha<i+2,i+1>(it)!=amap.template alpha<i+1,i+2>(it))
          res = false;
      }
      return res;
    }
  };
  // Specialization for i=CMap::dimension
  template <class CMap, unsigned int i>
  struct Is_removable_functor<CMap, i, 0>
  {
    static bool run(const CMap&, typename CMap::Dart_const_handle)
    { return true; }
  };
  // Specialization for i=CMap::dimension-1
  template <class CMap, unsigned int i>
  struct Is_removable_functor<CMap, i, 1>
  {
    static bool run(const CMap&, typename CMap::Dart_const_handle)
    { return true; }
  };
  /** Test if an i-cell can be removed.
   *  An i-cell can be removed if i==CMap::dimension or i==CMap::dimension-1,
   *     or if there are at most two (i+1)-cell incident to it.
   * @param adart a dart of the i-cell.
   * @return true iff the i-cell can be removed.
   */
  template < class CMap, unsigned int i >
  bool is_removable(const CMap& amap, typename CMap::Dart_const_handle adart)
  { return CGAL::Is_removable_functor<CMap, i>::run(amap,adart); }

  /** Remove an i-cell, 0<=i<dimension, and merge eventually both incident
   *  (i+1)-cells.
   *  @param amap the used generalized map.
   *  @param adart a dart of the i-cell to remove.
   *  @return the number of deleted darts.
   */
  template<class CMap, unsigned int i, unsigned int nmi>
  struct Remove_cell_functor
  {
    static size_t run(CMap& amap, typename CMap::Dart_handle adart)
    {
      CGAL_static_assertion ( i<CMap::dimension );
      CGAL_assertion( (is_removable<CMap,i>(amap, adart)) );

      size_t res = 0;

      typename CMap::Dart_handle d1, d2;
      typename CMap::Dart_handle dg1=amap.null_handle, dg2=amap.null_handle;

      int mark = amap.get_new_mark();
      int mark_modified_darts = amap.get_new_mark();

      std::deque<typename CMap::Dart_handle> to_erase;

      // First we store and mark all the darts of the i-cell to remove.
      for ( CGAL::GMap_dart_iterator_basic_of_cell<CMap,i> it(amap,adart,mark);
            it.cont(); ++it )
      {
        to_erase.push_back(it);
        if ( !amap.template is_free<i+1>(it) && dg1==amap.null_handle )
        { dg1=it; dg2=amap.template alpha<i+1>(it); }
        amap.mark(it, mark);
        ++res;
      }

      // We group the two (i+1)-cells incident if they exist.
      if ( dg1!=amap.null_handle )
        CGAL::internal::Group_attribute_functor_run<CMap, i+1>::
            run(&amap, dg1, dg2);

      // During the operation, we store in modified_darts the darts modified
      // to test after the loop the non void attributes that are split.
      std::deque<typename CMap::Dart_handle> modified_darts;

      // For each dart of the i-cell, we modify i-links of neighbors.
      typename std::deque<typename CMap::Dart_handle>::iterator it =
          to_erase.begin();
      for ( ; it!=to_erase.end(); ++it )
      {
        d1=amap.template alpha<i>(*it);

        if ( !amap.is_marked(d1, mark) )
        {
          d2=amap.template alpha<i+1,i>(*it);
          while ( amap.is_marked(d2, mark) )
          {
            d2=amap.template alpha<i+1,i>(d2);
          }

          if ( !amap.is_marked(d1, mark_modified_darts) )
          {
            CGAL_assertion( !amap.is_marked(d2, mark_modified_darts) );
            amap.template basic_link_alpha<i>(d1, d2);
            amap.mark(d1, mark_modified_darts);
            modified_darts.push_back(d1);
            // TODO push only one out of two dart ?
            if ( d2!=d1 )
            {
              modified_darts.push_back(d2);
              amap.mark(d2, mark_modified_darts);
            }
          }
        }
      }

      // We test the split of all the incident cells for all the non
      // void attributes.
      CMap::Helper::template Foreach_enabled_attributes_except
        <CGAL::internal::Test_split_attribute_functor<CMap,i>, i>::
        run(&amap, modified_darts, mark_modified_darts);

      // We remove all the darts of the i-cell.
      for ( it=to_erase.begin(); it!=to_erase.end(); ++it )
      { amap.erase_dart(*it); }

      CGAL_assertion( amap.is_whole_map_unmarked(mark) );
      amap.free_mark(mark);

      // If no attribute is enabled (or if only i-attributes are enabled),
      // the darts are not unmark by Foreach_enabled_attributes_except.
      // Thus we unmark them now.
      if ( !amap.is_whole_map_unmarked(mark_modified_darts) )
      {
        for ( it=modified_darts.begin();
              it!=modified_darts.end(); ++it )
          amap.unmark(*it, mark_modified_darts);
      }

      CGAL_assertion ( amap.is_whole_map_unmarked(mark_modified_darts) );
      amap.free_mark(mark_modified_darts);

#ifdef CGAL_GMAP_TEST_VALID_REMOVALS
      CGAL_assertion( amap.is_valid() );
#endif

      return res;
    }
  };

  /** Remove a d-cell, in a d-map (special case).
   *  @param amap the used combinatorial map.
   *  @param adart a dart of the volume to remove.
   *  @return the number of deleted darts.
   */
  template<class Gmap,unsigned int i>
  struct Remove_cell_functor<Gmap,i,0>
  {
    static size_t run(Gmap& amap, typename Gmap::Dart_handle adart)
    {
      int mark = amap.get_new_mark();
      std::deque<typename Gmap::Dart_handle> to_erase;
      size_t res = 0;

      std::deque<typename Gmap::Dart_handle> modified_darts;

      // We mark all the darts of the d-cell.
      for ( CGAL::GMap_dart_iterator_basic_of_cell<Gmap,Gmap::dimension>
            it(amap,adart,mark); it.cont(); ++it )
      {
        to_erase.push_back(it);
        amap.mark(it,mark);
        ++res;
      }

      // We unlink all the darts of the volume for alpha-d.
      typename std::deque<typename Gmap::Dart_handle>::iterator
        it = to_erase.begin();
      for ( it = to_erase.begin(); it != to_erase.end(); ++it )
      {
        if ( !amap.template is_free<Gmap::dimension>(*it) &&
             !amap.is_marked(amap.template alpha<Gmap::dimension>(*it), mark) )
        {
          modified_darts.push_back(amap.template alpha<Gmap::dimension>(*it));
          amap.template unlink_alpha<Gmap::dimension>(*it);
        }
      }

      // We test the split of all the incident cells for all the non
      // void attributes.
      Gmap::Helper::template Foreach_enabled_attributes_except
          <CGAL::internal::Test_split_attribute_functor<Gmap,i>,
          Gmap::dimension>::run(&amap, modified_darts);

      // We remove all the darts of the d-cell.
      for ( it = to_erase.begin(); it != to_erase.end(); ++it )
      { amap.erase_dart(*it); }

      CGAL_assertion( amap.is_whole_map_unmarked(mark) );
      amap.free_mark(mark);

#ifdef CGAL_GMAP_TEST_VALID_REMOVALS
      CGAL_assertion( amap.is_valid() );
#endif

      return res;
    }
  };

  /** Remove an i-cell, 0<=i<=dimension.
   * @param amap the used generalized map.
   * @param adart a dart of the i-cell to remove.
   * @return the number of deleted darts.
   */
  template< class Gmap, unsigned int i >
  size_t remove_cell(Gmap& amap, typename Gmap::Dart_handle adart)
  {
    return
        CGAL::Remove_cell_functor<Gmap,i,Gmap::dimension-i>::run(amap,adart);
  }

  /** Test if an i-cell can be contracted.
   *  An i-cell can be contracted if i==1
   *     or if there are at most two (i-1)-cell incident to it.
   * @param adart a dart of the i-cell.
   * @return true iff the i-cell can be contracted.
   */
  template <class Gmap, unsigned int i>
  struct Is_contractible_functor
  {
    static bool run(const Gmap& amap, typename Gmap::Dart_const_handle adart)
    {
      // TODO ? Optimisation possible to not test all the darts of the cell ?
      bool res = true;
      for ( CGAL::GMap_dart_const_iterator_of_cell<Gmap,i> it(amap, adart);
            res && it.cont(); ++it )
      {
        if (amap.template alpha<i-2,i-1>(it)!=amap.template alpha<i-1,i-2>(it))
          res = false;
      }
      return res;
    }
  };
  // Specialization for i=0
  template <class Gmap>
  struct Is_contractible_functor<Gmap, 0>
  {
    static bool run(const Gmap&, typename Gmap::Dart_const_handle)
    { return false; }
  };
  // Specialization for i=1
  template <class Gmap>
  struct Is_contractible_functor<Gmap, 1>
  {
    static bool run(const Gmap&, typename Gmap::Dart_const_handle)
    { return true; }
  };
  /** Test if an i-cell can be contracted.
   *  An i-cell can be contracted if i==1
   *     or if there are at most two (i-1)-cell incident to it.
   * @param adart a dart of the i-cell.
   * @return true iff the i-cell can be contracted.
   */
  template < class Gmap, unsigned int i >
  bool is_contractible(const Gmap& amap, typename Gmap::Dart_const_handle adart)
  { return CGAL::Is_contractible_functor<Gmap, i>::run(amap,adart); }

  /** Contract an i-cell, 1<=i<=dimension, and merge eventually both incident
   *  (i-1)-cells.
   * @param amap the used generalized map.
   * @param adart a dart of the i-cell to contract.
   * @return the number of deleted darts.
   */
  template<class Gmap, unsigned int i>
  struct Contract_cell_functor
  {
    static size_t run(Gmap& amap, typename Gmap::Dart_handle adart)
    {
      CGAL_static_assertion ( 1<=i && i<=Gmap::dimension );
      CGAL_assertion( (is_contractible<Gmap,i>(amap, adart)) );

      size_t res = 0;

      typename Gmap::Dart_handle d1, d2;
      typename Gmap::Dart_handle dg1=amap.null_handle, dg2=amap.null_handle;

      int mark = amap.get_new_mark();
      int mark_modified_darts = amap.get_new_mark();

      // First we store and mark all the darts of the i-cell to contract.
      std::deque<typename Gmap::Dart_handle> to_erase;
      for ( CGAL::GMap_dart_iterator_basic_of_cell<Gmap,i> it(amap,adart,mark);
            it.cont(); ++it )
      {
        to_erase.push_back(it);
        if ( !amap.template is_free<i-1>(it) && dg1==amap.null_handle )
        { dg1=it; dg2=amap.template alpha<i-1>(it); }
        amap.mark(it, mark);
        ++res;
      }

      // We group the two (i-1)-cells incident if they exist.
      if ( dg1!=amap.null_handle )
         CGAL::internal::Group_attribute_functor_run<Gmap,i-1>::
             run(&amap, dg1, dg2);

      // During the operation, we store in modified_darts the darts modified
      // to test after the loop the non void attributes that are split.
      std::deque<typename Gmap::Dart_handle> modified_darts;

      // For each dart of the i-cell, we modify i-links of neighbors.
      typename std::deque<typename Gmap::Dart_handle>::iterator it =
          to_erase.begin();
      for ( ; it!=to_erase.end(); ++it )
      {
        d1 = amap.template alpha<i>(*it);
        if ( !amap.is_marked(d1, mark) )
        {
          d2 = amap.template alpha<i-1,i>(*it);
          while ( amap.is_marked(d2, mark) )
          { d2 = amap.template alpha<i-1,i>(d2); }

          if ( !amap.is_marked(d1, mark_modified_darts) )
          {
            CGAL_assertion( !amap.is_marked(d2, mark_modified_darts) );
            amap.template basic_link_alpha<i>(d1, d2);
            amap.mark(d1, mark_modified_darts);
            modified_darts.push_back(d1);
            // TODO push only one out of two dart ?
            if ( d1!=d2 )
            {
              amap.mark(d2, mark_modified_darts);
              modified_darts.push_back(d2);
            }
          }
        }
      }

      // We test the split of all the incident cells for all the non
      // void attributes.
      Gmap::Helper::template Foreach_enabled_attributes_except
        <CGAL::internal::Test_split_attribute_functor<Gmap,i>, i>::
        run(&amap, modified_darts, mark_modified_darts);

      // We remove all the darts of the i-cell.
      for ( it=to_erase.begin(); it!=to_erase.end(); ++it )
      { amap.erase_dart(*it); }

      CGAL_assertion( amap.is_whole_map_unmarked(mark) );
      amap.free_mark(mark);

      // If no attribute is enabled (or if only i-attributes are enabled),
      // the darts are not unmark by Foreach_enabled_attributes_except.
      // Thus we unmark them now.
      if ( !amap.is_whole_map_unmarked(mark_modified_darts) )
      {
        for ( it=modified_darts.begin(); it!=modified_darts.end(); ++it )
          amap.unmark(*it, mark_modified_darts);
      }

      CGAL_assertion ( amap.is_whole_map_unmarked(mark_modified_darts) );
      amap.free_mark(mark_modified_darts);

#ifdef CGAL_GMAP_TEST_VALID_CONTRACTIONS
      CGAL_assertion( amap.is_valid() );
#endif

      return res;
    }
  };

  /** Contract an i-cell, 1<=i<=dimension.
   * @param amap the used generalized map.
   * @param adart a dart of the i-cell to remove.
   * @return the number of deleted darts.
   */
  template < class Gmap, unsigned int i >
  size_t contract_cell(Gmap& amap, typename Gmap::Dart_handle adart)
  { return CGAL::Contract_cell_functor<Gmap,i>::run(amap,adart); }

} // namespace CGAL

#endif // CGAL_GENERALIZED_MAP_OPERATIONS_H //
// EOF //
