// Copyright (c) 2016 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_GENERALIZED_MAP_OPERATIONS_H
#define CGAL_GENERALIZED_MAP_OPERATIONS_H 1

#include <CGAL/GMap_dart_const_iterators.h>
#include <CGAL/internal/Combinatorial_map_group_functors.h>
#include <CGAL/Combinatorial_map_basic_operations.h>

#include <deque>

namespace CGAL
{
  /** @file Generalized_map_operations.h
   * Some operations to modify a generalized map.
   */

  /** Test if an i-cell can be removed.
   *  An i-cell can be removed if i==GMap::dimension or i==GMap::dimension-1,
   *     or if there are at most two (i+1)-cell incident to it.
   * @param adart a dart of the i-cell.
   * @return true iff the i-cell can be removed.
   */
  template <class GMap, unsigned int i, unsigned int nmi=GMap::dimension-i>
  struct Is_removable_functor_gmap
  {
    static bool run(const GMap& amap, typename GMap::Dart_const_handle adart)
    {
      // TODO? Optimisation for dim-2, and to not test all the darts of the cell?
      bool res = true;
      for ( CGAL::GMap_dart_const_iterator_of_cell<GMap,i> it(amap, adart);
            res && it.cont(); ++it )
      {
        if (amap.template alpha<i+2,i+1>(it)!=amap.template alpha<i+1,i+2>(it))
          res = false;
      }
      return res;
    }
  };
  // Specialization for i=GMap::dimension
  template <class GMap, unsigned int i>
  struct Is_removable_functor_gmap<GMap, i, 0>
  {
    static bool run(const GMap&, typename GMap::Dart_const_handle)
    { return true; }
  };
  // Specialization for i=GMap::dimension-1
  template <class GMap, unsigned int i>
  struct Is_removable_functor_gmap<GMap, i, 1>
  {
    static bool run(const GMap&, typename GMap::Dart_const_handle)
    { return true; }
  };

  /** Remove an i-cell, 0<=i<dimension, and merge eventually both incident
   *  (i+1)-cells.
   *  @param amap the used generalized map.
   *  @param adart a dart of the i-cell to remove.
   *  @return the number of deleted darts.
   */
  template<class GMap, unsigned int i, unsigned int nmi>
  struct Remove_cell_functor_gmap
  {
    static size_t run(GMap& amap, typename GMap::Dart_handle adart,
                      bool update_attributes)
    {
      CGAL_static_assertion ( i<GMap::dimension );
      CGAL_assertion( (amap.template is_removable<i>(adart)) );

      size_t res = 0;

      typename GMap::Dart_handle d1, d2;
      typename GMap::Dart_handle dg1=amap.null_handle, dg2=amap.null_handle;

      typename GMap::size_type mark = amap.get_new_mark();
      typename GMap::size_type mark_modified_darts = amap.get_new_mark();

      std::deque<typename GMap::Dart_handle> to_erase;

      // First we store and mark all the darts of the i-cell to remove.
      for ( CGAL::GMap_dart_iterator_basic_of_cell<GMap,i> it(amap,adart,mark);
            it.cont(); ++it )
      {
        to_erase.push_back(it);
        if ( !amap.template is_free<i+1>(it) && dg1==amap.null_handle )
        { dg1=it; dg2=amap.template alpha<i+1>(it); }
        amap.mark(it, mark);
        ++res;
      }

      if (amap.are_attributes_automatically_managed())
      {
        // We group the two (i+1)-cells incident if they exist.
        if ( dg1!=amap.null_handle )
          CGAL::internal::GMap_group_attribute_functor_run<GMap, i+1>::
              run(amap, dg1, dg2);
      }

      // During the operation, we store in modified_darts the darts modified
      // to test after the loop the non void attributes that are split.
      std::deque<typename GMap::Dart_handle> modified_darts;

      // For each dart of the i-cell, we modify i-links of neighbors.
      typename std::deque<typename GMap::Dart_handle>::iterator it =
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

      if (amap.are_attributes_automatically_managed() && update_attributes)
      {
        // We test the split of all the incident cells for all the non
        // void attributes.
        GMap::Helper::template Foreach_enabled_attributes_except
          <CGAL::internal::GMap_test_split_attribute_functor<GMap,i>, i>::
          run(amap, modified_darts, mark_modified_darts);
      }

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
   *  @param amap the used generalized map.
   *  @param adart a dart of the volume to remove.
   *  @return the number of deleted darts.
   */
  template<class GMap,unsigned int i>
  struct Remove_cell_functor_gmap<GMap,i,0>
  {
    static size_t run(GMap& amap, typename GMap::Dart_handle adart,
                      bool update_attributes)
    {
      typename GMap::size_type mark = amap.get_new_mark();
      std::deque<typename GMap::Dart_handle> to_erase;
      size_t res = 0;

      std::deque<typename GMap::Dart_handle> modified_darts;

      // We mark all the darts of the d-cell.
      for ( CGAL::GMap_dart_iterator_basic_of_cell<GMap,GMap::dimension>
            it(amap,adart,mark); it.cont(); ++it )
      {
        to_erase.push_back(it);
        amap.mark(it,mark);
        ++res;
      }

      // We unlink all the darts of the volume for alpha-d.
      typename std::deque<typename GMap::Dart_handle>::iterator
        it = to_erase.begin();
      for ( it = to_erase.begin(); it != to_erase.end(); ++it )
      {
        if ( !amap.template is_free<GMap::dimension>(*it) &&
             !amap.is_marked(amap.template alpha<GMap::dimension>(*it), mark) )
        {
          if (amap.are_attributes_automatically_managed())
          {
            modified_darts.push_back(amap.template alpha<GMap::dimension>(*it));
          }
          amap.template unlink_alpha<GMap::dimension>(*it);
        }
      }

      if (amap.are_attributes_automatically_managed() && update_attributes)
      {
        // We test the split of all the incident cells for all the non
        // void attributes.
        GMap::Helper::template Foreach_enabled_attributes_except
            <CGAL::internal::GMap_test_split_attribute_functor<GMap,i>,
            GMap::dimension>::run(amap, modified_darts);
      }

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

  /** Test if an i-cell can be contracted.
   *  An i-cell can be contracted if i==1
   *     or if there are at most two (i-1)-cell incident to it.
   * @param adart a dart of the i-cell.
   * @return true iff the i-cell can be contracted.
   */
  template <class GMap, unsigned int i>
  struct Is_contractible_functor_gmap
  {
    static bool run(const GMap& amap, typename GMap::Dart_const_handle adart)
    {
      // TODO ? Optimisation possible to not test all the darts of the cell ?
      bool res = true;
      for ( CGAL::GMap_dart_const_iterator_of_cell<GMap,i> it(amap, adart);
            res && it.cont(); ++it )
      {
        if (amap.template alpha<i-2,i-1>(it)!=amap.template alpha<i-1,i-2>(it))
          res = false;
      }
      return res;
    }
  };
  // Specialization for i=0
  template <class GMap>
  struct Is_contractible_functor_gmap<GMap, 0>
  {
    static bool run(const GMap&, typename GMap::Dart_const_handle)
    { return false; }
  };
  // Specialization for i=1
  template <class GMap>
  struct Is_contractible_functor_gmap<GMap, 1>
  {
    static bool run(const GMap&, typename GMap::Dart_const_handle)
    { return true; }
  };

  /** Contract an i-cell, 1<=i<=dimension, and merge eventually both incident
   *  (i-1)-cells.
   * @param amap the used generalized map.
   * @param adart a dart of the i-cell to contract.
   * @return the number of deleted darts.
   */
  template<class GMap, unsigned int i>
  struct Contract_cell_functor_gmap
  {
    static size_t run(GMap& amap, typename GMap::Dart_handle adart,
                      bool update_attributes)
    {
      CGAL_static_assertion ( 1<=i && i<=GMap::dimension );
      CGAL_assertion( (amap.template is_contractible<i>(adart)) );

      size_t res = 0;

      typename GMap::Dart_handle d1, d2;
      typename GMap::Dart_handle dg1=amap.null_handle, dg2=amap.null_handle;

      typename GMap::size_type mark = amap.get_new_mark();
      typename GMap::size_type mark_modified_darts = amap.get_new_mark();

      // First we store and mark all the darts of the i-cell to contract.
      std::deque<typename GMap::Dart_handle> to_erase;
      for ( CGAL::GMap_dart_iterator_basic_of_cell<GMap,i> it(amap,adart,mark);
            it.cont(); ++it )
      {
        to_erase.push_back(it);
        if ( !amap.template is_free<i-1>(it) && dg1==amap.null_handle )
        { dg1=it; dg2=amap.template alpha<i-1>(it); }
        amap.mark(it, mark);
        ++res;
      }

      if ( amap.are_attributes_automatically_managed() )
      {
         // We group the two (i-1)-cells incident if they exist.
        if ( dg1!=amap.null_handle )
           CGAL::internal::GMap_group_attribute_functor_run<GMap,i-1>::
               run(amap, dg1, dg2);
      }

      // During the operation, we store in modified_darts the darts modified
      // to test after the loop the non void attributes that are split.
      std::deque<typename GMap::Dart_handle> modified_darts;

      // For each dart of the i-cell, we modify i-links of neighbors.
      typename std::deque<typename GMap::Dart_handle>::iterator it =
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

      if ( amap.are_attributes_automatically_managed() && update_attributes )
      {
        // We test the split of all the incident cells for all the non
        // void attributes.
        GMap::Helper::template Foreach_enabled_attributes_except
          <CGAL::internal::GMap_test_split_attribute_functor<GMap,i>, i>::
          run(amap, modified_darts, mark_modified_darts);
      }

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

} // namespace CGAL

#endif // CGAL_GENERALIZED_MAP_OPERATIONS_H //
// EOF //
