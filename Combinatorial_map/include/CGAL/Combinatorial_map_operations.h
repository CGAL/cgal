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
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_COMBINATORIAL_MAP_OPERATIONS_H
#define CGAL_COMBINATORIAL_MAP_OPERATIONS_H 1

#include <CGAL/Combinatorial_map_basic_operations.h>
#include <CGAL/Combinatorial_map_insertions.h>
#include <deque>
#include <stack>

namespace CGAL
{
  /** @file Combinatorial_map_operations.h
   * Some operations to modify a combinatorial map.
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
      for ( CGAL::CMap_dart_const_iterator_of_cell<CMap,i> it(amap, adart);
            res && it.cont(); ++it )
      {
        if ( amap.template beta<i+2,i+1>(it)!=
             amap.template beta<CGAL_BETAINV(i+1),i+2>(it) )
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

  /** Remove an i-cell, 0<i<dimension, and merge eventually both incident
   *  (i+1)-cells.
   *  @param amap the used combinatorial map.
   *  @param adart a dart of the i-cell to remove.
   *  @param update_attributes a boolean to update the enabled attributes
   *         (deprecated, now we use are_attributes_automatically_managed())
   *  @return the number of deleted darts.
   */
  template<class CMap, unsigned int i, unsigned int nmi>
  struct Remove_cell_functor
  {
    static size_t run(CMap& amap, typename CMap::Dart_handle adart, bool update_attributes)
    {
      CGAL_static_assertion ( 1<=i && i<CMap::dimension );
      CGAL_assertion( (is_removable<CMap,i>(amap, adart)) );

      size_t res = 0;

      typename CMap::Dart_handle d1, d2;
      typename CMap::Dart_handle dg1=amap.null_handle, dg2=amap.null_handle;

      int mark = amap.get_new_mark();
      int mark_modified_darts = amap.get_new_mark();

      std::deque<typename CMap::Dart_handle> to_erase;

      const int iinv = CGAL_BETAINV(i);

      // First we store and mark all the darts of the i-cell to remove.
      for ( CGAL::CMap_dart_iterator_basic_of_cell<CMap,i> it(amap,adart,mark);
            it.cont(); ++it )
      {
        to_erase.push_back(it);
        if ( !amap.template is_free<i+1>(it) && dg1==amap.null_handle )
        { dg1=it; dg2=amap.template beta<i+1>(it); }
        amap.mark(it, mark);
        ++res;
      }

      if (amap.are_attributes_automatically_managed() && update_attributes)
      {
        // We group the two (i+1)-cells incident if they exist.
        if ( dg1!=amap.null_handle )
          CGAL::internal::Group_attribute_functor_run<CMap, i+1>::
              run(&amap, dg1, dg2);
      }

      // During the operation, we store in modified_darts the darts modified
      // to test after the loop the non void attributes that are split.
      std::deque<typename CMap::Dart_handle> modified_darts;

      // If i==1, we modify beta1, thus in modified_darts we store all
      // the darts having beta0 modified, and in modified_darts2 all the
      // darts having beta1 modified. For i>1 all the modified darts are
      // stored in modified_darts.
      std::deque<typename CMap::Dart_handle> modified_darts2;
      std::deque<typename CMap::Dart_handle> &first_modified_darts=
          (i==1?modified_darts2:modified_darts);

      // For each dart of the i-cell, we modify i-links of neighbors.
      typename std::deque<typename CMap::Dart_handle>::iterator it =
          to_erase.begin();
      for ( ; it!=to_erase.end(); ++it )
      {
        d1=amap.template beta<iinv>(*it);
        while ( d1!=amap.null_dart_handle && amap.is_marked(d1, mark) )
        {
          d1=amap.template beta<i+1, iinv>(d1);
          if ( d1==amap. template beta<iinv>(*it) ) d1=amap.null_dart_handle;
        }

        if ( !amap.is_marked(d1, mark_modified_darts) )
        {
          d2=amap.template beta<i+1,i>(*it);
          while ( d2!=amap.null_dart_handle && amap.is_marked(d2, mark) )
          {
            d2=amap.template beta<i+1,i>(d2);
            if ( d2==amap.template beta<i+1,i>(*it) )
              d2=amap.null_dart_handle;
          }

          if ( !amap.is_marked(d2, mark_modified_darts) )
          {
            if ( d1!=amap.null_dart_handle )
            {
              if ( d2!=amap.null_dart_handle && d1!=d2 )
              {
                //d1->basic_link_beta(d2, i);
                amap.template basic_link_beta<i>(d1, d2);
                amap.mark(d1, mark_modified_darts);
                amap.mark(d2, mark_modified_darts);
                first_modified_darts.push_back(d1);
                modified_darts.push_back(d2);
                // TODO push only one out of two dart ?

                /*if ( i==1 )
                {
                  d2->basic_link_beta(d1, 0);
                  modified_darts.push_back(d2);
                }*/
                //            modified_darts2.push_back(d1);
              }
              else
              {
                if ( !amap.template is_free<i>(d1) )
                {
                  amap.template unlink_beta<i>(d1);
                  CGAL_assertion( !amap.is_marked(d1, mark_modified_darts) );
                  amap.mark(d1, mark_modified_darts);
                  first_modified_darts.push_back(d1);
                }
              }
            }
            else if ( d2!=amap.null_dart_handle )
            {
              if ( !amap.template is_free<iinv>(d2) )
              {
                amap.template unlink_beta<iinv>(d2);
                CGAL_assertion( !amap.is_marked(d2, mark_modified_darts) );
                amap.mark(d2, mark_modified_darts);
                modified_darts.push_back(d2);
              }
            }
          }
        }
        if ( amap.template is_free<i+1>(*it) &&
             !amap.template is_free<i>(*it) )
        {
          d1 = amap.template beta<i>(*it);
          if ( !amap.template is_free<iinv>(d1) )
          {
            amap.template unlink_beta<iinv>(d1);
            CGAL_assertion( !amap.is_marked(d1, mark_modified_darts) );
            amap.mark(d1, mark_modified_darts);
            modified_darts.push_back(d1);
          }
        }
      }

      if (amap.are_attributes_automatically_managed() && update_attributes)
      {
        // We test the split of all the incident cells for all the non
        // void attributes.
        if ( i==1 )
          CMap::Helper::template Foreach_enabled_attributes_except
              <CGAL::internal::Test_split_attribute_functor<CMap,i>, i>::
              run(&amap, modified_darts, modified_darts2,
                  mark_modified_darts);
        else
          CMap::Helper::template Foreach_enabled_attributes_except
              <CGAL::internal::Test_split_attribute_functor<CMap,i>, i>::
              run(&amap, modified_darts, mark_modified_darts);
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
        for ( typename std::deque<typename CMap::Dart_handle>::
                iterator it=modified_darts.begin();
              it!=modified_darts.end(); ++it )
          amap.unmark(*it, mark_modified_darts);
        if ( i==1 )
        {
          for ( typename std::deque<typename CMap::Dart_handle>::
                  iterator it=modified_darts2.begin();
                it!=modified_darts2.end(); ++it )
            amap.unmark(*it, mark_modified_darts);
        }
      }

      CGAL_assertion ( amap.is_whole_map_unmarked(mark_modified_darts) );
      amap.free_mark(mark_modified_darts);

#ifdef CGAL_CMAP_TEST_VALID_REMOVALS
      CGAL_assertion( amap.is_valid() );
#endif

      return res;
    }
  };

  /** Remove a d-cell, in a d-map (special case).
   *  @param amap the used combinatorial map.
   *  @param adart a dart of the volume to remove.
   *  @param update_attributes a boolean to update the enabled attributes
   *         (deprecated, now we use are_attributes_automatically_managed())
   *  @return the number of deleted darts.
   */
  template<class CMap,unsigned int i>
  struct Remove_cell_functor<CMap,i,0>
  {
    static size_t run(CMap& amap, typename CMap::Dart_handle adart, bool update_attributes)
    {
      int mark = amap.get_new_mark();
      std::deque<typename CMap::Dart_handle> to_erase;
      size_t res = 0;

      std::deque<typename CMap::Dart_handle> modified_darts;

      // We mark all the darts of the d-cell.
      for ( CGAL::CMap_dart_iterator_basic_of_cell<CMap,CMap::dimension>
            it(amap,adart,mark); it.cont(); ++it )
      {
        to_erase.push_back(it);
        amap.mark(it,mark);
        ++res;
      }

      // We unlink all the darts of the volume for beta-d.
      typename std::deque<typename CMap::Dart_handle>::iterator
        it = to_erase.begin();
      for ( it = to_erase.begin(); it != to_erase.end(); ++it )
      {
        if ( !amap.template is_free<CMap::dimension>(*it) &&
             !amap.is_marked(amap.template beta<CMap::dimension>(*it), mark) )
        {
          if (amap.are_attributes_automatically_managed() && update_attributes)
          {
            modified_darts.push_back(amap.template beta<CMap::dimension>(*it));
          }
          amap.template unlink_beta_for_involution<CMap::dimension>(*it);
        }
      }

      if (amap.are_attributes_automatically_managed() && update_attributes)
      {
        // We test the split of all the incident cells for all the non
        // void attributes.
        CMap::Helper::template Foreach_enabled_attributes_except
          <CGAL::internal::Test_split_attribute_functor<CMap,i>,
           CMap::dimension>::run(&amap, modified_darts);
      }

      // We remove all the darts of the d-cell.
      for ( it = to_erase.begin(); it != to_erase.end(); ++it )
      { amap.erase_dart(*it); }

      CGAL_assertion( amap.is_whole_map_unmarked(mark) );
      amap.free_mark(mark);

#ifdef CGAL_CMAP_TEST_VALID_REMOVALS
      CGAL_assertion( amap.is_valid() );
#endif

      return res;
    }
  };

  /** Remove a vertex, and merge eventually both incident edges.
   * @param amap the used combinatorial map.
   * @param adart a dart of the vertex to remove.
   * @param update_attributes a boolean to update the enabled attributes
   *        (deprecated, now we use are_attributes_automatically_managed())
   * @return the number of deleted darts.
   */
  template<class CMap,unsigned int nmi>
  struct Remove_cell_functor<CMap,0,nmi>
  {
    static size_t run(CMap& amap, typename CMap::Dart_handle adart,
                      bool update_attributes)
    {
      CGAL_assertion( (is_removable<CMap,0>(amap,adart)) );

      size_t res = 0;

      typename CMap::Dart_handle d1, d2;
      typename CMap::Dart_handle dg1=amap.null_handle, dg2=amap.null_handle;

      int mark = amap.get_new_mark();
//      int mark_modified_darts = amap.get_new_mark();

      // First we store and mark all the darts of the 0-cell to remove.
      std::deque<typename CMap::Dart_handle> to_erase;
      for ( CGAL::CMap_dart_iterator_basic_of_cell<CMap,0> it(amap,adart,mark);
            it.cont(); ++it )
      {
        to_erase.push_back(it);
        if ( !amap.template is_free<0>(it) && dg1==amap.null_handle )
        { dg1=it; dg2=amap.template beta<0>(it); }
        amap.mark(it, mark);
        ++res;
      }

      if (amap.are_attributes_automatically_managed() && update_attributes )
      {
        // We group the two edges incident if they exist.
        if ( dg1!=amap.null_handle )
          CGAL::internal::Group_attribute_functor_run<CMap, 1>::
              run(&amap, dg1, dg2);
      }

      // During the operation, we store in modified_darts the darts modified
      // by beta0 to test after the loop non void attributes that are split.
      std::deque<typename CMap::Dart_handle> modified_darts;
      // And we store in modified_darts2 all the darts having beta1 modified.
      std::deque<typename CMap::Dart_handle> modified_darts2;

      // For each dart of the vertex, we modify 0 and 1-links of neighbors.
      typename std::deque<typename CMap::Dart_handle>::iterator it =
          to_erase.begin();
      for ( ; it != to_erase.end(); ++it)
      {
        if ( !amap.template is_free<0>(*it) )
        {
          if ( !amap.template is_free<1>(*it) &&
               amap.template beta<0>(*it)!=(*it) )
          {
            amap.basic_link_beta_1(amap.template beta<0>(*it),
                                   amap.template beta<1>(*it));
            if (amap.are_attributes_automatically_managed() && update_attributes)
            {
              modified_darts.push_back(amap.template beta<0>(*it));
              modified_darts2.push_back(amap.template beta<1>(*it));
              // TODO push only one out of two dart ?
            }
          }
          else
          {
            amap.template dart_unlink_beta<1>(amap.template beta<0>(*it));
            if (amap.are_attributes_automatically_managed() && update_attributes)
            {
              modified_darts.push_back(amap.template beta<0>(*it));
            }
          }

          for ( unsigned int j=2; j<=CMap::dimension; ++j )
          {
            if ( !amap.is_free(*it,j) )
            {
              amap.basic_link_beta(amap.template beta<0>(*it),
                                   amap.beta(*it,j), j);
            //((*it)->beta(0))->basic_link_beta((*it)->beta(j),j);
            }
          }
        }
        else
        {
          if ( !amap.template is_free<1>(*it) )
          {
            amap.template dart_unlink_beta<0>(amap.template beta<1>(*it));
            if (amap.are_attributes_automatically_managed() && update_attributes)
            {
              modified_darts2.push_back(amap.template beta<1>(*it));
            }
          }

          for ( unsigned int j=2; j<=CMap::dimension; ++j )
          {
            if ( !amap.is_free(*it,j) )
            { amap.unlink_beta(*it, j); }
          }
        }
      }

      if (amap.are_attributes_automatically_managed() && update_attributes)
      {
        // We test the split of all the incident cells for all the non
        // void attributes.
        CMap::Helper::template Foreach_enabled_attributes_except
            <CGAL::internal::Test_split_attribute_functor<CMap,0>, 1>::
            run(&amap,modified_darts, modified_darts2);
      }

      // We remove all the darts of the 0-cell.
      for ( it=to_erase.begin(); it!=to_erase.end(); ++it )
      { amap.erase_dart(*it); }

      CGAL_assertion( amap.is_whole_map_unmarked(mark) );
      amap.free_mark(mark);

#ifdef CGAL_CMAP_TEST_VALID_REMOVALS
      CGAL_assertion( amap.is_valid() );
#endif

      return res;
    }
  };

  /** Remove an i-cell, 0<=i<=dimension.
   * @param amap the used combinatorial map.
   * @param adart a dart of the i-cell to remove.
   * @param update_attributes a boolean to update the enabled attributes
   * @return the number of deleted darts.
   */
  template < class CMap, unsigned int i >
  size_t remove_cell(CMap& amap, typename CMap::Dart_handle adart, bool update_attributes = true)
  {
    return
        CGAL::Remove_cell_functor<CMap,i,CMap::dimension-i>::run(amap,adart,update_attributes);
  }

  /** Test if an i-cell can be contracted.
   *  An i-cell can be contracted if i==1
   *     or if there are at most two (i-1)-cell incident to it.
   * @param adart a dart of the i-cell.
   * @return true iff the i-cell can be contracted.
   */
  template <class CMap, unsigned int i>
  struct Is_contractible_functor
  {
    static bool run(const CMap& amap, typename CMap::Dart_const_handle adart)
    {
      // TODO ? Optimisation possible to not test all the darts of the cell ?
      bool res = true;
      for ( CGAL::CMap_dart_const_iterator_of_cell<CMap,i> it(amap, adart);
            res && it.cont(); ++it )
      {
        if ( amap.template beta<i-2,i-1>(it)!=
             amap.template beta<i-1,CGAL_BETAINV(i-2)>(it) )
          res = false;
      }
      return res;
    }
  };
  // Specialization for i=0
  template <class CMap>
  struct Is_contractible_functor<CMap, 0>
  {
    static bool run(const CMap&, typename CMap::Dart_const_handle)
    { return false; }
  };
  // Specialization for i=1
  template <class CMap>
  struct Is_contractible_functor<CMap, 1>
  {
    static bool run(const CMap&, typename CMap::Dart_const_handle)
    { return true; }
  };
  /** Test if an i-cell can be contracted.
   *  An i-cell can be contracted if i==1
   *     or if there are at most two (i-1)-cell incident to it.
   * @param adart a dart of the i-cell.
   * @return true iff the i-cell can be contracted.
   */
  template < class CMap, unsigned int i >
  bool is_contractible(const CMap& amap, typename CMap::Dart_const_handle adart)
  { return CGAL::Is_contractible_functor<CMap, i>::run(amap,adart); }

  /** Contract an i-cell, 1<i<=dimension, and merge eventually both incident
   *  (i-1)-cells.
   * @param amap the used combinatorial map.
   * @param adart a dart of the i-cell to contract.
   * @return the number of deleted darts.
   */
  template<class CMap, unsigned int i>
  struct Contract_cell_functor
  {
    static size_t run(CMap& amap, typename CMap::Dart_handle adart)
    {
      CGAL_static_assertion ( 2<=i && i<=CMap::dimension );
      CGAL_assertion( (is_contractible<CMap,i>(amap, adart)) );

      size_t res = 0;

      typename CMap::Dart_handle d1, d2;
      typename CMap::Dart_handle dg1=amap.null_handle, dg2=amap.null_handle;

      int mark = amap.get_new_mark();
      int mark_modified_darts = amap.get_new_mark();

      const int imuinv = CGAL_BETAINV(i-1);

      // First we store and mark all the darts of the i-cell to contract.
      std::deque<typename CMap::Dart_handle> to_erase;
      for ( CGAL::CMap_dart_iterator_basic_of_cell<CMap,i> it(amap,adart,mark);
            it.cont(); ++it )
      {
        to_erase.push_back(it);
        if ( !amap.template is_free<i-1>(it) && dg1==amap.null_handle )
        { dg1=it; dg2=amap.template beta<i-1>(it); }
        amap.mark(it, mark);
        ++res;
      }

      if ( amap.are_attributes_automatically_managed() )
      {
        // We group the two (i+1)-cells incident if they exist.
        if ( dg1!=amap.null_handle )
          CGAL::internal::Group_attribute_functor_run<CMap,i-1>::
            run(&amap, dg1, dg2);
      }

      // During the operation, we store in modified_darts the darts modified
      // to test after the loop the non void attributes that are split.
      std::deque<typename CMap::Dart_handle> modified_darts;

      // For each dart of the i-cell, we modify i-links of neighbors.
      typename std::deque<typename CMap::Dart_handle>::iterator it =
          to_erase.begin();
      for ( ; it!=to_erase.end(); ++it )
      {
        d1 = amap.template beta<i>(*it);
        while ( d1!=amap.null_dart_handle && amap.is_marked(d1, mark) )
        {
          d1 = amap.template beta<imuinv,i>(d1);
          if (d1 == amap.template beta<i>(*it)) d1 = amap.null_dart_handle;
        }

        if ( !amap.is_marked(d1, mark_modified_darts) )
        {
          d2 = amap.template beta<i-1,i>(*it);
          while ( d2!=amap.null_dart_handle && amap.is_marked(d2, mark) )
          {
            d2 = amap.template beta<i-1,i>(d2);
            if ( d2==amap.template beta<i-1,i>(*it) )
              d2=amap.null_dart_handle;
          }

          if ( !amap.is_marked(d2, mark_modified_darts) )
          {
            if (d1 != amap.null_dart_handle)
            {
              if (d2 != amap.null_dart_handle && d1!=d2 )
              {
                amap.template basic_link_beta_for_involution<i>(d1, d2);
                amap.mark(d1, mark_modified_darts);
                amap.mark(d2, mark_modified_darts);
                modified_darts.push_back(d1);
                modified_darts.push_back(d2);
                // TODO push only one out of two dart ?
              }
              else
              {
                if ( !amap.template is_free<i>(d1) )
                {
                  amap.template unlink_beta<i>(d1);
                  CGAL_assertion( !amap.is_marked(d1, mark_modified_darts) );
                  amap.mark(d1, mark_modified_darts);
                  modified_darts.push_back(d1);
                }
              }
            }
            else if (d2 != amap.null_dart_handle)
            {
              if ( !amap.is_free(d2,i) )
              {
                amap.template unlink_beta<i>(d2);
                CGAL_assertion( !amap.is_marked(d2, mark_modified_darts) );
                amap.mark(d2, mark_modified_darts);
                modified_darts.push_back(d2);
              }
            }
          }
        }
        if (amap.is_free(*it,i-1) && !amap.is_free(*it,i))
        {
          d1 = amap.beta(*it,i);
          if ( !amap.is_free(d1,i) )
          {
            amap.template unlink_beta<i>(d1);
            CGAL_assertion( !amap.is_marked(d1, mark_modified_darts) );
            amap.mark(d1, mark_modified_darts);
            modified_darts.push_back(d1);
          }
        }
      }

      if ( amap.are_attributes_automatically_managed() )
      {
        // We test the split of all the incident cells for all the non
        // void attributes.
        CMap::Helper::template Foreach_enabled_attributes_except
          <CGAL::internal::Test_split_attribute_functor<CMap,i>, i>::
          run(&amap, modified_darts, mark_modified_darts);
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
        for ( typename std::deque<typename CMap::Dart_handle>::
              iterator it=modified_darts.begin();
              it!=modified_darts.end(); ++it )
          amap.unmark(*it, mark_modified_darts);
      }

      CGAL_assertion ( amap.is_whole_map_unmarked(mark_modified_darts) );
      amap.free_mark(mark_modified_darts);

#ifdef CGAL_CMAP_TEST_VALID_CONTRACTIONS
      CGAL_assertion( amap.is_valid() );
#endif

      return res;
    }
  };

  /** Contract an edge, and merge eventually both incident vertices.
   * @param amap the used combinatorial map.
   * @param adart a dart of the edge to contract.
   * @return the number of deleted darts.
   */
  template<class CMap>
  struct Contract_cell_functor<CMap,1>
  {
    static size_t run(CMap& amap, typename CMap::Dart_handle adart)
    {
      CGAL_assertion( (is_contractible<CMap,1>(amap,adart)) );

      size_t res = 0;

      typename CMap::Dart_handle d1, d2;
      typename CMap::Dart_handle dg1=amap.null_handle, dg2=amap.null_handle;

      int mark = amap.get_new_mark();
//      int mark_modified_darts = amap.get_new_mark();

      // First we store and mark all the darts of the 1-cell to contract.
      std::deque<typename CMap::Dart_handle> to_erase;
      for ( CGAL::CMap_dart_iterator_basic_of_cell<CMap,1> it(amap,adart,mark);
            it.cont(); ++it )
      {
        to_erase.push_back(it);
        if ( dg1==amap.null_handle && !amap.template is_free<0>(it) &&
             !amap.template is_free<1>(it) )
        { dg1=amap.template beta<0>(it); dg2=amap.template beta<1>(it); }
        amap.mark(it, mark);
        ++res;
      }

      if ( amap.are_attributes_automatically_managed() )
      {
        // We group the two vertices incident if they exist.
        if ( dg1!=amap.null_handle )
          CGAL::internal::Group_attribute_functor_run<CMap, 0, 1>::
            run(&amap, dg1, dg2);
      }

      // During the operation, we store in modified_darts the darts modified
      // by beta0 to test after the loop non void attributes that are split.
      std::deque<typename CMap::Dart_handle> modified_darts;
      // And we store in modified_darts2 all the darts having beta1 modified.
      std::deque<typename CMap::Dart_handle> modified_darts2;

      // For each dart of the cell, we modify link of neighbors.
      typename std::deque<typename CMap::Dart_handle>::iterator it =
        to_erase.begin();
      for ( ; it!=to_erase.end(); ++it )
      {
        if ( !amap.template is_free<0>(*it) )
        {
          if ( !amap.template is_free<1>(*it) )
          {
            if ( amap.template beta<1>(*it)!=*it )
            {
               /*modified_darts2.push_back((*it)->template beta<0>());
              if ( (*it)->beta(0)!=(*it)->beta(1) )*/
              if ( amap.are_attributes_automatically_managed() )
              {
                modified_darts.push_back(amap.template beta<1>(*it));
              }
              amap.basic_link_beta_1(amap.template beta<0>(*it),
                                     amap.template beta<1>(*it));
            }
          }
          else
          {
            if ( amap.are_attributes_automatically_managed() )
            {
              modified_darts2.push_back(amap.template beta<0>(*it));
            }
            amap.template dart_unlink_beta<1>(amap.template beta<0>(*it));
          }
        }
        else
        {
          if ( !amap.template is_free<1>(*it) )
          {
            if ( amap.are_attributes_automatically_managed() )
            {
              modified_darts.push_back(amap.template beta<1>(*it));
            }
            amap.template dart_unlink_beta<0>(amap.template beta<1>(*it));
          }
        }
      }

      // We remove all the darts of the cell.
      for ( it=to_erase.begin(); it!=to_erase.end(); ++it )
      { amap.erase_dart(*it); }

      CGAL_assertion( amap.is_whole_map_unmarked(mark) );
      amap.free_mark(mark);

      if ( amap.are_attributes_automatically_managed() )
      {
        // We test the split of all the incident cells for all the non
        // void attributes.
        CMap::Helper::template Foreach_enabled_attributes_except
          <CGAL::internal::Test_split_attribute_functor<CMap,1>, 1>::
          run(&amap, modified_darts, modified_darts2);
      }

#ifdef CGAL_CMAP_TEST_VALID_CONTRACTIONS
      CGAL_assertion( amap.is_valid() );
#endif

      return res;
    }
  };

  /** Contract an i-cell, 1<=i<=dimension.
   * @param amap the used combinatorial map.
   * @param adart a dart of the i-cell to remove.
   * @return the number of deleted darts.
   */
  template < class CMap, unsigned int i >
  size_t contract_cell(CMap& amap, typename CMap::Dart_handle adart)
  { return CGAL::Contract_cell_functor<CMap,i>::run(amap,adart); }

} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_OPERATIONS_H //
// EOF //
