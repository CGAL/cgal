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

namespace CGAL {

  /** @file Combinatorial_map_operations.h
   * Some operations to modify a combinatorial map.
   */

  /** Test if an i-cell can be removed.
   *  An i-cell can be removed if i==CMap::dimension or i==CMap::dimension-1,
   *     or if there are at most two (i+1)-cell incident to it.
   * @param adart a dart of the i-cell.
   * @return true iff the i-cell can be removed.
   */
  template < class CMap, unsigned int i >
  bool is_removable(const CMap& amap, typename CMap::Dart_const_handle adart)
  {
    CGAL_assertion( adart!=NULL );
    CGAL_static_assertion( 0<=i && i<=CMap::dimension );

    if ( i==CMap::dimension   ) return true;
    if ( i==CMap::dimension-1 ) return true;

    // TODO? Optimisation for dim-2, and to not test all the darts of the cell?
    bool res = true;
    for (CMap_dart_const_iterator_of_cell<CMap,i> it(amap, adart);
         res && it.cont(); ++it)
    {
      if (it->template beta<i+2>()->template beta<i+1>()!=
          it->template beta_inv<i+1>()->template beta<i+2>() )
        res = false;
    }
    return res;
  }

  /** Remove an i-cell, 0<i<dimension, and merge eventually both incident
   *  (i+1)-cells.
   *  @param amap the used combinatorial map.
   *  @param adart a dart of the i-cell to remove.
   *  @return the number of deleted darts.
   */
  template<class CMap, unsigned int i, unsigned int nmi>
  struct Remove_cell_functor
  {
    static size_t run(CMap& amap, typename CMap::Dart_handle adart)
    {
      CGAL_static_assertion ( 1<=i && i<CMap::dimension );
      CGAL_assertion( (is_removable<CMap,i>(amap, adart)) );

      size_t res = 0;

      typename CMap::Dart_handle d1, d2;
      typename CMap::Dart_handle dg1=NULL, dg2=NULL;

      int mark = amap.get_new_mark();
      int mark_modified_darts = amap.get_new_mark();

      std::deque<typename CMap::Dart_handle> to_erase;

      const int iinv = CGAL_BETAINV(i);

      // First we store and mark all the darts of the i-cell to remove.
      for ( CMap_dart_iterator_basic_of_cell<CMap,i> it(amap,adart,mark);
            it.cont(); ++it )
      {
        to_erase.push_back(it);
        if ( !it->template is_free<i+1>() && dg1==NULL )
        { dg1=it; dg2=it->template beta<i+1>(); }
        amap.mark(it, mark);
        ++res;
      }

      // We group the two (i+1)-cells incident if they exist.
      if ( dg1!=NULL )
        internal::Group_attribute_functor_run<CMap, i+1>::
            run(&amap, dg1, dg2);

      // Second we store all the incident cells that can be split by
      // the operation.
      typename std::deque<typename CMap::Dart_handle>::iterator it =
          to_erase.begin();
      //TODO remove ? or be sure that this is required
      for (; it != to_erase.end(); ++it)
        amap.update_dart_of_all_attributes(*it, mark);

      std::deque<typename CMap::Dart_handle> modified_darts;
      std::deque<typename CMap::Dart_handle> modified_darts2;

      // If i==1, we modify beta1, thus in modified_darts we store all
      // the darts having beta0 modified, and in modified_darts2 all the
      // darts having beta1 modified. Otherwise we store all the darts in
      // modified_darts.
      std::deque<typename CMap::Dart_handle> &first_modified_darts=
          (i==1?modified_darts2:modified_darts);

      // For each dart of the i-cell, we modify i-links of neighbors.
      for ( it=to_erase.begin(); it!=to_erase.end(); ++it )
      {
        d1=(*it)->template beta<iinv>();
        while ( d1!=CMap::null_dart_handle && amap.is_marked(d1, mark) )
        {
          d1=d1->template beta<i+1>()->template beta<iinv>();
          if ( d1==(*it)->template beta<iinv>() ) d1=CMap::null_dart_handle;
        }

        if ( !amap.is_marked(d1, mark_modified_darts) )
        {
          d2=(*it)->template beta<i+1>()->template beta<i>();
          while ( d2!=CMap::null_dart_handle && amap.is_marked(d2, mark) )
          {
            d2=d2->template beta<i+1>()->template beta<i>();
            if ( d2==(*it)->template beta<i+1>()->template beta<i>() )
              d2=CMap::null_dart_handle;
          }

          if ( !amap.is_marked(d2, mark_modified_darts) )
          {
            if ( d1!=CMap::null_dart_handle )
            {
              if ( d2!=CMap::null_dart_handle && d1!=d2 )
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
                if ( !d1->template is_free<i>() )
                {
                  d1->template unlink_beta<i>();
                  CGAL_assertion( !amap.is_marked(d1, mark_modified_darts) );
                  amap.mark(d1, mark_modified_darts);
                  first_modified_darts.push_back(d1);
                }
              }
            }
            else if ( d2!=CMap::null_dart_handle )
            {
              if ( !d2->template is_free<iinv>() )
              {
                d2->template unlink_beta<iinv>();
                CGAL_assertion( !amap.is_marked(d2, mark_modified_darts) );
                amap.mark(d2, mark_modified_darts);
                modified_darts.push_back(d2);
              }
            }
          }
        }
        if ( (*it)->template is_free<i+1>() &&
             !(*it)->template is_free<i>() )
        {
          d1 = (*it)->template beta<i>();
          if ( !d1->template is_free<iinv>() )
          {
            d1->template unlink_beta<iinv>();
            CGAL_assertion( !amap.is_marked(d1, mark_modified_darts) );
            amap.mark(d1, mark_modified_darts);
            modified_darts.push_back(d1);
          }
        }
      }

      // We test the split of all the incident cells for all the non
      // void attributes.
      if ( i==1 )
        CMap::Helper::template Foreach_enabled_attributes_except
            <internal::Test_split_attribute_functor<CMap,i>, i>::
            run(&amap, modified_darts, modified_darts2,
                mark_modified_darts);
      else
        CMap::Helper::template Foreach_enabled_attributes_except
            <internal::Test_split_attribute_functor<CMap,i>, i>::
            run(&amap, modified_darts, mark_modified_darts);

      // We remove all the darts of the i-cell.
      for (  it=to_erase.begin(); it!=to_erase.end(); ++it )
      { amap.erase_dart(*it); }

      CGAL_assertion( amap.is_whole_map_unmarked(mark) );
      amap.free_mark(mark);

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

      CGAL_expensive_postcondition( amap.is_valid() );
      assert( amap.is_valid() ); // TODO remove

      return res;
    }
  };

  /** Remove a d-cell, in a d-map (special case).
   *  @param amap the used combinatorial map.
   *  @param adart a dart of the volume to remove.
   *  @return the number of deleted darts.
   */
  template<class CMap,unsigned int i>
  struct Remove_cell_functor<CMap,i,0>
  {
    static size_t run(CMap& amap, typename CMap::Dart_handle adart)
    {
      CGAL_assertion( adart!=NULL );

      int mark = amap.get_new_mark();
      std::deque<typename CMap::Dart_handle> to_erase;
      size_t res = 0;

      std::deque<typename CMap::Dart_handle> modified_darts;

      // 1) We mark all the darts of the d-cell.
      for (CMap_dart_iterator_basic_of_cell<CMap,CMap::dimension>
           it(amap,adart,mark); it.cont(); ++it)
      {
        to_erase.push_back(it);
        amap.mark(it,mark);
        ++res;
      }

      typename std::deque<typename CMap::Dart_handle>::iterator
        it = to_erase.begin();

      // 3) We unlink all the darts of the volume for beta-d.
      for ( it = to_erase.begin(); it != to_erase.end(); ++it )
      {
        if ( !(*it)->template is_free<CMap::dimension>() &&
             !amap.is_marked((*it)->template beta<CMap::dimension>(), mark) )
        {
          modified_darts.push_back((*it)->template beta<CMap::dimension>());
          amap.template unlink_beta_for_involution<CMap::dimension>(*it);
        }
      }

      // 4) We test the split of all the incident cells for all the non
      // void attributes.
      CMap::Helper::template Foreach_enabled_attributes_except
          <internal::Test_split_attribute_functor<CMap,i>,
          CMap::dimension>::run(&amap, modified_darts);

      // 5) We remove all the darts of the d-cell.
      for ( it = to_erase.begin(); it != to_erase.end(); ++it )
      { amap.erase_dart(*it); }

      CGAL_assertion( amap.is_whole_map_unmarked(mark) );
      amap.free_mark(mark);

      CGAL_expensive_postcondition( amap.is_valid() );
      assert( amap.is_valid() ); // TO REMOVE

      return res;
    }
  };

  /** Remove a vertex, and merge eventually both incident edges.
   * @param amap the used combinatorial map.
   * @param adart a dart of the vertex to remove.
   * @return the number of deleted darts.
   */
  template<class CMap,unsigned int nmi>
  struct Remove_cell_functor<CMap,0,nmi>
  {
    static size_t run(CMap& amap, typename CMap::Dart_handle adart)
    {
      CGAL_assertion( (is_removable<CMap,0>(amap,adart)) );

      size_t res = 0;

      typename CMap::Dart_handle d1, d2;
      typename CMap::Dart_handle dg1=NULL, dg2=NULL;

      int mark = amap.get_new_mark();
//      int mark_modified_darts = amap.get_new_mark();

      std::deque<typename CMap::Dart_handle> to_erase;

      std::deque<typename CMap::Dart_handle> modified_darts;
      std::deque<typename CMap::Dart_handle> modified_darts2;

      // First we store and mark all the darts of the 0-cell to remove.
      for ( CMap_dart_iterator_basic_of_cell<CMap,0> it(amap,adart,mark);
            it.cont(); ++it )
      {
        to_erase.push_back(it);
        if ( !it->template is_free<0>() && dg1==NULL )
        { dg1=it; dg2=it->template beta<0>(); }
        amap.mark(it, mark);
        ++res;
      }

      // Second we store all the incident cells that can be split by
      // the operation.
      typename std::deque<typename CMap::Dart_handle>::iterator it =
        to_erase.begin();
      for (; it != to_erase.end(); ++it)
        amap.update_dart_of_all_attributes(*it, mark);

      // We group the two edges incident if they exist.
      if ( dg1!=NULL )
        internal::Group_attribute_functor_run<CMap, 1>::
            run(&amap, dg1, dg2);

      // For each dart of the vertex, we modify 0 and 1-links of neighbors.
      for ( it=to_erase.begin(); it != to_erase.end(); ++it)
      {
        if ( !(*it)->template is_free<0>() )
        {
          if ( !(*it)->template is_free<1>() &&
               (*it)->template beta<0>()!=(*it) )
          {
            amap.basic_link_beta_1((*it)->template beta<0>(),
                                   (*it)->template beta<1>());
            modified_darts.push_back((*it)->template beta<0>());
            modified_darts2.push_back((*it)->template beta<1>());
            // TODO push only one out of two dart ?
          }
          else
          {
            (*it)->template beta<0>()->template unlink_beta<1>();
            modified_darts.push_back((*it)->template beta<0>());
          }

          for ( unsigned int j=2; j<=CMap::dimension; ++j )
          {
            if ( !(*it)->is_free(j) )
            {
              // TODO push these darts in modified_darts ?
              // not sure this is required
              amap.basic_link_beta((*it)->template beta<0>(),
                                   (*it)->beta(j), j);
            //((*it)->beta(0))->basic_link_beta((*it)->beta(j),j);
            }
          }
        }
        else
        {
          if ( !(*it)->template is_free<1>() )
          {
            (*it)->template beta<1>()->template unlink_beta<0>();
            modified_darts2.push_back((*it)->template beta<1>());
          }

          for ( unsigned int j=2; j<=CMap::dimension; ++j )
          {
            if ( !(*it)->is_free(j) )
            {
              // TODO push these darts in modified_darts ?
              // not sure this is required
              amap.unlink_beta(*it, j);
            }
          }
        }
      }

      // We test the split of all the incident cells for all the non
      // void attributes.
      CMap::Helper::template Foreach_enabled_attributes_except
          <internal::Test_split_attribute_functor<CMap,0>, 1>::
                run(&amap,
                    modified_darts, modified_darts2);

      // We remove all the darts of the i-cell.
      for (  it=to_erase.begin(); it!=to_erase.end(); ++it )
      { amap.erase_dart(*it); }

      CGAL_assertion( amap.is_whole_map_unmarked(mark) );
      amap.free_mark(mark);

      CGAL_expensive_postcondition( amap.is_valid() );
      assert( amap.is_valid() ); // TO REMOVEE

      return res;
    }
  };

  /** Remove an i-cell, 0<=i<=dimension.
   * @param amap the used combinatorial map.
   * @param adart a dart of the i-cell to remove.
   * @return the number of deleted darts.
   */
  template < class CMap, unsigned int i >
  size_t remove_cell(CMap& amap, typename CMap::Dart_handle adart)
  { return Remove_cell_functor<CMap,i,CMap::dimension-i>::run(amap,adart); }

  /** Test if an i-cell can be contracted.
   *  An i-cell can be contracted if i==1
   *     or if there are at most two (i-1)-cell incident to it.
   * @param adart a dart of the i-cell.
   * @return true iff the i-cell can be contracted.
   */
  template < class CMap, unsigned int i >
  bool is_contractible(const CMap& amap, typename CMap::Dart_const_handle adart)
  {
    CGAL_assertion(adart != NULL);
    CGAL_static_assertion(0<=i && i<=CMap::dimension);

    if ( i==0 ) return false;
    if ( i==1 ) return true;

    // TODO ? Optimisation possible to not test all the darts of the cell ?
    bool res = true;
    for (CMap_dart_const_iterator_of_cell<CMap,i> it(amap, adart);
         res && it.cont(); ++it)
    {
      if ( it->template beta<i-2>()->template beta<i-1>()!=
           it->template beta<i-1>()->template beta_inv<i-2>() )
        res = false;
    }
    return res;
  }

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
      typename CMap::Dart_handle dg1=NULL, dg2=NULL;

      int mark = amap.get_new_mark();
      int mark_modified_darts = amap.get_new_mark();

      std::deque<typename CMap::Dart_handle> to_erase;

      const int imuinv = CGAL_BETAINV(i-1);

      // First we store and mark all the darts of the i-cell to contract.
      for ( CMap_dart_iterator_basic_of_cell<CMap,i> it(amap,adart,mark);
            it.cont(); ++it )
      {
        to_erase.push_back(it);
        if ( !it->template is_free<i-1>() && dg1==NULL )
        { dg1=it; dg2=it->template beta<i-1>(); }
        amap.mark(it, mark);
        ++res;
      }

      // We group the two (i+1)-cells incident if they exist.
      if ( dg1!=NULL )
         internal::Group_attribute_functor_run<CMap,i-1>::
             run(&amap, dg1, dg2);

      // Second we update the dart of the cell attributes on non marked darts.
      typename std::deque<typename CMap::Dart_handle>::iterator it =
          to_erase.begin();
      for (; it != to_erase.end(); ++it)
        amap.update_dart_of_all_attributes(*it, mark);

      std::deque<typename CMap::Dart_handle> modified_darts;

      // For each dart of the i-cell, we modify i-links of neighbors.
      for ( it=to_erase.begin(); it!=to_erase.end(); ++it )
      {
        d1 = (*it)->template beta<i>();
        while ( d1!=CMap::null_dart_handle && amap.is_marked(d1, mark) )
        {
          d1 = d1->template beta<imuinv>()->template beta<i>();
          if (d1 == (*it)->template beta<i>()) d1 = CMap::null_dart_handle;
        }

        if ( !amap.is_marked(d1, mark_modified_darts) )
        {
          d2 = (*it)->template beta<i-1>()->template beta<i>();
          while ( d2!=CMap::null_dart_handle && amap.is_marked(d2, mark) )
          {
            d2 = d2->template beta<i-1>()->template beta<i>();
            if ( d2==(*it)->template beta<i-1>()->template beta<i>() )
              d2=CMap::null_dart_handle;
          }

          if ( !amap.is_marked(d2, mark_modified_darts) )
          {
            if (d1 != CMap::null_dart_handle)
            {
              if (d2 != CMap::null_dart_handle && d1!=d2 )
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
                if ( !d1->template is_free<i>() )
                {
                  d1->template unlink_beta<i>();
                  CGAL_assertion( !amap.is_marked(d1, mark_modified_darts) );
                  amap.mark(d1, mark_modified_darts);
                  modified_darts.push_back(d1);
                }
              }
            }
            else if (d2 != CMap::null_dart_handle)
            {
              if ( !d2->is_free(i) )
              {
                d2->template unlink_beta<i>();
                CGAL_assertion( !amap.is_marked(d2, mark_modified_darts) );
                amap.mark(d2, mark_modified_darts);
                modified_darts.push_back(d2);
              }
            }
          }
        }
        if ((*it)->is_free(i-1) && !(*it)->is_free(i))
        {
          d1 = (*it)->beta(i);
          if ( !d1->is_free(i) )
          {
            d1->template unlink_beta<i>();
            CGAL_assertion( !amap.is_marked(d1, mark_modified_darts) );
            amap.mark(d1, mark_modified_darts);
            modified_darts.push_back(d1);
          }
        }
      }

      // We test the split of all the incident cells for all the non
      // void attributes.
      CMap::Helper::template Foreach_enabled_attributes_except
          <internal::Test_split_attribute_functor<CMap,i>, i>::
          run(&amap, modified_darts, mark_modified_darts);

      // We remove all the darts of the i-cell.
      for (  it=to_erase.begin(); it!=to_erase.end(); ++it )
      { amap.erase_dart(*it); }

      CGAL_assertion( amap.is_whole_map_unmarked(mark) );
      amap.free_mark(mark);

      if ( !amap.is_whole_map_unmarked(mark_modified_darts) )
      {
        for ( typename std::deque<typename CMap::Dart_handle>::
              iterator it=modified_darts.begin();
              it!=modified_darts.end(); ++it )
          amap.unmark(*it, mark_modified_darts);
      }

      // amap.display_darts(std::cout);

      CGAL_assertion ( amap.is_whole_map_unmarked(mark_modified_darts) );
      amap.free_mark(mark_modified_darts);

      CGAL_expensive_postcondition( amap.is_valid() );
      assert( amap.is_valid() );

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
      typename CMap::Dart_handle dg1=NULL, dg2=NULL;

      int mark = amap.get_new_mark();
//      int mark_modified_darts = amap.get_new_mark();

      std::deque<typename CMap::Dart_handle> to_erase;

      std::deque<typename CMap::Dart_handle> modified_darts;
      std::deque<typename CMap::Dart_handle> modified_darts2;

      // First we store and mark all the darts of the 1-cell to contract.
      for ( CMap_dart_iterator_basic_of_cell<CMap,1> it(amap,adart,mark);
            it.cont(); ++it )
      {
        to_erase.push_back(it);
        if ( dg1==NULL && !it->template is_free<0>() &&
             !it->template is_free<1>() )
        { dg1=it->template beta<0>(); dg2=it->template beta<1>(); }
        amap.mark(it, mark);
        ++res;
      }

      typename std::deque<typename CMap::Dart_handle>::iterator it =
        to_erase.begin();

      for (; it != to_erase.end(); ++it)
        amap.update_dart_of_all_attributes(*it, mark);

      // We group the two vertices incident if they exist.
      if ( dg1!=NULL )
         internal::Group_attribute_functor_run<CMap, 0, 1>::
             run(&amap, dg1, dg2);

      // 4) For each dart of the cell, we modify link of neighbors.
      for ( it=to_erase.begin(); it!=to_erase.end(); ++it )
      {
        if ( !(*it)->template is_free<0>() )
        {
          if ( !(*it)->template is_free<1>() )
          {
            if ( (*it)->template beta<1>()!=*it )
            {
               /*modified_darts2.push_back((*it)->template beta<0>());
              if ( (*it)->beta(0)!=(*it)->beta(1) )*/
              modified_darts.push_back((*it)->template beta<1>());
              amap.basic_link_beta_1((*it)->template beta<0>(),
                                     (*it)->template beta<1>());
            }
          }
          else
          {
            // TODO todegroup.push(Dart_pair((*it)->beta(0), *it));
            modified_darts2.push_back((*it)->template beta<0>());
            (*it)->template beta<0>()->template unlink_beta<1>();
          }
        }
        else
        {
          if ( !(*it)->template is_free<1>() )
          {
            // TODO todegroup.push(Dart_pair((*it)->beta(1), *it));
            modified_darts.push_back((*it)->template beta<1>());
            (*it)->template beta<1>()->template unlink_beta<0>();
          }
        }
      }

      // We test the split of all the incident cells for all the non
      // void attributes.
      CMap::Helper::template Foreach_enabled_attributes_except
          <internal::Test_split_attribute_functor<CMap,1>, 1>::
          run(&amap, modified_darts, modified_darts2);

      // 6) We remove all the darts of the cell.
      for (it = to_erase.begin(); it != to_erase.end(); ++it)
      { amap.erase_dart(*it); }

      CGAL_assertion( amap.is_whole_map_unmarked(mark) );
      amap.free_mark(mark);

      CGAL_expensive_postcondition( amap.is_valid() );
      assert( amap.is_valid() ); // TO REMOVE

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
  { return Contract_cell_functor<CMap,i>::run(amap,adart); }

} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_OPERATIONS_H //
// EOF //
