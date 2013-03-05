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
#ifndef CGAL_COMBINATORIAL_MAP_INSERTIONS_H
#define CGAL_COMBINATORIAL_MAP_INSERTIONS_H

namespace CGAL
{
/** @file Combinatorial_map_insertions.h
 * Insertion operations on combinatorial map.
 */

/** Insert a vertex in a given edge.
 * @param amap the used combinatorial map.
 * @param adart a dart of the edge (!=NULL && !=null_dart_handle).
 * @return a dart of the new vertex.
 */
template<class CMap>
typename CMap::Dart_handle
insert_cell_0_in_cell_1( CMap& amap, typename CMap::Dart_handle adart,
                         typename CMap::template
                         Attribute_handle<0>::type ah=NULL )
{
  CGAL_assertion(adart != NULL && adart!=CMap::null_dart_handle);
  typename CMap::Dart_handle d1, d2;
  int mark=amap.get_new_mark();

  // 1) We store all the darts of the edge.
  std::deque<typename CMap::Dart_handle> vect;
  int m=amap.get_new_mark();
  {
    for ( typename CMap::template Dart_of_cell_basic_range<1>::iterator
          it=amap.template darts_of_cell_basic<1>(adart, m).begin();
          it != amap.template darts_of_cell_basic<1>(adart, m).end(); ++it )
      vect.push_back(it);
  }

  // 2) For each dart of the cell, we modify link of neighbors.
  typename std::deque<typename CMap::Dart_handle>::iterator it = vect.begin();
  for (; it != vect.end(); ++it)
  {
    d1 = amap.create_dart();

    if (!(*it)->is_free(1))
    { amap.basic_link_beta_1(d1, (*it)->beta(1)); }

    for ( unsigned int dim=2; dim<=CMap::dimension; ++dim )
    {
      if (!(*it)->is_free(dim) && amap.is_marked((*it)->beta(dim), mark))
      {
        amap.basic_link_beta_for_involution((*it)->beta(dim), d1, dim);
        amap.basic_link_beta_for_involution(*it, (*it)->beta(dim)->
                                            template beta<1>(), dim);
      }
    }

    amap.basic_link_beta_1(*it, d1);

    // We copy all the attributes except for dim=0
    CMap::Helper::template Foreach_enabled_attributes_except
      <CGAL::internal::Group_attribute_functor_of_dart<CMap>, 0>::
      run(&amap,*it,d1);
    // We initialise the 0-atttrib to ah
    CGAL::internal::Set_i_attribute_of_dart_functor<CMap, 0>::
        run(&amap, d1, ah);
    amap.mark(*it, mark);
  }

  for (it = vect.begin(); it != vect.end(); ++it)
  {
    amap.unmark(*it, m);
    amap.unmark(*it, mark);
  }

  CGAL_assertion(amap.is_whole_map_unmarked(m));
  CGAL_assertion(amap.is_whole_map_unmarked(mark));

  amap.free_mark(m);
  amap.free_mark(mark);

  CGAL::internal::Degroup_attribute_functor_run<CMap, 1>::
      run(&amap, adart, adart->beta(1));

#ifdef CGAL_CMAP_TEST_VALID_INSERTIONS
  CGAL_assertion( amap.is_valid() );
#endif

  return adart->beta(1);
}

/** Insert a vertex in the given 2-cell which is splitted in triangles,
 * once for each inital edge of the facet.
 * @param amap the used combinatorial map.
 * @param adart a dart of the facet to triangulate.
 * @return A dart incident to the new vertex.
 */
template < class CMap >
typename CMap::Dart_handle
insert_cell_0_in_cell_2( CMap& amap, typename CMap::Dart_handle adart,
                         typename CMap::template
                         Attribute_handle<0>::type ah=NULL )
{
  CGAL_assertion(adart != NULL && adart!=CMap::null_dart_handle);
  typename CMap::Dart_handle first = adart, prev = NULL, cur = NULL;
  typename CMap::Dart_handle n1 = NULL, n2 = NULL;

  typename CMap::Dart_handle nn1 = NULL, nn2 = NULL;

  // If the facet is open, we search the dart 0-free
  while ( !first->template is_free<0>() && first->template beta<0>()!=adart )
    first = first->template beta<0>();

  // Mark used to mark darts already treated.
  int treated = amap.get_new_mark();

  // Stack of marked darts
  std::deque<typename CMap::Dart_handle> tounmark;

  // Now we run through the facet
  for ( CGAL::CMap_dart_iterator_basic_of_orbit<CMap, 1> it(amap, first);
        it.cont(); )
  {
    cur = it;
    ++it;
    amap.mark(cur, treated);
    tounmark.push_back(cur);

    if (!cur->template is_free<0>())
    {
      n1=amap.create_dart();
      amap.link_beta_0(cur, n1);
    }
    else n1 = NULL;

    if (!cur->template is_free<1>())
    {
      n2 = amap.create_dart();
      amap.link_beta_1(cur, n2);
    }
    else n2 = NULL;

    if ( n1!=NULL )
    {
      if ( n2!=NULL )
        amap.basic_link_beta_0(n1, n2);

      if ( prev!=NULL )
        amap.template basic_link_beta_for_involution<2>(prev, n1);

      CGAL::internal::Set_i_attribute_of_dart_functor<CMap, 0>::
          run(&amap, n1, ah);
    }

    for (unsigned int dim=3; dim<=CMap::dimension; ++dim)
    {
      if ( !adart->is_free(dim) )
      {
        if ( !amap.is_marked(cur->beta(dim), treated) )
        {
          if (n1!=NULL)
          {
            nn1=amap.create_dart();
            amap.link_beta_1(cur->beta(dim), nn1);
            amap.basic_link_beta_for_involution(n1, nn1, dim);
          }
          else nn1=NULL;

          if (n2!=NULL)
          {
            nn2=amap.create_dart();
            amap.link_beta_0(cur->beta(dim), nn2);
            amap.basic_link_beta_for_involution(n2, nn2, dim);
            CGAL::internal::Set_i_attribute_of_dart_functor<CMap, 0>::
                run(&amap, nn2, ah);
          }
          else nn2=NULL;

          if (nn1 != NULL && nn2 != NULL)
            amap.basic_link_beta_1(nn1, nn2);

          if (nn1 != NULL && prev != NULL)
            amap.template basic_link_beta_for_involution<2>
                (nn1, prev->beta(dim));

          amap.mark(cur->beta(dim), treated);
        }
        else
        {
          if ( n1!=NULL )
            amap.basic_link_beta_for_involution
                (n1, cur->beta(dim)->template beta<1>(), dim);
          if ( n2!=NULL )
            amap.basic_link_beta_for_involution
                (n2, cur->beta(dim)->template beta<0>(), dim);
        }
      }
    }

    prev = n2;
  }

  if (n2 != NULL)
  {
    amap.template basic_link_beta_for_involution<2>
        (first->template beta<0>(), n2);
    for (unsigned int dim=3; dim<=CMap::dimension; ++dim)
    {
      if ( !adart->is_free(dim) )
      {
        amap.template basic_link_beta_for_involution<2>
            (first->template beta<0>()->beta(dim), n2->beta(dim));
      }
    }
  }

  // Now we unmark all marked darts, and we degroup the new faces with the
  // initial one (if 2-attributes are non void).
  for ( typename std::deque<typename CMap::Dart_handle>::iterator
        itd=tounmark.begin(); itd!=tounmark.end(); ++itd )
  {
    amap.unmark(*itd, treated);
    for (unsigned int dim=3; dim<=CMap::dimension; ++dim)
    {
      if ( !(*itd)->is_free(dim) )
        amap.unmark((*itd)->beta(dim), treated);
    }
    if ( *itd!=adart )
      CGAL::internal::Degroup_attribute_functor_run<CMap, 2>::
          run(&amap, adart, *itd);
  }

  CGAL_assertion(amap.is_whole_map_unmarked(treated));
  amap.free_mark(treated);

#ifdef CGAL_CMAP_TEST_VALID_INSERTIONS
  CGAL_assertion( amap.is_valid() );
#endif

  return n1;
}
/** Insert a dangling edge in a 2-cell between given by a dart.
 * @param amap the used combinatorial map.
 * @param adart1 a first dart of the facet (!=NULL && !=null_dart_handle).
 * @return a dart of the new edge, not incident to the vertex of adart1.
 */
template<class CMap>
typename CMap::Dart_handle
insert_dangling_cell_1_in_cell_2( CMap& amap,
                                  typename CMap::Dart_handle adart1,
                                  typename CMap::template
                                  Attribute_handle<0>::type ah=NULL )
{
  CGAL_assertion(adart1!=NULL && adart1!=CMap::null_dart_handle);

  int mark1 = amap.get_new_mark();
  std::deque<typename CMap::Dart_handle> to_unmark;
  {
    for ( CMap_dart_iterator_basic_of_cell<CMap,0> it(amap,adart1,mark1);
          it.cont(); ++it )
    {
      to_unmark.push_back(it);
      amap.mark(it,mark1);
    }
  }

  typename CMap::Dart_handle d1 = NULL;
  typename CMap::Dart_handle d2 = NULL;
  unsigned int s1 = 0;

  int treated=amap.get_new_mark();

  CGAL::CMap_dart_iterator_basic_of_involution<CMap,1>
      it1(amap, adart1, treated);

  for ( ; it1.cont(); ++it1)
  {
    d1 = amap.create_dart();
    d2 = amap.create_dart();

    if ( amap.is_marked(it1, mark1) ) s1 = 0;
    else s1 = 1;

    if ( !it1->is_free(s1) )
    {
      if ( s1==0 )
        amap.link_beta_1(it1->template beta<0>(), d2);
      else
        amap.link_beta_0(it1->template beta<1>(), d2);
    }

    if (s1==0)
    {
      amap.link_beta_0(it1, d1);
      amap.link_beta_0(d1, d2);
    }
    else
    {
      amap.link_beta_1(it1, d1);
      amap.link_beta_1(d1, d2);
    }

    amap.template basic_link_beta_for_involution<2>(d1, d2);

    for ( unsigned int dim=3; dim<=CMap::dimension; ++dim)
    {
      if ( !it1->is_free(dim) &&
           amap.is_marked(it1->beta(dim), treated) )
      {
        amap.basic_link_beta_for_involution
            (it1->beta(dim)->beta_inv(s1), d1, dim);
        amap.basic_link_beta_for_involution
          (it1->beta(dim)->beta_inv(s1)->beta(2), d2, dim);
      }
    }
    CGAL::internal::Set_i_attribute_of_dart_functor<CMap, 0>::
        run(&amap, d1, ah);
    amap.mark(it1, treated);
  }

  amap.negate_mark(treated);
  for ( it1.rewind(); it1.cont(); ++it1 )
  { amap.mark(it1, treated); }

  CGAL_assertion( amap.is_whole_map_marked(treated) );
  amap.free_mark(treated);

  for ( typename std::deque<typename CMap::Dart_handle>::iterator
        it=to_unmark.begin(); it!=to_unmark.end(); ++it)
  { amap.unmark(*it, mark1); }

  CGAL_assertion( amap.is_whole_map_unmarked(mark1) );
  amap.free_mark(mark1);

#ifdef CGAL_CMAP_TEST_VALID_INSERTIONS
  CGAL_assertion( amap.is_valid() );
#endif

  return adart1->template beta<0>();
}

/** Test if an edge can be inserted onto a 2-cell between two given darts.
 * @param amap the used combinatorial map.
 * @param adart1 a first dart.
 * @param adart2 a second dart.
 * @return true iff an edge can be inserted between adart1 and adart2.
 */
template < class CMap >
bool is_insertable_cell_1_in_cell_2(const CMap& amap,
                                    typename CMap::Dart_const_handle adart1,
                                    typename CMap::Dart_const_handle adart2)
{
  CGAL_assertion(adart1 != NULL && adart2 != NULL);
  if ( adart1==adart2 ) return false;
  for ( CGAL::CMap_dart_const_iterator_of_orbit<CMap,1> it(amap,adart1);
        it.cont(); ++it )
  {
    if ( it==adart2 )  return true;
  }
  return false;
}

/** Insert an edge in a 2-cell between two given darts.
 * @param amap the used combinatorial map.
 * @param adart1 a first dart of the facet (!=NULL && !=null_dart_handle).
 * @param adart2 a second dart of the facet. If NULL insert a dangling edge.
 * @return a dart of the new edge, and not incident to the
 *         same vertex than adart1.
 */
template<class CMap>
typename CMap::Dart_handle
insert_cell_1_in_cell_2(CMap& amap,
                        typename CMap::Dart_handle adart1,
                        typename CMap::Dart_handle adart2)
{
  if ( adart2==NULL ) return insert_dangling_cell_1_in_cell_2(amap,adart1);

  CGAL_assertion(is_insertable_cell_1_in_cell_2<CMap>(amap, adart1, adart2));

  int m1=amap.get_new_mark();
  CGAL::CMap_dart_iterator_basic_of_involution<CMap,1> it1(amap, adart1, m1);

  int m2=amap.get_new_mark();
  CGAL::CMap_dart_iterator_basic_of_involution<CMap,1> it2(amap, adart2, m2);

  int mark1=amap.get_new_mark();
  std::deque<typename CMap::Dart_handle> to_unmark;
  {
    for ( CGAL::CMap_dart_iterator_basic_of_cell<CMap,0> it(amap,adart1,mark1);
          it.cont(); ++it )
    {
      to_unmark.push_back(it);
      amap.mark(it, mark1);
    }
  }

  typename CMap::Dart_handle d1=NULL;
  typename CMap::Dart_handle d2=NULL;
  unsigned int s1=0;

  int treated=amap.get_new_mark();

  for ( ; it1.cont(); ++it1, ++it2)
  {
    CGAL_assertion (it2.cont() );
    d1 = amap.create_dart();
    d2 = amap.create_dart();

    if ( amap.is_marked(it1, mark1) ) s1 = 0;
    else s1 = 1;

    if ( !it1->is_free(s1) )
    {
      if ( s1==0 ) amap.link_beta_1(it1->template beta<0>(), d2);
      else amap.link_beta_0(it1->template beta<1>(), d2);
    }

    if ( !it2->is_free(s1) )
    {
      if ( s1==0 ) amap.link_beta_1(it2->template beta<0>(), d1);
      else amap.link_beta_0(it2->template beta<1>(), d1);
    }

    if ( s1==0 )
    {
      amap.link_beta_0(it1, d1);
      amap.link_beta_0(it2, d2);
    }
    else
    {
      amap.link_beta_1(it1, d1);
      amap.link_beta_1(it2, d2);
    }
    amap.template basic_link_beta_for_involution<2>(d2, d1);

    for ( unsigned int dim=3; dim<=CMap::dimension; ++dim)
    {
      if ( !it1->is_free(dim) &&
           amap.is_marked(it1->beta(dim), treated) )
      {
        amap.basic_link_beta_for_involution
          (it1->beta(dim)->beta_inv(s1), d1, dim);
        amap.basic_link_beta_for_involution
          (it1->beta(dim)->beta_inv(s1)->beta(2), d2, dim);
      }
    }

    amap.mark(it1,treated);
  }

  CGAL::internal::Degroup_attribute_functor_run<CMap, 2>::run(&amap, d1, d2);

  amap.negate_mark(m1);
  amap.negate_mark(m2);
  it1.rewind(); it2.rewind();
  for ( ; it1.cont(); ++it1, ++it2)
  {
    amap.mark(it1,m1);
    amap.unmark(it1,treated);
    amap.mark(it2,m2);
  }
  amap.negate_mark(m1);
  amap.negate_mark(m2);
  CGAL_assertion( amap.is_whole_map_unmarked(m1) );
  CGAL_assertion( amap.is_whole_map_unmarked(m2) );
  CGAL_assertion( amap.is_whole_map_unmarked(treated) );
  amap.free_mark(m1);
  amap.free_mark(m2);
  amap.free_mark(treated);

  typename std::deque<typename CMap::Dart_handle>::iterator it =
    to_unmark.begin();
  for (; it != to_unmark.end(); ++it)
  { amap.unmark(*it, mark1); }
  CGAL_assertion( amap.is_whole_map_unmarked(mark1) );
  amap.free_mark(mark1);

#ifdef CGAL_CMAP_TEST_VALID_INSERTIONS
  CGAL_assertion( amap.is_valid() );
#endif

  return adart1->template beta<0>();
}

/** Test if a 2-cell can be inserted onto a given 3-cell along
 * a path of edges.
 * @param amap the used combinatorial map.
 * @param afirst iterator on the begining of the path.
 * @param alast  iterator on the end of the path.
 * @return true iff a 2-cell can be inserted along the path.
 */
template <class CMap, class InputIterator>
bool is_insertable_cell_2_in_cell_3(const CMap& amap,
                                    InputIterator afirst,
                                    InputIterator alast)
{
  CGAL_assertion( CMap::dimension>= 3 );

  // The path must have at least one dart.
  if (afirst==alast) return false;
  typename CMap::Dart_const_handle prec = NULL;
  typename CMap::Dart_const_handle od = NULL;

  for (InputIterator it(afirst); it!=alast; ++it)
  {
    // The path must contain only non empty darts.
    if (*it == NULL || *it==CMap::null_dart_handle) return false;

    // Two consecutive darts of the path must belong to two edges
    // incident to the same vertex of the same volume.
    if (prec != NULL)
    {
      od = prec->other_extremity();
      if ( od==CMap::null_dart_handle ) return false;

      // of and *it must belong to the same vertex of the same volume
      if ( !CGAL::belong_to_same_cell<CMap, 0, 2>(amap, od, *it) )
        return false;
    }
    prec = *it;
  }

  // The path must be closed.
  od = prec->other_extremity();
  if ( od==CMap::null_dart_handle ) return false;

  if (!CGAL::belong_to_same_cell<CMap, 0, 2>(amap, od, *afirst))
    return false;

  return true;
}

/** Insert a 2-cell in a given 3-cell along a path of darts.
 * @param amap the used combinatorial map.
 * @param afirst iterator on the begining of the path.
 * @param alast  iterator on the end of the path.
 * @return a dart of the new 2-cell.
 */
template<class CMap, class InputIterator>
typename CMap::Dart_handle
insert_cell_2_in_cell_3(CMap& amap, InputIterator afirst, InputIterator alast)
{
  CGAL_assertion(is_insertable_cell_2_in_cell_3(amap,afirst,alast));

  typename CMap::Dart_handle prec = NULL, d = NULL, dd = NULL, first = NULL;
  bool withBeta3 = false;

  {
    for (InputIterator it(afirst); !withBeta3 && it!=alast; ++it)
    {
      if (!(*it)->template is_free<2>()) withBeta3 = true;
    }
  }

  {
    for (InputIterator it(afirst); it!=alast; ++it)
    {
      d = amap.create_dart();
      if ( withBeta3 )
        dd = amap.create_dart();

      if (prec != NULL)
      {
        amap.basic_link_beta_0(prec, d);
        if (withBeta3)
          amap.basic_link_beta_1(prec->template beta<3>(), dd);
      }
      else first = d;

      if ( !(*it)->template is_free<2>() )
        amap.template basic_link_beta_for_involution<2>
            ((*it)->template beta<2>(), dd);

      amap.template link_beta_for_involution<2>(*it, d);
      if ( withBeta3 )
        amap.template link_beta_for_involution<3>(d, dd);

      prec = d;
    }
  }

  amap.basic_link_beta_0(prec, first);
  if ( withBeta3 )
  {
    amap.basic_link_beta_1(prec->template beta<3>(),
                           first->template beta<3>());
  }

  // Make copies of the new facet for dimension >=4
  for ( unsigned int dim=4; dim<=CMap::dimension; ++dim )
  {
    if ( !first->is_free(dim) )
    {
      typename CMap::Dart_handle first2 = NULL;
      prec = NULL;
      for ( CMap_dart_iterator_basic_of_orbit<CMap, 1> it(amap, first);
            it.cont(); ++it )
      {
        d = amap.create_dart();
        amap.basic_link_beta_for_involution(it->template beta<2>(), d, dim);
        if ( withBeta3 )
        {
          dd = amap.create_dart();
          amap.basic_link_beta_for_involution
              (it->template beta<2>()->template beta<3>(), dd, dim);
          amap.template basic_link_beta_for_involution<3>(d, dd);
        }
        if ( prec!=NULL )
        {
          amap.link_beta_0(prec, d);
          if ( withBeta3 )
          {
            amap.basic_link_beta_1(prec->template beta<3>(), dd);
          }
        }
        else first2 = prec;

        // We consider dim2=2 out of the loop to use link_beta instead of
        // basic _link_beta (to modify non void attributes only once).
        if ( !it->template is_free<2>() &&
             it->template beta<2>()->is_free(dim) )
          amap.template link_beta_for_involution<2>
              (it->template beta<2>()->beta(dim), d);
        if ( withBeta3 &&
             !it->template beta<3>()->template is_free<2>() &&
             it->template beta<3>()->template beta<2>()->is_free(dim) )
          amap.template link_beta_for_involution<2>
            (it->template beta<3>()->template beta<2>()->beta(dim), dd);

        for ( unsigned dim2=3; dim2<=CMap::dimension; ++dim2 )
        {
          if ( dim2+1!=dim && dim2!=dim && dim2!=dim+1 )
          {
            if ( !it->is_free(dim2) && it->beta(dim2)->is_free(dim) )
              amap.basic_link_beta_for_involution(it->beta(dim2)->beta(dim),
                                                  d, dim2);
            if ( withBeta3 && !it->template beta<3>()->is_free(dim2) &&
                 it->template beta<3>()->beta(dim2)->is_free(dim) )
              amap.basic_link_beta_for_involution
                (it->template beta<3>()->beta(dim2)->beta(dim), dd, dim2);
          }
        }
        prec = d;
      }
      amap.basic_link_beta_0( prec, first2 );
      if ( withBeta3 )
      {
        amap.basic_link_beta_1( prec->template beta<3>(),
                                first2->template beta<3>() );
      }
    }
  }

  // Degroup the attributes
  if ( withBeta3 )
  { // Here we cannot use Degroup_attribute_functor_run as new darts do not
    // have their 3-attribute
    CGAL::internal::Degroup_attribute_functor_run<CMap, 3>::
        run(&amap, first, first->template beta<3>());
  }

#ifdef CGAL_CMAP_TEST_VALID_INSERTIONS
  CGAL_assertion( amap.is_valid() );
#endif

  return first;
}

} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_INSERTIONS_H
