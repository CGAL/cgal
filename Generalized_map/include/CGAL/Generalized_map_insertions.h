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
#ifndef CGAL_GENERALIZED_MAP_INSERTIONS_H
#define CGAL_GENERALIZED_MAP_INSERTIONS_H

namespace CGAL
{
/** @file Generalized_map_insertions.h
 * Insertion operations on generalized map.
 */

/** Insert a vertex in a given edge.
 * @param amap the used generalized map.
 * @param adart a dart of the edge (!=NULL).
 * @return a dart of the new vertex.
 */
template<class CMap>
typename CMap::Dart_handle
insert_cell_0_in_cell_1( CMap& amap, typename CMap::Dart_handle adart,
                         typename CMap::template
                         Attribute_handle<0>::type ah=CMap::null_handle )
{
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

    if (!amap.template is_free<0>(*it) &&
        amap.is_marked(amap.template alpha<0>(*it), mark))
      amap.template basic_link_alpha<1>(d1, amap.template alpha<0,0>(*it));

    amap.template basic_link_alpha<0>(*it, d1);
    amap.mark(*it, mark);

    for ( unsigned int dim=2; dim<=CMap::dimension; ++dim )
    {
      if (!amap.is_free(*it, dim) && amap.is_marked(amap.alpha(*it, dim), mark))
      {
        amap.basic_link_alpha(amap.beta(*it, dim, 0), d1, dim);
      }
    }

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

  if ( !amap.template is_free<1>(amap.template alpha<0>(adart)) )
  {
    CGAL::internal::Degroup_attribute_functor_run<CMap, 1>::
      run(&amap, adart, amap.template alpha<0,1>(adart));
  }

#ifdef CGAL_CMAP_TEST_VALID_INSERTIONS
  CGAL_assertion( amap.is_valid() );
#endif

  return amap.template alpha<0>(adart);
}

/** Insert a vertex in the given 2-cell which is splitted in triangles,
 * once for each inital edge of the facet.
 * @param amap the used generalized map.
 * @param adart a dart of the facet to triangulate.
 * @return A dart incident to the new vertex.
 */
template < class CMap >
typename CMap::Dart_handle
insert_cell_0_in_cell_2( CMap& amap, typename CMap::Dart_handle adart,
                         typename CMap::template
                         Attribute_handle<0>::type ah=CMap::null_handle )
{
  CGAL_assertion(adart!=amap.null_handle);

  typename CMap::Dart_handle d1=amap.null_handle, d2=amap.null_handle;

  // Mark used to mark darts already treated.
  int treated = amap.get_new_mark();
  int m = amap.get_new_mark();

  // Stack of darts of the face
  std::deque<typename CMap::Dart_handle> vect;
  {
    for ( typename CMap::template Dart_of_cell_basic_range<2>::iterator
            it=amap.template darts_of_cell_basic<2>(adart,m).begin();
          it != amap.template darts_of_cell_basic<2>(adart,m).end(); ++it )
      vect.push_back(it);
  }

  // Stack of darts to degroup (one dart per edge of the face)
  std::deque<typename CMap::Dart_handle> todegroup;
  {
    for ( typename CMap::template Dart_of_cell_basic_range<2,2>::iterator
            it=amap.template darts_of_cell_basic<2,2>(adart).begin();
          it != amap.template darts_of_cell_basic<2,2>(adart).end(); ++it )
      if ( it!=adart && it!=amap.template alpha<0>(adart) )
        todegroup.push_back(it);
  }

  // 2) For each dart of the cell, we modify link of neighbors.
  typename std::deque<typename CMap::Dart_handle>::iterator it = vect.begin();
  for (; it != vect.end(); ++it)
  {
    d1 = amap.create_dart();
    d2 = amap.create_dart();
    amap.template basic_link_alpha<0>(d1, d2);
    amap.mark(*it, treated);

    amap.template basic_link_alpha<1>(*it, d1);

    if (!amap.template is_free<0>(*it) &&
        amap.is_marked(amap.template alpha<0>(*it), treated))
      amap.template basic_link_alpha<1>(d2, amap.template alpha<0,1,0>(*it));

    for ( unsigned int dim=3; dim<=CMap::dimension; ++dim )
    {
      if (!amap.is_free(*it, dim) && amap.is_marked(amap.alpha(*it, dim), treated))
      {
        amap.basic_link_alpha(amap.beta(*it, dim, 1), d1, dim);
        amap.basic_link_alpha(amap.beta(*it, dim, 1, 0), d2, dim);
      }
    }

    // We copy all the attributes except for dim=1
    CMap::Helper::template Foreach_enabled_attributes_except
      <CGAL::internal::Group_attribute_functor_of_dart<CMap>, 1>::
      run(&amap,*it,d1);
    // We initialise the 0-atttrib to ah
    CGAL::internal::Set_i_attribute_of_dart_functor<CMap, 0>::
        run(&amap, d2, ah);
  }

  for (it = vect.begin(); it != vect.end(); ++it)
  {
    amap.unmark(*it, m);
    amap.unmark(*it, treated);
  }

  CGAL_assertion(amap.is_whole_map_unmarked(m));
  CGAL_assertion(amap.is_whole_map_unmarked(treated));
  amap.free_mark(m);
  amap.free_mark(treated);

  for (it = todegroup.begin(); it != todegroup.end(); ++it)
  {
    CGAL::internal::Degroup_attribute_functor_run<CMap, 2>::
      run(&amap, adart, *it);
  }

#ifdef CGAL_CMAP_TEST_VALID_INSERTIONS
  CGAL_assertion( amap.is_valid() );
#endif

  return amap.template alpha<1,0>(adart);
}

/** Test if an edge can be inserted onto a 2-cell between two given darts.
 * @param amap the used generalized map.
 * @param adart1 a first dart.
 * @param adart2 a second dart.
 * @return true iff an edge can be inserted between adart1 and adart2.
 */
template < class CMap >
bool is_insertable_cell_1_in_cell_2(const CMap& amap,
                                    typename CMap::Dart_const_handle adart1,
                                    typename CMap::Dart_const_handle adart2)
{
  if ( adart1==adart2 || adart1==amap.template alpha<0>(adart2) ) return false;
  for ( CGAL::CMap_dart_const_iterator_of_orbit<CMap,0,1> it(amap,adart1);
        it.cont(); ++it )
  {
    if ( it==adart2 )  return true;
  }
  return false;
}

/** Insert an edge in a 2-cell between two given darts.
 * @param amap the used generalized map.
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
  if ( adart2!=amap.null_handle)
  {
    CGAL_assertion(is_insertable_cell_1_in_cell_2<CMap>(amap, adart1, adart2));
  }

  int m1=amap.get_new_mark();
  CGAL::CMap_dart_iterator_basic_of_involution<CMap,1> it1(amap, adart1, m1);

  int m2=amap.get_new_mark();
  CGAL::CMap_dart_iterator_basic_of_involution<CMap,1> it2(amap, adart2, m2);

  typename CMap::Dart_handle d1=amap.null_handle;
  typename CMap::Dart_handle d2=amap.null_handle;

  int treated=amap.get_new_mark();

  for ( ; it1.cont(); ++it1)
  {
    CGAL_assertion (it2.cont() );
    d1 = amap.create_dart();
    d2 = amap.create_dart();
    amap.template basic_link_alpha<0>(d1, d2);
    amap.mark(it1,treated);

    if ( !amap.template is_free<1>(it1) &&
         amap.is_marked(amap.template alpha<1>(it1), treated) )
    {
      amap.template basic_link_alpha<2>(amap.template alpha<1,1>(it1), d1);
      amap.template basic_link_alpha<2>(amap.template alpha<2,0>(d1), d2);
    }

    for ( unsigned int dim=3; dim<=CMap::dimension; ++dim)
    {
      if ( !amap.is_free(it1, dim) &&
           amap.is_marked(amap.alpha(it1, dim), treated) )
      {
        amap.basic_link_alpha(amap.alpha(it1, dim, 1), d1, dim);
        amap.basic_link_alpha(amap.alpha(d1, dim, 0), d2, dim);
      }
    }

    amap.template link_alpha<1>(it1, d1);
    if ( adart2!=amap.null_handle )
    {
      amap.template link_alpha<1>(it2, d2);
      ++it2;
    }
  }

  if ( !amap.template is_free<2>(d1) && d2!=amap.null_handle )
    CGAL::internal::Degroup_attribute_functor_run<CMap, 2>::
      run(&amap, d1, amap.template alpha<2>(d1));

  amap.negate_mark(m1);
  it1.rewind();

  if ( adart2!=amap.null_handle )
  { it2.rewind(); amap.negate_mark(m2); }

  for ( ; it1.cont(); ++it1, ++it2)
  {
    amap.mark(it1,m1);
    amap.unmark(it1,treated);
    if ( adart2!=amap.null_handle ) amap.mark(it2,m2);
  }
  amap.negate_mark(m1);
  CGAL_assertion( amap.is_whole_map_unmarked(m1) );
  CGAL_assertion( amap.is_whole_map_unmarked(treated) );
  amap.free_mark(m1);
  amap.free_mark(treated);

  if ( adart2!=amap.null_handle )
  {
    amap.negate_mark(m2);
    CGAL_assertion( amap.is_whole_map_unmarked(m2) );
  }
  amap.free_mark(m2);

#ifdef CGAL_CMAP_TEST_VALID_INSERTIONS
  CGAL_assertion( amap.is_valid() );
#endif

  return amap.template alpha<1>(adart1);
}

/** Test if a 2-cell can be inserted onto a given 3-cell along
 * a path of edges.
 * @param amap the used generalized map.
 * @param afirst iterator on the begining of the path.
 * @param alast  iterator on the end of the path.
 * @return true iff a 2-cell can be inserted along the path.
 * the path is a sequence of dartd, one per edge
 * where the face will be inserted.
 */
template <class CMap, class InputIterator>
bool is_insertable_cell_2_in_cell_3(const CMap& amap,
                                    InputIterator afirst,
                                    InputIterator alast)
{
  CGAL_assertion( CMap::dimension>= 3 );

  // The path must have at least one dart.
  if (afirst==alast) return false;
  typename CMap::Dart_const_handle prec = amap.null_handle;
  typename CMap::Dart_const_handle od = amap.null_handle;

  for (InputIterator it(afirst); it!=alast; ++it)
  {
    // The path must contain only non empty darts.
    if (*it == amap.null_handle) return false;

    // Two consecutive darts of the path must belong to two edges
    // incident to the same vertex of the same volume.
    if (prec != amap.null_handle)
    {
      if ( amap.template is_free<0>(prec) ) return false;

      // alpha0(prec) and *it must belong to the same vertex of the same volume
      if ( !CGAL::belong_to_same_cell<CMap, 0, 2>
           (amap,  amap.template alpha<0>(prec), *it) )
        return false;
    }
    prec = *it;
  }

  // The path must be closed.
  if ( amap.template is_free<0>(prec) ) return false;
  if (!CGAL::belong_to_same_cell<CMap, 0, 2>
      (amap, amap.template alpha<0>(prec), *afirst))
    return false;

  return true;
}

/** Insert a 2-cell in a given 3-cell along a path of darts.
 * @param amap the used generalized map.
 * @param afirst iterator on the begining of the path.
 * @param alast  iterator on the end of the path.
 * the path is a sequence of dartd, one per edge
 * where the face will be inserted.
 * @return a dart of the new 2-cell.
 */
template<class CMap, class InputIterator>
typename CMap::Dart_handle
insert_cell_2_in_cell_3(CMap& amap, InputIterator afirst, InputIterator alast)
{
  CGAL_assertion(is_insertable_cell_2_in_cell_3(amap,afirst,alast));

  typename CMap::Dart_handle prec = amap.null_handle, d = amap.null_handle,
    dd = amap.null_handle, first = amap.null_handle;
  bool withAlpha3 = false;

  int treated = amap.get_new_mark();

  {
    for (InputIterator it(afirst); !withAlpha3 && it!=alast; ++it)
    {
      if (!amap.template is_free<2>(*it)) withAlpha3 = true;
    }
  }

  {
    for (InputIterator it(afirst); it!=alast; ++it)
    {
      d = amap.create_dart();
      amap.template basic_link_alpha<0>(d, amap.create_dart());

      if ( withAlpha3 )
      {
        dd = amap.create_dart();
        amap.template basic_link_alpha<0>(dd, amap.create_dart());
      }

      if ( prec!=amap.null_handle )
      {
        amap.template basic_link_alpha<1>(prec, d);
        if (withAlpha3)
          amap.template basic_link_alpha<1>(amap.template alpha<3>(prec), dd);
      }
      else first = amap.template alpha<0>(d);

      if ( !amap.template is_free<2>(*it) )
      {
        amap.template link_alpha<2>(amap.template alpha<2>(*it), dd);
      }

      amap.template link_alpha<2>(*it, d);
      if (withAlpha3) amap.template basic_link_alpha<3>(d, dd);

      prec = amap.template alpha<0>(d);
    }
  }

  amap.template basic_link_alpha<1>(prec, first);
  if ( withAlpha3 )
  {
    amap.template basic_link_alpha<1>(amap.template alpha<3>(prec),
                                      amap.template alpha<3>(first));
  }

  // Make copies of the new facet for dimension >=4
  /*  for ( unsigned int dim=4; dim<=CMap::dimension; ++dim )
  {
    if ( !amap.is_free(*it, dim) )
    {
          ddd = amap.create_dart();
          amap.template basic_link_alpha<0>(ddd, amap.create_dart());
          amap.basic_link_alpha(d, ddd, dim);
          amap.basic_link_alpha(amap.template alpha<0>(d),
                                amap.template alpha<0>(ddd), dim);

          if ( withAlpha3 )
          {
            dddd = amap.create_dart();
            amap.template basic_link_alpha<0>(dddd, amap.create_dart());

            amap.basic_link_alpha(dd, dddd, dim);
            amap.basic_link_alpha(amap.template alpha<0>(dd),
                                  amap.template alpha<0>(dddd), dim);
          }



        }
        }*/

  // Make copies of the new facet for dimension >=4
  for ( unsigned int dim=4; dim<=CMap::dimension; ++dim )
  {
    if ( !amap.is_free(first, dim) )
    {
      typename CMap::Dart_handle first2 = amap.null_handle;
      prec = amap.null_handle;
      for ( CMap_dart_iterator_basic_of_orbit<CMap,0,1> it(amap, first);
            it.cont(); ++it )
      {
        d = amap.create_dart();
        amap.basic_link_alpha(amap.template alpha<2>(it), d, dim);
        if ( withAlpha3 )
        {
          dd = amap.create_dart();
          amap.basic_link_alpha_for_involution
            (amap.template alpha<2,3>(it), dd, dim);
          amap.template basic_link_alpha_for_involution<3>(d, dd);
        }
        if ( prec!=amap.null_handle )
        {
          amap.link_alpha_0(prec, d);
          if ( withAlpha3 )
          {
            amap.basic_link_alpha_1(amap.template alpha<3>(prec), dd);
          }
        }
        else first2 = prec;

        // We consider dim2=2 out of the loop to use link_alpha instead of
        // basic _link_alpha (to modify non void attributes only once).
        if ( !amap.template is_free<2>(it) &&
             amap.is_free(amap.template alpha<2>(it), dim) )
          amap.template link_alpha_for_involution<2>
            (amap.alpha(it,2,dim), d);
        if ( withAlpha3 &&
             !amap.template is_free<2>(amap.template alpha<3>(it)) &&
             amap.is_free(amap.template alpha<3,2>(it), dim) )
          amap.template link_alpha_for_involution<2>(amap.alpha(it,3,2,dim), dd);

        for ( unsigned int dim2=3; dim2<=CMap::dimension; ++dim2 )
        {
          if ( dim2+1!=dim && dim2!=dim && dim2!=dim+1 )
          {
            if ( !amap.is_free(it, dim2) && amap.is_free(amap.alpha(it, dim2), dim) )
              amap.basic_link_alpha_for_involution(amap.alpha(it, dim2, dim),
                                                  d, dim2);
            if ( withAlpha3 && !amap.is_free(amap.template alpha<3>(it), dim2) &&
                 amap.is_free(amap.alpha(it, 3, dim2), dim) )
              amap.basic_link_alpha_for_involution
                (amap.alpha(it, 3, dim2, dim), dd, dim2);
          }
        }
        prec = d;
      }
      amap.basic_link_alpha_0( prec, first2 );
      if ( withAlpha3 )
      {
        amap.basic_link_alpha_1( amap.template alpha<3>(prec),
                                amap.template alpha<3>(first2) );
      }
    }
  }

  // Degroup the attributes
  if ( withAlpha3 )
  { // Here we cannot use Degroup_attribute_functor_run as new darts do not
    // have their 3-attribute
    CGAL::internal::Degroup_attribute_functor_run<CMap, 3>::
        run(&amap, first, amap.template alpha<3>(first));
  }

#ifdef CGAL_CMAP_TEST_VALID_INSERTIONS
  CGAL_assertion( amap.is_valid() );
#endif

  return first;
}

} // namespace CGAL

#endif // CGAL_GENERALIZED_MAP_INSERTIONS_H
