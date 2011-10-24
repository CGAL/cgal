// Copyright (c) 2010 CNRS, LIRIS, http://liris.cnrs.fr/, All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
#ifndef CGAL_LINEAR_CELL_COMPLEX_OPERATIONS_H
#define CGAL_LINEAR_CELL_COMPLEX_OPERATIONS_H 1

#include <CGAL/Cell_iterators.h>
#include <CGAL/Combinatorial_map_operations.h>
#include <vector>

namespace CGAL {

/** @file Combinatorial_map_with_embedding_operations.h
 * Basic operators to modify an embedded combinatorial map.
 * basic operators to modify an embedded combinatorial map: insert a
 * vertex in a facet,  insertion or bursting of a cell.
 */


/** Compute the normal of the given facet.
 * @param amap the used combinatorial map.
 * @param adart a dart incident to the facet.
 * @return the normal of the facet.
 */
template <class LCC>
typename LCC::Vector compute_normal_of_cell_2
(const LCC& amap, typename LCC::Dart_const_handle adart)
{
  // TODO Better approximation by using Newell's method
  // Nx += (Vy - V'y) * (Vz + V'z);
  // Ny += (Vz - V'z) * (Vx + V'x);
  // Nz += (Vx - V'x) * (Vy + V'y);
  // But problem with functor since this is not the sum of normal vectors.
  
  typedef typename LCC::Point Point;
  typedef typename LCC::Vector Vector;
  typename LCC::Dart_const_handle start=adart;
  Vector normal(CGAL::NULL_VECTOR);

  while ( !start->is_free(0) && start->beta(0)!=adart )
    start = start->beta(0);

  if ( start->is_free(1) || start->beta(1)->other_extremity()==NULL )
    return normal;

  unsigned int nb = 0;
  adart = start->beta(1);
  
  const Point* prev = &LCC::point(start);
  const Point* curr = &LCC::point(adart);
  for ( ; adart!=start && adart->other_extremity()!=NULL; adart=adart->beta(1) )
  {
    const Point* next = &LCC::point(adart->other_extremity());
    if ( !typename LCC::Traits::Collinear_3()(*prev, *curr, *next) )
    {
      normal = typename LCC::Traits::Construct_sum_of_vectors()
        (normal, typename LCC::Traits::Construct_normal_3()(*prev, *curr, *next));
      prev = curr;
      ++nb;
    }
    curr = next;
  }
  
  if ( nb<2 ) return normal;       
  return (typename LCC::Traits::Construct_scaled_vector()(normal, 1.0/nb));  
  //  return normal / std::sqrt(normal * normal);
}

/** Compute the normal of the given vertex.
 * @param amap the used combinatorial map.
 * @param adart a dart incident to the vertex.
 * @return the normal of the vertex.
 */
template <class LCC>
typename LCC::Vector compute_normal_of_cell_0
(const LCC& amap, typename LCC::Dart_const_handle adart)
{
  typedef typename LCC::Point Point;
  typedef typename LCC::Vector Vector;
  Vector normal(CGAL::NULL_VECTOR);
  unsigned int nb = 0;
  
  for ( CMap_one_dart_per_incident_cell_const_iterator<LCC,2,0> it(amap, adart); 
        it.cont(); ++it)
  {
    normal = typename LCC::Traits::Construct_sum_of_vectors()
      (normal, CGAL::compute_normal_of_cell_2(amap,it));
    ++nb;
  }
  
  if ( nb<2 ) return normal;       
  return (typename LCC::Traits::Construct_scaled_vector()(normal, 1.0/nb));
  //  return normal / std::sqrt(normal * normal);
}

/** Compute the dual of a combinatorial map.
 * @param amap1 the initial map.
 * @param amap2 the map in which we build the dual of amap1.
 * @param adart a dart of the initial map, NULL by default.
 * @return adart of the dual map, the dual of adart if adart!=NULL.
 */
template<class Map>
typename Map::Dart_handle dual(Map& amap1, Map& amap2, 
                               typename Map::Dart_handle adart=NULL)
{
  CGAL_assertion( amap1.is_without_boundary(Map::dimension) );

  typedef typename Map::Dart_handle Dart_handle;
  typedef typename Map::Dart_range::iterator Dart_iterator;

  std::map< Dart_handle, Dart_handle > dual;
  Dart_handle d, d2, res = NULL;
  
  // We clear the amap2. TODO return a new amap ? (but we need to make
  // a copy contructor and =operator...)
  amap2.clear();
  
  // We create a copy of all the dart of the map.
  for (Dart_iterator it=amap1.darts().begin(); it!=amap1.darts().end(); ++it)
    {
      dual[it] = amap2.create_dart();

      if ( it==adart && res==NULL ) res = dual[it];
    }
  
  // Then we link the darts by using the dual formula :
  // G(B,b1,b2,...,bn-1,bn) => dual(G)=(B, b(n-1)obn, b(n-2)obn,...,b1obn, bn)
  // We suppose darts are run in the same order for both maps.
  Dart_iterator it2=amap2.darts().begin();
  for (Dart_iterator it=amap1.darts().begin(); it!=amap1.darts().end(); 
       ++it, ++it2)
  {
    d = it2; // The supposition on the order allows to avoid d=dual[it];
    CGAL_assertion(it2 == dual[it]);

    // First case outside the loop since we need to use link_beta1
    if ( it->beta(Map::dimension)->beta(Map::dimension-1)!=Map::null_dart_handle )
      amap2. template
        link_beta<1>(d, 
                     dual[it->beta(Map::dimension)->beta(Map::dimension-1)]);

    // and during the loop we use link_beta(d1,d2,i)
    for (unsigned int i=Map::dimension-2; i>=1; --i)
    {
      if ( it->beta(Map::dimension)->beta(i)!=Map::null_dart_handle )
        amap2.link_beta(d, dual[it->beta(Map::dimension)->beta(i)],
                        Map::dimension-i);
    }
    CGAL_assertion ( !it->is_free(Map::dimension) );
    amap2.link_beta(d, dual[it->beta(Map::dimension)],Map::dimension);
  }
  
  // Now the map amap is topologically correct, we just need to add
  // its geometry to each vertex (the barycenter of the corresponding
  // volume in the initial map).
  it2 = amap2.darts().begin();
  for (Dart_iterator it(amap1.darts().begin()); it!=amap1.darts().end();
       ++it, ++it2)
  {
    if (Map::vertex_attribute(it2) == NULL)
    {
      amap2.set_vertex_attribute(it2,
                                 amap2.create_vertex_attribute
                                 (amap1.barycenter<Map::dimension>(it)));
    }
  }

  //  CGAL_postcondition(amap2.is_valid());

  if ( res==NULL ) res = amap2.darts().begin();
  return res;
}

} // namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_OPERATIONS_H //
// EOF //
