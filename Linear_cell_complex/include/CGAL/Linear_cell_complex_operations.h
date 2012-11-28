// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_LINEAR_CELL_COMPLEX_OPERATIONS_H
#define CGAL_LINEAR_CELL_COMPLEX_OPERATIONS_H 1

#include <CGAL/Cell_iterators.h>
#include <CGAL/Combinatorial_map_operations.h>
#include <vector>

namespace CGAL {

  /** @file Linear_cell_complex_operations.h
   * Basic operators on  a linear cell complex.
   */


  /** Compute the normal of the given facet.
   * @param amap the used linear cell complex.
   * @param adart a dart incident to the facet.
   * @return the normal of the facet.
   */
  template <class LCC>
  typename LCC::Vector compute_normal_of_cell_2
  (const LCC& /*amap*/, typename LCC::Dart_const_handle adart)
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
    for ( ; adart!=start && adart->other_extremity()!=NULL;
          adart=adart->beta(1) )
    {
      const Point* next = &LCC::point(adart->other_extremity());
      if ( !typename LCC::Traits::Collinear_3()(*prev, *curr, *next) )
      {
        normal = typename LCC::Traits::Construct_sum_of_vectors()
          (normal, typename LCC::Traits::Construct_normal_3()
           (*prev, *curr, *next));
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
   * @param amap the used linear cell complex.
   * @param adart a dart incident to the vertex.
   * @return the normal of the vertex.
   */
  template <class LCC>
  typename LCC::Vector compute_normal_of_cell_0
  (const LCC& amap, typename LCC::Dart_const_handle adart)
  {
    typedef typename LCC::Vector Vector;
    Vector normal(CGAL::NULL_VECTOR);
    unsigned int nb = 0;
  
    for ( CMap_one_dart_per_incident_cell_const_iterator<LCC,2,0>
            it(amap, adart); it.cont(); ++it )
    {
      normal = typename LCC::Traits::Construct_sum_of_vectors()
        (normal, CGAL::compute_normal_of_cell_2(amap,it));
      ++nb;
    }
  
    if ( nb<2 ) return normal;       
    return (typename LCC::Traits::Construct_scaled_vector()(normal, 1.0/nb));
  }

} // namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_OPERATIONS_H //
// EOF //
