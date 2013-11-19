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

    while ( !amap.template is_free<0>(start) &&
            amap.template beta<0>(start)!=adart )
      start = amap.template beta<0>(start);

    if ( amap.template is_free<1>(start) ||
         amap.other_extremity(amap.template beta<1>(start))==LCC::null_handle )
      return normal;

    unsigned int nb = 0;
    adart = amap.template beta<1>(start);

    const Point* prev = &amap.point(start);
    const Point* curr = &amap.point(adart);
    for ( ; adart!=start && amap.other_extremity(adart)!=LCC::null_handle;
          adart=amap.template beta<1>(adart) )
    {
      const Point* next = &amap.point( amap.other_extremity(adart));
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
  // Compute the barycenter of a given i-cell
  // General case, 1<i<=dimension
  template<class LCC, unsigned int i, unsigned int dim=LCC::dimension>
  struct Barycenter_functor
  {
    static typename LCC::Point run(const LCC& amap,
                                   typename LCC::Dart_const_handle adart)
    {
      CGAL_static_assertion(0<i && i<=LCC::dimension);
      CGAL_assertion(adart != LCC::null_handle);

      typename LCC::Vector vec
        (typename LCC::Traits::Construct_vector()(CGAL::ORIGIN,
                                                  amap.point(adart)));
      unsigned int nb = 1;

      CGAL::CMap_one_dart_per_incident_cell_const_iterator<LCC,0,i,i>
          it(amap, adart);
      for ( ++it; it.cont(); ++it)
      {
        vec = typename LCC::Traits::Construct_sum_of_vectors()
          (vec, typename LCC::Traits::Construct_vector()(CGAL::ORIGIN,
                                                    amap.point(it) ));
        ++nb;
      }

      return typename LCC::Traits::Construct_translated_point()
        (CGAL::ORIGIN, typename LCC::Traits::Construct_scaled_vector()
         (vec, 1.0/nb));
    }
  };

  // Compute the barycenter of a given 1-cell
  template<class LCC, unsigned int dim>
  struct Barycenter_functor<LCC, 1, dim>
  {
    static typename LCC::Point run(const LCC& amap,
                                   typename LCC::Dart_const_handle adart)
    {
      CGAL_static_assertion(1<=LCC::dimension);
      CGAL_assertion(adart != LCC::null_handle);
      typename LCC::Dart_const_handle d2=amap.other_extremity(adart);
      if (d2==amap.null_handle) return amap.point(adart);
      return typename LCC::Traits::Construct_midpoint()
        (amap.point(adart),
         amap.point(d2));
    }
  };

  // Compute the barycenter of a given 2-cell
  template<class LCC, unsigned int dim>
  struct Barycenter_functor<LCC, 2, dim>
  {
    static typename LCC::Point run(const LCC& amap,
                                   typename LCC::Dart_const_handle adart)
    {
      CGAL_static_assertion(2<=LCC::dimension);
      CGAL_assertion(adart != LCC::null_handle);

      typename LCC::Vector vec
        (typename LCC::Traits::Construct_vector()(CGAL::ORIGIN,
                                                  amap.point(adart)));
      unsigned int nb = 1;

        typename LCC::template Dart_of_cell_range<2,2>::const_iterator
          vhit  = amap.template darts_of_cell<2,2>(adart).begin(),
          vhend = amap.template darts_of_cell<2,2>(adart).end();
      for( ++vhit; vhit!=vhend; ++vhit )
      {
        vec = typename LCC::Traits::Construct_sum_of_vectors()
          (vec, typename LCC::Traits::Construct_vector()(CGAL::ORIGIN,
                                                         amap.point(vhit) ));
        ++nb;
      }
      return typename LCC::Traits::Construct_translated_point()
        (CGAL::ORIGIN, typename LCC::Traits::Construct_scaled_vector()
         (vec, 1.0/nb));
    }
  };

} // namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_OPERATIONS_H //
// EOF //
