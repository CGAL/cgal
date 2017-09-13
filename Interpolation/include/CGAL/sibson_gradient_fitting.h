// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
//
// Author(s)     : Julia Floetotto

#ifndef CGAL_SIBSON_GRADIENT_FITTING_H
#define CGAL_SIBSON_GRADIENT_FITTING_H

#include <CGAL/license/Interpolation.h>

#include <CGAL/Origin.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/regular_neighbor_coordinates_2.h>

#include <iterator>
#include <utility>

namespace CGAL {

template < class ForwardIterator, class Functor, class Traits>
typename Traits::Vector_d
sibson_gradient_fitting(ForwardIterator first, ForwardIterator beyond,
                        const typename
                        std::iterator_traits<ForwardIterator>::
                          value_type::second_type& norm,
                        const typename
                        std::iterator_traits<ForwardIterator>::value_type
                          ::first_type& p,
                        Functor function_value,
                        const Traits& traits)
{
  CGAL_precondition( first!=beyond && norm!=0);
  typedef typename Traits::Aff_transformation_d Aff_transformation;
  typedef typename Traits::FT                   Coord_type;

  const typename Traits::Point_d& bare_p = traits.construct_point_d_object()(p);

  typename Functor::result_type fn = function_value(bare_p);
  CGAL_assertion(fn.second); //function value of p is valid

  typename Traits::Vector_d pn =
      traits.construct_vector_d_object()(NULL_VECTOR);
  Aff_transformation scaling, m,
      Hn(traits.construct_null_matrix_d_object()());

  for(;first!=beyond; ++first)
  {
    const typename Traits::Point_d& bare_f = traits.construct_point_d_object()(first->first);

    Coord_type square_dist = traits.compute_squared_distance_d_object()(bare_f, bare_p);
    CGAL_assertion(square_dist != 0);
    Coord_type scale = first->second / (norm*square_dist);
    typename Traits::Vector_d d = traits.construct_vector_d_object()(bare_p, bare_f);

    //compute the vector pn:
    typename Functor::result_type f = function_value(bare_f);
    CGAL_assertion(f.second);//function value of first->first is valid
    pn = pn + traits.construct_scaled_vector_d_object()
         (d,scale * (f.first - fn.first));

    //compute the matrix Hn:
    m = traits.construct_outer_product_d_object()(d);
    scaling = traits.construct_scaling_matrix_d_object()(scale);

    Hn = traits.construct_sum_matrix_d_object()(Hn, scaling * m);
  }

  return Hn.inverse().transform(pn);
}

template < class Triangul, class OutputIterator, class Functor,
           class CoordFunctor, class Traits>
OutputIterator
sibson_gradient_fitting(const Triangul& tr,
                        OutputIterator out,
                        Functor function_value,
                        CoordFunctor compute_coordinates,
                        const Traits& traits)
{
  typedef typename Triangul::Point                        Point;
  typedef typename Traits::FT                             Coord_type;

  std::vector< std::pair< Point, Coord_type > > coords;
  Coord_type norm;

  typename Triangul::Finite_vertices_iterator vit = tr.finite_vertices_begin();
  for(; vit != tr.finite_vertices_end(); ++vit){
    //test if vit is a convex hull vertex, otherwise do nothing
    if (!tr.is_edge(vit, tr.infinite_vertex()))
    {
      norm = compute_coordinates(tr, vit, std::back_inserter(coords)).second;
      *out++ = std::make_pair(traits.construct_point_d_object()(vit->point()),
                              sibson_gradient_fitting(coords.begin(),
                                                      coords.end(),
                                                      norm, vit->point(),
                                                      function_value,
                                                      traits));
      coords.clear();
    }
  }
  return out;
}

//the following functions allow to fit the gradients for all points in
// a triangulation except the convex hull points.
// -> _nn2: natural_neighbor_coordinates_2
// -> _rn2: regular_neighbor_coordinates_2
// -> _sn2_3: surface_neighbor_coordinates_2_3
template < class Dt, class OutputIterator, class Functor, class Traits>
OutputIterator
sibson_gradient_fitting_nn_2(const Dt& dt,
                             OutputIterator out,
                             Functor function_value,
                             const Traits& traits)
{
  typedef typename std::back_insert_iterator<
                     std::vector<
                       std::pair< typename Traits::Point_d,
                                  typename  Traits::FT > > >   CoordInserter;

  return sibson_gradient_fitting(dt, out, function_value,
                                 natural_neighbor_coordinates_2_object< Dt,
                                                                        CoordInserter >(),
                                 traits);
}

template < class Rt, class OutputIterator, class Functor, class Traits>
OutputIterator
sibson_gradient_fitting_rn_2(const Rt& rt,
                             OutputIterator out,
                             Functor function_value,
                             const Traits& traits)
{
  typedef typename std::back_insert_iterator<
                     std::vector<
                       std::pair< typename Traits::Weighted_point_d,
                                  typename Traits::FT > > >   CoordInserter;

  return sibson_gradient_fitting(rt, out, function_value,
                                 regular_neighbor_coordinates_2_object< Rt,
                                                                        CoordInserter >(),
                                 traits);
}

} //namespace CGAL

#endif // CGAL_SIBSON_GRADIENT_FITTING_H
