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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Julia Floetotto

#ifndef CGAL_SIBSON_GRADIENT_FITTING_H
#define CGAL_SIBSON_GRADIENT_FITTING_H

#include <CGAL/license/Interpolation.h>

#include <CGAL/Origin.h>
#include <CGAL/interpolation_functions.h> // AF: for V2P
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/regular_neighbor_coordinates_2.h>

#include <iterator>
#include <utility>

namespace CGAL {

template < class ForwardIterator, class Functor, class Traits>
typename Traits::Vector_d
sibson_gradient_fitting(ForwardIterator first,
                        ForwardIterator beyond,
                        const typename
                        std::iterator_traits<ForwardIterator>::
                          value_type::second_type& norm,
                        const typename Traits::Point_d& bare_p,
                        const typename Functor::result_type::first_type fn,
                        Functor function_value,
                        const Traits& traits)
{
  internal::interpolation::V2P<Traits> v2p(traits);
  CGAL_precondition( first!=beyond && norm!=0);
  typedef typename Traits::Aff_transformation_d Aff_transformation;
  typedef typename Traits::FT                   Coord_type;

  typename Traits::Vector_d pn =
      traits.construct_vector_d_object()(NULL_VECTOR);
  Aff_transformation scaling, m,
      Hn(traits.construct_null_matrix_d_object()());

  for(;first!=beyond; ++first)
  {
    const typename Traits::Point_d& bare_f = v2p(first->first);
    Coord_type square_dist = traits.compute_squared_distance_d_object()(bare_f, bare_p);
    CGAL_assertion(square_dist != 0);
    Coord_type scale = first->second / (norm*square_dist);
    typename Traits::Vector_d d = traits.construct_vector_d_object()(bare_p, bare_f);

    //compute the vector pn:
    typename Functor::result_type f = function_value(first->first);
    CGAL_assertion(f.second);//function value of first->first is valid
    pn = pn + traits.construct_scaled_vector_d_object()
         (d,scale * (f.first - fn));

    //compute the matrix Hn:
    m = traits.construct_outer_product_d_object()(d);
    scaling = traits.construct_scaling_matrix_d_object()(scale);

    Hn = traits.construct_sum_matrix_d_object()(Hn, scaling * m);
  }
  return Hn.inverse().transform(pn);
}

  
template < class Triangul, class ForwardIterator, class Functor, class Traits, class VH>
typename Traits::Vector_d
sibson_gradient_fitting(const Triangul& tr,
                        ForwardIterator first,
                        ForwardIterator beyond,
                        const typename
                        std::iterator_traits<ForwardIterator>::
                          value_type::second_type& norm,
                        VH vh,
                        Functor function_value,
                        const Traits& traits,
                        const typename Traits::Point_d& ignored)
{
  const typename Traits::Point_d& bare_p = traits.construct_point_d_object()(vh->point());
  typename Functor::result_type fn = function_value(bare_p);
 
  CGAL_assertion(fn.second);
                                                         
  return sibson_gradient_fitting(first,
                                 beyond,
                                 norm,
                                 bare_p,
                                 fn.first,
                                 function_value,
                                 traits);
}

template < class Triangul, class ForwardIterator, class Functor, class Traits, class VH>
typename Traits::Vector_d
sibson_gradient_fitting(const Triangul& tr,
                        ForwardIterator first,
                        ForwardIterator beyond,
                        const typename
                        std::iterator_traits<ForwardIterator>::
                          value_type::second_type& norm,
                        VH vh,
                        Functor function_value,
                        const Traits& traits,
                        typename Triangul::Vertex_handle ignored)
{
  const typename Traits::Point_d& bare_p = traits.construct_point_d_object()(vh->point());
  typename Functor::result_type fn = function_value(vh);
   CGAL_assertion(fn.second);
                                                         
  return sibson_gradient_fitting(first,
                                 beyond,
                                 norm,
                                 bare_p,
                                 fn.first,
                                 function_value,
                                 traits);
}

  
template < class Triangul, class OutputIterator, class Functor,
           class CoordFunctor, class OIF, class Traits>
OutputIterator
sibson_gradient_fitting_internal(const Triangul& tr,
                        OutputIterator out,
                        Functor function_value,
                        CoordFunctor compute_coordinates,
                        OIF fct,
                        const Traits& traits)
{
  typedef typename Traits::FT                             Coord_type;

  std::vector< std::pair< typename Functor::argument_type, Coord_type> > coords;
  Coord_type norm;

  typedef typename internal::interpolation::Output_iterator_functor_selector<Triangul, Traits,
                                                                             typename Functor::argument_type,
                                                                             Coord_type>::type     Coord_OIF;   
  
  typename Triangul::Finite_vertices_iterator vit = tr.finite_vertices_begin();

  for(; vit != tr.finite_vertices_end(); ++vit){
    //test if vit is a convex hull vertex, otherwise do nothing
    if (!tr.is_edge(vit, tr.infinite_vertex()))
    {
      norm = compute_coordinates(tr, vit, std::back_inserter(coords), Coord_OIF()).second;

       *out++ =  fct(std::make_pair(vit,
                                  sibson_gradient_fitting(tr,
                                                          coords.begin(),
                                                          coords.end(),
                                                          norm,
                                                          vit,
                                                          function_value,
                                                          traits,
                                                          typename Functor::argument_type())));
      
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
template < class Dt, class OutputIterator, class Functor, class OIF, class Traits>
OutputIterator
sibson_gradient_fitting_nn_2(const Dt& dt,
                             OutputIterator out,
                             Functor function_value,
                             OIF fct,
                             const Traits& traits)
{
  typedef typename std::back_insert_iterator<
                     std::vector<
                       std::pair< typename Functor::argument_type,
                                  typename  Traits::FT > > >   CoordInserter;

  typedef typename internal::interpolation::Output_iterator_functor_selector<Dt, Traits,
                                                                             typename Functor::argument_type,
                                                                             typename Traits::FT>::type     Coord_OIF;

  return sibson_gradient_fitting_internal(dt, out, function_value,
                                 natural_neighbor_coordinates_2_object<Dt,
                                                                       CoordInserter,
                                                                       Coord_OIF>(),
                                 fct,
                                 traits);
}


template < class Dt, class OutputIterator, class Functor, class Traits>
OutputIterator
sibson_gradient_fitting_nn_2(const Dt& dt,
                             OutputIterator out,
                             Functor function_value,    
                             const Traits& traits)
{
  return sibson_gradient_fitting_nn_2(dt, out, function_value,
                                      internal::interpolation::Vertex2Point<Dt, typename Traits::Vector_d>(), traits);
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
