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

#include <CGAL/Interpolation/internal/helpers.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/regular_neighbor_coordinates_2.h>

#include <CGAL/Origin.h>

#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>

#include <iterator>
#include <utility>
#include <vector>

namespace CGAL {

template < class ForwardIterator, class ValueFunctor, class Traits >
typename Traits::Vector_d
sibson_gradient_fitting(ForwardIterator first,
                        ForwardIterator beyond,
                        const typename std::iterator_traits<ForwardIterator>::value_type::second_type& norm,
                        const typename Traits::Point_d& bare_p,
                        const typename ValueFunctor::result_type::first_type fn,
                        ValueFunctor value_function,
                        const Traits& traits)
{
  CGAL_precondition( first != beyond && norm != 0);

  typedef typename Traits::Aff_transformation_d Aff_transformation;
  typedef typename Traits::FT                   Coord_type;

  typename Traits::Vector_d pn = traits.construct_vector_d_object()(NULL_VECTOR);
  Aff_transformation scaling, m, Hn(traits.construct_null_matrix_d_object()());
  Interpolation::internal::Extract_bare_point<Traits> cp(traits);

  for(; first!=beyond; ++first)
  {
    const typename Traits::Point_d& bare_f = cp(first->first);
    Coord_type square_dist = traits.compute_squared_distance_d_object()(bare_f, bare_p);
    CGAL_assertion(square_dist != 0);
    Coord_type scale = first->second / (norm * square_dist);
    typename Traits::Vector_d d = traits.construct_vector_d_object()(bare_p, bare_f);

    // compute the vector pn:
    typename ValueFunctor::result_type f = value_function(first->first);
    CGAL_assertion(f.second); // function value of first->first is valid
    pn = pn + traits.construct_scaled_vector_d_object()(d, scale * (f.first - fn));

    // compute the matrix Hn:
    m = traits.construct_outer_product_d_object()(d);
    scaling = traits.construct_scaling_matrix_d_object()(scale);
    Hn = traits.construct_sum_matrix_d_object()(Hn, scaling * m);
  }
  return Hn.inverse().transform(pn);
}


// The next three functions are used to call the value functor for different
// types of arguments and pass a final (bare) point + value to the function above.
template < class ForwardIterator, class ValueFunctor, class Traits, class VH >
typename Traits::Vector_d
sibson_gradient_fitting_internal(ForwardIterator first,
                                 ForwardIterator beyond,
                                 const typename std::iterator_traits<
                                                  ForwardIterator>::value_type::second_type& norm,
                                 VH vh,
                                 ValueFunctor value_function,
                                 const Traits& traits,
                                 const typename Traits::Point_d& /*dummy*/)
{
  const typename Traits::Point_d& bare_p = traits.construct_point_d_object()(vh->point());
  typename ValueFunctor::result_type fn = value_function(bare_p);
  CGAL_assertion(fn.second);

  return sibson_gradient_fitting(first, beyond, norm, bare_p, fn.first, value_function, traits);
}


template < class ForwardIterator, class ValueFunctor, class Traits, class VH >
typename Traits::Vector_d
sibson_gradient_fitting_internal(ForwardIterator first,
                                 ForwardIterator beyond,
                                 const typename std::iterator_traits<
                                                  ForwardIterator>::value_type::second_type& norm,
                                 VH vh,
                                 ValueFunctor value_function,
                                 const Traits& traits,
                                 const typename Traits::Weighted_point_d& /*dummy*/)
{
  const typename Traits::Point_d& bare_p = traits.construct_point_d_object()(vh->point());
  typename ValueFunctor::result_type fn = value_function(vh->point());
  CGAL_assertion(fn.second);

  return sibson_gradient_fitting(first, beyond, norm, bare_p, fn.first, value_function, traits);
}


template < class ForwardIterator, class ValueFunctor, class Traits, class VH >
typename Traits::Vector_d
sibson_gradient_fitting_internal(ForwardIterator first,
                                 ForwardIterator beyond,
                                 const typename std::iterator_traits<
                                                  ForwardIterator>::value_type::second_type& norm,
                                 VH vh,
                                 ValueFunctor value_function,
                                 const Traits& traits,
                                 VH /*dummy*/)
{
  const typename Traits::Point_d& bare_p = traits.construct_point_d_object()(vh->point());
  typename ValueFunctor::result_type fn = value_function(vh);
  CGAL_assertion(fn.second);

  return sibson_gradient_fitting(first, beyond, norm, bare_p, fn.first, value_function, traits);
}


template < class Tr, class OutputIterator, class OutputFunctor,
           class ValueFunctor, class CoordFunctor, class Traits >
OutputIterator
sibson_gradient_fitting_internal(const Tr& tr,
                                 OutputIterator out,
                                 OutputFunctor fct,
                                 ValueFunctor value_function,
                                 CoordFunctor compute_coordinates,
                                 const Traits& traits)
{
  typedef typename Traits::FT                             Coord_type;
  typedef typename CoordFunctor::Function                 Coord_OutputFunctor;
  typedef typename Tr::Vertex_handle                      Vertex_handle;

  Coord_type norm;
  std::vector<std::pair<typename ValueFunctor::argument_type, Coord_type> > coords;

  typename Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin();
  for(; vit != tr.finite_vertices_end(); ++vit)
  {
    // test if vit is a convex hull vertex, otherwise do nothing
    if(!tr.is_edge(vit, tr.infinite_vertex()))
    {
      norm = compute_coordinates(tr, vit, std::back_inserter(coords), Coord_OutputFunctor()).second;

      *out++ = fct(std::make_pair(vit,
                                  sibson_gradient_fitting_internal(coords.begin(),
                                                                   coords.end(),
                                                                   norm,
                                                                   Vertex_handle(vit),
                                                                   value_function,
                                                                   traits,
                                                                   typename ValueFunctor::argument_type())));

      coords.clear();
    }
  }

  return out;
}


// The following functions allow to fit the gradients for all points in
// a triangulation except the convex hull points.
// -> _nn2: natural_neighbor_coordinates_2
// -> _rn2: regular_neighbor_coordinates_2
// -> _sn2_3: surface_neighbor_coordinates_2_3
template < class Dt, class OutputIterator, class OutputFunctor, class ValueFunctor, class Traits >
OutputIterator
sibson_gradient_fitting_nn_2(const Dt& dt,
                             OutputIterator out,
                             OutputFunctor fct,
                             ValueFunctor value_function,
                             const Traits& traits)
{
  typedef typename Traits::FT                                        FT;
  typedef typename ValueFunctor::argument_type                       VF_arg_type;
  typedef typename std::back_insert_iterator<std::vector<
                     std::pair<VF_arg_type, FT> > >                  CoordInserter;

  // If the functor evaluates at points (and not vertices), then we must convert
  // the output of the coordinates computations - a pair<vertex_handle, FT> -
  // to a pair<point, FT>
  typedef typename boost::mpl::if_<
    boost::is_same<VF_arg_type, typename Dt::Point>,
    Interpolation::internal::Extract_point_in_pair<Dt, FT>,
    CGAL::Identity<std::pair<VF_arg_type, FT> >
  >::type                                                            Coord_OutputFunctor;

  return sibson_gradient_fitting_internal(dt, out, fct, value_function,
                                          natural_neighbor_coordinates_2_object<Dt,
                                                                                CoordInserter,
                                                                                Coord_OutputFunctor>(),
                                          traits);
}


// Same as above but without OutputFunctor. Default to extracting the point, for backward compatibility.
template < class Dt, class OutputIterator, class ValueFunctor, class Traits >
OutputIterator
sibson_gradient_fitting_nn_2(const Dt& dt,
                             OutputIterator out,
                             ValueFunctor value_function,
                             const Traits& traits)
{
  typedef typename Traits::Vector_d                                     Vector_d;
  typedef Interpolation::internal::Extract_point_in_pair<Dt, Vector_d>  OutputFunctor;

  return sibson_gradient_fitting_nn_2(dt, out, OutputFunctor(), value_function, traits);
}


template < class Rt, class OutputIterator, class OutputFunctor, class ValueFunctor, class Traits >
OutputIterator
sibson_gradient_fitting_rn_2(const Rt& rt,
                             OutputIterator out,
                             OutputFunctor fct,
                             ValueFunctor value_function,
                             const Traits& traits)
{
  typedef typename Traits::FT                                        FT;
  typedef typename ValueFunctor::argument_type                       VF_arg_type;
  typedef typename std::back_insert_iterator<std::vector<
                       std::pair<VF_arg_type, FT> > >                CoordInserter;

  // If the functor evaluates at points (and not vertices), then we must convert
  // the output of the coordinates computations - a pair<vertex_handle, FT> -
  // to a pair<point, FT>
  typedef typename boost::mpl::if_<
    boost::is_same<VF_arg_type, typename Rt::Weighted_point>,
    Interpolation::internal::Extract_point_in_pair<Rt, FT>,
    CGAL::Identity<std::pair<VF_arg_type, FT> >
  >::type                                                            Coord_OutputFunctor;

  return sibson_gradient_fitting_internal(rt, out, fct, value_function,
                                          regular_neighbor_coordinates_2_object<Rt,
                                                                                CoordInserter,
                                                                                Coord_OutputFunctor>(),
                                          traits);
}


// Same as above but without OutputFunctor. Default to extracting the point, for backward compatibility.
template < class Rt, class OutputIterator, class ValueFunctor, class Traits >
OutputIterator
sibson_gradient_fitting_rn_2(const Rt& rt,
                             OutputIterator out,
                             ValueFunctor value_function,
                             const Traits& traits)
{
  typedef typename Traits::Vector_d                                     Vector_d;
  typedef Interpolation::internal::Extract_point_in_pair<Rt, Vector_d>  OutputFunctor;

  return sibson_gradient_fitting_rn_2(rt, out, OutputFunctor(), value_function, traits);
}

} //namespace CGAL

#endif // CGAL_SIBSON_GRADIENT_FITTING_H
