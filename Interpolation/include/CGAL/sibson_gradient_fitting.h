// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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
#include <functional>

#include <boost/any.hpp>
#include <boost/mpl/if.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/utility/result_of.hpp>

#include <iterator>
#include <utility>
#include <vector>

#include <type_traits>
#include <functional>

namespace CGAL {

template < class ForwardIterator, class ValueFunctor, class Traits, class Point >
typename Traits::Vector_d
sibson_gradient_fitting(ForwardIterator first,
                        ForwardIterator beyond,
                        const typename std::iterator_traits<ForwardIterator>::value_type::second_type& norm,
                        const Point& p,
                        const typename boost::result_of<
                                         ValueFunctor(typename std::iterator_traits<ForwardIterator>::value_type::first_type)>
                                         ::type::first_type fn,
                        ValueFunctor value_function,
                        const Traits& traits)
{
  CGAL_precondition( first != beyond && norm != 0);

  typedef typename std::iterator_traits<ForwardIterator>::value_type::first_type  arg_type;
  typedef typename boost::result_of<ValueFunctor(arg_type)>::type                 value_functor_result_type;

  typedef typename Traits::Aff_transformation_d Aff_transformation;
  typedef typename Traits::FT                   Coord_type;
  typedef typename Traits::Point_d              Bare_point;

  typename Traits::Vector_d pn = traits.construct_vector_d_object()(NULL_VECTOR);
  Aff_transformation scaling, m, Hn(traits.construct_null_matrix_d_object()());
  Interpolation::internal::Extract_bare_point<Traits> cp(traits);
  const Bare_point& bp = cp(p);

  for(; first!=beyond; ++first)
  {
    const typename Traits::Point_d& bare_f = cp(first->first);
    Coord_type square_dist = traits.compute_squared_distance_d_object()(bare_f, bp);
    CGAL_assertion(square_dist != 0);
    Coord_type scale = first->second / (norm * square_dist);
    typename Traits::Vector_d d = traits.construct_vector_d_object()(bp, bare_f);

    // compute the vector pn:
    value_functor_result_type f = value_function(first->first);
    CGAL_assertion(f.second); // function value of first->first is valid
    pn = pn + traits.construct_scaled_vector_d_object()(d, scale * (f.first - fn));

    // compute the matrix Hn:
    m = traits.construct_outer_product_d_object()(d);
    scaling = traits.construct_scaling_matrix_d_object()(scale);
    Hn = traits.construct_sum_matrix_d_object()(Hn, scaling * m);
  }
  return Hn.inverse().transform(pn);
}

// for backward compatibility
template < class ForwardIterator, class ValueFunctor, class Traits >
typename Traits::Vector_d
sibson_gradient_fitting(ForwardIterator first,
                        ForwardIterator beyond,
                        const typename std::iterator_traits<ForwardIterator>::value_type::second_type& norm,
                        const typename std::iterator_traits<ForwardIterator>::value_type::first_type& p,
                        ValueFunctor value_function,
                        const Traits& traits)
{
  typedef typename std::iterator_traits<ForwardIterator>::value_type::first_type  arg_type;
  typedef typename boost::result_of<ValueFunctor(arg_type)>::type                 value_functor_result_type;

  value_functor_result_type fn = value_function(p);
  CGAL_assertion(fn.second);

  return sibson_gradient_fitting(first, beyond, norm, p, fn.first, value_function, traits);
}

// The next three functions are used to call the value functor for different
// types of arguments and pass a final (bare) point + value to the function above.
template < class ForwardIterator, class ValueFunctor, class Traits, class VH >
typename Traits::Vector_d
sibson_gradient_fitting_internal_with_dummy(ForwardIterator first,
                                            ForwardIterator beyond,
                                            const typename std::iterator_traits<
                                            ForwardIterator>::value_type::second_type& norm,
                                            VH vh,
                                            ValueFunctor value_function,
                                            const Traits& traits,
                                            const typename Traits::Point_d& /*dummy*/)
{
  typedef typename std::iterator_traits<ForwardIterator>::value_type::first_type  arg_type;
  typedef typename boost::result_of<ValueFunctor(arg_type)>::type                 value_functor_result_type;

  const typename Traits::Point_d& bare_p = traits.construct_point_d_object()(vh->point());
  value_functor_result_type fn = value_function(bare_p);
  CGAL_assertion(fn.second);

  return sibson_gradient_fitting(first, beyond, norm, bare_p, fn.first, value_function, traits);
}


template < class ForwardIterator, class ValueFunctor, class Traits, class VH >
typename Traits::Vector_d
sibson_gradient_fitting_internal_with_dummy(ForwardIterator first,
                                            ForwardIterator beyond,
                                            const typename std::iterator_traits<
                                            ForwardIterator>::value_type::second_type& norm,
                                            VH vh,
                                            ValueFunctor value_function,
                                            const Traits& traits,
                                            const typename Traits::Weighted_point_d& /*dummy*/)
{
  typedef typename std::iterator_traits<ForwardIterator>::value_type::first_type  arg_type;
  typedef typename boost::result_of<ValueFunctor(arg_type)>::type                 value_functor_result_type;

  value_functor_result_type fn = value_function(vh->point());
  CGAL_assertion(fn.second);

  return sibson_gradient_fitting(first, beyond, norm, vh->point(), fn.first, value_function, traits);
}


template < class ForwardIterator, class ValueFunctor, class Traits, class VH >
typename Traits::Vector_d
sibson_gradient_fitting_internal_with_dummy(ForwardIterator first,
                                            ForwardIterator beyond,
                                            const typename std::iterator_traits<
                                            ForwardIterator>::value_type::second_type& norm,
                                            VH vh,
                                            ValueFunctor value_function,
                                            const Traits& traits,
                                            VH /*dummy*/)
{
  typedef typename std::iterator_traits<ForwardIterator>::value_type::first_type  arg_type;
  typedef typename boost::result_of<ValueFunctor(arg_type)>::type                 value_functor_result_type;

  const typename Traits::Point_d& bare_p = traits.construct_point_d_object()(vh->point());
  value_functor_result_type fn = value_function(vh);
  CGAL_assertion(fn.second);

  return sibson_gradient_fitting(first, beyond, norm, bare_p, fn.first, value_function, traits);
}


template < class ValueFunctorArgType,
           class Tr, class OutputIterator, class OutputFunctor,
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
  std::vector<std::pair<ValueFunctorArgType, Coord_type> > coords;

  typename Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin();
  for(; vit != tr.finite_vertices_end(); ++vit)
  {
    // test if vit is a convex hull vertex, otherwise do nothing
    if(!tr.is_edge(vit, tr.infinite_vertex()))
    {
      norm = compute_coordinates(tr, vit, std::back_inserter(coords), Coord_OutputFunctor()).second;

      *out++ = fct(std::make_pair(vit,
                                  sibson_gradient_fitting_internal_with_dummy(coords.begin(),
                                                                              coords.end(),
                                                                              norm,
                                                                              Vertex_handle(vit),
                                                                              value_function,
                                                                              traits,
                                                                              ValueFunctorArgType())));

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

// The ugly distinction below is needed to make it work with lambdas for C++11 because std::is_constructible
// is used, which is C++11 (there is a boost equivalent, but it is said (by boost) to be relying on C++11 features
// to properly work...)
template < class Dt, class OutputIterator, class OutputFunctor, class ValueFunctor, class Traits >
OutputIterator
sibson_gradient_fitting_nn_2(const Dt& dt,
                             OutputIterator out,
                             OutputFunctor fct,
                             ValueFunctor value_function,
                             const Traits& traits,
                             // Some SFINAE to distinguish whether the argument type
                             // of the value functor is 'DT::Point' or 'DT::Vertex_handle'
                             typename boost::enable_if_c<
                               std::is_constructible<
                                 std::function<boost::any(typename Dt::Point)>,
                                 ValueFunctor
                             >::value>::type* = nullptr)
{
  typedef typename Traits::FT                                        FT;
  typedef typename Dt::Point                                         VF_arg_type;
  typedef typename std::back_insert_iterator<std::vector<
                     std::pair<VF_arg_type, FT> > >                  CoordInserter;
  typedef Interpolation::internal::Extract_point_in_pair<Dt, FT>     Coord_OutputFunctor;

  return sibson_gradient_fitting_internal<VF_arg_type>(dt, out, fct, value_function,
                                                       natural_neighbor_coordinates_2_object<Dt,
                                                                                             CoordInserter,
                                                                                             Coord_OutputFunctor>(),
                                                       traits);
}

template < class Dt, class OutputIterator, class OutputFunctor, class ValueFunctor, class Traits >
OutputIterator
sibson_gradient_fitting_nn_2(const Dt& dt,
                             OutputIterator out,
                             OutputFunctor fct,
                             ValueFunctor value_function,
                             const Traits& traits,
                             typename boost::enable_if_c<
                               std::is_constructible<
                                 std::function<boost::any(typename Dt::Vertex_handle)>,
                                 ValueFunctor
                             >::value>::type* = nullptr)
{
  typedef typename Traits::FT                                        FT;
  typedef typename Dt::Vertex_handle                                 VF_arg_type;
  typedef typename std::back_insert_iterator<std::vector<
                     std::pair<VF_arg_type, FT> > >                  CoordInserter;
  typedef CGAL::Identity<std::pair<VF_arg_type, FT> >                Coord_OutputFunctor;

  return sibson_gradient_fitting_internal<VF_arg_type>(dt, out, fct, value_function,
                                                       natural_neighbor_coordinates_2_object<Dt,
                                                                                             CoordInserter,
                                                                                             Coord_OutputFunctor>(),
                                                       traits);
}

// Same as above but without OutputFunctor.
// Defaults to extracting the point, for backward compatibility.
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

// See above for the explanation.
template < class Rt, class OutputIterator, class OutputFunctor, class ValueFunctor, class Traits >
OutputIterator
sibson_gradient_fitting_rn_2(const Rt& rt,
                             OutputIterator out,
                             OutputFunctor fct,
                             ValueFunctor value_function,
                             const Traits& traits,
                             // Some SFINAE to distinguish whether the argument type
                             // of the value functor is 'Rt::Point' (weighted point) or 'Rt::Vertex_handle'
                             typename boost::enable_if_c<
                               std::is_constructible<
                                 std::function<boost::any(typename Rt::Point)>,
                                 ValueFunctor
                             >::value>::type* = nullptr)
{
  typedef typename Traits::FT                                        FT;
  typedef typename Rt::Point                                         VF_arg_type;
  typedef typename std::back_insert_iterator<std::vector<
                       std::pair<VF_arg_type, FT> > >                CoordInserter;
  typedef Interpolation::internal::Extract_point_in_pair<Rt, FT>     Coord_OutputFunctor;

  return sibson_gradient_fitting_internal<VF_arg_type>(rt, out, fct, value_function,
                                                       regular_neighbor_coordinates_2_object<Rt,
                                                                                             CoordInserter,
                                                                                             Coord_OutputFunctor>(),
                                                       traits);
}

template < class Rt, class OutputIterator, class OutputFunctor, class ValueFunctor, class Traits >
OutputIterator
sibson_gradient_fitting_rn_2(const Rt& rt,
                             OutputIterator out,
                             OutputFunctor fct,
                             ValueFunctor value_function,
                             const Traits& traits,
                             typename boost::enable_if_c<
                               std::is_constructible<
                                 std::function<boost::any(typename Rt::Vertex_handle)>,
                                 ValueFunctor
                             >::value>::type* = nullptr)
{
  typedef typename Traits::FT                                        FT;
  typedef typename Rt::Vertex_handle                                 VF_arg_type;
  typedef typename std::back_insert_iterator<std::vector<
                       std::pair<VF_arg_type, FT> > >                CoordInserter;
  typedef CGAL::Identity<std::pair<VF_arg_type, FT> >                Coord_OutputFunctor;

  return sibson_gradient_fitting_internal<VF_arg_type>(rt, out, fct, value_function,
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
