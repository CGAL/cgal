// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julia Floetotto

#ifndef CGAL_INTERPOLATION_FUNCTIONS_H
#define CGAL_INTERPOLATION_FUNCTIONS_H

#include <CGAL/license/Interpolation.h>

#include <CGAL/Interpolation/internal/helpers.h>
#include <CGAL/double.h>
#include <CGAL/use.h>

#include <boost/utility/result_of.hpp>
#include <iterator>
#include <utility>
#include <vector>

namespace CGAL {

// Functor class for accessing the function values/gradients
template< class Map >
struct Data_access
  : public CGAL::cpp98::unary_function<typename Map::key_type,
                                       std::pair<typename Map::mapped_type, bool> >
{
  typedef typename Map::mapped_type   Data_type;
  typedef typename Map::key_type      Key_type;

  Data_access(const Map& m): map(m){}

  std::pair< Data_type, bool>
  operator()(const Key_type& p) const
  {
    typename Map::const_iterator mit = map.find(p);
    if(mit!= map.end())
      return std::make_pair(mit->second, true);
    return std::make_pair(Data_type(), false);
  }

  const Map& map;
};

//the interpolation functions:
  template < class ForwardIterator, class ValueFunctor>
typename boost::result_of<
           ValueFunctor(typename std::iterator_traits<ForwardIterator>::value_type::first_type)>
             ::type::first_type
linear_interpolation(ForwardIterator first, ForwardIterator beyond,
                     const typename std::iterator_traits<ForwardIterator>::value_type::second_type& norm,
                     ValueFunctor value_function)
{
  CGAL_precondition(first != beyond);
  CGAL_precondition(norm > 0);

  typedef typename std::iterator_traits<ForwardIterator>::value_type::first_type  arg_type;
  typedef typename boost::result_of<ValueFunctor(arg_type)>::type                 result_type;
  typedef typename result_type::first_type                                        Value_type;

  result_type val = value_function(first->first);
  CGAL_assertion(val.second);
  Value_type result = (first->second / norm) * val.first;
  ++first;
  for(; first!=beyond; ++first)
  {
    val = value_function(first->first);
    result = result + (first->second / norm) * val.first;
  }
  return result;
}

template < class ForwardIterator, class ValueFunctor, class GradFunctor, class Traits, class Point >
std::pair<
  typename boost::result_of<
             ValueFunctor(typename std::iterator_traits<ForwardIterator>::value_type::first_type)>
               ::type::first_type,
  bool>
quadratic_interpolation(ForwardIterator first, ForwardIterator beyond,
                        const typename std::iterator_traits<ForwardIterator>::value_type::second_type& norm,
                        const Point& p,
                        ValueFunctor value_function,
                        GradFunctor gradient_function,
                        const Traits& traits)
{
  CGAL_precondition(first != beyond);
  CGAL_precondition(norm > 0);

  typedef typename std::iterator_traits<ForwardIterator>::value_type::first_type  arg_type;
  typedef typename boost::result_of<ValueFunctor(arg_type)>::type                 value_functor_result_type;
  typedef typename boost::result_of<GradFunctor(arg_type)>::type                  gradient_functor_result_type;
  typedef typename value_functor_result_type::first_type                          Value_type;

  typedef typename Traits::Point_d                                                Bare_point;

  Interpolation::internal::Extract_bare_point<Traits> cp(traits);
  const Bare_point& bp = cp(p);

  Value_type result(0);
  value_functor_result_type f;
  gradient_functor_result_type grad;

  for(; first!=beyond; ++first)
  {
    f = value_function(first->first);
    grad = gradient_function(first->first);

    //test if value and gradient are correctly retrieved:
    CGAL_assertion(f.second);
    if(!grad.second)
      return std::make_pair(Value_type(0), false);

    result += (first->second/norm) * (f.first + grad.first *
                 traits.construct_scaled_vector_d_object()
                 (traits.construct_vector_d_object()(cp(first->first), bp), 0.5));
  }

  return std::make_pair(result, true);
}


template < class ForwardIterator, class ValueFunctor, class GradFunctor, class Traits, class Point >
std::pair<
  typename boost::result_of<
             ValueFunctor(typename std::iterator_traits<ForwardIterator>::value_type::first_type)>
               ::type::first_type,
  bool>
sibson_c1_interpolation(ForwardIterator first, ForwardIterator beyond,
                        const typename std::iterator_traits<ForwardIterator>::value_type::second_type& norm,
                        const Point& p,
                        ValueFunctor value_function,
                        GradFunctor gradient_function,
                        const Traits& traits)
{
  CGAL_precondition(first != beyond);
  CGAL_precondition(norm >0);

  typedef typename std::iterator_traits<ForwardIterator>::value_type::first_type  arg_type;
  typedef typename boost::result_of<ValueFunctor(arg_type)>::type                 value_functor_result_type;
  typedef typename boost::result_of<GradFunctor(arg_type)>::type                  gradient_functor_result_type;
  typedef typename value_functor_result_type::first_type                          Value_type;

  typedef typename Traits::FT                                                     Coord_type;
  typedef typename Traits::Point_d                                                Bare_point;

  Interpolation::internal::Extract_bare_point<Traits> cp(traits);
  const Bare_point& bp = cp(p);

  Coord_type term1(0), term2(term1), term3(term1), term4(term1);
  Value_type linear_int(0), gradient_int(0);
  value_functor_result_type f;
  gradient_functor_result_type grad;

  for(; first!=beyond; ++first)
  {
    f = value_function(first->first);
    grad = gradient_function(first->first);

    CGAL_assertion(f.second);
    if(!grad.second)
      return std::make_pair(Value_type(0), false); //the values are not correct

    Coord_type coeff = first->second/norm;
    Coord_type squared_dist = traits.compute_squared_distance_d_object()(cp(first->first), bp);
    Coord_type dist = CGAL_NTS sqrt(squared_dist);

    if(squared_dist == 0)
    {
      ForwardIterator it = first;
      CGAL_USE(it);
      CGAL_assertion(++it == beyond);
      return std::make_pair(f.first, true);
    }

    //three different terms to mix linear and gradient
    //interpolation
    term1 += coeff / dist;
    term2 += coeff * squared_dist;
    term3 += coeff * dist;

    linear_int += coeff * f.first;

    gradient_int += (coeff/dist) * (f.first + grad.first *
                       traits.construct_vector_d_object()(cp(first->first), bp));
  }

  term4 = term3 / term1;
  gradient_int = gradient_int / term1;

  return std::make_pair((term4 * linear_int + term2 * gradient_int) /
                        (term4 + term2), true);
}

// This method works with rational number types:
// modification of Sibson's interpolant without sqrt
// following a proposition by Gunther Rote:
//
// The general scheme:
//  Coord_type inv_weight = f(dist); //i.e. dist^2
//           term1 +=  coeff/inv_weight;
//    term2 +=  coeff * squared_dist;
//    term3 +=  coeff*(squared_dist/inv_weight);
//           gradient_int += (coeff/inv_weight) * (vh->get_value()+ vh->get_gradient() * (p - vh->point()));

template < class ForwardIterator, class ValueFunctor, class GradFunctor, class Traits, class Point >
std::pair<
  typename boost::result_of<
             ValueFunctor(typename std::iterator_traits<ForwardIterator>::value_type::first_type)>
               ::type::first_type,
  bool>
sibson_c1_interpolation_square(ForwardIterator first, ForwardIterator beyond,
                               const typename std::iterator_traits<ForwardIterator>::value_type::second_type& norm,
                               const Point& p,
                               ValueFunctor value_function,
                               GradFunctor gradient_function,
                               const Traits& traits)
{
  CGAL_precondition(first != beyond);
  CGAL_precondition(norm > 0);

  typedef typename std::iterator_traits<ForwardIterator>::value_type::first_type  arg_type;
  typedef typename boost::result_of<ValueFunctor(arg_type)>::type                 value_functor_result_type;
  typedef typename boost::result_of<GradFunctor(arg_type)>::type                  gradient_functor_result_type;
  typedef typename value_functor_result_type::first_type                          Value_type;

  typedef typename Traits::FT                                                     Coord_type;
  typedef typename Traits::Point_d                                                Bare_point;

  Interpolation::internal::Extract_bare_point<Traits> cp(traits);
  const Bare_point& bp = cp(p);

  Coord_type term1(0), term2(term1), term3(term1), term4(term1);
  Value_type linear_int(0), gradient_int(0);
  value_functor_result_type f;
  gradient_functor_result_type grad;

  for(; first!=beyond; ++first)
  {
    f = value_function(first->first);
    grad = gradient_function(first->first);
    CGAL_assertion(f.second);

    if(!grad.second)
      return std::make_pair(Value_type(0), false); // the gradient is not known

    Coord_type coeff = first->second/norm;
    Coord_type squared_dist = traits.compute_squared_distance_d_object()(cp(first->first), bp);

    if(squared_dist ==0)
    {
      ForwardIterator it = first;
      CGAL_USE(it);
      CGAL_assertion(++it == beyond);
      return std::make_pair(f.first, true);
    }

    //three different terms to mix linear and gradient
    //interpolation
    term1 += coeff / squared_dist;
    term2 += coeff * squared_dist;
    term3 += coeff;

    linear_int += coeff * f.first;

    gradient_int += (coeff/squared_dist) * (f.first + grad.first *
                       traits.construct_vector_d_object()(cp(first->first), bp));
  }

  term4 = term3/ term1;
  gradient_int = gradient_int / term1;

  return std::make_pair((term4 * linear_int + term2 * gradient_int)/
                        (term4 + term2), true);
}


template < class RandomAccessIterator, class ValueFunctor, class GradFunctor,
           class Traits, class Point_>
std::pair<
  typename boost::result_of<
             ValueFunctor(typename std::iterator_traits<RandomAccessIterator>::value_type::first_type)>
               ::type::first_type,
  bool>
farin_c1_interpolation(RandomAccessIterator first,
                       RandomAccessIterator beyond,
                       const typename std::iterator_traits<RandomAccessIterator>::value_type::second_type& norm,
                       const Point_& /*p*/,
                       ValueFunctor value_function,
                       GradFunctor gradient_function,
                       const Traits& traits)
{
  CGAL_precondition(first != beyond);
  CGAL_precondition(norm >0);

  typedef typename std::iterator_traits<RandomAccessIterator>::value_type::first_type  arg_type;
  typedef typename boost::result_of<ValueFunctor(arg_type)>::type                 value_functor_result_type;
  typedef typename boost::result_of<GradFunctor(arg_type)>::type                  gradient_functor_result_type;
  typedef typename value_functor_result_type::first_type                          Value_type;

  typedef typename Traits::FT                                                     Coord_type;

  Interpolation::internal::Extract_bare_point<Traits> cp(traits);
  value_functor_result_type f;
  gradient_functor_result_type grad;

  int n = static_cast<int>(beyond - first);
  if(n == 1)
  {
    f = value_function(first->first);
    CGAL_assertion(f.second);
    return std::make_pair(f.first, true);
  }

  //there must be one or at least three NN-neighbors:
  CGAL_assertion(n > 2);

  RandomAccessIterator it2, it;

  Value_type result(0);
  const Coord_type fac3(3);

  std::vector< std::vector<Value_type> > ordinates(n, std::vector<Value_type>(n, Value_type(0)));

  for(int i=0; i<n; ++i)
  {
    it = first + i;
    Coord_type coord_i_square = CGAL_NTS square(it->second);

    // for later: the function value of it->first:
    f = value_function(it->first);
    CGAL_assertion(f.second);
    ordinates[i][i] = f.first;

    // control point = data point
    result += coord_i_square * it->second* ordinates[i][i];

    // compute tangent plane control point (one 2, one 1 entry)
    Value_type res_i(0);
    for(int j=0; j<n; ++j)
    {
      if(i!=j)
      {
        it2 = first + j;

        grad = gradient_function(it->first);
        if(!grad.second)
          return std::make_pair(Value_type(0), false); // the gradient is not known

        //ordinates[i][j] = (p_j - p_i) * g_i
        ordinates[i][j] = grad.first *
          traits.construct_vector_d_object()(cp(it->first),cp(it2->first));

        // a point in the tangent plane:
        // 3( f(p_i) + (1/3)(p_j - p_i) * g_i)
        // => 3*f(p_i) + (p_j - p_i) * g_i
        res_i += (fac3 * ordinates[i][i] + ordinates[i][j])* it2->second;
      }
    }
    // res_i already multiplied by three
    result += coord_i_square *res_i;
  }

  // the third type of control points: three 1 entries i,j,k
  for(int i=0; i< n; ++i)
  {
    for(int j=i+1; j< n; ++j)
    {
      for(int k=j+1; k<n; ++k)
      {
        // add 6* (u_i*u_j*u_k) * b_ijk
        //  b_ijk = 1.5 * a - 0.5*c
        // where
        // c : average of the three data control points
        // a : 1.5*a = 1/12 * (ord[i][j] + ord[i][k] + ord[j][i] +
        //             ord[j][k] + ord[k][i]+ ord[k][j])
        // =>  6 * b_ijk = 3*(f_i + f_j + f_k) + 0.5*a
        result += (Coord_type(2.0)*(ordinates[i][i]+ ordinates[j][j]+
                                    ordinates[k][k])
                   + Coord_type(0.5)*(ordinates[i][j] + ordinates[i][k]
                                      + ordinates[j][i] +
                                      ordinates[j][k] + ordinates[k][i]+
                                      ordinates[k][j]))
                  * (first+i)->second *(first+j)->second *(first+k)->second ;
      }
    }
  }

  return std::make_pair(result/(CGAL_NTS square(norm)*norm), true);
}

} // namespace CGAL

#endif // CGAL_INTERPOLATION_FUNCTIONS_H
