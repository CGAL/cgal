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
// Author(s)     : Julia Floetotto

#ifndef CGAL_INTERPOLATION_GRADIENT_FITTING_TRAITS_2_H
#define CGAL_INTERPOLATION_GRADIENT_FITTING_TRAITS_2_H

#include <CGAL/license/Interpolation.h>

#include <CGAL/aff_transformation_tags.h>

namespace CGAL {

//-----------------------------------------------------------------------//
//                          Interpolation_gradient_fitting_traits_2
//-----------------------------------------------------------------------//
// The class meets the requirement of the concept InterpolationTraits2
// and GradientFittingTraits2: it defines the geometric
// operations used by the interpolation methods and the gradient
// fitting function.

template <class Aff_2>
class Construct_sum_matrix_2
{
public:
  typedef Aff_2 Aff_transformation_2;

  Aff_transformation_2
  operator()(const Aff_transformation_2& tr1,
             const Aff_transformation_2& tr2) const
  {
    return Aff_transformation_2(tr1.m(0,0) + tr2.m(0,0),
                                tr1.m(0,1) + tr2.m(0,1),
                                tr1.m(0,2) + tr2.m(0,2),
                                tr1.m(1,0) + tr2.m(1,0),
                                tr1.m(1,1) + tr2.m(1,1),
                                tr1.m(1,2) + tr2.m(1,2));
  }
};

template <class Aff_2>
class Construct_null_matrix_2
{
public:
  typedef Aff_2 Aff_transformation_2;

  Aff_transformation_2
  operator()() const
  {
    return Aff_transformation_2(0,0,0,0,0,0);
  }
};

template <class Aff_2>
class Construct_scaling_matrix_2
{
public:
  typedef Aff_2 Aff_transformation_2;

  Aff_transformation_2
  operator()(const typename Aff_transformation_2::R::FT & scale) const
  {
    return Aff_transformation_2(SCALING, scale);
  }
};

template < class R >
class Construct_outer_product_2
{
public:
  typedef typename R::Aff_transformation_2       Aff_transformation_2;
  typedef typename R::Vector_2                   Vector_2;

  Aff_transformation_2
  operator()(const Vector_2& v) const
  {
    return Aff_transformation_2(v.x()*v.x(),
                                v.x()*v.y(), v.x()*v.y(),
                                v.y()*v.y());
  }
};

template <class R>
class Interpolation_gradient_fitting_traits_2
{
public:
  typedef R                                          Rep;

  typedef typename Rep::FT                           FT;
  typedef typename Rep::Point_2                      Point_d;
  typedef typename Rep::Weighted_point_2             Weighted_point_d;
  typedef typename Rep::Vector_2                     Vector_d;

  typedef typename Rep::Construct_point_2            Construct_point_d;
  typedef typename Rep::Construct_vector_2           Construct_vector_d;
  typedef typename Rep::Construct_scaled_vector_2    Construct_scaled_vector_d;
  //only one not needed by gradient fitting:
  typedef typename Rep::Compute_squared_distance_2   Compute_squared_distance_d;


  //additional types for gradient computation:
  typedef typename Rep::Aff_transformation_2         Aff_transformation_d;

  typedef Construct_null_matrix_2<Aff_transformation_d>
      Construct_null_matrix_d;
  typedef Construct_scaling_matrix_2<Aff_transformation_d>
      Construct_scaling_matrix_d;
  typedef Construct_sum_matrix_2<Aff_transformation_d> Construct_sum_matrix_d;
  typedef Construct_outer_product_2<Rep>             Construct_outer_product_d;


  Construct_outer_product_d
  construct_outer_product_d_object() const
  {return Construct_outer_product_d();}

  Construct_sum_matrix_d
  construct_sum_matrix_d_object() const
  {return Construct_sum_matrix_d();}

  Construct_scaling_matrix_d
  construct_scaling_matrix_d_object() const
  {return Construct_scaling_matrix_d();}

  Construct_null_matrix_d
  construct_null_matrix_d_object() const
  {return Construct_null_matrix_d();}

  //also in the traits without gradient computation:
  Construct_scaled_vector_d
  construct_scaled_vector_d_object()const
  {return Construct_scaled_vector_d();}

  Construct_point_d
  construct_point_d_object()const
  {return Construct_point_d();}

  Construct_vector_d
  construct_vector_d_object()const
  {return Construct_vector_d();}

  Compute_squared_distance_d
  compute_squared_distance_d_object()const
  {return Compute_squared_distance_d();}
};

} //namespace CGAL

#endif // CGAL_INTERPOLATION_GRADIENT_FITTING_TRAITS_2_H
