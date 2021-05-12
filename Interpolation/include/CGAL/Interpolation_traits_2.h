// Copyright (c) 2003   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julia Floetotto

#ifndef CGAL_INTERPOLATION_TRAITS_2_H
#define CGAL_INTERPOLATION_TRAITS_2_H

#include <CGAL/license/Interpolation.h>

namespace CGAL {

//-----------------------------------------------------------------------//
//                          Interpolation_traits_2
//-----------------------------------------------------------------------//
// The class Interpolation_traits_2 is needed to define the geometric
// operations used by the interpolation methods.

template <class R>
class Interpolation_traits_2
{
public:
  typedef typename R::FT                     FT;
  typedef typename R::Point_2                Point_d;
  typedef typename R::Weighted_point_2       Weighted_point_d;
  typedef typename R::Vector_2               Vector_d;

  typedef typename R::Construct_point_2      Construct_point_d;
  typedef typename R::Construct_vector_2     Construct_vector_d;
  typedef typename R::Construct_scaled_vector_2
                                             Construct_scaled_vector_d;
  typedef typename R::Compute_squared_distance_2
                                             Compute_squared_distance_d;

  Construct_point_d
  construct_point_d_object()const
  {return Construct_point_d();}

  Construct_scaled_vector_d
  construct_scaled_vector_d_object()const
    {return Construct_scaled_vector_d();}

  Construct_vector_d
  construct_vector_d_object()const
    {return Construct_vector_d();}

  Compute_squared_distance_d
  compute_squared_distance_d_object()const
    {return Compute_squared_distance_d();}
};

} //namespace CGAL

#endif // CGAL_INTERPOLATION_TRAITS_2_H
